from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
import tempfile
from pathlib import Path
from typing import Iterable

from .clustering import cluster_pose_hypotheses
from .compound_retrieval import run_compound_retrieval
from .export import write_complex_pdb
from .gnina import run_gnina_ensemble
from .homolog_search import DEFAULT_API_URL, run_homolog_search
from .interactions import summarize_interactions
from .ligand import LigandRequest, prepare_ligand, read_manifest
from .models import LigandResult
from .pocket import build_pocket_definitions
from .receptor import prepare_receptor
from .uncertainty import assess_uncertainty
from .validity import filter_poses_with_posebusters
from .run_report import write_run_report


def _top_score(result: LigandResult) -> float:
    if result.top_score is None:
        return float("-inf")
    return result.top_score


def _output_name(
    *,
    compound_rank: int,
    ligand_name: str,
    pocket_id: str,
    hypothesis_rank: int,
    multi_pocket: bool,
) -> str:
    if multi_pocket:
        return (
            f"{compound_rank:02d}__{ligand_name}__{pocket_id}__"
            f"hypothesis_{hypothesis_rank:02d}.pdb"
        )
    return f"{compound_rank:02d}__{ligand_name}__hypothesis_{hypothesis_rank:02d}.pdb"



def _best_cnn_score_text(records: list[object]) -> str:
    scores: list[float] = []

    for record in records:
        value = getattr(record, "cnn_score", None)
        if isinstance(value, (int, float)):
            scores.append(float(value))

    if not scores:
        return ""

    return f"{max(scores):.9f}"


def _write_docking_attempt_summary(
    output_dir: Path,
    rows: list[dict[str, object]],
) -> Path | None:
    if not rows:
        return None

    output_path = output_dir / "docking_attempt_summary.csv"
    fieldnames = [
        "compound",
        "pocket",
        "raw_poses",
        "accepted_poses",
        "rejected_poses",
        "status",
        "best_raw_cnn_score",
        "best_accepted_cnn_score",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    return output_path


def run_pipeline(
    *,
    receptor_pdb: Path,
    ligand_requests: Iterable[LigandRequest],
    data_root: Path,
    output_dir: Path,
    seeds: list[int],
    center_x: float | None,
    center_y: float | None,
    center_z: float | None,
    size_x: float | None,
    size_y: float | None,
    size_z: float | None,
    autobox_ligand: Path | None,
    fpocket_padding: float,
    fpocket_pocket: int | None,
    fpocket_top_n: int,
    max_hypotheses: int,
    cluster_threshold: float,
    exhaustiveness: int,
    num_modes: int,
    cnn_scoring: str,
    ph: float,
    gnina_bin: str,
    fpocket_bin: str,
    obabel_bin: str,
    pdb2pqr_bin: str,
    meeko_receptor_bin: str,
    meeko_ligand_bin: str,
    posebusters_bin: str,
    skip_validity: bool,
    keep_workdir: bool,
    overwrite: bool,
    cpu: int | None,
    device: int | None,
    fasta_path: Path | None = None,
    homolog_api_url: str = DEFAULT_API_URL,
    homolog_timeout_seconds: int = 7200,
    auto_retrieve_ligands: bool = False,
    auto_retrieve_max_candidates: int = 20,
    auto_retrieve_fetch_structures: bool = True,
    auto_retrieve_pubchem_timeout_seconds: int = 60,
) -> list[Path]:
    cache_root = data_root / "cache"
    work_root = data_root / "work"

    cache_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    ligand_requests = list(ligand_requests)

    if auto_retrieve_ligands:
        if fasta_path is None:
            raise ValueError("--auto-retrieve-ligands requires --fasta")
        if ligand_requests:
            raise ValueError(
                "--auto-retrieve-ligands should not be combined with manual ligand inputs"
            )

    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Output directory is not empty: {output_dir}. Use --overwrite."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    if overwrite:
        for path in output_dir.glob("*.pdb"):
            path.unlink()
        for path in output_dir.glob("homolog_search_*.json"):
            path.unlink()

    homology_executor: ThreadPoolExecutor | None = None
    homology_future = None

    if fasta_path is not None and auto_retrieve_ligands:
        print(f"[HOMOLOGY] Running CPU homolog search before ligand retrieval: {fasta_path}")
        homology_result = run_homolog_search(
            fasta_path=Path(fasta_path),
            output_dir=output_dir,
            api_url=homolog_api_url,
            timeout_seconds=homolog_timeout_seconds,
        )

        if homology_result.get("status") != "ok":
            raise RuntimeError(
                "Auto ligand retrieval requires successful target evidence. "
                f"Homology error: {homology_result.get('error')}"
            )

        target_evidence_output = homology_result.get("target_evidence")
        if not target_evidence_output:
            raise RuntimeError("Homology completed but did not return target_evidence path")

        homolog_summary_output = homology_result.get("summary_output")
        retrieval_dir = output_dir / "stage4a_compound_retrieval"

        print(f"[COMPOUND_RETRIEVAL] Running Stage 4A into: {retrieval_dir}")
        retrieval_outputs = run_compound_retrieval(
            target_evidence_path=Path(target_evidence_output),
            homolog_summary_path=Path(homolog_summary_output) if homolog_summary_output else None,
            output_dir=retrieval_dir,
            max_candidates=auto_retrieve_max_candidates,
            fetch_structures=auto_retrieve_fetch_structures,
            pubchem_timeout_seconds=auto_retrieve_pubchem_timeout_seconds,
        )

        generated_manifest = retrieval_outputs["docking_manifest"]
        ligand_requests = read_manifest(generated_manifest)

        if not ligand_requests:
            raise RuntimeError(
                "Auto ligand retrieval produced no dockable ligand requests. "
                f"Check {generated_manifest}"
            )

        print(f"[COMPOUND_RETRIEVAL] Docking manifest: {generated_manifest}")
        print(f"[COMPOUND_RETRIEVAL] Dockable ligands: {len(ligand_requests)}")

    elif fasta_path is not None:
        print(f"[HOMOLOGY] Starting CPU homolog search in background: {fasta_path}")
        homology_executor = ThreadPoolExecutor(
            max_workers=1,
            thread_name_prefix="compoundrank-homology",
        )
        homology_future = homology_executor.submit(
            run_homolog_search,
            fasta_path=Path(fasta_path),
            output_dir=output_dir,
            api_url=homolog_api_url,
            timeout_seconds=homolog_timeout_seconds,
        )

    temporary: tempfile.TemporaryDirectory[str] | None = None

    if keep_workdir:
        work_dir = Path(tempfile.mkdtemp(prefix="compoundrank-", dir=work_root))
        cleanup = False
    else:
        temporary = tempfile.TemporaryDirectory(prefix="compoundrank-", dir=work_root)
        work_dir = Path(temporary.name)
        cleanup = True

    print(f"[PATHS] Work directory: {work_dir}")
    print(f"[PATHS] Cache directory: {cache_root}")
    print(f"[PATHS] Final results: {output_dir}")

    try:
        receptor = prepare_receptor(
            receptor_pdb,
            cache_root,
            ph=ph,
            pdb2pqr_bin=pdb2pqr_bin,
            meeko_receptor_bin=meeko_receptor_bin,
        )

        pockets = build_pocket_definitions(
            receptor_pdb=receptor.display_pdb,
            work_dir=work_dir / "pocket",
            explicit_values=(center_x, center_y, center_z, size_x, size_y, size_z),
            autobox_ligand=autobox_ligand,
            fpocket_padding=fpocket_padding,
            fpocket_pocket=fpocket_pocket,
            fpocket_top_n=fpocket_top_n,
            fpocket_bin=fpocket_bin,
        )

        print(f"[POCKET] Testing {len(pockets)} pocket definition(s)")
        for pocket in pockets:
            print(f"[POCKET] {pocket.pocket_id}: {pocket.source or pocket.mode}")

        ligand_results: list[LigandResult] = []
        docking_attempt_rows: list[dict[str, object]] = []

        for request in ligand_requests:
            print(f"\n[LIGAND] Preparing {request.name}")
            ligand = prepare_ligand(
                request,
                cache_root,
                ph=ph,
                obabel_bin=obabel_bin,
                meeko_ligand_bin=meeko_ligand_bin,
            )

            all_valid_records = []
            total_failures = 0

            for pocket in pockets:
                print(f"\n[POCKET RUN] {ligand.name}: {pocket.pocket_id}")

                raw_records = run_gnina_ensemble(
                    receptor,
                    ligand,
                    pocket,
                    seeds,
                    work_dir / "docking",
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes,
                    cnn_scoring=cnn_scoring,
                    gnina_bin=gnina_bin,
                    cpu=cpu,
                    device=device,
                )

                try:
                    valid_records, failures = filter_poses_with_posebusters(
                        raw_records,
                        receptor.display_pdb,
                        work_dir / "validity" / ligand.name / pocket.pocket_id,
                        posebusters_bin=posebusters_bin,
                        skip=skip_validity,
                    )
                except RuntimeError as error:
                    if "Every GNINA pose failed" not in str(error):
                        raise

                    valid_records = []
                    failures = list(raw_records)

                    print(
                        f"[VALIDITY] {ligand.name} {pocket.pocket_id}: "
                        "accepted 0/"
                        f"{len(raw_records)}; rejected {len(raw_records)}; "
                        "skipping this exploratory pocket"
                    )

                total_failures += len(failures)
                all_valid_records.extend(valid_records)

                print(
                    f"[VALIDITY] {ligand.name} {pocket.pocket_id}: "
                    f"accepted {len(valid_records)}/{len(raw_records)}; "
                    f"rejected {len(failures)}"
                )

                docking_attempt_rows.append(
                    {
                        "compound": ligand.name,
                        "pocket": pocket.pocket_id,
                        "raw_poses": len(raw_records),
                        "accepted_poses": len(valid_records),
                        "rejected_poses": len(failures),
                        "status": "accepted" if valid_records else "failed_posebusters",
                        "best_raw_cnn_score": _best_cnn_score_text(raw_records),
                        "best_accepted_cnn_score": _best_cnn_score_text(valid_records),
                    }
                )

            print(f"[VALIDITY] Failure records: {total_failures}")

            clusters = cluster_pose_hypotheses(
                all_valid_records,
                rmsd_threshold=cluster_threshold,
            )

            uncertainty, reasons = assess_uncertainty(clusters, len(seeds))

            top_score = clusters[0].representative.cnn_score if clusters else None

            ligand_results.append(
                LigandResult(
                    ligand=ligand,
                    clusters=clusters,
                    uncertainty=uncertainty,
                    uncertainty_reasons=reasons,
                    top_score=top_score,
                )
            )

            print(
                f"[HYPOTHESES] {ligand.name}: {len(clusters)} clusters; "
                f"uncertainty={uncertainty}"
            )

        ligand_results.sort(key=_top_score, reverse=True)

        written: list[Path] = []
        multi_pocket = len(pockets) > 1

        for compound_rank, result in enumerate(ligand_results, start=1):
            selected_clusters = result.clusters[:max_hypotheses]

            for hypothesis_rank, cluster in enumerate(selected_clusters, start=1):
                evidence = summarize_interactions(
                    receptor.display_pdb,
                    cluster.representative.molecule,
                )

                output_path = output_dir / _output_name(
                    compound_rank=compound_rank,
                    ligand_name=result.ligand.name,
                    pocket_id=cluster.representative.pocket_id,
                    hypothesis_rank=hypothesis_rank,
                    multi_pocket=multi_pocket,
                )

                write_complex_pdb(
                    output_path,
                    receptor.display_pdb,
                    result,
                    cluster,
                    evidence,
                    compound_priority_rank=compound_rank,
                    hypothesis_rank=hypothesis_rank,
                    hypothesis_count=len(selected_clusters),
                )

                written.append(output_path)

        print("\n=== COMPOUND PRIORITY BY GNINA CNN SCORE ===")

        for rank, result in enumerate(ligand_results, start=1):
            print(
                f"{rank:2d}. {result.ligand.name}: "
                f"CNNscore={result.top_score if result.top_score is not None else 'N/A'}; "
                f"uncertainty={result.uncertainty}; "
                f"clusters={len(result.clusters)}"
            )

        print("\nFinal PDB files:")
        for path in written:
            print(path)

        attempt_summary_path = _write_docking_attempt_summary(
            output_dir,
            docking_attempt_rows,
        )
        if attempt_summary_path is not None:
            print(f"[REPORT] Docking attempt summary: {attempt_summary_path}")

        if homology_future is not None:
            print("\n[HOMOLOGY] Waiting for CPU homolog search to finish")
            homology_result = homology_future.result()

            if homology_result.get("status") == "ok":
                print("[HOMOLOGY] Completed successfully")
                print(f"[HOMOLOGY] Summary: {homology_result.get('summary_output')}")

                target_evidence_output = homology_result.get("target_evidence")
                if target_evidence_output:
                    print(f"[HOMOLOGY] Target evidence: {target_evidence_output}")

                target_evidence_report = homology_result.get("target_evidence_report")
                if target_evidence_report:
                    print(f"[HOMOLOGY] Target evidence report: {target_evidence_report}")

                print(f"[HOMOLOGY] Counts: {homology_result.get('result_counts')}")
            else:
                print("[HOMOLOGY] Failed; docking outputs were still preserved")
                print(f"[HOMOLOGY] Error file: {homology_result.get('error_output')}")
                print(f"[HOMOLOGY] Error: {homology_result.get('error')}")

        run_report_path = write_run_report(
            output_dir=output_dir,
        )
        print(f"\n[REPORT] Run report: {run_report_path}")

        return written

    finally:
        if homology_executor is not None:
            homology_executor.shutdown(
                wait=False,
                cancel_futures=True,
            )

        if cleanup and temporary is not None:
            temporary.cleanup()
        else:
            print(f"[DEBUG] Work directory retained: {work_dir}")
