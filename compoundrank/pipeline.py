from __future__ import annotations

import csv
import json
from concurrent.futures import ThreadPoolExecutor
import tempfile
from pathlib import Path
from typing import Iterable

from rdkit import Chem

from .clustering import cluster_pose_hypotheses
from .compound_retrieval import run_compound_retrieval
from .export import write_complex_pdb
from .gnina import run_gnina_ensemble
from .homolog_search import DEFAULT_API_URL, run_homolog_search
from .interactions import summarize_interactions
from .ligand import (
    LigandEligibilityError,
    LigandRequest,
    assess_ligand_file,
    prepare_ligand,
    read_manifest,
)
from .models import LigandResult, PoseRecord
from .pocket import (
    build_pocket_definitions,
    write_pocket_definitions,
)
from .pocket_evidence import (
    load_pocket_evidence,
    score_pocket_biological_evidence,
    write_pocket_biological_evidence,
)
from .pocket_selection import (
    rank_pocket_attempts,
    summarize_pocket_attempt,
    write_pocket_selection_summary,
)
from .pose_recovery import (
    evaluate_scored_pose_sdf,
    write_scored_pose_outputs,
)
from .receptor import prepare_receptor
from .reference_evidence_workflow import (
    run_reference_evidence_workflow,
)
from .uniprot_acquisition import (
    fetch_uniprot_entry,
)
from .uncertainty import assess_uncertainty
from .validity import filter_poses_with_posebusters
from .run_report import write_run_report
from .subprocess_utils import CommandTimeoutError


def _top_score(result: LigandResult) -> float:
    if result.top_score is None:
        return float("-inf")
    return result.top_score


def _accepted_records_for_selected_pocket(
    *,
    valid_records_by_pocket: dict[
        str,
        list[PoseRecord],
    ],
    selected_pocket_id: str,
) -> list[PoseRecord]:
    """Return only accepted poses from the selected pocket."""

    return list(
        valid_records_by_pocket.get(
            selected_pocket_id,
            [],
        )
    )


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


def _preserve_posebusters_artifacts(
    *,
    validity_dir: Path,
    output_dir: Path,
    ligand_name: str,
    pocket_id: str,
) -> list[Path]:
    """Copy PoseBusters audit inputs and reports to final results."""

    artifact_names = (
        "posebusters_input.sdf",
        "posebusters_report.csv",
    )

    destination_dir = (
        output_dir
        / "posebusters_reports"
        / ligand_name
        / pocket_id
    )

    copied_paths: list[Path] = []

    for artifact_name in artifact_names:
        source_path = (
            validity_dir / artifact_name
        )

        if not source_path.is_file():
            continue

        destination_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        destination_path = (
            destination_dir / artifact_name
        )

        destination_path.write_bytes(
            source_path.read_bytes()
        )

        copied_paths.append(destination_path)

    return copied_paths


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



def _write_pose_records_sdf(
    records: list[PoseRecord],
    output_path: Path,
) -> Path:
    """Persist selected-pocket poses outside the temporary work tree."""
    if not records:
        raise ValueError(
            "No selected-pocket pose records were supplied."
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    writer = Chem.SDWriter(str(output_path))
    written = 0

    try:
        for record in records:
            molecule = Chem.Mol(record.molecule)

            molecule.SetProp(
                "_Name",
                (
                    f"{record.ligand_name}_"
                    f"{record.pocket_id}_"
                    f"seed_{record.seed}_"
                    f"pose_{record.pose_number}"
                ),
            )

            if record.cnn_score is not None:
                molecule.SetDoubleProp(
                    "CNNscore",
                    float(record.cnn_score),
                )

            if record.cnn_affinity is not None:
                molecule.SetDoubleProp(
                    "CNNaffinity",
                    float(record.cnn_affinity),
                )

            if record.minimized_affinity is not None:
                molecule.SetDoubleProp(
                    "minimizedAffinity",
                    float(record.minimized_affinity),
                )

            molecule.SetIntProp(
                "seed",
                int(record.seed),
            )
            molecule.SetIntProp(
                "pose_number",
                int(record.pose_number),
            )
            molecule.SetProp(
                "pocket_id",
                str(record.pocket_id),
            )
            molecule.SetIntProp(
                "pocket_rank",
                int(record.pocket_rank),
            )

            if record.fpocket_score is not None:
                molecule.SetDoubleProp(
                    "fpocket_score",
                    float(record.fpocket_score),
                )

            writer.write(molecule)
            written += 1
    finally:
        writer.close()

    if (
        written == 0
        or not output_path.is_file()
        or output_path.stat().st_size == 0
    ):
        raise RuntimeError(
            "Selected-pocket poses could not be persisted."
        )

    return output_path


def _run_selected_pocket_pose_recovery(
    *,
    reference_ligand: Path,
    records: list[PoseRecord],
    output_dir: Path,
    ligand_name: str,
    pocket_id: str,
    rmsd_threshold: float,
    autobox_ligand: Path | None = None,
) -> tuple[dict[str, object], dict[str, Path]]:
    """Evaluate the selected pocket only after ordinary selection."""
    pose_set_path = (
        output_dir
        / "pose_recovery_selected_pocket_poses.sdf"
    )

    _write_pose_records_sdf(
        records,
        pose_set_path,
    )

    summary = evaluate_scored_pose_sdf(
        reference_ligand=reference_ligand,
        poses_sdf=pose_set_path,
        rmsd_threshold=rmsd_threshold,
    )

    same_file_as_autobox = False

    if autobox_ligand is not None:
        try:
            same_file_as_autobox = (
                reference_ligand.resolve()
                == Path(autobox_ligand).resolve()
            )
        except OSError:
            same_file_as_autobox = (
                str(reference_ligand)
                == str(autobox_ligand)
            )

    summary.update(
        {
            "evaluated_compound": ligand_name,
            "evaluated_pocket_id": pocket_id,
            "evaluation_stage": (
                "after normal GNINA scoring, "
                "PoseBusters filtering, and "
                "pocket selection"
            ),
            "reference_ligand_used_for_posthoc_evaluation": True,
            "reference_ligand_used_for_pocket_selection": False,
            "reference_ligand_used_for_box_definition": (
                same_file_as_autobox
            ),
            "reference_ligand_used_for_docking": (
                same_file_as_autobox
            ),
            "reference_ligand_also_supplied_as_autobox_ligand": (
                same_file_as_autobox
            ),
            "autobox_ligand": (
                str(autobox_ligand)
                if autobox_ligand is not None
                else None
            ),
            "evaluated_pose_source": (
                "all raw GNINA poses from the normally "
                "selected pocket across configured seeds"
            ),
        }
    )

    outputs = write_scored_pose_outputs(
        summary,
        output_dir,
    )

    return summary, outputs


def _write_ligand_eligibility_report(
    output_dir: Path,
    rows: list[dict[str, object]],
) -> tuple[Path, Path]:
    """Write ligand eligibility decisions in CSV and JSON formats."""

    csv_path = output_dir / "ligand_eligibility_report.csv"
    json_path = output_dir / "ligand_eligibility_report.json"

    status_counts: dict[str, int] = {}

    for row in rows:
        status = str(
            row.get("eligibility_status")
            or "unknown"
        )
        status_counts[status] = (
            status_counts.get(status, 0) + 1
        )

    eligible_count = sum(
        bool(row.get("eligible"))
        for row in rows
    )

    payload = {
        "summary": {
            "evaluated_count": len(rows),
            "eligible_count": eligible_count,
            "excluded_count": (
                len(rows) - eligible_count
            ),
            "status_counts": status_counts,
        },
        "eligibility_configuration": {
            "workflow": (
                "standard_small_molecule_gnina"
            ),
            "decision_source": (
                "RDKit structure-based assessment"
            ),
            "note": (
                "Excluded compounds remain documented but "
                "are not sent to Open Babel, Meeko, GNINA, "
                "or PoseBusters."
            ),
        },
        "ligands": rows,
    }

    json_path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    fieldnames = [
        "compound",
        "source_type",
        "source_value",
        "eligible",
        "eligibility_status",
        "recommended_workflow",
        "molecular_weight",
        "heavy_atom_count",
        "rotatable_bond_count",
        "formal_charge",
        "amide_bond_count",
        "fragment_count",
        "eligibility_reasons",
        "error",
    ]

    with csv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in rows:
            csv_row = dict(row)

            reasons = csv_row.get(
                "eligibility_reasons",
                [],
            )

            if isinstance(reasons, list):
                csv_row["eligibility_reasons"] = (
                    "; ".join(
                        str(reason)
                        for reason in reasons
                    )
                )

            writer.writerow(csv_row)

    return csv_path, json_path


def _resolve_pocket_evidence_json(
    *,
    pocket_evidence_json: Path | None,
    auto_reference_evidence: bool,
    reference_uniprot_accession: str | None,
    reference_uniprot_json: Path | None,
    reference_pdb_id: str | None,
    reference_chain_id: str | None,
    receptor_chain_id: str | None,
    reference_pdb: Path | None,
    reference_evidence_timeout_seconds: float,
    fasta_path: Path | None,
    receptor_pdb: Path,
    output_dir: Path,
) -> Path | None:
    if not auto_reference_evidence:
        return (
            Path(pocket_evidence_json)
            if pocket_evidence_json is not None
            else None
        )

    if pocket_evidence_json is not None:
        raise ValueError(
            "Automatic reference evidence cannot be "
            "combined with manual pocket evidence."
        )

    if fasta_path is None:
        raise ValueError(
            "Automatic reference evidence requires "
            "a submitted FASTA file."
        )

    source_count = sum(
        value is not None
        for value in (
            reference_uniprot_accession,
            reference_uniprot_json,
        )
    )

    if source_count != 1:
        raise ValueError(
            "Automatic reference evidence requires "
            "exactly one UniProt accession or "
            "UniProt JSON file."
        )

    if reference_evidence_timeout_seconds <= 0:
        raise ValueError(
            "Reference-evidence timeout must be "
            "greater than zero."
        )

    if reference_uniprot_json is not None:
        source_path = Path(
            reference_uniprot_json
        )

        payload = json.loads(
            source_path.read_text(
                encoding="utf-8"
            )
        )

        response_metadata = {
            "source": "input_json",
            "input_path": str(
                source_path
            ),
        }
    else:
        payload, response_metadata = (
            fetch_uniprot_entry(
                str(
                    reference_uniprot_accession
                ),
                timeout_seconds=(
                    reference_evidence_timeout_seconds
                ),
            )
        )

    workflow_output = (
        output_dir
        / "automatic_reference_evidence"
    )

    print(
        "[REFERENCE EVIDENCE] Generating "
        "automatic biological pocket evidence"
    )

    result = run_reference_evidence_workflow(
        payload=payload,
        submitted_fasta=Path(
            fasta_path
        ),
        receptor_pdb=Path(
            receptor_pdb
        ),
        output_dir=workflow_output,
        response_metadata=(
            response_metadata
        ),
        requested_pdb_id=(
            reference_pdb_id
        ),
        reference_chain_id=(
            reference_chain_id
        ),
        receptor_chain_id=(
            receptor_chain_id
        ),
        reference_pdb=(
            Path(reference_pdb)
            if reference_pdb is not None
            else None
        ),
        timeout_seconds=(
            reference_evidence_timeout_seconds
        ),
    )

    evidence_path = Path(
        result[
            "pocket_evidence_path"
        ]
    )

    print(
        "[REFERENCE EVIDENCE] "
        f"Generated: {evidence_path}"
    )

    return evidence_path


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
    fpocket_merge_nearby: bool = False,
    fpocket_merge_distance: float = 4.0,
    pocket_evidence_json: Path | None = None,
    auto_reference_evidence: bool = False,
    reference_uniprot_accession: str | None = None,
    reference_uniprot_json: Path | None = None,
    reference_pdb_id: str | None = None,
    reference_chain_id: str | None = None,
    receptor_chain_id: str | None = None,
    reference_pdb: Path | None = None,
    reference_evidence_timeout_seconds: float = 60.0,
    gnina_timeout_seconds: int | None = 3600,
    fasta_path: Path | None = None,
    homolog_api_url: str = DEFAULT_API_URL,
    homolog_timeout_seconds: int = 7200,
    auto_retrieve_ligands: bool = False,
    auto_retrieve_mode: str = "rules-only",
    auto_retrieve_max_candidates: int = 20,
    auto_retrieve_fetch_structures: bool = True,
    auto_retrieve_pubchem_timeout_seconds: int = 60,
    reference_ligand: Path | None = None,
    pose_recovery_rmsd_threshold: float = 2.0,
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

        for filename in (
            "pocket_biological_evidence.json",
            "pose_recovery_selected_pocket_poses.sdf",
            "pose_set_recovery_summary.json",
            "pose_set_recovery_metrics.csv",
            "pose_set_recovery_report.md",
        ):
            stale_path = output_dir / filename

            if stale_path.exists():
                stale_path.unlink()

    effective_pocket_evidence_json = (
        _resolve_pocket_evidence_json(
            pocket_evidence_json=(
                pocket_evidence_json
            ),
            auto_reference_evidence=(
                auto_reference_evidence
            ),
            reference_uniprot_accession=(
                reference_uniprot_accession
            ),
            reference_uniprot_json=(
                reference_uniprot_json
            ),
            reference_pdb_id=(
                reference_pdb_id
            ),
            reference_chain_id=(
                reference_chain_id
            ),
            receptor_chain_id=(
                receptor_chain_id
            ),
            reference_pdb=(
                reference_pdb
            ),
            reference_evidence_timeout_seconds=(
                reference_evidence_timeout_seconds
            ),
            fasta_path=fasta_path,
            receptor_pdb=receptor_pdb,
            output_dir=output_dir,
        )
    )

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
        print(
            f"[COMPOUND_RETRIEVAL] Retrieval mode: "
            f"{auto_retrieve_mode}"
        )
        retrieval_outputs = run_compound_retrieval(
            target_evidence_path=Path(target_evidence_output),
            fasta_path=Path(fasta_path),
            homolog_summary_path=Path(homolog_summary_output) if homolog_summary_output else None,
            output_dir=retrieval_dir,
            max_candidates=auto_retrieve_max_candidates,
            fetch_structures=auto_retrieve_fetch_structures,
            pubchem_timeout_seconds=auto_retrieve_pubchem_timeout_seconds,
            retrieval_mode=auto_retrieve_mode,
        )

        generated_manifest = retrieval_outputs["docking_manifest"]
        ligand_requests = read_manifest(generated_manifest)

        if not ligand_requests:
            if auto_retrieve_mode == "generic-strict":
                skip_reason = (
                    "No externally supported candidates met the "
                    "automatic docking threshold. Exploratory evidence "
                    "remains available for manual review."
                )
            else:
                skip_reason = (
                    "Automatic ligand retrieval produced no candidates "
                    "that were approved for docking and associated with "
                    "a usable molecular structure."
                )

            skip_payload = {
                "stage": "docking",
                "status": "skipped",
                "pipeline_outcome": (
                    "completed_without_docking"
                ),
                "reason_code": (
                    "no_dockable_ligands"
                ),
                "reason": skip_reason,
                "retrieval_mode": (
                    auto_retrieve_mode
                ),
                "dockable_ligand_count": 0,
                "docking_manifest": str(
                    generated_manifest
                ),
                "retrieval_outputs": {
                    str(name): str(value)
                    for name, value
                    in retrieval_outputs.items()
                },
                "downstream_stages": {
                    "receptor_preparation": (
                        "skipped"
                    ),
                    "pocket_detection": (
                        "skipped"
                    ),
                    "ligand_preparation": (
                        "skipped"
                    ),
                    "gnina_docking": "skipped",
                    "pose_validation": (
                        "skipped"
                    ),
                    "pose_ranking": "skipped",
                },
            }

            skip_path = (
                output_dir
                / "docking_skipped.json"
            )

            skip_path.write_text(
                json.dumps(
                    skip_payload,
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )

            print(
                "[DOCKING SKIPPED] No ligand "
                "passed the automatic evidence "
                "and structure gates."
            )
            print(
                "[DOCKING SKIPPED] Reason: "
                f"{skip_reason}"
            )
            print(
                "[DOCKING SKIPPED] Manifest: "
                f"{generated_manifest}"
            )
            print(
                "[DOCKING SKIPPED] Status: "
                f"{skip_path}"
            )

            run_report_path = write_run_report(
                output_dir=output_dir,
            )

            print(
                "\n[REPORT] Run report: "
                f"{run_report_path}"
            )

            return []

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

    if pose_recovery_rmsd_threshold <= 0:
        raise ValueError(
            "pose_recovery_rmsd_threshold must be greater than zero"
        )

    if reference_ligand is not None:
        reference_ligand = Path(reference_ligand)

        if (
            not reference_ligand.is_file()
            or reference_ligand.stat().st_size == 0
        ):
            raise FileNotFoundError(
                "Reference ligand was not found or is empty: "
                f"{reference_ligand}"
            )

        if len(ligand_requests) != 1:
            raise ValueError(
                "--reference-ligand currently requires exactly one "
                "docked ligand"
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

        # Detect geometric pockets on the original protein coordinates.
        # The protonated display structure is retained for validation/output,
        # while GNINA uses the prepared PDBQT receptor.
        pockets = build_pocket_definitions(
            receptor_pdb=receptor.source_pdb,
            work_dir=work_dir / "pocket",
            explicit_values=(center_x, center_y, center_z, size_x, size_y, size_z),
            autobox_ligand=autobox_ligand,
            fpocket_padding=fpocket_padding,
            fpocket_pocket=fpocket_pocket,
            fpocket_top_n=fpocket_top_n,
            fpocket_bin=fpocket_bin,
            fpocket_merge_nearby=(
                fpocket_merge_nearby
            ),
            fpocket_merge_distance=(
                fpocket_merge_distance
            ),
        )

        print(f"[POCKET] Testing {len(pockets)} pocket definition(s)")
        for pocket in pockets:
            print(f"[POCKET] {pocket.pocket_id}: {pocket.source or pocket.mode}")

        pocket_definitions_path = write_pocket_definitions(
            output_dir / "pocket_definitions.json",
            pockets,
        )
        print(
            f"[POCKET] Definitions: "
            f"{pocket_definitions_path}"
        )

        pocket_biological_scores: dict[
            str,
            dict[str, object],
        ] = {}

        if effective_pocket_evidence_json is not None:
            pocket_evidence = (
                load_pocket_evidence(
                    Path(
                        effective_pocket_evidence_json
                    )
                )
            )

            (
                pocket_biological_scores,
                pocket_biological_report,
            ) = score_pocket_biological_evidence(
                pockets,
                fpocket_output_dir=(
                    work_dir
                    / "pocket"
                    / "fpocket_receptor_out"
                ),
                evidence=pocket_evidence,
            )

            biological_report_path = (
                write_pocket_biological_evidence(
                    output_dir
                    / "pocket_biological_evidence.json",
                    pocket_biological_report,
                )
            )

            print(
                "[POCKET EVIDENCE] "
                f"Mode={pocket_evidence['selection_mode']}; "
                f"origin={pocket_evidence['evidence_origin']}; "
                f"residues="
                f"{len(pocket_evidence['residues'])}"
            )
            print(
                "[POCKET EVIDENCE] Report: "
                f"{biological_report_path}"
            )

        ligand_results: list[LigandResult] = []
        docking_attempt_rows: list[dict[str, object]] = []
        pocket_selection_rows: list[dict[str, object]] = []
        ligand_eligibility_rows: list[dict[str, object]] = []

        for request in ligand_requests:
            print(f"\n[LIGAND] Preparing {request.name}")

            try:
                ligand = prepare_ligand(
                    request,
                    cache_root,
                    ph=ph,
                    obabel_bin=obabel_bin,
                    meeko_ligand_bin=meeko_ligand_bin,
                )

            except LigandEligibilityError as error:
                assessment = dict(error.assessment)

                ligand_eligibility_rows.append(
                    {
                        "compound": request.name,
                        "source_type": request.source_type,
                        "source_value": request.value,
                        **assessment,
                        "error": str(error),
                    }
                )

                _write_ligand_eligibility_report(
                    output_dir,
                    ligand_eligibility_rows,
                )

                print(
                    "[LIGAND EXCLUDED] "
                    f"{request.name}: "
                    f"{assessment.get(
                        'eligibility_status',
                        'excluded',
                    )}"
                )

                for reason in assessment.get(
                    "eligibility_reasons",
                    [],
                ):
                    print(
                        "[LIGAND EXCLUDED] Reason: "
                        f"{reason}"
                    )

                print(
                    "[LIGAND EXCLUDED] "
                    "Recommended workflow: "
                    f"{assessment.get(
                        'recommended_workflow',
                        'manual_review',
                    )}"
                )

                continue

            assessment = assess_ligand_file(
                ligand.source_sdf
            )

            ligand_eligibility_rows.append(
                {
                    "compound": request.name,
                    "source_type": request.source_type,
                    "source_value": request.value,
                    **assessment,
                    "error": "",
                }
            )

            # Update the report throughout long runs so it
            # survives later failures or manual termination.
            _write_ligand_eligibility_report(
                output_dir,
                ligand_eligibility_rows,
            )

            all_valid_records = []
            total_failures = 0
            ligand_pocket_rows: list[dict[str, object]] = []
            raw_records_by_pocket: dict[
                str,
                list[PoseRecord],
            ] = {}

            valid_records_by_pocket: dict[
                str,
                list[PoseRecord],
            ] = {}

            records_for_hypotheses: list[
                PoseRecord
            ] = []

            for pocket in pockets:
                print(f"\n[POCKET RUN] {ligand.name}: {pocket.pocket_id}")

                try:
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
                        timeout_seconds=gnina_timeout_seconds,
                    )

                except CommandTimeoutError as error:
                    print(
                        "[GNINA TIMEOUT] "
                        f"{ligand.name} "
                        f"{pocket.pocket_id}: exceeded "
                        f"{error.timeout_seconds:g} seconds"
                    )

                    docking_attempt_rows.append(
                        {
                            "compound": ligand.name,
                            "pocket": pocket.pocket_id,
                            "raw_poses": 0,
                            "accepted_poses": 0,
                            "rejected_poses": 0,
                            "status": "timed_out",
                            "best_raw_cnn_score": "",
                            "best_accepted_cnn_score": "",
                        }
                    )

                    print(
                        "[GNINA TIMEOUT] Continuing to "
                        "the next pocket."
                    )

                    continue

                raw_records_by_pocket[pocket.pocket_id] = list(
                    raw_records
                )

                validity_output_dir = (
                    work_dir
                    / "validity"
                    / ligand.name
                    / pocket.pocket_id
                )

                try:
                    valid_records, failures = filter_poses_with_posebusters(
                        raw_records,
                        receptor.display_pdb,
                        validity_output_dir,
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

                preserved_validity_artifacts = (
                    _preserve_posebusters_artifacts(
                        validity_dir=validity_output_dir,
                        output_dir=output_dir,
                        ligand_name=ligand.name,
                        pocket_id=pocket.pocket_id,
                    )
                )

                for artifact_path in (
                    preserved_validity_artifacts
                ):
                    print(
                        "[VALIDITY] Preserved audit "
                        f"artifact: {artifact_path}"
                    )

                total_failures += len(failures)
                all_valid_records.extend(valid_records)

                valid_records_by_pocket[
                    pocket.pocket_id
                ] = list(valid_records)

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

                pocket_attempt_row = (
                    summarize_pocket_attempt(
                        ligand_name=ligand.name,
                        pocket=pocket,
                        raw_records=raw_records,
                        accepted_records=valid_records,
                        rejected_pose_count=len(failures),
                    )
                )

                pocket_attempt_row.update(
                    pocket_biological_scores.get(
                        pocket.pocket_id,
                        {},
                    )
                )

                ligand_pocket_rows.append(
                    pocket_attempt_row
                )

            ranked_pocket_rows = rank_pocket_attempts(
                ligand_pocket_rows
            )
            pocket_selection_rows.extend(
                ranked_pocket_rows
            )

            if ranked_pocket_rows:
                selected_pocket = ranked_pocket_rows[0]

                selected_pocket_id = str(
                    selected_pocket["pocket_id"]
                )

                records_for_hypotheses = (
                    _accepted_records_for_selected_pocket(
                        valid_records_by_pocket=(
                            valid_records_by_pocket
                        ),
                        selected_pocket_id=(
                            selected_pocket_id
                        ),
                    )
                )

                print(
                    "[POCKET SELECTION] "
                    f"{ligand.name}: "
                    f"{selected_pocket['pocket_id']}; "
                    f"fpocket rank="
                    f"{selected_pocket['pocket_rank']}; "
                    f"CNNscore="
                    f"{selected_pocket['top_cnn_score']}; "
                    f"score source="
                    f"{selected_pocket['score_source']}"
                )

                if reference_ligand is not None:
                    selected_records = (
                        raw_records_by_pocket.get(
                            selected_pocket_id,
                            [],
                        )
                    )

                    if not selected_records:
                        raise RuntimeError(
                            "The selected pocket's raw GNINA "
                            "records could not be located."
                        )

                    pose_summary, pose_outputs = (
                        _run_selected_pocket_pose_recovery(
                            reference_ligand=reference_ligand,
                            records=selected_records,
                            output_dir=output_dir,
                            ligand_name=ligand.name,
                            pocket_id=selected_pocket_id,
                            rmsd_threshold=(
                                pose_recovery_rmsd_threshold
                            ),
                            autobox_ligand=(
                                autobox_ligand
                            ),
                        )
                    )

                    print(
                        "\n[POSE_RECOVERY] "
                        f"Evaluated compound: {ligand.name}"
                    )
                    print(
                        "[POSE_RECOVERY] Normally selected pocket: "
                        f"{selected_pocket_id}"
                    )
                    print(
                        "[POSE_RECOVERY] Top CNN pose RMSD: "
                        f"{pose_summary['top_cnn_pose']['heavy_atom_rmsd']:.3f} Å"
                    )
                    print(
                        "[POSE_RECOVERY] Best sampled RMSD: "
                        f"{pose_summary['best_sampled_pose']['heavy_atom_rmsd']:.3f} Å"
                    )
                    print(
                        "[POSE_RECOVERY] Sampling pass: "
                        f"{pose_summary['sampling_pass']}"
                    )
                    print(
                        "[POSE_RECOVERY] Ranking pass: "
                        f"{pose_summary['ranking_pass']}"
                    )

                    for label, pose_path in pose_outputs.items():
                        print(
                            f"[POSE_RECOVERY] {label}: "
                            f"{pose_path}"
                        )

            print(f"[VALIDITY] Failure records: {total_failures}")

            clusters = cluster_pose_hypotheses(
                records_for_hypotheses,
                rmsd_threshold=cluster_threshold,
            )

            uncertainty, reasons = assess_uncertainty(
                clusters,
                len(seeds),
            )

            top_score = (
                clusters[0].representative.cnn_score
                if clusters
                else None
            )

            if clusters:
                ligand_results.append(
                    LigandResult(
                        ligand=ligand,
                        clusters=clusters,
                        uncertainty=uncertainty,
                        uncertainty_reasons=reasons,
                        top_score=top_score,
                    )
                )
            else:
                print(
                    "[RANKING EXCLUDED] "
                    f"{ligand.name}: no valid pose clusters "
                    "were produced."
                )

            print(
                f"[HYPOTHESES] {ligand.name}: "
                f"{len(clusters)} clusters; "
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

        pocket_selection_paths = (
            write_pocket_selection_summary(
                output_dir,
                pocket_selection_rows,
            )
        )

        if pocket_selection_paths is not None:
            selection_csv, selection_json = (
                pocket_selection_paths
            )
            print(
                "[REPORT] Pocket selection CSV: "
                f"{selection_csv}"
            )
            print(
                "[REPORT] Pocket selection JSON: "
                f"{selection_json}"
            )

        eligibility_csv, eligibility_json = (
            _write_ligand_eligibility_report(
                output_dir,
                ligand_eligibility_rows,
            )
        )

        print(
            "[REPORT] Ligand eligibility CSV: "
            f"{eligibility_csv}"
        )
        print(
            "[REPORT] Ligand eligibility JSON: "
            f"{eligibility_json}"
        )

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
