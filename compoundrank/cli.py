from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from .ligand import combine_requests, read_manifest
from .paths import (
    require_absolute_external_dir,
    require_absolute_external_file,
    sanitize_name,
)
from .aligned_receptor_ensemble import (
    validate_receptor_ensemble_options,
)
from .homolog_search import DEFAULT_API_URL
from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "CompoundRank: GNINA-scored docking hypotheses with physical-validity "
            "filtering, clustering, interaction evidence, and uncertainty."
        )
    )
    parser.add_argument(
        "--receptor",
        default=None,
        help=(
            "Optional absolute receptor PDB path. If omitted with --fasta, "
            "CompoundRank enters pure-FASTA mode and attempts receptor "
            "structure prediction/acquisition before docking."
        ),
    )
    parser.add_argument(
        "--skip-structure-prediction",
        action="store_true",
        help=(
            "Pure-FASTA mode only: do not run a structure predictor. "
            "Write no_receptor_structure_available artifacts instead."
        ),
    )
    parser.add_argument(
        "--structure-predictor-bin",
        default=os.environ.get(
            "COMPOUNDRANK_STRUCTURE_PREDICTOR_BIN",
            "colabfold_batch",
        ),
        help=(
            "Executable used for pure-FASTA receptor prediction. "
            "Default: colabfold_batch or COMPOUNDRANK_STRUCTURE_PREDICTOR_BIN."
        ),
    )
    parser.add_argument(
        "--structure-predictor-extra-arg",
        action="append",
        default=[],
        help=(
            "Extra argument passed to the structure predictor. "
            "Repeat this option for multiple arguments."
        ),
    )
    parser.add_argument(
        "--receptor-ensemble-json",
        default=None,
        help=(
            "Optional absolute path to a "
            "structure_ensemble.v0.1 "
            "manifest. The current "
            "integration validates and "
            "records the ensemble but does "
            "not change docking behavior."
        ),
    )
    parser.add_argument(
        "--aligned-receptor-ensemble-json",
        default=None,
        help=(
            "Optional absolute path to an "
            "aligned_receptor_ensemble.v0.1 "
            "manifest. Accepted conformers are "
            "validated, receptor-prepared, and "
            "docked independently using shared "
            "aligned pocket coordinates."
        ),
    )
    parser.add_argument("--data-root", required=True, help="Absolute external data root")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--run-name", default=None)

    parser.add_argument("--ligand-file", action="append", default=[])
    parser.add_argument("--ligand-cid", action="append", default=[])
    parser.add_argument(
        "--ligand-smiles",
        action="append",
        default=[],
        help='Repeatable NAME=SMILES value',
    )
    parser.add_argument("--ligand-manifest", default=None)
    parser.add_argument(
        "--auto-retrieve-ligands",
        action="store_true",
        help=(
            "Run CPU target evidence first, then Stage 4A compound retrieval, "
            "then dock the generated ligand manifest. Requires --fasta."
        ),
    )
    parser.add_argument("--auto-retrieve-max-candidates", type=int, default=6)
    parser.add_argument(
        "--auto-retrieve-mode",
        choices=("rules-only", "hybrid", "generic-strict"),
        default="rules-only",
        help=(
            "Select Stage 4A retrieval behavior. generic-strict disables "
            "all local target-specific seed compounds."
        ),
    )
    parser.add_argument(
        "--auto-retrieve-no-fetch-structures",
        action="store_true",
        help="Run Stage 4A without PubChem structure fetching. Usually not suitable for docking.",
    )
    parser.add_argument("--auto-retrieve-pubchem-timeout-seconds", type=int, default=60)

    parser.add_argument(
        "--fasta",
        type=Path,
        default=None,
        help=(
            "Optional protein FASTA to send to the CPU annotation/homolog API. "
            "Runs concurrently with docking and writes homolog_search_*.json "
            "into the run output directory."
        ),
    )
    parser.add_argument(
        "--homolog-api-url",
        default=DEFAULT_API_URL,
        help="CPU homolog/annotation API endpoint for --fasta.",
    )
    parser.add_argument(
        "--homolog-timeout-seconds",
        type=int,
        default=7200,
        help="Maximum time to wait for the CPU homolog/annotation API request.",
    )


    parser.add_argument("--center-x", type=float)
    parser.add_argument("--center-y", type=float)
    parser.add_argument("--center-z", type=float)
    parser.add_argument("--size-x", type=float)
    parser.add_argument("--size-y", type=float)
    parser.add_argument("--size-z", type=float)
    parser.add_argument("--box-json", default=None, help="Path to a reference_box.json file containing center_x/center_y/center_z/size_x/size_y/size_z.")
    parser.add_argument("--autobox-ligand", default=None)
    parser.add_argument(
        "--reference-ligand",
        default=None,
        help=(
            "Optional cognate reference-ligand SDF used only after "
            "normal docking and pocket selection to calculate pose-recovery RMSD. "
            "The reference ligand does not influence pocket detection, docking, "
            "scoring, filtering, clustering, or selection."
        ),
    )
    parser.add_argument(
        "--pose-recovery-rmsd-threshold",
        type=float,
        default=2.0,
        help=(
            "Heavy-atom RMSD threshold in angstroms for an optional "
            "cognate pose-recovery benchmark."
        ),
    )
    parser.add_argument("--fpocket-padding", type=float, default=4.0)
    parser.add_argument("--fpocket-pocket", type=int, default=None)
    parser.add_argument(
        "--fpocket-top-n",
        type=int,
        default=1,
        help=(
            "When no explicit box or autobox ligand is provided, dock the top N "
            "fpocket pockets instead of only the highest-scoring pocket."
        ),
    )
    parser.add_argument(
        "--fpocket-merge-nearby",
        action="store_true",
        help=(
            "Append pairwise merged boxes for selected "
            "fpocket pockets whose closest alpha-sphere "
            "vertices are within the merge threshold. "
            "Independent top-N pockets remain available."
        ),
    )
    parser.add_argument(
        "--fpocket-merge-distance",
        type=float,
        default=4.0,
        help=(
            "Maximum alpha-sphere vertex distance in "
            "angstroms for constructing a merged "
            "fpocket candidate."
        ),
    )

    parser.add_argument(
        "--pocket-evidence-json",
        default=None,
        help=(
            "Optional pocket_evidence.v0.1 JSON containing "
            "curated, homolog-transferred, or expert-supplied "
            "functional residues. The file controls whether "
            "evidence is report-only or may prioritize "
            "biologically supported pockets."
        ),
    )

    parser.add_argument(
        "--auto-reference-evidence",
        action="store_true",
        help=(
            "Automatically acquire UniProt functional-site "
            "annotations, select a linked PDB reference, "
            "transfer the residues to the submitted receptor, "
            "and use the generated pocket evidence."
        ),
    )

    parser.add_argument(
        "--reference-uniprot-accession",
        default=None,
        help=(
            "Optional UniProt accession override for "
            "--auto-reference-evidence. When omitted, "
            "EXORCIST attempts conservative discovery from "
            "the submitted FASTA header or filename."
        ),
    )

    parser.add_argument(
        "--reference-uniprot-json",
        default=None,
        help=(
            "Optional previously downloaded UniProt JSON "
            "used instead of accession discovery or live "
            "UniProt retrieval."
        ),
    )

    parser.add_argument(
        "--reference-pdb-id",
        default=None,
        help=(
            "Optional override for the automatically ranked "
            "UniProt-linked PDB structure."
        ),
    )

    parser.add_argument(
        "--reference-chain",
        default=None,
        help=(
            "Optional chain override for the reference PDB."
        ),
    )

    parser.add_argument(
        "--receptor-chain",
        default=None,
        help=(
            "Optional receptor chain onto which functional "
            "residues should be transferred."
        ),
    )

    parser.add_argument(
        "--reference-pdb",
        default=None,
        help=(
            "Optional local PDB file for the selected "
            "reference structure instead of downloading it."
        ),
    )

    parser.add_argument(
        "--reference-evidence-timeout-seconds",
        type=float,
        default=60.0,
        help=(
            "Timeout for UniProt and PDB retrieval during "
            "automatic reference-evidence generation."
        ),
    )

    parser.add_argument("--seeds", nargs="+", type=int, default=[2026, 3101, 4202])
    parser.add_argument("--exhaustiveness", type=int, default=32)
    parser.add_argument("--num-modes", type=int, default=20)
    parser.add_argument(
        "--cnn-scoring",
        choices=["none", "rescore", "refinement", "all"],
        default="refinement",
    )
    parser.add_argument("--cluster-threshold", type=float, default=2.0)
    parser.add_argument("--max-hypotheses", type=int, default=3)
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--cpu", type=int, default=None)
    parser.add_argument("--device", type=int, default=None)

    parser.add_argument(
        "--gnina-timeout-seconds",
        type=int,
        default=3600,
        help=(
            "Maximum runtime for each GNINA seed/pocket job. "
            "Use 0 to disable the timeout."
        ),
    )
    parser.add_argument("--gnina-bin", default="gnina")
    parser.add_argument("--fpocket-bin", default="fpocket")
    parser.add_argument("--obabel-bin", default="obabel")
    parser.add_argument("--pdb2pqr-bin", default="pdb2pqr")
    parser.add_argument("--meeko-receptor-bin", default="mk_prepare_receptor.py")
    parser.add_argument("--meeko-ligand-bin", default="mk_prepare_ligand.py")
    parser.add_argument(
        "--posebusters-bin",
        default=os.path.expanduser("~/.venvs/posebusters/bin/bust"),
    )
    parser.add_argument("--skip-validity", action="store_true")
    parser.add_argument("--keep-workdir", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--reference-sequence-search-email",
        default=None,
        help=(
            "Contact email required by EMBL-EBI "
            "when sequence-based UniProt discovery "
            "is needed because the FASTA contains "
            "no accession."
        ),
    )
    parser.add_argument(
        "--reference-sequence-search-timeout-seconds",
        type=float,
        default=600.0,
        help=(
            "Maximum total seconds allowed for the "
            "remote UniProtKB sequence search."
        ),
    )

    return parser




def load_box_json(path):
    """Load explicit docking box values from a reference_box.json-style file."""
    box_path = Path(path).expanduser()
    if not box_path.exists():
        raise FileNotFoundError(f"Box JSON not found: {box_path}")

    with box_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    required = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Box JSON missing required keys: {', '.join(missing)}")

    return {key: float(data[key]) for key in required}


def _candidate_predicted_receptor_pdbs(folding_dir: Path) -> list[Path]:
    """Return predicted receptor PDB candidates sorted by conservative preference."""
    if not folding_dir.exists():
        return []

    candidates = [
        path
        for path in folding_dir.rglob("*.pdb")
        if path.is_file()
    ]

    def score(path: Path) -> tuple[int, int, int, str]:
        name = path.name.lower()

        if "rank_001" in name or "rank_1" in name:
            rank_score = 0
        elif "model_1" in name:
            rank_score = 1
        else:
            rank_score = 2

        # Prefer relaxed models if present, but do not confuse "unrelaxed"
        # with "relaxed".
        if "relaxed" in name and "unrelaxed" not in name:
            relaxation_score = 0
        elif "unrelaxed" in name:
            relaxation_score = 1
        else:
            relaxation_score = 2

        return (rank_score, relaxation_score, len(name), name)

    return sorted(candidates, key=score)


def _write_pure_fasta_no_receptor_artifacts(
    *,
    output_dir: Path,
    fasta_path: Path,
    folding_dir: Path,
    reason: str,
    attempted_prediction: bool,
    command: list[str] | None = None,
    exit_code: int | None = None,
) -> None:
    """Write a clean abstention artifact when pure FASTA cannot produce a receptor."""
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "stage": "structure_acquisition",
        "status": "skipped",
        "pipeline_outcome": "completed_without_docking",
        "reason_code": "no_receptor_structure_available",
        "reason": reason,
        "fasta": str(fasta_path),
        "structure_prediction": {
            "attempted": attempted_prediction,
            "folding_dir": str(folding_dir),
            "command": command,
            "exit_code": exit_code,
        },
        "downstream_stages": {
            "structure_validation": "skipped",
            "target_evidence": "skipped",
            "compound_retrieval": "skipped",
            "receptor_preparation": "skipped",
            "pocket_detection": "skipped",
            "ligand_preparation": "skipped",
            "gnina_docking": "skipped",
            "pose_validation": "skipped",
            "pose_ranking": "skipped",
        },
    }

    for filename in (
        "pure_fasta_structure_status.json",
        "docking_skipped.json",
    ):
        (output_dir / filename).write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    report = output_dir / "compoundrank_run_report.md"
    report.write_text(
        "\n".join(
            [
                "# CompoundRank Pure-FASTA Run Report",
                "",
                "## Outcome",
                "",
                "Docking was skipped.",
                "",
                f"- Reason code: `{payload['reason_code']}`",
                f"- Reason: {reason}",
                f"- FASTA: `{fasta_path}`",
                f"- Structure prediction attempted: {attempted_prediction}",
                f"- Folding directory: `{folding_dir}`",
                "",
                "## Interpretation",
                "",
                "No receptor structure was available for docking. This is a clean",
                "abstention result, not a candidate-discovery result.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _resolve_receptor_from_pure_fasta(
    *,
    fasta_path: Path,
    output_dir: Path,
    structure_predictor_bin: str,
    structure_predictor_extra_args: list[str],
    skip_structure_prediction: bool,
) -> Path | None:
    """Resolve a receptor PDB from FASTA by reusing or running structure prediction."""
    folding_dir = output_dir / "folding" / "colabfold"
    folding_dir.mkdir(parents=True, exist_ok=True)

    existing = _candidate_predicted_receptor_pdbs(folding_dir)
    if existing:
        receptor = existing[0]
        (output_dir / "receptor_pdb.txt").write_text(
            str(receptor) + "\n",
            encoding="utf-8",
        )
        print(f"[PURE FASTA] Reusing predicted receptor: {receptor}")
        return receptor

    if skip_structure_prediction:
        reason = (
            "Pure-FASTA mode was requested without a receptor, but structure "
            "prediction was disabled by --skip-structure-prediction."
        )
        _write_pure_fasta_no_receptor_artifacts(
            output_dir=output_dir,
            fasta_path=fasta_path,
            folding_dir=folding_dir,
            reason=reason,
            attempted_prediction=False,
        )
        print("[PURE FASTA] Docking skipped: no receptor structure available.")
        return None

    predictor = shutil.which(structure_predictor_bin)
    if predictor is None:
        reason = (
            "Pure-FASTA mode was requested without a receptor, but the structure "
            f"predictor executable was not found: {structure_predictor_bin}."
        )
        _write_pure_fasta_no_receptor_artifacts(
            output_dir=output_dir,
            fasta_path=fasta_path,
            folding_dir=folding_dir,
            reason=reason,
            attempted_prediction=False,
        )
        print("[PURE FASTA] Docking skipped: structure predictor not found.")
        return None

    command = [
        predictor,
        str(fasta_path),
        str(folding_dir),
        *structure_predictor_extra_args,
    ]

    print("[PURE FASTA] Running structure prediction:")
    print("[PURE FASTA] " + " ".join(command))

    stdout_path = folding_dir / "structure_prediction_stdout.log"
    stderr_path = folding_dir / "structure_prediction_stderr.log"

    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w",
        encoding="utf-8",
    ) as stderr:
        completed = subprocess.run(
            command,
            stdout=stdout,
            stderr=stderr,
            text=True,
            check=False,
        )

    candidates = _candidate_predicted_receptor_pdbs(folding_dir)
    if not candidates:
        reason = (
            "Structure prediction finished but no receptor PDB was found in "
            "the expected folding output directory."
        )
        _write_pure_fasta_no_receptor_artifacts(
            output_dir=output_dir,
            fasta_path=fasta_path,
            folding_dir=folding_dir,
            reason=reason,
            attempted_prediction=True,
            command=command,
            exit_code=completed.returncode,
        )
        print("[PURE FASTA] Docking skipped: no predicted receptor PDB found.")
        return None

    receptor = candidates[0]

    manifest = {
        "schema_version": "pure_fasta_structure_prediction.v0.1",
        "status": "complete",
        "fasta": str(fasta_path),
        "folding_dir": str(folding_dir),
        "command": command,
        "exit_code": completed.returncode,
        "selected_receptor_pdb": str(receptor),
        "candidate_receptor_pdbs": [str(path) for path in candidates],
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
    }

    (output_dir / "structure_prediction_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_dir / "receptor_pdb.txt").write_text(
        str(receptor) + "\n",
        encoding="utf-8",
    )

    print(f"[PURE FASTA] Selected receptor: {receptor}")
    return receptor


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    pure_fasta_mode = args.receptor is None and args.fasta is not None
    manual_ligand_inputs = any(
        (
            args.ligand_file,
            args.ligand_cid,
            args.ligand_smiles,
            args.ligand_manifest,
        )
    )

    if (
        pure_fasta_mode
        and not manual_ligand_inputs
        and not args.auto_retrieve_ligands
    ):
        args.auto_retrieve_ligands = True
        args.auto_retrieve_mode = "generic-strict"
        print(
            "[PURE FASTA] No ligands supplied; enabling "
            "generic-strict ligand retrieval by default."
        )
    receptor = (
        require_absolute_external_file(args.receptor, "Receptor PDB")
        if args.receptor
        else None
    )
    receptor_ensemble_json = (
        require_absolute_external_file(
            args.receptor_ensemble_json,
            "Receptor ensemble manifest",
        )
        if args.receptor_ensemble_json
        else None
    )
    aligned_receptor_ensemble_json = (
        require_absolute_external_file(
            args.aligned_receptor_ensemble_json,
            "Aligned receptor ensemble manifest",
        )
        if args.aligned_receptor_ensemble_json
        else None
    )

    validate_receptor_ensemble_options(
        receptor_ensemble_json,
        aligned_receptor_ensemble_json,
    )
    if receptor is not None and receptor.suffix.lower() != ".pdb":
        raise ValueError("--receptor must be a PDB file")
    data_root = require_absolute_external_dir(args.data_root, "Data root", create=True)

    # Early pure-FASTA structure acquisition/abstention gate.
    #
    # This must happen before older ligand/reference/manifest handling because
    # several legacy CLI paths assume a receptor-backed run. In pure-FASTA mode,
    # no receptor exists yet.
    if pure_fasta_mode:
        if args.fasta is None:
            raise ValueError(
                "Pure-FASTA mode requires --fasta when --receptor is omitted."
            )

        if receptor_ensemble_json is not None or aligned_receptor_ensemble_json is not None:
            raise ValueError(
                "Receptor ensemble options require --receptor until "
                "pure-FASTA ensemble generation is implemented."
            )

        fasta_for_structure = require_absolute_external_file(
            args.fasta,
            "FASTA",
        )

        if args.output_dir:
            output_dir = Path(args.output_dir)
        else:
            run_label = sanitize_name(
                args.run_name or f"pure_fasta_{fasta_for_structure.stem}"
            )
            output_dir = data_root / "results" / run_label

        if not output_dir.is_absolute():
            raise ValueError("--output-dir must be an absolute path in pure-FASTA mode")

        # Preserve the resolved output directory for the later legacy path if
        # structure prediction succeeds and the normal pipeline continues.
        args.output_dir = str(output_dir)

        receptor = _resolve_receptor_from_pure_fasta(
            fasta_path=fasta_for_structure,
            output_dir=output_dir,
            structure_predictor_bin=args.structure_predictor_bin,
            structure_predictor_extra_args=args.structure_predictor_extra_arg,
            skip_structure_prediction=args.skip_structure_prediction,
        )

        if receptor is None:
            return 0

    automatic_reference_options = (
        args.reference_uniprot_accession,
        args.reference_uniprot_json,
        args.reference_pdb_id,
        args.reference_chain,
        args.receptor_chain,
        args.reference_pdb,
    )

    if args.auto_reference_evidence:
        if args.pocket_evidence_json:
            raise ValueError(
                "--auto-reference-evidence cannot be combined "
                "with --pocket-evidence-json"
            )

        if args.fasta is None:
            raise ValueError(
                "--auto-reference-evidence requires --fasta"
            )

        source_count = sum(
            value is not None
            for value in (
                args.reference_uniprot_accession,
                args.reference_uniprot_json,
            )
        )

        if source_count > 1:
            raise ValueError(
                "--auto-reference-evidence accepts at most "
                "one of --reference-uniprot-accession or "
                "--reference-uniprot-json"
            )

        if args.reference_evidence_timeout_seconds <= 0:
            raise ValueError(
                "--reference-evidence-timeout-seconds must "
                "be greater than zero"
            )

    elif any(
        value is not None
        for value in automatic_reference_options
    ):
        raise ValueError(
            "Automatic reference-evidence options require "
            "--auto-reference-evidence"
        )

    if args.output_dir:
        output_dir = require_absolute_external_dir(
            args.output_dir,
            "Output directory",
            create=True,
        )
    else:
        run_name = sanitize_name(args.run_name or receptor.stem)
        output_dir = require_absolute_external_dir(
            data_root / "results" / run_name,
            "Output directory",
            create=True,
        )

    if args.auto_retrieve_ligands:
        if args.fasta is None:
            raise ValueError("--auto-retrieve-ligands requires --fasta")
        if args.ligand_manifest or args.ligand_file or args.ligand_cid or args.ligand_smiles:
            raise ValueError(
                "--auto-retrieve-ligands should not be combined with manual ligand inputs"
            )
        requests = []
    else:
        manifest_requests = []
        if args.ligand_manifest:
            manifest_path = require_absolute_external_file(
                args.ligand_manifest,
                "Ligand manifest",
            )
            manifest_requests = read_manifest(manifest_path)

        # Validate file ligand paths as external before passing them onward.
        ligand_files = [
            str(require_absolute_external_file(value, "Ligand file"))
            for value in args.ligand_file
        ]
        requests = combine_requests(
            ligand_files,
            args.ligand_cid,
            args.ligand_smiles,
            manifest_requests,
        )
        normalized_requests = []
        for request in requests:
            if request.source_type == "file":
                external_path = require_absolute_external_file(
                    request.value,
                    f"Ligand file ({request.name})",
                )
                normalized_requests.append(
                    type(request)(request.name, request.source_type, str(external_path))
                )
            else:
                normalized_requests.append(request)
        requests = normalized_requests

    autobox_ligand = None
    if args.autobox_ligand:
        autobox_ligand = require_absolute_external_file(
            args.autobox_ligand,
            "Autobox ligand",
        )

    reference_ligand = None
    if args.reference_ligand:
        reference_ligand = require_absolute_external_file(
            args.reference_ligand,
            "Reference ligand",
        )

        if reference_ligand.suffix.lower() not in {
            ".sdf",
            ".sd",
        }:
            raise ValueError(
                "--reference-ligand must be an SDF file"
            )

    if args.pose_recovery_rmsd_threshold <= 0:
        raise ValueError(
            "--pose-recovery-rmsd-threshold must be greater than zero"
        )

    reference_uniprot_json = None

    if args.reference_uniprot_json:
        reference_uniprot_json = (
            require_absolute_external_file(
                args.reference_uniprot_json,
                "Reference UniProt JSON",
            )
        )

    reference_pdb = None

    if args.reference_pdb:
        reference_pdb = (
            require_absolute_external_file(
                args.reference_pdb,
                "Reference PDB",
            )
        )

        if reference_pdb.suffix.lower() != ".pdb":
            raise ValueError(
                "--reference-pdb must be a PDB file"
            )

    center_x = args.center_x
    center_y = args.center_y
    center_z = args.center_z
    size_x = args.size_x
    size_y = args.size_y
    size_z = args.size_z

    if args.box_json:
        if autobox_ligand is not None:
            raise ValueError("Use either --box-json or --autobox-ligand, not both.")
        explicit_cli_box_values = (center_x, center_y, center_z, size_x, size_y, size_z)
        if any(value is not None for value in explicit_cli_box_values):
            raise ValueError("Use either --box-json or explicit --center/--size values, not both.")

        box_values = load_box_json(args.box_json)
        center_x = box_values["center_x"]
        center_y = box_values["center_y"]
        center_z = box_values["center_z"]
        size_x = box_values["size_x"]
        size_y = box_values["size_y"]
        size_z = box_values["size_z"]

    run_pipeline(
        receptor_pdb=receptor,
        receptor_ensemble_json=(
            receptor_ensemble_json
        ),
        aligned_receptor_ensemble_json=(
            aligned_receptor_ensemble_json
        ),
        ligand_requests=requests,
        data_root=data_root,
        output_dir=output_dir,
        fasta_path=args.fasta,
        homolog_api_url=args.homolog_api_url,
        homolog_timeout_seconds=args.homolog_timeout_seconds,
        auto_retrieve_ligands=args.auto_retrieve_ligands,
        auto_retrieve_mode=args.auto_retrieve_mode,
        auto_retrieve_max_candidates=args.auto_retrieve_max_candidates,
        auto_retrieve_fetch_structures=not args.auto_retrieve_no_fetch_structures,
        auto_retrieve_pubchem_timeout_seconds=args.auto_retrieve_pubchem_timeout_seconds,
        seeds=args.seeds,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        autobox_ligand=autobox_ligand,
        reference_ligand=reference_ligand,
        pose_recovery_rmsd_threshold=(
            args.pose_recovery_rmsd_threshold
        ),
        fpocket_padding=args.fpocket_padding,
        fpocket_pocket=args.fpocket_pocket,
        fpocket_top_n=args.fpocket_top_n,
        fpocket_merge_nearby=(
            args.fpocket_merge_nearby
        ),
        fpocket_merge_distance=(
            args.fpocket_merge_distance
        ),
        pocket_evidence_json=(
            Path(
                args.pocket_evidence_json
            ).expanduser().resolve()
            if args.pocket_evidence_json
            else None
        ),
        auto_reference_evidence=(
            args.auto_reference_evidence
        ),
        reference_uniprot_accession=(
            args.reference_uniprot_accession
        ),
        reference_uniprot_json=(
            reference_uniprot_json
        ),
        reference_pdb_id=(
            args.reference_pdb_id
        ),
        reference_chain_id=(
            args.reference_chain
        ),
        receptor_chain_id=(
            args.receptor_chain
        ),
        reference_pdb=reference_pdb,
        reference_evidence_timeout_seconds=(
            args.reference_evidence_timeout_seconds
        ),
        reference_sequence_search_email=(
            args.reference_sequence_search_email
        ),
        reference_sequence_search_timeout_seconds=(
            args.reference_sequence_search_timeout_seconds
        ),
        max_hypotheses=args.max_hypotheses,
        cluster_threshold=args.cluster_threshold,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        cnn_scoring=args.cnn_scoring,
        ph=args.ph,
        gnina_bin=args.gnina_bin,
        gnina_timeout_seconds=(
            args.gnina_timeout_seconds
            if args.gnina_timeout_seconds > 0
            else None
        ),
        fpocket_bin=args.fpocket_bin,
        obabel_bin=args.obabel_bin,
        pdb2pqr_bin=args.pdb2pqr_bin,
        meeko_receptor_bin=args.meeko_receptor_bin,
        meeko_ligand_bin=args.meeko_ligand_bin,
        posebusters_bin=args.posebusters_bin,
        skip_validity=args.skip_validity,
        keep_workdir=args.keep_workdir,
        overwrite=args.overwrite,
        cpu=args.cpu,
        device=args.device,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
