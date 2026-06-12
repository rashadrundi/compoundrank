from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from rdkit import Chem

from pose_validation import validate_poses


SCORE_PROPERTIES = (
    "CNNscore",
    "CNNaffinity",
    "minimizedAffinity",
    "minimizedRMSD",
)


def require_file(path: Path, label: str) -> Path:
    path = path.resolve()

    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")

    if path.stat().st_size == 0:
        raise RuntimeError(f"{label} is empty: {path}")

    return path


def resolve_gnina(executable: str) -> str:
    resolved = shutil.which(executable)

    if resolved is None:
        candidate = Path(executable).expanduser()

        if candidate.is_file():
            return str(candidate.resolve())

        raise RuntimeError(
            f"GNINA executable was not found: {executable}"
        )

    return resolved


def get_float_property(
    molecule: Chem.Mol,
    property_name: str,
) -> float | None:
    if not molecule.HasProp(property_name):
        return None

    try:
        return float(molecule.GetProp(property_name))
    except ValueError:
        return None


def extract_pose_scores(sdf_path: Path) -> list[dict[str, Any]]:
    supplier = Chem.SDMolSupplier(
        str(sdf_path),
        removeHs=False,
    )

    rows: list[dict[str, Any]] = []

    for pose_number, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue

        row: dict[str, Any] = {
            "pose": pose_number,
            "atom_count": molecule.GetNumAtoms(),
            "heavy_atom_count": sum(
                atom.GetAtomicNum() > 1
                for atom in molecule.GetAtoms()
            ),
        }

        for property_name in SCORE_PROPERTIES:
            row[property_name] = get_float_property(
                molecule,
                property_name,
            )

        rows.append(row)

    if not rows:
        raise RuntimeError(
            f"GNINA completed, but no readable poses were found in "
            f"{sdf_path}"
        )

    return rows


def validate_box_arguments(args: argparse.Namespace) -> str:
    explicit_values = (
        args.center_x,
        args.center_y,
        args.center_z,
        args.size_x,
        args.size_y,
        args.size_z,
    )

    explicit_count = sum(
        value is not None
        for value in explicit_values
    )

    has_autobox = args.autobox_ligand is not None

    if has_autobox and explicit_count:
        raise ValueError(
            "Use either --autobox-ligand or all six explicit box "
            "coordinates, not both."
        )

    if has_autobox:
        return "autobox"

    if explicit_count == 6:
        return "explicit"

    if explicit_count:
        raise ValueError(
            "Explicit-box mode requires --center-x, --center-y, "
            "--center-z, --size-x, --size-y and --size-z."
        )

    raise ValueError(
        "A search box is required. Use --autobox-ligand or all six "
        "explicit box values."
    )


def build_command(
    gnina_executable: str,
    receptor: Path,
    ligand: Path,
    poses_output: Path,
    args: argparse.Namespace,
    box_mode: str,
) -> list[str]:
    command = [
        gnina_executable,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
    ]

    if box_mode == "autobox":
        autobox_ligand = require_file(
            Path(args.autobox_ligand),
            "Autobox ligand",
        )

        command += [
            "--autobox_ligand",
            str(autobox_ligand),
            "--autobox_add",
            str(args.autobox_add),
        ]
    else:
        command += [
            "--center_x",
            str(args.center_x),
            "--center_y",
            str(args.center_y),
            "--center_z",
            str(args.center_z),
            "--size_x",
            str(args.size_x),
            "--size_y",
            str(args.size_y),
            "--size_z",
            str(args.size_z),
        ]

    command += [
        "--exhaustiveness",
        str(args.exhaustiveness),
        "--num_modes",
        str(args.num_modes),
        "--seed",
        str(args.seed),
        "--cnn_scoring",
        args.cnn_scoring,
        "--out",
        str(poses_output),
    ]

    if args.cpu is not None:
        command += ["--cpu", str(args.cpu)]

    if args.device is not None:
        command += ["--device", str(args.device)]

    if args.flexres is not None:
        command += ["--flexres", args.flexres]

    return command


def run_command_with_live_log(
    command: list[str],
    log_path: Path,
) -> tuple[int, float]:
    started = time.monotonic()

    with log_path.open(
        "w",
        encoding="utf-8",
        errors="replace",
    ) as log_handle:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        if process.stdout is None:
            raise RuntimeError("Could not capture GNINA output")

        for line in process.stdout:
            print(line, end="", flush=True)
            log_handle.write(line)
            log_handle.flush()

        return_code = process.wait()

    elapsed = time.monotonic() - started
    return return_code, elapsed


def save_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run GNINA docking, save structured outputs, and "
            "optionally validate poses against a reference ligand."
        )
    )

    parser.add_argument(
        "--receptor",
        required=True,
    )

    parser.add_argument(
        "--ligand",
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        required=True,
    )

    parser.add_argument(
        "--prefix",
        default="docking",
    )

    parser.add_argument(
        "--gnina-bin",
        default="gnina",
    )

    parser.add_argument("--center-x", type=float)
    parser.add_argument("--center-y", type=float)
    parser.add_argument("--center-z", type=float)
    parser.add_argument("--size-x", type=float)
    parser.add_argument("--size-y", type=float)
    parser.add_argument("--size-z", type=float)

    parser.add_argument(
        "--autobox-ligand",
        default=None,
    )

    parser.add_argument(
        "--autobox-add",
        type=float,
        default=4.0,
    )

    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=32,
    )

    parser.add_argument(
        "--num-modes",
        type=int,
        default=20,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
    )

    parser.add_argument(
        "--cnn-scoring",
        choices=("none", "rescore", "refinement", "all"),
        default="rescore",
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=None,
    )

    parser.add_argument(
        "--device",
        type=int,
        default=None,
    )

    parser.add_argument(
        "--flexres",
        default=None,
        help=(
            "Comma-separated flexible receptor residues in GNINA "
            "chain:residue format, for example A:32,B:32,A:50."
        ),
    )

    parser.add_argument(
        "--reference-ligand",
        default=None,
        help=(
            "Optional experimental ligand SDF. When provided, "
            "pose_validation.py is run automatically."
        ),
    )

    parser.add_argument(
        "--validation-threshold",
        type=float,
        default=2.0,
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    args = parser.parse_args()

    receptor = require_file(
        Path(args.receptor),
        "Receptor",
    )
    ligand = require_file(
        Path(args.ligand),
        "Ligand",
    )

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    poses_output = output_dir / f"{args.prefix}_poses.sdf"
    log_output = output_dir / f"{args.prefix}_gnina.log"
    scores_output = output_dir / f"{args.prefix}_scores.json"
    summary_output = output_dir / f"{args.prefix}_summary.json"

    generated_files = (
        poses_output,
        log_output,
        scores_output,
        summary_output,
    )

    existing = [
        path for path in generated_files
        if path.exists()
    ]

    if existing and not args.overwrite:
        names = "\n".join(f"  - {path}" for path in existing)
        raise FileExistsError(
            "Output files already exist. Use --overwrite to replace "
            f"them:\n{names}"
        )

    for path in existing:
        path.unlink()

    box_mode = validate_box_arguments(args)
    gnina_executable = resolve_gnina(args.gnina_bin)

    command = build_command(
        gnina_executable=gnina_executable,
        receptor=receptor,
        ligand=ligand,
        poses_output=poses_output,
        args=args,
        box_mode=box_mode,
    )

    print("Running GNINA:")
    print(" ".join(command))
    print()

    return_code, elapsed_seconds = run_command_with_live_log(
        command=command,
        log_path=log_output,
    )

    if return_code != 0:
        raise RuntimeError(
            f"GNINA failed with exit code {return_code}. "
            f"See log: {log_output}"
        )

    if not poses_output.is_file():
        raise RuntimeError(
            f"GNINA exited successfully but did not create: "
            f"{poses_output}"
        )

    pose_scores = extract_pose_scores(poses_output)

    save_json(
        scores_output,
        {
            "pose_count": len(pose_scores),
            "poses": pose_scores,
        },
    )

    validation_summary = None

    if args.reference_ligand is not None:
        reference_ligand = require_file(
            Path(args.reference_ligand),
            "Reference ligand",
        )

        validation_summary = validate_poses(
            poses_path=poses_output,
            reference_path=reference_ligand,
            output_dir=output_dir,
            prefix=f"{args.prefix}_validation",
            pass_threshold=args.validation_threshold,
        )

    summary: dict[str, Any] = {
        "status": "complete",
        "command": command,
        "return_code": return_code,
        "elapsed_seconds": round(elapsed_seconds, 3),
        "inputs": {
            "receptor": str(receptor),
            "ligand": str(ligand),
            "box_mode": box_mode,
        },
        "settings": {
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "seed": args.seed,
            "cnn_scoring": args.cnn_scoring,
            "flexres": args.flexres,
            "autobox_add": (
                args.autobox_add
                if box_mode == "autobox"
                else None
            ),
            "center": (
                {
                    "x": args.center_x,
                    "y": args.center_y,
                    "z": args.center_z,
                }
                if box_mode == "explicit"
                else None
            ),
            "size": (
                {
                    "x": args.size_x,
                    "y": args.size_y,
                    "z": args.size_z,
                }
                if box_mode == "explicit"
                else None
            ),
        },
        "pose_count": len(pose_scores),
        "pose_scores": pose_scores,
        "validation": validation_summary,
        "files": {
            "poses_sdf": str(poses_output),
            "gnina_log": str(log_output),
            "scores_json": str(scores_output),
            "summary_json": str(summary_output),
        },
    }

    save_json(summary_output, summary)

    print("\nDocking complete.")
    print(f"Poses: {len(pose_scores)}")
    print(f"Elapsed: {elapsed_seconds:.1f} seconds")
    print(f"Pose file: {poses_output}")
    print(f"Log: {log_output}")
    print(f"Summary: {summary_output}")

    if validation_summary is not None:
        best = validation_summary["best_pose"]

        print("\nValidation:")
        print(f"Best pose: {best['pose']}")
        print(f"Best RMSD: {best['rmsd']:.3f} Å")
        print(
            "Result:",
            "PASS" if validation_summary["passed"] else "FAIL",
        )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
