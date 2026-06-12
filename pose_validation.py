from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem


SCORE_PROPERTIES = (
    "CNNscore",
    "CNNaffinity",
    "minimizedAffinity",
)


def load_first_molecule(path: Path) -> Chem.Mol:
    if not path.is_file():
        raise FileNotFoundError(f"Molecule file does not exist: {path}")

    supplier = Chem.SDMolSupplier(str(path), removeHs=False)

    molecule = next(
        (mol for mol in supplier if mol is not None),
        None,
    )

    if molecule is None:
        raise RuntimeError(f"Could not read a molecule from: {path}")

    return molecule


def load_pose_molecules(path: Path) -> list[Chem.Mol]:
    if not path.is_file():
        raise FileNotFoundError(f"Pose file does not exist: {path}")

    poses = [
        molecule
        for molecule in Chem.SDMolSupplier(
            str(path),
            removeHs=False,
        )
        if molecule is not None
    ]

    if not poses:
        raise RuntimeError(f"No valid poses were found in: {path}")

    return poses


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


def calculate_symmetry_aware_rmsd(
    docked_molecule: Chem.Mol,
    reference_molecule: Chem.Mol,
) -> tuple[float, int]:
    docked_heavy = Chem.RemoveHs(docked_molecule)
    reference_heavy = Chem.RemoveHs(reference_molecule)

    if docked_heavy.GetNumAtoms() != reference_heavy.GetNumAtoms():
        raise RuntimeError(
            "Docked and reference molecules have different heavy-atom "
            f"counts: {docked_heavy.GetNumAtoms()} versus "
            f"{reference_heavy.GetNumAtoms()}"
        )

    mappings = reference_heavy.GetSubstructMatches(
        docked_heavy,
        uniquify=False,
        useChirality=True,
        maxMatches=100000,
    )

    if not mappings:
        mappings = reference_heavy.GetSubstructMatches(
            docked_heavy,
            uniquify=False,
            useChirality=False,
            maxMatches=100000,
        )

    if not mappings:
        raise RuntimeError(
            "Could not establish an atom mapping between the docked "
            "and reference molecules."
        )

    docked_conformer = docked_heavy.GetConformer()
    reference_conformer = reference_heavy.GetConformer()

    best_rmsd = float("inf")

    for mapping in mappings:
        squared_distances = []

        for docked_index, reference_index in enumerate(mapping):
            docked_point = docked_conformer.GetAtomPosition(
                docked_index
            )
            reference_point = reference_conformer.GetAtomPosition(
                reference_index
            )

            squared_distances.append(
                (docked_point.x - reference_point.x) ** 2
                + (docked_point.y - reference_point.y) ** 2
                + (docked_point.z - reference_point.z) ** 2
            )

        rmsd = float(np.sqrt(np.mean(squared_distances)))
        best_rmsd = min(best_rmsd, rmsd)

    return best_rmsd, len(mappings)


def validate_poses(
    poses_path: Path,
    reference_path: Path,
    output_dir: Path,
    prefix: str,
    pass_threshold: float = 2.0,
) -> dict[str, Any]:
    reference = load_first_molecule(reference_path)
    poses = load_pose_molecules(poses_path)

    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []

    for pose_number, pose in enumerate(poses, start=1):
        rmsd, mapping_count = calculate_symmetry_aware_rmsd(
            docked_molecule=pose,
            reference_molecule=reference,
        )

        row: dict[str, Any] = {
            "pose": pose_number,
            "rmsd": rmsd,
            "atom_mappings_tested": mapping_count,
        }

        for property_name in SCORE_PROPERTIES:
            row[property_name] = get_float_property(
                pose,
                property_name,
            )

        rows.append(row)

    rows.sort(key=lambda row: row["rmsd"])

    best = rows[0]
    best_pose_number = int(best["pose"])
    best_pose = poses[best_pose_number - 1]

    pose_one = next(
        row for row in rows if row["pose"] == 1
    )

    csv_path = output_dir / f"{prefix}_rmsd.csv"
    json_path = output_dir / f"{prefix}_summary.json"
    best_pose_path = output_dir / f"{prefix}_best_pose.sdf"

    with csv_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pose",
                "rmsd",
                "atom_mappings_tested",
                "CNNscore",
                "CNNaffinity",
                "minimizedAffinity",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    sdf_writer = Chem.SDWriter(str(best_pose_path))
    sdf_writer.write(best_pose)
    sdf_writer.close()

    passed = best["rmsd"] <= pass_threshold

    summary = {
        "poses_file": str(poses_path),
        "reference_ligand": str(reference_path),
        "pose_count": len(poses),
        "pass_threshold_angstrom": pass_threshold,
        "passed": passed,
        "best_pose": best,
        "pose_1": pose_one,
        "top_10_by_rmsd": rows[:10],
        "files": {
            "rmsd_csv": str(csv_path),
            "summary_json": str(json_path),
            "best_pose_sdf": str(best_pose_path),
        },
    }

    json_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary


def format_number(value: float | None) -> str:
    if value is None:
        return "N/A"

    return f"{value:.3f}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare docked ligand poses with an experimental "
            "reference ligand using symmetry-aware heavy-atom RMSD."
        )
    )

    parser.add_argument(
        "--poses",
        required=True,
        help="SDF containing one or more docked poses.",
    )

    parser.add_argument(
        "--reference",
        required=True,
        help="Experimental reference ligand SDF.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
    )

    parser.add_argument(
        "--prefix",
        default="pose_validation",
    )

    parser.add_argument(
        "--pass-threshold",
        type=float,
        default=2.0,
        help="Maximum RMSD in Angstroms for a passing result.",
    )

    args = parser.parse_args()

    summary = validate_poses(
        poses_path=Path(args.poses),
        reference_path=Path(args.reference),
        output_dir=Path(args.output_dir),
        prefix=args.prefix,
        pass_threshold=args.pass_threshold,
    )

    print("\nTop 10 poses by crystallographic RMSD:\n")

    for row in summary["top_10_by_rmsd"]:
        print(
            f"Pose {row['pose']:2d}: "
            f"RMSD={row['rmsd']:.3f} Å  "
            f"CNNscore={format_number(row['CNNscore'])}  "
            f"CNNaffinity={format_number(row['CNNaffinity'])}  "
            f"Affinity={format_number(row['minimizedAffinity'])}"
        )

    best = summary["best_pose"]
    pose_one = summary["pose_1"]

    print("\n=== POSE VALIDATION RESULT ===")
    print(f"Best pose: {best['pose']}")
    print(f"Best RMSD: {best['rmsd']:.3f} Å")
    print(f"GNINA pose 1 RMSD: {pose_one['rmsd']:.3f} Å")

    if summary["passed"]:
        print("PASS: A crystallographic-like pose was recovered.")
    else:
        print("FAIL: No pose passed the RMSD threshold.")

    print("\nSaved files:")

    for name, path in summary["files"].items():
        print(f"  {name}: {path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
