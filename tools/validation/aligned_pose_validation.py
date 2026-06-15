from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

from robust_pose_rmsd import calculate_symmetry_aware_rmsd


SCORE_PROPERTIES = (
    "CNNscore",
    "CNNaffinity",
    "minimizedAffinity",
)


def require_file(path: Path, label: str) -> Path:
    path = path.resolve()

    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")

    return path


def read_ca_atoms(path: Path) -> dict[tuple[str, int], np.ndarray]:
    atoms: dict[tuple[str, int], np.ndarray] = {}

    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue

        if line[12:16].strip() != "CA":
            continue

        if line[16].strip() not in {"", "A"}:
            continue

        chain = line[21].strip()

        try:
            residue_number = int(line[22:26])
            xyz = np.array(
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ],
                dtype=float,
            )
        except ValueError:
            continue

        atoms[(chain, residue_number)] = xyz

    if not atoms:
        raise RuntimeError(f"No CA atoms found in: {path}")

    return atoms


def calculate_fit(
    mobile: np.ndarray,
    reference: np.ndarray,
) -> dict[str, Any]:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)

    mobile_zero = mobile - mobile_center
    reference_zero = reference - reference_center

    covariance = mobile_zero.T @ reference_zero
    u, _, vt = np.linalg.svd(covariance)

    rotation = vt.T @ u.T

    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T

    fitted = mobile_zero @ rotation.T + reference_center

    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (fitted - reference) ** 2,
                    axis=1,
                )
            )
        )
    )

    return {
        "rotation": rotation,
        "mobile_center": mobile_center,
        "reference_center": reference_center,
        "rmsd": rmsd,
    }


def choose_receptor_alignment(
    mobile_path: Path,
    reference_path: Path,
) -> dict[str, Any]:
    mobile_atoms = read_ca_atoms(mobile_path)
    reference_atoms = read_ca_atoms(reference_path)

    mappings = (
        {"A": "A", "B": "B"},
        {"A": "B", "B": "A"},
    )

    fits = []

    for chain_mapping in mappings:
        mobile_points = []
        reference_points = []

        for mobile_chain, reference_chain in chain_mapping.items():
            residue_numbers = sorted(
                residue_number
                for chain, residue_number in mobile_atoms
                if chain == mobile_chain
            )

            for residue_number in residue_numbers:
                mobile_key = (mobile_chain, residue_number)
                reference_key = (reference_chain, residue_number)

                if reference_key not in reference_atoms:
                    continue

                mobile_points.append(mobile_atoms[mobile_key])
                reference_points.append(reference_atoms[reference_key])

        if len(mobile_points) < 3:
            continue

        fit = calculate_fit(
            np.asarray(mobile_points),
            np.asarray(reference_points),
        )

        fit["chain_mapping"] = chain_mapping
        fit["matched_ca_atoms"] = len(mobile_points)
        fits.append(fit)

    if not fits:
        raise RuntimeError(
            "Could not establish a receptor alignment."
        )

    return min(fits, key=lambda item: item["rmsd"])


def transform_molecule(
    molecule: Chem.Mol,
    alignment: dict[str, Any],
) -> Chem.Mol:
    transformed = Chem.Mol(molecule)
    conformer = transformed.GetConformer()

    rotation = alignment["rotation"]
    mobile_center = alignment["mobile_center"]
    reference_center = alignment["reference_center"]

    for atom_index in range(transformed.GetNumAtoms()):
        point = conformer.GetAtomPosition(atom_index)

        coordinate = np.array(
            [point.x, point.y, point.z],
            dtype=float,
        )

        new_coordinate = (
            (coordinate - mobile_center) @ rotation.T
            + reference_center
        )

        conformer.SetAtomPosition(
            atom_index,
            Point3D(
                float(new_coordinate[0]),
                float(new_coordinate[1]),
                float(new_coordinate[2]),
            ),
        )

    return transformed


def get_float_property(
    molecule: Chem.Mol,
    name: str,
) -> float | None:
    if not molecule.HasProp(name):
        return None

    try:
        return float(molecule.GetProp(name))
    except ValueError:
        return None


def load_reference_ligand(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
    )

    molecule = next(
        (mol for mol in supplier if mol is not None),
        None,
    )

    if molecule is None:
        raise RuntimeError(
            f"Could not read reference ligand: {path}"
        )

    return molecule


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Align a docked receptor to an experimental receptor, "
            "transform all ligand poses, and calculate heavy-atom RMSD."
        )
    )

    parser.add_argument("--poses", required=True)
    parser.add_argument("--mobile-receptor", required=True)
    parser.add_argument("--reference-receptor", required=True)
    parser.add_argument("--reference-ligand", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--prefix",
        default="aligned_pose_validation",
    )
    parser.add_argument(
        "--pass-threshold",
        type=float,
        default=2.0,
    )

    args = parser.parse_args()

    poses_path = require_file(
        Path(args.poses),
        "Docked poses",
    )
    mobile_receptor = require_file(
        Path(args.mobile_receptor),
        "Mobile receptor",
    )
    reference_receptor = require_file(
        Path(args.reference_receptor),
        "Reference receptor",
    )
    reference_ligand_path = require_file(
        Path(args.reference_ligand),
        "Reference ligand",
    )

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    alignment = choose_receptor_alignment(
        mobile_path=mobile_receptor,
        reference_path=reference_receptor,
    )

    reference_ligand = load_reference_ligand(
        reference_ligand_path
    )

    supplier = Chem.SDMolSupplier(
        str(poses_path),
        removeHs=False,
    )

    rows: list[dict[str, Any]] = []
    transformed_poses: dict[int, Chem.Mol] = {}
    skipped_pose_numbers: list[int] = []

    for pose_number, molecule in enumerate(
        supplier,
        start=1,
    ):
        if molecule is None:
            skipped_pose_numbers.append(pose_number)
            continue

        transformed = transform_molecule(
            molecule,
            alignment,
        )

        rmsd, mapping_count = calculate_symmetry_aware_rmsd(
            docked_molecule=transformed,
            reference_molecule=reference_ligand,
        )

        row: dict[str, Any] = {
            "pose": pose_number,
            "rmsd": rmsd,
            "atom_mappings_tested": mapping_count,
        }

        for property_name in SCORE_PROPERTIES:
            row[property_name] = get_float_property(
                molecule,
                property_name,
            )

        rows.append(row)
        transformed_poses[pose_number] = transformed

    if not rows:
        raise RuntimeError(
            "No readable docked poses could be validated."
        )

    rows.sort(key=lambda row: row["rmsd"])

    best = rows[0]
    best_pose_number = int(best["pose"])

    csv_path = output_dir / f"{args.prefix}_rmsd.csv"
    json_path = output_dir / f"{args.prefix}_summary.json"
    best_pose_path = (
        output_dir / f"{args.prefix}_best_pose.sdf"
    )

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

    writer = Chem.SDWriter(str(best_pose_path))
    writer.write(transformed_poses[best_pose_number])
    writer.close()

    passed = best["rmsd"] <= args.pass_threshold

    pose_one = next(
        (
            row
            for row in rows
            if row["pose"] == 1
        ),
        None,
    )

    summary = {
        "status": "complete",
        "passed": passed,
        "pass_threshold_angstrom": args.pass_threshold,
        "receptor_alignment": {
            "chain_mapping": alignment["chain_mapping"],
            "matched_ca_atoms": alignment["matched_ca_atoms"],
            "ca_rmsd": alignment["rmsd"],
        },
        "sdf_record_count": len(supplier),
        "validated_pose_count": len(rows),
        "skipped_pose_numbers": skipped_pose_numbers,
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

    print("\nReceptor alignment:")
    print(
        "Chain mapping:",
        alignment["chain_mapping"],
    )
    print(
        f"Matched CA atoms: "
        f"{alignment['matched_ca_atoms']}"
    )
    print(
        f"Receptor CA RMSD: "
        f"{alignment['rmsd']:.3f} Å"
    )

    if skipped_pose_numbers:
        print(
            "Unreadable SDF pose records:",
            skipped_pose_numbers,
        )

    print("\nTop 10 poses by experimental RMSD:\n")

    for row in rows[:10]:
        print(
            f"Pose {row['pose']:2d}: "
            f"RMSD={row['rmsd']:.3f} Å  "
            f"CNNscore={row['CNNscore']:.3f}  "
            f"CNNaffinity={row['CNNaffinity']:.3f}  "
            f"Affinity={row['minimizedAffinity']:.3f}"
        )

    print("\n=== PREDICTED-RECEPTOR VALIDATION ===")
    print(f"Best pose: {best_pose_number}")
    print(f"Best RMSD: {best['rmsd']:.3f} Å")

    if pose_one is not None:
        print(
            f"GNINA pose 1 RMSD: "
            f"{pose_one['rmsd']:.3f} Å"
        )

    print(
        "Result:",
        "PASS" if passed else "FAIL",
    )

    print("\nSaved:")
    print(csv_path)
    print(json_path)
    print(best_pose_path)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
