from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from aligned_pose_validation import (
    choose_receptor_alignment,
    get_float_property,
    require_file,
    transform_molecule,
)


SCORE_PROPERTIES = (
    "CNNscore",
    "CNNaffinity",
    "minimizedAffinity",
)


def load_first_sdf(path: Path, sanitize: bool = True) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
        sanitize=sanitize,
    )

    molecule = next(
        (mol for mol in supplier if mol is not None),
        None,
    )

    if molecule is None:
        raise RuntimeError(f"Could not read molecule: {path}")

    return molecule


def parse_meeko_index_pairs(pdbqt_path: Path) -> list[tuple[int, int]]:
    numbers: list[int] = []

    for line in pdbqt_path.read_text(errors="replace").splitlines():
        if not line.startswith("REMARK INDEX MAP"):
            continue

        fields = line.split()[3:]

        for field in fields:
            numbers.append(int(field))

    if not numbers:
        raise RuntimeError(
            f"No REMARK INDEX MAP records found in {pdbqt_path}"
        )

    if len(numbers) % 2 != 0:
        raise RuntimeError(
            "The Meeko index map contains an odd number of integers."
        )

    return [
        (numbers[index], numbers[index + 1])
        for index in range(0, len(numbers), 2)
    ]


def choose_pdbqt_to_source_mapping(
    pairs: list[tuple[int, int]],
    source: Chem.Mol,
    example_pose: Chem.Mol,
) -> dict[int, int]:
    # Meeko usually writes:
    # original-SDF-index, PDBQT-index
    #
    # Test both possible orientations and choose the one whose
    # atom elements agree with the GNINA output pose.
    candidates = [
        {
            second - 1: first - 1
            for first, second in pairs
        },
        {
            first - 1: second - 1
            for first, second in pairs
        },
    ]

    scored = []

    for mapping in candidates:
        compared = 0
        matching_elements = 0

        for pose_index, pose_atom in enumerate(example_pose.GetAtoms()):
            if pose_index not in mapping:
                continue

            source_index = mapping[pose_index]

            if source_index < 0 or source_index >= source.GetNumAtoms():
                continue

            source_atom = source.GetAtomWithIdx(source_index)
            compared += 1

            if pose_atom.GetAtomicNum() == source_atom.GetAtomicNum():
                matching_elements += 1

        scored.append(
            (
                matching_elements,
                compared,
                mapping,
            )
        )

    matching_elements, compared, best_mapping = max(
        scored,
        key=lambda item: (
            item[0],
            item[1],
        ),
    )

    if compared != example_pose.GetNumAtoms():
        raise RuntimeError(
            "Meeko index map did not cover every GNINA output atom: "
            f"{compared}/{example_pose.GetNumAtoms()}"
        )

    if matching_elements != compared:
        raise RuntimeError(
            "GNINA output atom order does not agree with the Meeko "
            f"index map: {matching_elements}/{compared} elements matched."
        )

    return best_mapping


def source_to_reference_mappings(
    source: Chem.Mol,
    reference: Chem.Mol,
) -> list[tuple[int, ...]]:
    source_heavy = Chem.RemoveHs(source)
    reference_heavy = Chem.RemoveHs(reference)

    if source_heavy.GetNumAtoms() != reference_heavy.GetNumAtoms():
        raise RuntimeError(
            "Source and reference ligands have different heavy-atom "
            f"counts: {source_heavy.GetNumAtoms()} versus "
            f"{reference_heavy.GetNumAtoms()}"
        )

    mappings = reference_heavy.GetSubstructMatches(
        source_heavy,
        uniquify=False,
        useChirality=True,
        maxMatches=100000,
    )

    if not mappings:
        mappings = reference_heavy.GetSubstructMatches(
            source_heavy,
            uniquify=False,
            useChirality=False,
            maxMatches=100000,
        )

    if not mappings:
        raise RuntimeError(
            "Could not map the original prepared ligand to the "
            "crystallographic reference ligand."
        )

    return list(mappings)


def build_original_to_heavy_index(
    molecule: Chem.Mol,
) -> dict[int, int]:
    mapping: dict[int, int] = {}
    heavy_index = 0

    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue

        mapping[atom.GetIdx()] = heavy_index
        heavy_index += 1

    return mapping


def calculate_indexed_rmsd(
    transformed_pose: Chem.Mol,
    source: Chem.Mol,
    reference: Chem.Mol,
    pdbqt_to_source: dict[int, int],
    source_original_to_heavy: dict[int, int],
    source_reference_mappings: list[tuple[int, ...]],
) -> tuple[float, int]:
    pose_conformer = transformed_pose.GetConformer()
    reference_heavy = Chem.RemoveHs(reference)
    reference_conformer = reference_heavy.GetConformer()

    pose_heavy_records = []

    for pose_index, pose_atom in enumerate(
        transformed_pose.GetAtoms()
    ):
        if pose_atom.GetAtomicNum() <= 1:
            continue

        if pose_index not in pdbqt_to_source:
            raise RuntimeError(
                f"No source-ligand index for pose atom {pose_index + 1}"
            )

        source_original_index = pdbqt_to_source[pose_index]

        if source_original_index not in source_original_to_heavy:
            raise RuntimeError(
                "A GNINA heavy atom mapped to a source hydrogen: "
                f"pose atom {pose_index + 1}"
            )

        source_heavy_index = source_original_to_heavy[
            source_original_index
        ]

        source_atom = source.GetAtomWithIdx(source_original_index)

        if source_atom.GetAtomicNum() != pose_atom.GetAtomicNum():
            raise RuntimeError(
                "Element mismatch between GNINA pose and source "
                f"ligand at pose atom {pose_index + 1}"
            )

        pose_heavy_records.append(
            (
                pose_index,
                source_heavy_index,
            )
        )

    expected_heavy_atoms = len(source_original_to_heavy)

    if len(pose_heavy_records) != expected_heavy_atoms:
        raise RuntimeError(
            "The GNINA pose does not contain the expected number of "
            f"heavy atoms: {len(pose_heavy_records)} versus "
            f"{expected_heavy_atoms}"
        )

    best_rmsd = float("inf")

    for source_to_reference in source_reference_mappings:
        squared_distances = []

        for pose_index, source_heavy_index in pose_heavy_records:
            reference_index = source_to_reference[
                source_heavy_index
            ]

            pose_point = pose_conformer.GetAtomPosition(pose_index)
            reference_point = reference_conformer.GetAtomPosition(
                reference_index
            )

            squared_distances.append(
                (pose_point.x - reference_point.x) ** 2
                + (pose_point.y - reference_point.y) ** 2
                + (pose_point.z - reference_point.z) ** 2
            )

        rmsd = float(
            np.sqrt(
                np.mean(squared_distances)
            )
        )

        best_rmsd = min(best_rmsd, rmsd)

    return best_rmsd, len(source_reference_mappings)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate GNINA poses after receptor alignment using "
            "Meeko's atom index map instead of SDF bond perception."
        )
    )

    parser.add_argument("--poses", required=True)
    parser.add_argument("--mobile-receptor", required=True)
    parser.add_argument("--reference-receptor", required=True)
    parser.add_argument("--reference-ligand", required=True)
    parser.add_argument("--source-ligand", required=True)
    parser.add_argument("--prepared-ligand-pdbqt", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--prefix",
        default="indexed_aligned_validation",
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
    source_ligand_path = require_file(
        Path(args.source_ligand),
        "Original source ligand",
    )
    prepared_ligand_pdbqt = require_file(
        Path(args.prepared_ligand_pdbqt),
        "Prepared ligand PDBQT",
    )

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    alignment = choose_receptor_alignment(
        mobile_path=mobile_receptor,
        reference_path=reference_receptor,
    )

    source_ligand = load_first_sdf(
        source_ligand_path,
        sanitize=True,
    )
    reference_ligand = load_first_sdf(
        reference_ligand_path,
        sanitize=True,
    )

    # sanitize=False keeps GNINA records readable even when its SDF
    # bond-order assignment creates valence warnings.
    supplier = Chem.SDMolSupplier(
        str(poses_path),
        removeHs=False,
        sanitize=False,
    )

    pose_records = list(supplier)

    example_pose = next(
        (pose for pose in pose_records if pose is not None),
        None,
    )

    if example_pose is None:
        raise RuntimeError("No GNINA poses could be read")

    index_pairs = parse_meeko_index_pairs(
        prepared_ligand_pdbqt
    )

    pdbqt_to_source = choose_pdbqt_to_source_mapping(
        pairs=index_pairs,
        source=source_ligand,
        example_pose=example_pose,
    )

    source_reference = source_to_reference_mappings(
        source=source_ligand,
        reference=reference_ligand,
    )

    source_original_to_heavy = build_original_to_heavy_index(
        source_ligand
    )

    rows: list[dict[str, Any]] = []
    transformed_poses: dict[int, Chem.Mol] = {}
    skipped: dict[str, str] = {}

    for pose_number, molecule in enumerate(
        pose_records,
        start=1,
    ):
        if molecule is None:
            skipped[str(pose_number)] = "RDKit could not read record"
            continue

        try:
            transformed = transform_molecule(
                molecule,
                alignment,
            )

            rmsd, mapping_count = calculate_indexed_rmsd(
                transformed_pose=transformed,
                source=source_ligand,
                reference=reference_ligand,
                pdbqt_to_source=pdbqt_to_source,
                source_original_to_heavy=source_original_to_heavy,
                source_reference_mappings=source_reference,
            )

        except Exception as error:
            skipped[str(pose_number)] = str(error)
            continue

        row: dict[str, Any] = {
            "pose": pose_number,
            "rmsd": rmsd,
            "mapped_heavy_atoms": len(
                source_original_to_heavy
            ),
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
            "No poses could be validated with the index map."
        )

    rows.sort(key=lambda row: row["rmsd"])

    best = rows[0]
    best_pose_number = int(best["pose"])
    pose_one = next(
        (
            row
            for row in rows
            if row["pose"] == 1
        ),
        None,
    )

    csv_path = output_dir / f"{args.prefix}_rmsd.csv"
    json_path = output_dir / f"{args.prefix}_summary.json"
    best_pose_path = output_dir / f"{args.prefix}_best_pose.sdf"

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
                "mapped_heavy_atoms",
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

    summary = {
        "status": "complete",
        "passed": passed,
        "pass_threshold_angstrom": args.pass_threshold,
        "receptor_alignment": {
            "chain_mapping": alignment["chain_mapping"],
            "matched_ca_atoms": alignment["matched_ca_atoms"],
            "ca_rmsd": alignment["rmsd"],
        },
        "sdf_record_count": len(pose_records),
        "validated_pose_count": len(rows),
        "skipped_poses": skipped,
        "mapped_heavy_atoms": len(source_original_to_heavy),
        "reference_atom_mappings": len(source_reference),
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
    print("Chain mapping:", alignment["chain_mapping"])
    print(
        f"Matched CA atoms: {alignment['matched_ca_atoms']}"
    )
    print(
        f"Receptor CA RMSD: {alignment['rmsd']:.3f} Å"
    )

    print(
        f"\nFull heavy-atom mapping: "
        f"{len(source_original_to_heavy)} atoms"
    )
    print(f"Validated poses: {len(rows)}/{len(pose_records)}")

    if skipped:
        print("Skipped poses:", skipped)

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

    print("Result:", "PASS" if passed else "FAIL")

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
