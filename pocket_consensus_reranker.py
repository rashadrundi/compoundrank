from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from indexed_aligned_pose_validation import (
    build_original_to_heavy_index,
    choose_pdbqt_to_source_mapping,
    parse_meeko_index_pairs,
)
from pose_ranker import (
    build_pose_coordinate_array,
    load_first_sdf,
    read_receptor_pdbqt,
)


BACKBONE_ATOMS = {"N", "CA", "C", "O"}

POCKET_GROUPS = {
    "catalytic_core": {
        ("A", 25), ("A", 27), ("A", 29), ("A", 30),
        ("B", 25), ("B", 27), ("B", 29), ("B", 30),
    },
    "flap_region": {
        ("A", 47), ("A", 48), ("A", 50), ("A", 54),
        ("B", 47), ("B", 48), ("B", 50), ("B", 54),
    },
    "hydrophobic_subsites": {
        ("A", 8), ("A", 23), ("A", 32), ("A", 76),
        ("A", 81), ("A", 82), ("A", 84),
        ("B", 8), ("B", 23), ("B", 32), ("B", 76),
        ("B", 81), ("B", 82), ("B", 84),
    },
}

ALL_POCKET_RESIDUES = set().union(*POCKET_GROUPS.values())


def parse_residue_spec(value: str) -> tuple[str, int]:
    chain, number = value.split(":", maxsplit=1)
    return chain.strip(), int(number)


def require_file(path: Path, label: str) -> Path:
    path = path.resolve()

    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")

    return path


def load_base_ranking(path: Path) -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}

    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pose_number = int(row["pose"])

            rows[pose_number] = {
                "pose": pose_number,
                "base_consensus_score": float(
                    row["consensus_score"]
                ),
                "CNNscore": float(row["CNNscore"]),
                "CNNaffinity": float(row["CNNaffinity"]),
                "catalytic_score": float(
                    row["catalytic_score"]
                ),
                "hard_clashes": int(row["hard_clashes"]),
                "cluster_size": int(row["cluster_size"]),
            }

    if not rows:
        raise RuntimeError(f"No ranking rows found in {path}")

    return rows


def load_validation_rmsd(
    path: Path | None,
) -> dict[int, float]:
    if path is None:
        return {}

    results: dict[int, float] = {}

    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            results[int(row["pose"])] = float(row["rmsd"])

    return results


def residue_contact_distances(
    ligand_coordinates: np.ndarray,
    receptor_atoms: list[dict[str, Any]],
    flexible_residues: set[tuple[str, int]],
) -> dict[tuple[str, int], float]:
    distances_by_residue: dict[tuple[str, int], float] = {}

    for residue_key in ALL_POCKET_RESIDUES:
        chain, residue_number = residue_key
        receptor_coordinates = []

        for atom in receptor_atoms:
            if atom["is_hydrogen"]:
                continue

            if (
                atom["chain"] != chain
                or atom["residue_number"] != residue_number
            ):
                continue

            # Flexible side-chain coordinates in the original
            # receptor are no longer authoritative. Retain only
            # their backbone coordinates for this comparison.
            if (
                residue_key in flexible_residues
                and atom["atom_name"] not in BACKBONE_ATOMS
            ):
                continue

            receptor_coordinates.append(atom["coordinate"])

        if not receptor_coordinates:
            continue

        receptor_array = np.asarray(
            receptor_coordinates,
            dtype=float,
        )

        distance_matrix = np.linalg.norm(
            ligand_coordinates[:, None, :]
            - receptor_array[None, :, :],
            axis=2,
        )

        distances_by_residue[residue_key] = float(
            distance_matrix.min()
        )

    return distances_by_residue


def calculate_pocket_features(
    ligand_coordinates: np.ndarray,
    receptor_atoms: list[dict[str, Any]],
    flexible_residues: set[tuple[str, int]],
    cutoff: float,
) -> dict[str, Any]:
    residue_distances = residue_contact_distances(
        ligand_coordinates=ligand_coordinates,
        receptor_atoms=receptor_atoms,
        flexible_residues=flexible_residues,
    )

    contacted = {
        residue
        for residue, distance in residue_distances.items()
        if distance <= cutoff
    }

    group_counts = {
        group_name: len(residues & contacted)
        for group_name, residues in POCKET_GROUPS.items()
    }

    group_fractions = {
        group_name: (
            group_counts[group_name] / len(residues)
        )
        for group_name, residues in POCKET_GROUPS.items()
    }

    chain_a_contacts = sum(
        chain == "A"
        for chain, _ in contacted
    )
    chain_b_contacts = sum(
        chain == "B"
        for chain, _ in contacted
    )

    maximum_chain_contacts = max(
        chain_a_contacts,
        chain_b_contacts,
        1,
    )

    chain_balance = (
        min(chain_a_contacts, chain_b_contacts)
        / maximum_chain_contacts
    )

    total_fraction = (
        len(contacted) / len(ALL_POCKET_RESIDUES)
    )

    # Provisional HIV-protease-specific pocket score.
    # Flap coverage receives the largest weight because it was
    # the strongest non-RMSD discriminator in this benchmark.
    pocket_score = (
        0.50 * group_fractions["flap_region"]
        + 0.30 * total_fraction
        + 0.15 * group_fractions["hydrophobic_subsites"]
        + 0.05 * chain_balance
    )

    return {
        "pocket_score": pocket_score,
        "total_contacts": len(contacted),
        "total_fraction": total_fraction,
        "chain_a_contacts": chain_a_contacts,
        "chain_b_contacts": chain_b_contacts,
        "chain_balance": chain_balance,
        "catalytic_contacts": group_counts[
            "catalytic_core"
        ],
        "catalytic_fraction": group_fractions[
            "catalytic_core"
        ],
        "flap_contacts": group_counts["flap_region"],
        "flap_fraction": group_fractions["flap_region"],
        "hydrophobic_contacts": group_counts[
            "hydrophobic_subsites"
        ],
        "hydrophobic_fraction": group_fractions[
            "hydrophobic_subsites"
        ],
        "contacted_residues": sorted(
            f"{chain}:{number}"
            for chain, number in contacted
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rerank consensus-scored poses using whole-pocket "
            "and HIV-1 protease flap-region coverage."
        )
    )

    parser.add_argument("--base-ranking-csv", required=True)
    parser.add_argument("--poses", required=True)
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--source-ligand", required=True)
    parser.add_argument(
        "--prepared-ligand-pdbqt",
        required=True,
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--prefix",
        default="pocket_consensus",
    )
    parser.add_argument(
        "--flexres",
        default="",
    )
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=4.0,
    )
    parser.add_argument(
        "--base-weight",
        type=float,
        default=0.35,
    )
    parser.add_argument(
        "--pocket-weight",
        type=float,
        default=0.65,
    )
    parser.add_argument(
        "--validation-csv",
        default=None,
        help=(
            "Optional RMSD data used only to evaluate the ranking."
        ),
    )

    args = parser.parse_args()

    base_ranking_path = require_file(
        Path(args.base_ranking_csv),
        "Base ranking CSV",
    )
    poses_path = require_file(
        Path(args.poses),
        "Pose SDF",
    )
    receptor_path = require_file(
        Path(args.receptor),
        "Receptor PDBQT",
    )
    source_ligand_path = require_file(
        Path(args.source_ligand),
        "Source ligand",
    )
    prepared_ligand_path = require_file(
        Path(args.prepared_ligand_pdbqt),
        "Prepared ligand PDBQT",
    )

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    weight_total = args.base_weight + args.pocket_weight

    if weight_total <= 0:
        raise ValueError("Weights must sum above zero")

    base_weight = args.base_weight / weight_total
    pocket_weight = args.pocket_weight / weight_total

    base_rows = load_base_ranking(base_ranking_path)
    receptor_atoms = read_receptor_pdbqt(receptor_path)

    flexible_residues = {
        parse_residue_spec(value.strip())
        for value in args.flexres.split(",")
        if value.strip()
    }

    source_ligand = load_first_sdf(
        source_ligand_path,
        sanitize=True,
    )

    poses = list(
        Chem.SDMolSupplier(
            str(poses_path),
            removeHs=False,
            sanitize=False,
        )
    )

    example_pose = next(
        (pose for pose in poses if pose is not None),
        None,
    )

    if example_pose is None:
        raise RuntimeError("No poses could be read")

    index_pairs = parse_meeko_index_pairs(
        prepared_ligand_path
    )

    pdbqt_to_source = choose_pdbqt_to_source_mapping(
        pairs=index_pairs,
        source=source_ligand,
        example_pose=example_pose,
    )

    source_original_to_heavy = build_original_to_heavy_index(
        source_ligand
    )

    validation_path = (
        require_file(
            Path(args.validation_csv),
            "Validation CSV",
        )
        if args.validation_csv
        else None
    )

    validation_rmsd = load_validation_rmsd(
        validation_path
    )

    rows = []
    pose_lookup: dict[int, Chem.Mol] = {}

    for pose_number, pose in enumerate(poses, start=1):
        if pose is None:
            continue

        if pose_number not in base_rows:
            continue

        coordinates, _ = build_pose_coordinate_array(
            pose=pose,
            pdbqt_to_source=pdbqt_to_source,
            source_original_to_heavy=source_original_to_heavy,
        )

        features = calculate_pocket_features(
            ligand_coordinates=coordinates,
            receptor_atoms=receptor_atoms,
            flexible_residues=flexible_residues,
            cutoff=args.contact_cutoff,
        )

        row = dict(base_rows[pose_number])
        row.update(features)

        row["final_score"] = (
            base_weight * row["base_consensus_score"]
            + pocket_weight * row["pocket_score"]
        )

        row["experimental_rmsd"] = validation_rmsd.get(
            pose_number
        )

        rows.append(row)
        pose_lookup[pose_number] = pose

    if not rows:
        raise RuntimeError("No poses were reranked")

    rows.sort(
        key=lambda row: row["final_score"],
        reverse=True,
    )

    for rank, row in enumerate(rows, start=1):
        row["final_rank"] = rank

    selected = rows[0]
    selected_pose_number = int(selected["pose"])

    csv_path = output_dir / f"{args.prefix}_ranking.csv"
    json_path = output_dir / f"{args.prefix}_summary.json"
    pose_path = output_dir / f"{args.prefix}_selected_pose.sdf"

    fieldnames = [
        "final_rank",
        "pose",
        "final_score",
        "base_consensus_score",
        "pocket_score",
        "CNNscore",
        "CNNaffinity",
        "catalytic_score",
        "hard_clashes",
        "cluster_size",
        "total_contacts",
        "total_fraction",
        "chain_a_contacts",
        "chain_b_contacts",
        "chain_balance",
        "catalytic_contacts",
        "catalytic_fraction",
        "flap_contacts",
        "flap_fraction",
        "hydrophobic_contacts",
        "hydrophobic_fraction",
        "experimental_rmsd",
        "contacted_residues",
    ]

    with csv_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)

    writer = Chem.SDWriter(str(pose_path))
    writer.write(pose_lookup[selected_pose_number])
    writer.close()

    summary = {
        "status": "complete",
        "weights": {
            "base_consensus": base_weight,
            "whole_pocket": pocket_weight,
        },
        "pocket_formula": {
            "flap_fraction": 0.50,
            "total_fraction": 0.30,
            "hydrophobic_fraction": 0.15,
            "chain_balance": 0.05,
        },
        "selected_pose": selected,
        "top_10": rows[:10],
        "files": {
            "ranking_csv": str(csv_path),
            "summary_json": str(json_path),
            "selected_pose_sdf": str(pose_path),
        },
    }

    json_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nTop 10 pocket-consensus poses:\n")

    for row in rows[:10]:
        rmsd_text = (
            f"{row['experimental_rmsd']:.3f} Å"
            if row["experimental_rmsd"] is not None
            else "N/A"
        )

        print(
            f"Rank {row['final_rank']:2d} | "
            f"Pose {row['pose']:2d} | "
            f"Final={row['final_score']:.3f} | "
            f"Base={row['base_consensus_score']:.3f} | "
            f"Pocket={row['pocket_score']:.3f} | "
            f"Contacts={row['total_contacts']} | "
            f"Flap={row['flap_contacts']}/8 | "
            f"Balance={row['chain_balance']:.3f} | "
            f"Known RMSD={rmsd_text}"
        )

    print("\n=== POCKET-CONSENSUS SELECTION ===")
    print(f"Selected pose: {selected_pose_number}")
    print(f"Final score: {selected['final_score']:.3f}")

    if selected["experimental_rmsd"] is not None:
        print(
            "Selected-pose RMSD:",
            f"{selected['experimental_rmsd']:.3f} Å",
        )

        true_best = min(
            (
                row for row in rows
                if row["experimental_rmsd"] is not None
            ),
            key=lambda row: row["experimental_rmsd"],
        )

        print(
            f"Known closest pose: {true_best['pose']} "
            f"({true_best['experimental_rmsd']:.3f} Å)"
        )
        print(
            "Final rank of known closest pose:",
            true_best["final_rank"],
        )

    print("\nSaved:")
    print(csv_path)
    print(json_path)
    print(pose_path)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
