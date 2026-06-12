from __future__ import annotations

import argparse
import csv
import json
import math
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


DEFAULT_WEIGHTS = {
    "cnn_score": 0.25,
    "cnn_affinity": 0.15,
    "catalytic_geometry": 0.35,
    "steric_validity": 0.15,
    "cluster_support": 0.10,
}


def require_file(path: Path, label: str) -> Path:
    path = path.resolve()

    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")

    if path.stat().st_size == 0:
        raise RuntimeError(f"{label} is empty: {path}")

    return path


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


def parse_residue_spec(value: str) -> tuple[str, int]:
    fields = value.split(":", maxsplit=1)

    if len(fields) != 2:
        raise ValueError(
            f"Invalid residue specification: {value}. "
            "Expected CHAIN:NUMBER, such as A:25."
        )

    chain = fields[0].strip()

    if not chain:
        raise ValueError(f"Missing chain in residue specification: {value}")

    try:
        number = int(fields[1])
    except ValueError as error:
        raise ValueError(
            f"Invalid residue number in: {value}"
        ) from error

    return chain, number


def read_receptor_pdbqt(
    path: Path,
) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []

    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue

        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain = line[21].strip()

        try:
            residue_number = int(line[22:26])
            coordinate = np.array(
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ],
                dtype=float,
            )
        except ValueError:
            continue

        fields = line.split()
        atom_type = fields[-1] if fields else ""

        # AutoDock hydrogen atom types begin with H.
        is_hydrogen = atom_type.upper().startswith("H")

        atoms.append(
            {
                "atom_name": atom_name,
                "residue_name": residue_name,
                "chain": chain,
                "residue_number": residue_number,
                "coordinate": coordinate,
                "atom_type": atom_type,
                "is_hydrogen": is_hydrogen,
            }
        )

    if not atoms:
        raise RuntimeError(f"No receptor atoms found in: {path}")

    return atoms


def collect_catalytic_coordinates(
    receptor_atoms: list[dict[str, Any]],
    residue_specs: list[tuple[str, int]],
    catalytic_atom_names: set[str],
) -> list[np.ndarray]:
    groups: list[np.ndarray] = []

    for chain, residue_number in residue_specs:
        coordinates = [
            atom["coordinate"]
            for atom in receptor_atoms
            if atom["chain"] == chain
            and atom["residue_number"] == residue_number
            and atom["atom_name"] in catalytic_atom_names
        ]

        if not coordinates:
            raise RuntimeError(
                "No catalytic atoms found for "
                f"{chain}:{residue_number}; expected one of "
                f"{sorted(catalytic_atom_names)}"
            )

        groups.append(np.asarray(coordinates, dtype=float))

    return groups


def build_pose_coordinate_array(
    pose: Chem.Mol,
    pdbqt_to_source: dict[int, int],
    source_original_to_heavy: dict[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    heavy_count = len(source_original_to_heavy)

    coordinates = np.full(
        (heavy_count, 3),
        np.nan,
        dtype=float,
    )
    atomic_numbers = np.zeros(heavy_count, dtype=int)

    conformer = pose.GetConformer()

    for pose_index, atom in enumerate(pose.GetAtoms()):
        if atom.GetAtomicNum() <= 1:
            continue

        if pose_index not in pdbqt_to_source:
            raise RuntimeError(
                f"No source index for pose atom {pose_index + 1}"
            )

        source_original_index = pdbqt_to_source[pose_index]

        if source_original_index not in source_original_to_heavy:
            raise RuntimeError(
                "A pose heavy atom mapped to a source hydrogen."
            )

        heavy_index = source_original_to_heavy[
            source_original_index
        ]

        point = conformer.GetAtomPosition(pose_index)

        coordinates[heavy_index] = [
            point.x,
            point.y,
            point.z,
        ]
        atomic_numbers[heavy_index] = atom.GetAtomicNum()

    if np.isnan(coordinates).any():
        raise RuntimeError(
            "The atom-index map did not cover every heavy atom."
        )

    return coordinates, atomic_numbers


def gaussian_score(
    distance: float,
    ideal: float,
    sigma: float,
) -> float:
    return math.exp(
        -0.5 * ((distance - ideal) / sigma) ** 2
    )


def calculate_catalytic_geometry(
    ligand_coordinates: np.ndarray,
    atomic_numbers: np.ndarray,
    catalytic_groups: list[np.ndarray],
    ideal_distance: float,
    sigma: float,
    maximum_contact_distance: float,
) -> dict[str, Any]:
    hetero_indices = np.where(
        np.isin(atomic_numbers, [7, 8, 16])
    )[0]

    if not len(hetero_indices):
        return {
            "score": 0.0,
            "best_heteroatom": None,
            "distances": [],
            "maximum_distance": None,
            "passes_contact_gate": False,
        }

    best: dict[str, Any] | None = None

    for heavy_index in hetero_indices:
        coordinate = ligand_coordinates[heavy_index]

        distances = []

        for catalytic_coordinates in catalytic_groups:
            atom_distances = np.linalg.norm(
                catalytic_coordinates - coordinate,
                axis=1,
            )
            distances.append(float(atom_distances.min()))

        individual_scores = [
            gaussian_score(
                distance=distance,
                ideal=ideal_distance,
                sigma=sigma,
            )
            for distance in distances
        ]

        # Geometric mean requires good contact with every catalytic
        # residue rather than excellent contact with only one.
        combined_score = float(
            np.prod(individual_scores)
            ** (1.0 / len(individual_scores))
        )

        maximum_distance = max(distances)
        passes_gate = (
            maximum_distance <= maximum_contact_distance
        )

        if not passes_gate:
            excess = (
                maximum_distance - maximum_contact_distance
            )
            combined_score *= math.exp(
                -0.5 * (excess / sigma) ** 2
            )

        candidate = {
            "score": combined_score,
            "best_heteroatom": int(heavy_index + 1),
            "element": int(atomic_numbers[heavy_index]),
            "distances": distances,
            "maximum_distance": maximum_distance,
            "passes_contact_gate": passes_gate,
        }

        if best is None or candidate["score"] > best["score"]:
            best = candidate

    if best is None:
        raise RuntimeError("Could not score catalytic geometry")

    return best


def calculate_steric_validity(
    ligand_coordinates: np.ndarray,
    receptor_atoms: list[dict[str, Any]],
    excluded_residues: set[tuple[str, int]],
) -> dict[str, Any]:
    receptor_coordinates = np.asarray(
        [
            atom["coordinate"]
            for atom in receptor_atoms
            if not atom["is_hydrogen"]
            and (
                atom["chain"],
                atom["residue_number"],
            )
            not in excluded_residues
        ],
        dtype=float,
    )

    if not len(receptor_coordinates):
        raise RuntimeError(
            "No receptor heavy atoms remained for clash checking."
        )

    difference = (
        ligand_coordinates[:, None, :]
        - receptor_coordinates[None, :, :]
    )
    distances = np.linalg.norm(difference, axis=2)

    hard_clashes = int(np.sum(distances < 1.55))
    soft_clashes = int(
        np.sum(
            (distances >= 1.55)
            & (distances < 1.90)
        )
    )

    minimum_distance = float(distances.min())

    score = math.exp(
        -2.0 * hard_clashes
        -0.20 * soft_clashes
    )

    return {
        "score": score,
        "hard_clashes": hard_clashes,
        "soft_clashes": soft_clashes,
        "minimum_distance": minimum_distance,
    }


def direct_pose_rmsd(
    first: np.ndarray,
    second: np.ndarray,
) -> float:
    return float(
        np.sqrt(
            np.mean(
                np.sum(
                    (first - second) ** 2,
                    axis=1,
                )
            )
        )
    )


def calculate_cluster_support(
    pose_coordinates: list[np.ndarray],
    threshold: float,
) -> tuple[list[int], list[float]]:
    pose_count = len(pose_coordinates)
    cluster_sizes = []

    for first_index in range(pose_count):
        neighbors = 0

        for second_index in range(pose_count):
            rmsd = direct_pose_rmsd(
                pose_coordinates[first_index],
                pose_coordinates[second_index],
            )

            if rmsd <= threshold:
                neighbors += 1

        cluster_sizes.append(neighbors)

    maximum_size = max(cluster_sizes)

    normalized = [
        size / maximum_size
        for size in cluster_sizes
    ]

    return cluster_sizes, normalized


def percentile_scores(
    values: list[float | None],
) -> list[float]:
    valid_indices = [
        index
        for index, value in enumerate(values)
        if value is not None and np.isfinite(value)
    ]

    result = [0.0] * len(values)

    if not valid_indices:
        return result

    if len(valid_indices) == 1:
        result[valid_indices[0]] = 1.0
        return result

    sorted_indices = sorted(
        valid_indices,
        key=lambda index: float(values[index]),
    )

    for rank, original_index in enumerate(sorted_indices):
        result[original_index] = rank / (
            len(sorted_indices) - 1
        )

    return result


def load_validation_rmsd(
    path: Path | None,
) -> dict[int, float]:
    if path is None:
        return {}

    mapping: dict[int, float] = {}

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            mapping[int(row["pose"])] = float(row["rmsd"])

    return mapping


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rank GNINA poses using CNN scores, catalytic geometry, "
            "steric validity, and pose-cluster support."
        )
    )

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
        default="consensus_pose_ranking",
    )

    parser.add_argument(
        "--catalytic-residue",
        action="append",
        required=True,
        help=(
            "Catalytic residue in CHAIN:NUMBER form. "
            "Repeat for every residue in the catalytic group."
        ),
    )

    parser.add_argument(
        "--catalytic-atoms",
        default="OD1,OD2",
    )

    parser.add_argument(
        "--flexres",
        default="",
        help=(
            "Comma-separated flexible residues excluded from rigid "
            "clash checking, such as A:32,B:32,A:50."
        ),
    )

    parser.add_argument(
        "--ideal-contact-distance",
        type=float,
        default=2.9,
    )
    parser.add_argument(
        "--contact-sigma",
        type=float,
        default=0.8,
    )
    parser.add_argument(
        "--maximum-catalytic-distance",
        type=float,
        default=4.0,
    )
    parser.add_argument(
        "--cluster-threshold",
        type=float,
        default=2.0,
    )

    parser.add_argument(
        "--weight-cnn-score",
        type=float,
        default=DEFAULT_WEIGHTS["cnn_score"],
    )
    parser.add_argument(
        "--weight-cnn-affinity",
        type=float,
        default=DEFAULT_WEIGHTS["cnn_affinity"],
    )
    parser.add_argument(
        "--weight-catalytic",
        type=float,
        default=DEFAULT_WEIGHTS[
            "catalytic_geometry"
        ],
    )
    parser.add_argument(
        "--weight-steric",
        type=float,
        default=DEFAULT_WEIGHTS["steric_validity"],
    )
    parser.add_argument(
        "--weight-cluster",
        type=float,
        default=DEFAULT_WEIGHTS["cluster_support"],
    )

    parser.add_argument(
        "--validation-csv",
        default=None,
        help=(
            "Optional experimental RMSD CSV used only to evaluate "
            "the ranking. It is never included in the score."
        ),
    )

    args = parser.parse_args()

    poses_path = require_file(
        Path(args.poses),
        "Pose SDF",
    )
    receptor_path = require_file(
        Path(args.receptor),
        "Receptor",
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

    weights = {
        "cnn_score": args.weight_cnn_score,
        "cnn_affinity": args.weight_cnn_affinity,
        "catalytic_geometry": args.weight_catalytic,
        "steric_validity": args.weight_steric,
        "cluster_support": args.weight_cluster,
    }

    total_weight = sum(weights.values())

    if total_weight <= 0:
        raise ValueError("Ranking weights must sum above zero")

    weights = {
        name: value / total_weight
        for name, value in weights.items()
    }

    receptor_atoms = read_receptor_pdbqt(receptor_path)

    catalytic_specs = [
        parse_residue_spec(value)
        for value in args.catalytic_residue
    ]

    catalytic_atom_names = {
        name.strip()
        for name in args.catalytic_atoms.split(",")
        if name.strip()
    }

    catalytic_groups = collect_catalytic_coordinates(
        receptor_atoms=receptor_atoms,
        residue_specs=catalytic_specs,
        catalytic_atom_names=catalytic_atom_names,
    )

    excluded_residues = {
        parse_residue_spec(value.strip())
        for value in args.flexres.split(",")
        if value.strip()
    }

    source_ligand = load_first_sdf(
        source_ligand_path,
        sanitize=True,
    )

    supplier = Chem.SDMolSupplier(
        str(poses_path),
        removeHs=False,
        sanitize=False,
    )
    poses = list(supplier)

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

    records: list[dict[str, Any]] = []
    valid_poses: list[Chem.Mol] = []
    pose_coordinates: list[np.ndarray] = []

    for pose_number, pose in enumerate(poses, start=1):
        if pose is None:
            continue

        coordinates, atomic_numbers = build_pose_coordinate_array(
            pose=pose,
            pdbqt_to_source=pdbqt_to_source,
            source_original_to_heavy=source_original_to_heavy,
        )

        catalytic = calculate_catalytic_geometry(
            ligand_coordinates=coordinates,
            atomic_numbers=atomic_numbers,
            catalytic_groups=catalytic_groups,
            ideal_distance=args.ideal_contact_distance,
            sigma=args.contact_sigma,
            maximum_contact_distance=(
                args.maximum_catalytic_distance
            ),
        )

        steric = calculate_steric_validity(
            ligand_coordinates=coordinates,
            receptor_atoms=receptor_atoms,
            excluded_residues=excluded_residues,
        )

        records.append(
            {
                "pose": pose_number,
                "CNNscore": get_float_property(
                    pose,
                    "CNNscore",
                ),
                "CNNaffinity": get_float_property(
                    pose,
                    "CNNaffinity",
                ),
                "minimizedAffinity": get_float_property(
                    pose,
                    "minimizedAffinity",
                ),
                "catalytic_score": catalytic["score"],
                "catalytic_heteroatom": catalytic[
                    "best_heteroatom"
                ],
                "catalytic_distances": catalytic[
                    "distances"
                ],
                "catalytic_maximum_distance": catalytic[
                    "maximum_distance"
                ],
                "catalytic_gate_passed": catalytic[
                    "passes_contact_gate"
                ],
                "steric_score": steric["score"],
                "hard_clashes": steric["hard_clashes"],
                "soft_clashes": steric["soft_clashes"],
                "minimum_receptor_distance": steric[
                    "minimum_distance"
                ],
            }
        )

        valid_poses.append(pose)
        pose_coordinates.append(coordinates)

    if not records:
        raise RuntimeError("No poses could be scored")

    cluster_sizes, cluster_scores = (
        calculate_cluster_support(
            pose_coordinates=pose_coordinates,
            threshold=args.cluster_threshold,
        )
    )

    cnn_score_percentiles = percentile_scores(
        [record["CNNscore"] for record in records]
    )
    cnn_affinity_percentiles = percentile_scores(
        [record["CNNaffinity"] for record in records]
    )

    for index, record in enumerate(records):
        record["cnn_score_percentile"] = (
            cnn_score_percentiles[index]
        )
        record["cnn_affinity_percentile"] = (
            cnn_affinity_percentiles[index]
        )
        record["cluster_size"] = cluster_sizes[index]
        record["cluster_score"] = cluster_scores[index]

        consensus_score = (
            weights["cnn_score"]
            * record["cnn_score_percentile"]
            + weights["cnn_affinity"]
            * record["cnn_affinity_percentile"]
            + weights["catalytic_geometry"]
            * record["catalytic_score"]
            + weights["steric_validity"]
            * record["steric_score"]
            + weights["cluster_support"]
            * record["cluster_score"]
        )

        # Strong target-specific penalty: HIV-1 protease poses
        # should place a ligand heteroatom near both Asp25 residues.
        if not record["catalytic_gate_passed"]:
            consensus_score *= 0.50

        record["consensus_score"] = consensus_score

    validation_path = (
        Path(args.validation_csv).resolve()
        if args.validation_csv
        else None
    )
    validation_rmsd = load_validation_rmsd(
        validation_path
    )

    for record in records:
        record["experimental_rmsd"] = validation_rmsd.get(
            int(record["pose"])
        )

    records.sort(
        key=lambda record: record["consensus_score"],
        reverse=True,
    )

    for rank, record in enumerate(records, start=1):
        record["consensus_rank"] = rank

    best = records[0]
    best_pose_number = int(best["pose"])

    pose_lookup = {
        int(record["pose"]): pose
        for record, pose in zip(records, valid_poses)
    }

    # The records were sorted, so rebuild the lookup from original
    # pose numbering instead of relying on list position.
    pose_lookup = {}

    for pose_number, pose in enumerate(poses, start=1):
        if pose is not None:
            pose_lookup[pose_number] = pose

    csv_path = output_dir / f"{args.prefix}_ranking.csv"
    json_path = output_dir / f"{args.prefix}_summary.json"
    best_pose_path = output_dir / f"{args.prefix}_selected_pose.sdf"

    fieldnames = [
        "consensus_rank",
        "pose",
        "consensus_score",
        "CNNscore",
        "CNNaffinity",
        "minimizedAffinity",
        "cnn_score_percentile",
        "cnn_affinity_percentile",
        "catalytic_score",
        "catalytic_heteroatom",
        "catalytic_distances",
        "catalytic_maximum_distance",
        "catalytic_gate_passed",
        "steric_score",
        "hard_clashes",
        "soft_clashes",
        "minimum_receptor_distance",
        "cluster_size",
        "cluster_score",
        "experimental_rmsd",
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
        writer.writerows(records)

    writer = Chem.SDWriter(str(best_pose_path))
    writer.write(pose_lookup[best_pose_number])
    writer.close()

    summary = {
        "status": "complete",
        "pose_count": len(records),
        "weights": weights,
        "settings": {
            "catalytic_residues": args.catalytic_residue,
            "catalytic_atoms": sorted(
                catalytic_atom_names
            ),
            "ideal_contact_distance": (
                args.ideal_contact_distance
            ),
            "maximum_catalytic_distance": (
                args.maximum_catalytic_distance
            ),
            "cluster_threshold": args.cluster_threshold,
            "flexres_excluded_from_clash_check": sorted(
                f"{chain}:{number}"
                for chain, number in excluded_residues
            ),
            "minimized_affinity_used": False,
        },
        "selected_pose": best,
        "top_10": records[:10],
        "files": {
            "ranking_csv": str(csv_path),
            "summary_json": str(json_path),
            "selected_pose_sdf": str(best_pose_path),
        },
    }

    json_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nTop 10 consensus-ranked poses:\n")

    for record in records[:10]:
        rmsd_text = (
            f"{record['experimental_rmsd']:.3f} Å"
            if record["experimental_rmsd"] is not None
            else "N/A"
        )

        distances = ", ".join(
            f"{distance:.2f}"
            for distance in record["catalytic_distances"]
        )

        print(
            f"Rank {record['consensus_rank']:2d} | "
            f"Pose {record['pose']:2d} | "
            f"Score={record['consensus_score']:.3f} | "
            f"CNN={record['CNNscore']:.3f} | "
            f"Dyad={record['catalytic_score']:.3f} "
            f"({distances} Å) | "
            f"Clashes={record['hard_clashes']} | "
            f"Cluster={record['cluster_size']} | "
            f"Known RMSD={rmsd_text}"
        )

    print("\n=== CONSENSUS SELECTION ===")
    print(f"Selected pose: {best_pose_number}")
    print(
        f"Consensus score: "
        f"{best['consensus_score']:.3f}"
    )

    if best["experimental_rmsd"] is not None:
        print(
            f"Selected-pose experimental RMSD: "
            f"{best['experimental_rmsd']:.3f} Å"
        )

        known_rows = [
            record
            for record in records
            if record["experimental_rmsd"] is not None
        ]

        true_best = min(
            known_rows,
            key=lambda record: record[
                "experimental_rmsd"
            ],
        )

        print(
            f"Known closest pose: {true_best['pose']} "
            f"({true_best['experimental_rmsd']:.3f} Å)"
        )
        print(
            f"Consensus rank of known closest pose: "
            f"{true_best['consensus_rank']}"
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
