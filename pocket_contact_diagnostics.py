from __future__ import annotations

from pathlib import Path

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


POSES_PATH = Path(
    "output/hiv1_protease_dimer/docking/predicted_flexible/"
    "darunavir_predicted_flexible_poses.sdf"
)

RECEPTOR_PATH = Path(
    "output/hiv1_protease_dimer/docking/receptor/"
    "hiv1_protease_dimer_prepared.pdbqt"
)

SOURCE_LIGAND_PATH = Path(
    "output/hiv1_protease_dimer/docking/ligand/"
    "darunavir_pH7_4.sdf"
)

PREPARED_LIGAND_PATH = Path(
    "output/hiv1_protease_dimer/docking/ligand/"
    "darunavir_prepared.pdbqt"
)

POSE_NUMBERS = (1, 2, 7, 14, 18)

CONTACT_CUTOFF = 4.0
HBOND_MINIMUM = 2.2
HBOND_MAXIMUM = 3.6

FLEXIBLE_RESIDUES = {
    ("A", 32),
    ("B", 32),
    ("A", 50),
    ("B", 50),
    ("A", 82),
    ("B", 82),
}

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


def residue_contact_distances(
    ligand_coordinates: np.ndarray,
    receptor_atoms: list[dict],
) -> dict[tuple[str, int], float]:
    result: dict[tuple[str, int], float] = {}

    for chain, residue_number in ALL_POCKET_RESIDUES:
        residue_key = (chain, residue_number)

        atom_coordinates = []

        for atom in receptor_atoms:
            if atom["is_hydrogen"]:
                continue

            if (
                atom["chain"] != chain
                or atom["residue_number"] != residue_number
            ):
                continue

            # Flexible side-chain coordinates in the original receptor
            # are no longer authoritative. Use their fixed backbone only.
            if (
                residue_key in FLEXIBLE_RESIDUES
                and atom["atom_name"] not in BACKBONE_ATOMS
            ):
                continue

            atom_coordinates.append(atom["coordinate"])

        if not atom_coordinates:
            continue

        receptor_coordinates = np.asarray(atom_coordinates)

        distances = np.linalg.norm(
            ligand_coordinates[:, None, :]
            - receptor_coordinates[None, :, :],
            axis=2,
        )

        result[residue_key] = float(distances.min())

    return result


def backbone_hbond_candidates(
    ligand_coordinates: np.ndarray,
    atomic_numbers: np.ndarray,
    receptor_atoms: list[dict],
) -> list[dict]:
    ligand_hetero_indices = np.where(
        np.isin(atomic_numbers, [7, 8])
    )[0]

    contacts = []

    for atom in receptor_atoms:
        if atom["is_hydrogen"]:
            continue

        if atom["atom_name"] not in {"N", "O"}:
            continue

        residue_key = (
            atom["chain"],
            atom["residue_number"],
        )

        if residue_key not in ALL_POCKET_RESIDUES:
            continue

        for heavy_index in ligand_hetero_indices:
            distance = float(
                np.linalg.norm(
                    ligand_coordinates[heavy_index]
                    - atom["coordinate"]
                )
            )

            if HBOND_MINIMUM <= distance <= HBOND_MAXIMUM:
                contacts.append({
                    "chain": atom["chain"],
                    "residue_number": atom["residue_number"],
                    "residue_name": atom["residue_name"],
                    "receptor_atom": atom["atom_name"],
                    "ligand_heavy_atom": int(heavy_index + 1),
                    "distance": distance,
                })

    contacts.sort(key=lambda item: item["distance"])
    return contacts


source_ligand = load_first_sdf(
    SOURCE_LIGAND_PATH,
    sanitize=True,
)

poses = list(
    Chem.SDMolSupplier(
        str(POSES_PATH),
        removeHs=False,
        sanitize=False,
    )
)

example_pose = next(
    pose for pose in poses
    if pose is not None
)

index_pairs = parse_meeko_index_pairs(
    PREPARED_LIGAND_PATH
)

pdbqt_to_source = choose_pdbqt_to_source_mapping(
    pairs=index_pairs,
    source=source_ligand,
    example_pose=example_pose,
)

source_original_to_heavy = build_original_to_heavy_index(
    source_ligand
)

receptor_atoms = read_receptor_pdbqt(
    RECEPTOR_PATH
)

for pose_number in POSE_NUMBERS:
    pose = poses[pose_number - 1]

    if pose is None:
        print(f"\nPose {pose_number}: unreadable")
        continue

    coordinates, atomic_numbers = build_pose_coordinate_array(
        pose=pose,
        pdbqt_to_source=pdbqt_to_source,
        source_original_to_heavy=source_original_to_heavy,
    )

    residue_distances = residue_contact_distances(
        ligand_coordinates=coordinates,
        receptor_atoms=receptor_atoms,
    )

    hbond_candidates = backbone_hbond_candidates(
        ligand_coordinates=coordinates,
        atomic_numbers=atomic_numbers,
        receptor_atoms=receptor_atoms,
    )

    contacted_residues = {
        residue
        for residue, distance in residue_distances.items()
        if distance <= CONTACT_CUTOFF
    }

    chain_a_contacts = sum(
        chain == "A"
        for chain, _ in contacted_residues
    )
    chain_b_contacts = sum(
        chain == "B"
        for chain, _ in contacted_residues
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

    print(f"\n=== Pose {pose_number} ===")
    print(
        f"Total pocket contacts: "
        f"{len(contacted_residues)}/{len(ALL_POCKET_RESIDUES)}"
    )
    print(
        f"Chain contacts: A={chain_a_contacts}, "
        f"B={chain_b_contacts}, "
        f"balance={chain_balance:.3f}"
    )

    for group_name, residues in POCKET_GROUPS.items():
        contacted = sorted(
            residue
            for residue in residues
            if residue in contacted_residues
        )

        print(
            f"{group_name}: "
            f"{len(contacted)}/{len(residues)} "
            f"{contacted}"
        )

    print(
        f"Backbone N/O contact candidates: "
        f"{len(hbond_candidates)}"
    )

    for contact in hbond_candidates[:12]:
        print(
            f"  {contact['chain']}:"
            f"{contact['residue_name']}"
            f"{contact['residue_number']} "
            f"{contact['receptor_atom']} -> "
            f"ligand heavy #{contact['ligand_heavy_atom']} "
            f"{contact['distance']:.2f} Å"
        )
