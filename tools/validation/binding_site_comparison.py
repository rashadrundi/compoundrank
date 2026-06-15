from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem

from aligned_pose_validation import choose_receptor_alignment


def read_pdb_atoms(path: Path) -> dict[tuple[str, int, str], dict[str, np.ndarray]]:
    residues: dict[
        tuple[str, int, str],
        dict[str, np.ndarray],
    ] = {}

    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue

        altloc = line[16].strip()

        if altloc not in {"", "A"}:
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

        if atom_name.startswith("H"):
            continue

        key = (chain, residue_number, residue_name)
        residues.setdefault(key, {})[atom_name] = coordinate

    return residues


def load_ligand_coordinates(path: Path) -> np.ndarray:
    molecule = next(
        (
            mol
            for mol in Chem.SDMolSupplier(
                str(path),
                removeHs=False,
            )
            if mol is not None
        ),
        None,
    )

    if molecule is None:
        raise RuntimeError(f"Could not read ligand: {path}")

    conformer = molecule.GetConformer()
    coordinates = []

    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue

        point = conformer.GetAtomPosition(atom.GetIdx())
        coordinates.append([point.x, point.y, point.z])

    return np.asarray(coordinates, dtype=float)


def transform_coordinate(
    coordinate: np.ndarray,
    alignment: dict[str, Any],
) -> np.ndarray:
    return (
        (coordinate - alignment["mobile_center"])
        @ alignment["rotation"].T
        + alignment["reference_center"]
    )


def rmsd(
    mobile: list[np.ndarray],
    reference: list[np.ndarray],
) -> float:
    mobile_array = np.asarray(mobile)
    reference_array = np.asarray(reference)

    return float(
        np.sqrt(
            np.mean(
                np.sum(
                    (mobile_array - reference_array) ** 2,
                    axis=1,
                )
            )
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare a predicted receptor binding site against an "
            "experimental receptor after CA alignment."
        )
    )

    parser.add_argument("--predicted-receptor", required=True)
    parser.add_argument("--reference-receptor", required=True)
    parser.add_argument("--reference-ligand", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cutoff", type=float, default=5.0)

    args = parser.parse_args()

    predicted_path = Path(args.predicted_receptor)
    reference_path = Path(args.reference_receptor)
    ligand_path = Path(args.reference_ligand)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    alignment = choose_receptor_alignment(
        mobile_path=predicted_path,
        reference_path=reference_path,
    )

    predicted = read_pdb_atoms(predicted_path)
    reference = read_pdb_atoms(reference_path)
    ligand_coordinates = load_ligand_coordinates(ligand_path)

    contact_residues = []

    for residue_key, atoms in reference.items():
        minimum_distance = float("inf")

        for coordinate in atoms.values():
            distances = np.linalg.norm(
                ligand_coordinates - coordinate,
                axis=1,
            )
            minimum_distance = min(
                minimum_distance,
                float(distances.min()),
            )

        if minimum_distance <= args.cutoff:
            contact_residues.append(
                (residue_key, minimum_distance)
            )

    chain_mapping = alignment["chain_mapping"]
    reverse_mapping = {
        reference_chain: predicted_chain
        for predicted_chain, reference_chain
        in chain_mapping.items()
    }

    rows = []
    all_predicted_coordinates = []
    all_reference_coordinates = []

    for reference_key, contact_distance in contact_residues:
        reference_chain, residue_number, residue_name = reference_key

        predicted_chain = reverse_mapping.get(reference_chain)

        if predicted_chain is None:
            continue

        predicted_candidates = [
            key
            for key in predicted
            if key[0] == predicted_chain
            and key[1] == residue_number
        ]

        if not predicted_candidates:
            continue

        predicted_key = predicted_candidates[0]

        reference_atoms = reference[reference_key]
        predicted_atoms = predicted[predicted_key]

        common_atoms = sorted(
            set(reference_atoms) & set(predicted_atoms)
        )

        if not common_atoms:
            continue

        mobile_coordinates = []
        reference_coordinates = []

        for atom_name in common_atoms:
            transformed = transform_coordinate(
                predicted_atoms[atom_name],
                alignment,
            )

            mobile_coordinates.append(transformed)
            reference_coordinates.append(
                reference_atoms[atom_name]
            )

            all_predicted_coordinates.append(transformed)
            all_reference_coordinates.append(
                reference_atoms[atom_name]
            )

        residue_rmsd = rmsd(
            mobile_coordinates,
            reference_coordinates,
        )

        ca_distance = None

        if "CA" in common_atoms:
            transformed_ca = transform_coordinate(
                predicted_atoms["CA"],
                alignment,
            )

            ca_distance = float(
                np.linalg.norm(
                    transformed_ca - reference_atoms["CA"]
                )
            )

        rows.append({
            "reference_chain": reference_chain,
            "predicted_chain": predicted_chain,
            "residue_number": residue_number,
            "residue_name": residue_name,
            "ligand_contact_distance": contact_distance,
            "matched_heavy_atoms": len(common_atoms),
            "residue_heavy_atom_rmsd": residue_rmsd,
            "ca_distance": ca_distance,
            "matched_atom_names": ",".join(common_atoms),
        })

    rows.sort(
        key=lambda row: row["residue_heavy_atom_rmsd"],
        reverse=True,
    )

    overall_pocket_rmsd = rmsd(
        all_predicted_coordinates,
        all_reference_coordinates,
    )

    csv_path = output_dir / "binding_site_comparison.csv"
    json_path = output_dir / "binding_site_comparison.json"

    with csv_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "reference_chain",
                "predicted_chain",
                "residue_number",
                "residue_name",
                "ligand_contact_distance",
                "matched_heavy_atoms",
                "residue_heavy_atom_rmsd",
                "ca_distance",
                "matched_atom_names",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "contact_cutoff_angstrom": args.cutoff,
        "receptor_alignment": {
            "chain_mapping": alignment["chain_mapping"],
            "matched_ca_atoms": alignment["matched_ca_atoms"],
            "ca_rmsd": alignment["rmsd"],
        },
        "contact_residue_count": len(rows),
        "overall_binding_site_heavy_atom_rmsd": overall_pocket_rmsd,
        "residues_ranked_by_difference": rows,
        "files": {
            "csv": str(csv_path),
            "json": str(json_path),
        },
    }

    json_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nReceptor alignment:")
    print(
        f"Overall CA RMSD: {alignment['rmsd']:.3f} Å"
    )
    print(
        f"Binding-site heavy-atom RMSD: "
        f"{overall_pocket_rmsd:.3f} Å"
    )
    print(f"Contact residues: {len(rows)}")

    print("\nResidues with the largest differences:\n")

    for row in rows[:15]:
        ca_text = (
            f"{row['ca_distance']:.3f}"
            if row["ca_distance"] is not None
            else "N/A"
        )

        print(
            f"{row['reference_chain']}:"
            f"{row['residue_name']}{row['residue_number']:02d}  "
            f"heavy-atom RMSD={row['residue_heavy_atom_rmsd']:.3f} Å  "
            f"CA shift={ca_text} Å  "
            f"ligand distance={row['ligand_contact_distance']:.3f} Å"
        )

    print("\nSaved:")
    print(csv_path)
    print(json_path)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
