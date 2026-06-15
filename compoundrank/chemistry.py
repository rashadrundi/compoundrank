from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D


def load_first_sdf(path: Path, sanitize: bool = True) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
        sanitize=sanitize,
        strictParsing=False,
    )
    molecule = next((mol for mol in supplier if mol is not None), None)
    if molecule is None:
        raise RuntimeError(f"Could not read molecule: {path}")
    return molecule


def load_sdf_records(path: Path, sanitize: bool = False) -> list[Chem.Mol]:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
        sanitize=sanitize,
        strictParsing=False,
    )
    return [mol for mol in supplier if mol is not None]


def parse_meeko_index_pairs(pdbqt_path: Path) -> list[tuple[int, int]]:
    numbers: list[int] = []
    for line in pdbqt_path.read_text(errors="replace").splitlines():
        if not line.startswith("REMARK INDEX MAP"):
            continue
        for field in line.split()[3:]:
            numbers.append(int(field))
    if not numbers:
        raise RuntimeError(f"No REMARK INDEX MAP in {pdbqt_path}")
    if len(numbers) % 2:
        raise RuntimeError("Meeko index map contains an odd number of integers")
    return [
        (numbers[index], numbers[index + 1])
        for index in range(0, len(numbers), 2)
    ]


def choose_pose_to_source_mapping(
    pairs: list[tuple[int, int]],
    source: Chem.Mol,
    example_pose: Chem.Mol,
) -> dict[int, int]:
    candidates = [
        {second - 1: first - 1 for first, second in pairs},
        {first - 1: second - 1 for first, second in pairs},
    ]
    scored: list[tuple[int, int, dict[int, int]]] = []
    for mapping in candidates:
        compared = 0
        matches = 0
        for pose_index, pose_atom in enumerate(example_pose.GetAtoms()):
            source_index = mapping.get(pose_index)
            if source_index is None:
                continue
            if source_index < 0 or source_index >= source.GetNumAtoms():
                continue
            compared += 1
            if source.GetAtomWithIdx(source_index).GetAtomicNum() == pose_atom.GetAtomicNum():
                matches += 1
        scored.append((matches, compared, mapping))
    matches, compared, mapping = max(scored, key=lambda item: (item[0], item[1]))
    if compared != example_pose.GetNumAtoms() or matches != compared:
        raise RuntimeError(
            "Could not align GNINA atom order with the Meeko index map: "
            f"elements={matches}/{compared}, pose_atoms={example_pose.GetNumAtoms()}"
        )
    return mapping


def original_to_heavy_index(source: Chem.Mol) -> dict[int, int]:
    mapping: dict[int, int] = {}
    heavy_index = 0
    for atom in source.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        mapping[atom.GetIdx()] = heavy_index
        heavy_index += 1
    return mapping


def reconstruct_heavy_pose(
    source: Chem.Mol,
    docked_pose: Chem.Mol,
    pose_to_source: dict[int, int],
) -> Chem.Mol:
    source_heavy = Chem.RemoveHs(source, sanitize=False)
    original_to_heavy = original_to_heavy_index(source)
    result = Chem.Mol(source_heavy)
    result.RemoveAllConformers()
    conformer = Chem.Conformer(result.GetNumAtoms())
    conformer.Set3D(True)
    docked_conformer = docked_pose.GetConformer()
    assigned: set[int] = set()

    for pose_index, pose_atom in enumerate(docked_pose.GetAtoms()):
        if pose_atom.GetAtomicNum() <= 1:
            continue
        source_original_index = pose_to_source.get(pose_index)
        if source_original_index is None:
            raise RuntimeError(f"No source mapping for pose atom {pose_index + 1}")
        heavy_index = original_to_heavy.get(source_original_index)
        if heavy_index is None:
            raise RuntimeError(
                f"Pose heavy atom {pose_index + 1} mapped to a source hydrogen"
            )
        source_atom = result.GetAtomWithIdx(heavy_index)
        if source_atom.GetAtomicNum() != pose_atom.GetAtomicNum():
            raise RuntimeError(
                f"Element mismatch at pose atom {pose_index + 1}: "
                f"{pose_atom.GetSymbol()} vs {source_atom.GetSymbol()}"
            )
        point = docked_conformer.GetAtomPosition(pose_index)
        conformer.SetAtomPosition(
            heavy_index,
            Point3D(float(point.x), float(point.y), float(point.z)),
        )
        assigned.add(heavy_index)

    if len(assigned) != result.GetNumAtoms():
        missing = sorted(set(range(result.GetNumAtoms())) - assigned)
        raise RuntimeError(
            f"Only reconstructed {len(assigned)}/{result.GetNumAtoms()} heavy atoms; "
            f"missing={missing}"
        )

    result.AddConformer(conformer, assignId=True)
    Chem.SanitizeMol(result)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    return result


def heavy_coordinates(molecule: Chem.Mol) -> np.ndarray:
    conformer = molecule.GetConformer()
    coordinates = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        point = conformer.GetAtomPosition(atom.GetIdx())
        coordinates.append((point.x, point.y, point.z))
    return np.asarray(coordinates, dtype=float)


def direct_heavy_rmsd(first: Chem.Mol, second: Chem.Mol) -> float:
    first_coordinates = heavy_coordinates(first)
    second_coordinates = heavy_coordinates(second)
    if first_coordinates.shape != second_coordinates.shape:
        raise ValueError(
            f"Pose coordinate shapes differ: {first_coordinates.shape} vs "
            f"{second_coordinates.shape}"
        )
    return float(
        np.sqrt(
            np.mean(np.sum((first_coordinates - second_coordinates) ** 2, axis=1))
        )
    )


def write_sdf(path: Path, molecules: Iterable[Chem.Mol]) -> int:
    writer = Chem.SDWriter(str(path))
    count = 0
    try:
        for molecule in molecules:
            writer.write(molecule)
            count += 1
    finally:
        writer.close()
    return count
