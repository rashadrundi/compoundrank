from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem

from .models import InteractionEvidence


HYDROPHOBIC_RESIDUES = {
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PHE",
    "TRP",
    "TYR",
    "PRO",
}
POLAR_ELEMENTS = {7, 8, 16}


@dataclass(frozen=True)
class ProteinAtom:
    chain: str
    residue_name: str
    residue_number: int
    atom_name: str
    element: str
    coordinate: np.ndarray

    @property
    def residue_label(self) -> str:
        chain = self.chain or "_"
        return f"{self.residue_name}{self.residue_number}:{chain}"


def _infer_element(line: str, atom_name: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    stripped = "".join(character for character in atom_name if character.isalpha())
    if not stripped:
        return ""
    if len(stripped) >= 2 and stripped[:2].upper() in {"CL", "BR", "FE", "ZN", "MG"}:
        return stripped[:2].upper()
    return stripped[0].upper()


def read_protein_atoms(path: Path) -> list[ProteinAtom]:
    atoms: list[ProteinAtom] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name.startswith("H"):
            continue
        try:
            coordinate = np.asarray(
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ],
                dtype=float,
            )
            residue_number = int(line[22:26])
        except ValueError:
            continue
        atoms.append(
            ProteinAtom(
                chain=line[21].strip(),
                residue_name=line[17:20].strip(),
                residue_number=residue_number,
                atom_name=atom_name,
                element=_infer_element(line, atom_name),
                coordinate=coordinate,
            )
        )
    if not atoms:
        raise RuntimeError(f"No protein atoms found in receptor: {path}")
    return atoms


def summarize_interactions(
    receptor_pdb: Path,
    ligand: Chem.Mol,
    *,
    contact_cutoff: float = 4.0,
    polar_cutoff: float = 3.5,
    max_residues: int = 12,
) -> InteractionEvidence:
    protein_atoms = read_protein_atoms(receptor_pdb)
    conformer = ligand.GetConformer()
    ligand_atoms = []
    for atom in ligand.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        point = conformer.GetAtomPosition(atom.GetIdx())
        ligand_atoms.append(
            (
                atom,
                np.asarray([point.x, point.y, point.z], dtype=float),
            )
        )

    residue_minimum: dict[str, float] = {}
    polar_residues: set[str] = set()
    hydrophobic_residues: set[str] = set()

    for protein_atom in protein_atoms:
        minimum = float("inf")
        for ligand_atom, ligand_coordinate in ligand_atoms:
            distance = float(np.linalg.norm(protein_atom.coordinate - ligand_coordinate))
            minimum = min(minimum, distance)
            protein_is_polar = protein_atom.element in {"N", "O", "S"}
            ligand_is_polar = ligand_atom.GetAtomicNum() in POLAR_ELEMENTS
            if protein_is_polar and ligand_is_polar and distance <= polar_cutoff:
                polar_residues.add(protein_atom.residue_label)
            if (
                protein_atom.residue_name in HYDROPHOBIC_RESIDUES
                and protein_atom.element in {"C", "S"}
                and ligand_atom.GetAtomicNum() in {6, 16, 17, 35, 53}
                and distance <= contact_cutoff
            ):
                hydrophobic_residues.add(protein_atom.residue_label)
        if minimum <= contact_cutoff:
            current = residue_minimum.get(protein_atom.residue_label, float("inf"))
            residue_minimum[protein_atom.residue_label] = min(current, minimum)

    ordered_contacts = sorted(residue_minimum, key=lambda key: residue_minimum[key])
    closest = min(residue_minimum.values()) if residue_minimum else None
    return InteractionEvidence(
        contact_residues=ordered_contacts[:max_residues],
        polar_contact_candidates=sorted(polar_residues)[:max_residues],
        hydrophobic_contact_residues=sorted(hydrophobic_residues)[:max_residues],
        closest_residue_distance=closest,
    )
