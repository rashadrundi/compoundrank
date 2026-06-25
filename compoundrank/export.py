from __future__ import annotations

from pathlib import Path
from textwrap import wrap

from rdkit import Chem

from .models import InteractionEvidence, LigandResult, PoseCluster


def _remark_lines(text: str) -> list[str]:
    lines: list[str] = []
    for chunk in wrap(text, width=68, break_long_words=False, break_on_hyphens=False):
        lines.append(f"REMARK 900 {chunk}")
    return lines or ["REMARK 900"]


def _format_optional(value: float | None, digits: int = 4) -> str:
    return "N/A" if value is None else f"{value:.{digits}f}"


def _ligand_pdb_records(molecule: Chem.Mol) -> list[str]:
    block = Chem.MolToPDBBlock(molecule)
    lines: list[str] = []

    for raw_line in block.splitlines():
        if raw_line.startswith(("ATOM", "HETATM")):
            line = "HETATM" + raw_line[6:]
            if len(line) < 80:
                line = line.ljust(80)
            line = line[:17] + "LIG" + line[20:21] + "Z" + f"{1:4d}" + line[26:]
            lines.append(line.rstrip())

        elif raw_line.startswith("CONECT"):
            lines.append(raw_line)

    return lines


def write_complex_pdb(
    output_path: Path,
    receptor_pdb: Path,
    ligand_result: LigandResult,
    cluster: PoseCluster,
    interactions: InteractionEvidence,
    *,
    compound_priority_rank: int,
    hypothesis_rank: int,
    hypothesis_count: int,
) -> None:
    representative = cluster.representative
    remarks: list[str] = []

    metadata = [
        "COMPOUNDRANK COMPUTATIONAL BINDING HYPOTHESIS",
        f"COMPOUND {ligand_result.ligand.name}",
        f"COMPOUND PRIORITY RANK {compound_priority_rank}",
        f"HYPOTHESIS {hypothesis_rank} OF {hypothesis_count}",
        "ORDERING SOURCE GNINA CNN SCORE ONLY",
        f"GNINA CNN SCORE {_format_optional(representative.cnn_score)}",
        f"GNINA CNN AFFINITY {_format_optional(representative.cnn_affinity)}",
        f"GNINA MINIMIZED AFFINITY {_format_optional(representative.minimized_affinity)}",
        f"CLUSTER MEMBERS {cluster.member_count}",
        f"SEEDS REPRESENTED {len(cluster.seeds)}",
        f"POSE CONFIDENCE {ligand_result.uncertainty.upper()}",
        f"POCKET ID {representative.pocket_id}",
        f"POCKET SOURCE {representative.pocket_source or 'unknown'}",
        (
            "RECEPTOR CONFORMER "
            f"{representative.receptor_conformer_id}"
        ),
        "PHYSICAL VALIDITY POSEBUSTERS PASS",
    ]

    if representative.fpocket_score is not None:
        metadata.append(f"FPOCKET SCORE {representative.fpocket_score:.4f}")

    if interactions.closest_residue_distance is not None:
        metadata.append(
            f"CLOSEST PROTEIN CONTACT {interactions.closest_residue_distance:.2f} ANGSTROM"
        )

    if interactions.contact_residues:
        metadata.append("CONTACT RESIDUES " + ", ".join(interactions.contact_residues))

    if interactions.polar_contact_candidates:
        metadata.append(
            "POLAR CONTACT CANDIDATES "
            + ", ".join(interactions.polar_contact_candidates)
        )

    if interactions.hydrophobic_contact_residues:
        metadata.append(
            "HYDROPHOBIC CONTACT RESIDUES "
            + ", ".join(interactions.hydrophobic_contact_residues)
        )

    for reason in ligand_result.uncertainty_reasons:
        metadata.append("UNCERTAINTY NOTE " + reason)

    metadata.append(
        "INTERACTIONS ARE DISTANCE-BASED EVIDENCE, NOT PROOF OF BINDING OR EFFICACY"
    )

    for item in metadata:
        remarks.extend(_remark_lines(item))

    receptor_lines: list[str] = []
    for line in receptor_pdb.read_text(errors="replace").splitlines():
        if line.startswith(("END", "CONECT")):
            continue
        receptor_lines.append(line)

    ligand_lines = _ligand_pdb_records(representative.molecule)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "\n".join([*remarks, *receptor_lines, "TER", *ligand_lines, "END"]) + "\n",
        encoding="utf-8",
    )
