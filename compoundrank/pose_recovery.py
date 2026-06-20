"""Generic pose-recovery metrics for protein-ligand docking benchmarks.

This module is target-agnostic. It compares a known reference ligand pose
against a docked hypothesis pose in the same protein coordinate frame.

It intentionally does not hard-code any compound, virus, target, or PDB ID.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdFMCS


WATER_RESNAMES = {"HOH", "WAT", "DOD"}


@dataclass(frozen=True)
class LigandGroup:
    resname: str
    chain: str
    resseq: str
    icode: str
    altloc: str


@dataclass
class PoseRecoveryMetrics:
    reference_ligand: str
    docked_pose: str
    reference_atom_count: int
    docked_atom_count: int
    reference_center_x: float
    reference_center_y: float
    reference_center_z: float
    docked_center_x: float
    docked_center_y: float
    docked_center_z: float
    center_distance: float
    ordered_coordinate_rmsd: float | None
    symmetric_nearest_neighbor_rmsd: float
    openbabel_rmsd: float | None
    openbabel_minimized_rmsd: float | None
    interpretation: str
    limitations: list[str]


def parse_atom_line(line: str) -> dict[str, Any] | None:
    record = line[0:6].strip()
    if record not in {"ATOM", "HETATM"}:
        return None

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None

    element = line[76:78].strip() if len(line) >= 78 else ""
    atom_name = line[12:16].strip()
    if not element:
        element = "".join(ch for ch in atom_name if ch.isalpha())[:1].upper() or "?"

    return {
        "record": record,
        "atom_name": atom_name,
        "altloc": line[16].strip() or "-",
        "resname": line[17:20].strip(),
        "chain": line[21].strip() or "-",
        "resseq": line[22:26].strip(),
        "icode": line[26].strip() or "-",
        "x": x,
        "y": y,
        "z": z,
        "element": element.upper(),
        "raw": line.rstrip("\n"),
    }


def read_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for line in path.read_text(errors="replace").splitlines():
        atom = parse_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def group_key(atom: dict[str, Any]) -> LigandGroup:
    return LigandGroup(
        resname=str(atom["resname"]),
        chain=str(atom["chain"]),
        resseq=str(atom["resseq"]),
        icode=str(atom["icode"]),
        altloc=str(atom["altloc"]),
    )


def is_candidate_ligand_atom(atom: dict[str, Any]) -> bool:
    if atom["record"] != "HETATM":
        return False
    if str(atom["resname"]).upper() in WATER_RESNAMES:
        return False
    return True


def select_ligand_atoms(
    path: Path,
    *,
    resname: str | None = None,
    chain: str | None = None,
    resseq: str | None = None,
) -> list[dict[str, Any]]:
    atoms = [atom for atom in read_atoms(path) if is_candidate_ligand_atom(atom)]

    if resname is not None:
        atoms = [atom for atom in atoms if atom["resname"] == resname]
    if chain is not None:
        atoms = [atom for atom in atoms if atom["chain"] == chain]
    if resseq is not None:
        atoms = [atom for atom in atoms if atom["resseq"] == resseq]

    if not atoms:
        raise ValueError(f"No non-water HETATM ligand atoms found in {path}")

    grouped: dict[LigandGroup, list[dict[str, Any]]] = defaultdict(list)
    for atom in atoms:
        grouped[group_key(atom)].append(atom)

    if len(grouped) == 1:
        return atoms

    # If the file contains multiple non-water HETATM groups and the caller did
    # not fully disambiguate, choose the largest group as the most likely ligand.
    # The report will disclose this limitation.
    largest_key = max(grouped, key=lambda key: len(grouped[key]))
    return grouped[largest_key]


def heavy_atoms(atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [atom for atom in atoms if atom["element"] != "H"]


def center(atoms: list[dict[str, Any]]) -> tuple[float, float, float]:
    return (
        sum(float(atom["x"]) for atom in atoms) / len(atoms),
        sum(float(atom["y"]) for atom in atoms) / len(atoms),
        sum(float(atom["z"]) for atom in atoms) / len(atoms),
    )


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


def ordered_rmsd(
    reference_atoms: list[dict[str, Any]],
    docked_atoms: list[dict[str, Any]],
) -> float | None:
    if len(reference_atoms) != len(docked_atoms):
        return None

    squared = []
    for ref, dock in zip(reference_atoms, docked_atoms):
        squared.append(
            (float(ref["x"]) - float(dock["x"])) ** 2
            + (float(ref["y"]) - float(dock["y"])) ** 2
            + (float(ref["z"]) - float(dock["z"])) ** 2
        )

    return math.sqrt(sum(squared) / len(squared))


def one_way_nearest_neighbor_rmsd(
    source_atoms: list[dict[str, Any]],
    target_atoms: list[dict[str, Any]],
) -> float:
    squared_min_distances = []

    for src in source_atoms:
        src_xyz = (float(src["x"]), float(src["y"]), float(src["z"]))
        min_dist = min(
            distance(src_xyz, (float(tgt["x"]), float(tgt["y"]), float(tgt["z"])))
            for tgt in target_atoms
        )
        squared_min_distances.append(min_dist**2)

    return math.sqrt(sum(squared_min_distances) / len(squared_min_distances))


def symmetric_nearest_neighbor_rmsd(
    reference_atoms: list[dict[str, Any]],
    docked_atoms: list[dict[str, Any]],
) -> float:
    forward = one_way_nearest_neighbor_rmsd(reference_atoms, docked_atoms)
    reverse = one_way_nearest_neighbor_rmsd(docked_atoms, reference_atoms)
    return (forward + reverse) / 2.0


def parse_obrms_output(text: str) -> float | None:
    """Parse Open Babel obrms output.

    Expected example:
        RMSD reference.pdb:test.pdb 1.64412
    """
    for line in text.splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        try:
            return float(parts[-1])
        except ValueError:
            continue
    return None


def run_obrms(
    reference_ligand: Path,
    docked_pose: Path,
    *,
    minimize: bool = False,
    obrms_bin: str = "obrms",
) -> float | None:
    """Run Open Babel obrms if available.

    Returns None when obrms is unavailable or fails. This keeps the core
    benchmark portable while using a chemically mapped RMSD when possible.
    """
    if shutil.which(obrms_bin) is None:
        return None

    cmd = [obrms_bin]
    if minimize:
        cmd.append("-m")
    cmd.extend([str(reference_ligand), str(docked_pose)])

    try:
        result = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except Exception:
        return None

    combined = "\n".join(part for part in [result.stdout, result.stderr] if part)
    if result.returncode != 0 and not combined.strip():
        return None

    return parse_obrms_output(combined)


def interpret_pose_recovery(
    center_distance_value: float,
    nearest_neighbor_value: float,
    openbabel_value: float | None = None,
) -> str:
    """Interpret pose recovery using the strongest available metric.

    If chemically mapped Open Babel RMSD is available, use that as the primary
    shape/pose metric. Otherwise fall back to same-coordinate-frame geometric
    nearest-neighbor overlap.
    """
    pose_metric = openbabel_value if openbabel_value is not None else nearest_neighbor_value

    if center_distance_value <= 2.0 and pose_metric <= 2.5:
        return "strong_pose_recovery"
    if center_distance_value <= 4.0 and pose_metric <= 4.5:
        return "partial_pose_recovery"
    return "shifted_or_failed_pose_recovery"


def compare_pose_recovery(
    *,
    reference_ligand: Path,
    docked_pose: Path,
    reference_resname: str | None = None,
    reference_chain: str | None = None,
    reference_resseq: str | None = None,
    docked_resname: str | None = None,
    docked_chain: str | None = None,
    docked_resseq: str | None = None,
    use_openbabel: bool = True,
    obrms_bin: str = "obrms",
) -> PoseRecoveryMetrics:
    ref_atoms_all = select_ligand_atoms(
        reference_ligand,
        resname=reference_resname,
        chain=reference_chain,
        resseq=reference_resseq,
    )
    dock_atoms_all = select_ligand_atoms(
        docked_pose,
        resname=docked_resname,
        chain=docked_chain,
        resseq=docked_resseq,
    )

    ref_atoms = heavy_atoms(ref_atoms_all)
    dock_atoms = heavy_atoms(dock_atoms_all)

    ref_center = center(ref_atoms)
    dock_center = center(dock_atoms)
    center_distance_value = distance(ref_center, dock_center)

    ordered = ordered_rmsd(ref_atoms, dock_atoms)
    nearest_neighbor = symmetric_nearest_neighbor_rmsd(ref_atoms, dock_atoms)
    openbabel_rmsd = (
        run_obrms(reference_ligand, docked_pose, minimize=False, obrms_bin=obrms_bin)
        if use_openbabel
        else None
    )
    openbabel_minimized_rmsd = (
        run_obrms(reference_ligand, docked_pose, minimize=True, obrms_bin=obrms_bin)
        if use_openbabel
        else None
    )


    limitations = [
        "This comparison assumes reference and docked pose are already in the same protein coordinate frame.",
        "Nearest-neighbor RMSD is a geometric overlap metric, not a chemically atom-mapped RMSD.",
        "Open Babel RMSD is included when obrms is available and should be preferred over ordered-coordinate RMSD.",
    ]
    if ordered is None:
        limitations.append(
            "Ordered coordinate RMSD was not calculated because reference and docked heavy-atom counts differ."
        )
    else:
        limitations.append(
            "Ordered coordinate RMSD assumes atom ordering is comparable between the two ligand files."
        )

    return PoseRecoveryMetrics(
        reference_ligand=str(reference_ligand),
        docked_pose=str(docked_pose),
        reference_atom_count=len(ref_atoms),
        docked_atom_count=len(dock_atoms),
        reference_center_x=ref_center[0],
        reference_center_y=ref_center[1],
        reference_center_z=ref_center[2],
        docked_center_x=dock_center[0],
        docked_center_y=dock_center[1],
        docked_center_z=dock_center[2],
        center_distance=center_distance_value,
        ordered_coordinate_rmsd=ordered,
        symmetric_nearest_neighbor_rmsd=nearest_neighbor,
        openbabel_rmsd=openbabel_rmsd,
        openbabel_minimized_rmsd=openbabel_minimized_rmsd,
        interpretation=interpret_pose_recovery(center_distance_value, nearest_neighbor, openbabel_minimized_rmsd or openbabel_rmsd),
        limitations=limitations,
    )


def write_outputs(metrics: PoseRecoveryMetrics, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    json_path = output_dir / "pose_recovery_metrics.json"
    csv_path = output_dir / "pose_recovery_metrics.csv"
    report_path = output_dir / "pose_recovery_report.md"

    data = asdict(metrics)
    json_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data.keys()))
        writer.writeheader()
        writer.writerow(data)

    ordered = (
        f"{metrics.ordered_coordinate_rmsd:.3f}"
        if metrics.ordered_coordinate_rmsd is not None
        else "not calculated"
    )
    ob_rmsd = (
        f"{metrics.openbabel_rmsd:.3f}"
        if metrics.openbabel_rmsd is not None
        else "not calculated"
    )
    ob_min = (
        f"{metrics.openbabel_minimized_rmsd:.3f}"
        if metrics.openbabel_minimized_rmsd is not None
        else "not calculated"
    )

    report = f"""# Pose Recovery Report

## Inputs

- Reference ligand: `{metrics.reference_ligand}`
- Docked pose: `{metrics.docked_pose}`

## Metrics

| Metric | Value |
|---|---:|
| Reference heavy atoms | {metrics.reference_atom_count} |
| Docked heavy atoms | {metrics.docked_atom_count} |
| Center distance | {metrics.center_distance:.3f} Å |
| Ordered coordinate RMSD | {ordered} Å |
| Symmetric nearest-neighbor RMSD | {metrics.symmetric_nearest_neighbor_rmsd:.3f} Å |
| Open Babel RMSD | {ob_rmsd} Å |
| Open Babel minimized/symmetry RMSD | {ob_min} Å |

## Interpretation

`{metrics.interpretation}`

## Limitations

"""
    for item in metrics.limitations:
        report += f"- {item}\n"

    report_path.write_text(report, encoding="utf-8")

    return {
        "pose_recovery_metrics_json": json_path,
        "pose_recovery_metrics_csv": csv_path,
        "pose_recovery_report": report_path,
    }


def _load_first_sdf_molecule(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
    )

    molecule = next(
        (
            mol
            for mol in supplier
            if mol is not None
        ),
        None,
    )

    if molecule is None:
        raise ValueError(
            f"Could not read an SDF molecule from {path}"
        )

    return molecule


def _heavy_rdkit_molecule(
    molecule: Chem.Mol,
) -> Chem.Mol:
    result = Chem.RemoveHs(
        Chem.Mol(molecule),
        sanitize=True,
    )
    Chem.SanitizeMol(result)
    return result


def _molecule_property_float(
    molecule: Chem.Mol,
    names: tuple[str, ...],
) -> float | None:
    available = {
        name.casefold(): name
        for name in molecule.GetPropNames()
    }

    for requested in names:
        actual = available.get(
            requested.casefold()
        )

        if actual is None:
            continue

        try:
            return float(
                molecule.GetProp(actual)
            )
        except (TypeError, ValueError):
            continue

    return None


def _coordinate_rmsd_from_map(
    probe: Chem.Mol,
    reference: Chem.Mol,
    atom_map: list[tuple[int, int]],
) -> float:
    """Calculate mapped RMSD without translating or rotating coordinates."""
    if not atom_map:
        raise ValueError(
            "An empty atom map cannot be used for RMSD."
        )

    probe_conf = probe.GetConformer()
    reference_conf = reference.GetConformer()

    squared_distance = 0.0

    for probe_index, reference_index in atom_map:
        probe_position = probe_conf.GetAtomPosition(
            probe_index
        )
        reference_position = (
            reference_conf.GetAtomPosition(
                reference_index
            )
        )

        dx = probe_position.x - reference_position.x
        dy = probe_position.y - reference_position.y
        dz = probe_position.z - reference_position.z

        squared_distance += (
            dx * dx
            + dy * dy
            + dz * dz
        )

    return math.sqrt(
        squared_distance / len(atom_map)
    )


def symmetry_aware_nofit_rmsd(
    probe: Chem.Mol,
    reference: Chem.Mol,
) -> tuple[float, int]:
    """Calculate complete, element-mapped RMSD in the receptor frame.

    Bond order is deliberately ignored during mapping because equivalent
    ChEMBL, RCSB, RDKit, Open Babel, and GNINA files may encode aromaticity
    or tautomeric bond orders differently. Atom elements, connectivity,
    ring membership, complete atom coverage, and coordinates are retained.
    """
    probe = _heavy_rdkit_molecule(
        probe
    )
    reference = _heavy_rdkit_molecule(
        reference
    )

    if (
        probe.GetNumAtoms()
        != reference.GetNumAtoms()
    ):
        raise ValueError(
            "Heavy-atom counts differ: "
            f"probe={probe.GetNumAtoms()}, "
            f"reference={reference.GetNumAtoms()}"
        )

    if (
        probe.GetNumBonds()
        != reference.GetNumBonds()
    ):
        raise ValueError(
            "Heavy-atom bond counts differ: "
            f"probe={probe.GetNumBonds()}, "
            f"reference={reference.GetNumBonds()}"
        )

    mcs = rdFMCS.FindMCS(
        [probe, reference],
        atomCompare=(
            rdFMCS.AtomCompare.CompareElements
        ),
        bondCompare=(
            rdFMCS.BondCompare.CompareAny
        ),
        matchValences=False,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        matchChiralTag=False,
        timeout=30,
    )

    if mcs.canceled:
        raise RuntimeError(
            "Complete atom mapping timed out."
        )

    if (
        mcs.numAtoms != probe.GetNumAtoms()
        or mcs.numBonds != probe.GetNumBonds()
    ):
        raise ValueError(
            "Could not establish a complete "
            "element-and-connectivity mapping: "
            f"MCS atoms={mcs.numAtoms}/"
            f"{probe.GetNumAtoms()}, "
            f"MCS bonds={mcs.numBonds}/"
            f"{probe.GetNumBonds()}"
        )

    query = Chem.MolFromSmarts(
        mcs.smartsString
    )

    if query is None:
        raise RuntimeError(
            "Could not construct the complete MCS query."
        )

    probe_matches = probe.GetSubstructMatches(
        query,
        uniquify=True,
        useChirality=False,
        maxMatches=512,
    )

    reference_matches = (
        reference.GetSubstructMatches(
            query,
            uniquify=True,
            useChirality=False,
            maxMatches=512,
        )
    )

    if (
        not probe_matches
        or not reference_matches
    ):
        raise RuntimeError(
            "Complete MCS mapping produced no atom matches."
        )

    best_rmsd = math.inf
    mapping_count = 0

    for probe_match in probe_matches:
        for reference_match in reference_matches:
            atom_map = list(
                zip(
                    probe_match,
                    reference_match,
                )
            )

            rmsd = _coordinate_rmsd_from_map(
                probe,
                reference,
                atom_map,
            )

            mapping_count += 1

            if rmsd < best_rmsd:
                best_rmsd = rmsd

    return best_rmsd, mapping_count


def evaluate_scored_pose_sdf(
    *,
    reference_ligand: Path,
    poses_sdf: Path,
    rmsd_threshold: float = 2.0,
) -> dict[str, Any]:
    """Evaluate all scored GNINA poses against a cognate SDF pose."""
    reference = _load_first_sdf_molecule(
        reference_ligand
    )

    supplier = Chem.SDMolSupplier(
        str(poses_sdf),
        removeHs=False,
    )

    pose_records: list[dict[str, Any]] = []
    mapping_failures: list[dict[str, Any]] = []

    for pose_index, molecule in enumerate(
        supplier,
        start=1,
    ):
        if molecule is None:
            mapping_failures.append(
                {
                    "pose_index": pose_index,
                    "error": "RDKit could not read the pose.",
                }
            )
            continue

        cnnscore = _molecule_property_float(
            molecule,
            (
                "CNNscore",
                "CNN_score",
            ),
        )

        try:
            rmsd, mapping_count = (
                symmetry_aware_nofit_rmsd(
                    molecule,
                    reference,
                )
            )
        except Exception as exc:
            mapping_failures.append(
                {
                    "pose_index": pose_index,
                    "error": str(exc),
                }
            )
            continue

        pose_records.append(
            {
                "pose_index": pose_index,
                "cnnscore": cnnscore,
                "cnnaffinity": (
                    _molecule_property_float(
                        molecule,
                        (
                            "CNNaffinity",
                            "CNN_affinity",
                        ),
                    )
                ),
                "affinity": (
                    _molecule_property_float(
                        molecule,
                        (
                            "minimizedAffinity",
                            "minimized_affinity",
                            "affinity",
                        ),
                    )
                ),
                "heavy_atom_rmsd": rmsd,
                "mapping_count": mapping_count,
            }
        )

    if not pose_records:
        raise RuntimeError(
            "No poses could be chemically mapped "
            "to the reference ligand."
        )

    scored_records = [
        record
        for record in pose_records
        if record["cnnscore"] is not None
    ]

    if not scored_records:
        raise RuntimeError(
            "No chemically mapped poses contained "
            "a GNINA CNNscore."
        )

    by_score = sorted(
        scored_records,
        key=lambda record: (
            -float(record["cnnscore"]),
            float(record["heavy_atom_rmsd"]),
        ),
    )

    by_rmsd = sorted(
        scored_records,
        key=lambda record: (
            float(record["heavy_atom_rmsd"]),
            -float(record["cnnscore"]),
        ),
    )

    top_cnn_pose = by_score[0]
    best_sampled_pose = by_rmsd[0]

    sampling_pass = (
        float(
            best_sampled_pose[
                "heavy_atom_rmsd"
            ]
        )
        <= rmsd_threshold
    )

    ranking_pass = (
        float(
            top_cnn_pose[
                "heavy_atom_rmsd"
            ]
        )
        <= rmsd_threshold
    )

    if sampling_pass and ranking_pass:
        overall = (
            "cognate_pose_recovery_and_ranking_pass"
        )
    elif sampling_pass:
        overall = (
            "sampling_pass_ranking_failure"
        )
    else:
        overall = (
            "pose_recovery_failure"
        )

    return {
        "reference_ligand": str(
            reference_ligand
        ),
        "poses_sdf": str(
            poses_sdf
        ),
        "rmsd_method": (
            "symmetry-aware complete heavy-atom "
            "coordinate RMSD without translation "
            "or rotation"
        ),
        "bond_order_mapping": (
            "bond order ignored; atom elements, "
            "connectivity, ring membership, and "
            "complete atom coverage required"
        ),
        "rmsd_threshold_angstrom": (
            rmsd_threshold
        ),
        "mapped_pose_count": len(
            pose_records
        ),
        "mapping_failure_count": len(
            mapping_failures
        ),
        "mapping_failures": (
            mapping_failures
        ),
        "top_cnn_pose": top_cnn_pose,
        "best_sampled_pose": (
            best_sampled_pose
        ),
        "sampling_pass": sampling_pass,
        "ranking_pass": ranking_pass,
        "overall": overall,
        "poses_by_cnnscore": by_score,
    }


def write_scored_pose_outputs(
    summary: dict[str, Any],
    output_dir: Path,
) -> dict[str, Path]:
    """Write batch pose-recovery JSON, CSV, and Markdown outputs."""
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    json_path = (
        output_dir
        / "pose_set_recovery_summary.json"
    )
    csv_path = (
        output_dir
        / "pose_set_recovery_metrics.csv"
    )
    report_path = (
        output_dir
        / "pose_set_recovery_report.md"
    )

    json_path.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    records = (
        summary.get("poses_by_cnnscore")
        or []
    )

    if records:
        with csv_path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(
                    records[0].keys()
                ),
            )
            writer.writeheader()
            writer.writerows(records)
    else:
        csv_path.write_text(
            "",
            encoding="utf-8",
        )

    top_pose = summary["top_cnn_pose"]
    best_pose = summary[
        "best_sampled_pose"
    ]

    report = f"""# Scored Pose-Recovery Report

## Inputs

- Reference ligand: `{summary['reference_ligand']}`
- GNINA pose file: `{summary['poses_sdf']}`
- RMSD method: {summary['rmsd_method']}
- Pass threshold: {summary['rmsd_threshold_angstrom']:.3f} Å

## Benchmark Result

| Metric | Value |
|---|---:|
| Chemically mapped poses | {summary['mapped_pose_count']} |
| Mapping failures | {summary['mapping_failure_count']} |
| Top CNN pose index | {top_pose['pose_index']} |
| Top CNN score | {top_pose['cnnscore']:.6f} |
| Top CNN pose RMSD | {top_pose['heavy_atom_rmsd']:.3f} Å |
| Best sampled pose index | {best_pose['pose_index']} |
| Best sampled RMSD | {best_pose['heavy_atom_rmsd']:.3f} Å |
| Sampling pass | {'yes' if summary['sampling_pass'] else 'no'} |
| Ranking pass | {'yes' if summary['ranking_pass'] else 'no'} |

## Overall

`{summary['overall']}`

## Interpretation Limits

- RMSD is calculated in the existing receptor coordinate frame without translating or rotating the ligand.
- Complete element-and-connectivity mapping is required.
- Bond order is ignored only to tolerate equivalent aromaticity or tautomer encodings across molecular file sources.
- A cognate redocking pass validates pose recovery under benchmark conditions; it does not prove biological inhibition.
"""

    report_path.write_text(
        report,
        encoding="utf-8",
    )

    return {
        "pose_set_recovery_summary": (
            json_path
        ),
        "pose_set_recovery_metrics": (
            csv_path
        ),
        "pose_set_recovery_report": (
            report_path
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare docked ligand poses against "
            "a reference ligand pose."
        )
    )

    parser.add_argument(
        "--reference-ligand",
        required=True,
        type=Path,
    )

    pose_group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    pose_group.add_argument(
        "--docked-pose",
        type=Path,
        help=(
            "Single PDB pose for the legacy "
            "coordinate-overlap workflow."
        ),
    )

    pose_group.add_argument(
        "--poses-sdf",
        type=Path,
        help=(
            "Scored GNINA multi-pose SDF for "
            "chemically mapped cognate redocking."
        ),
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--rmsd-threshold",
        type=float,
        default=2.0,
    )

    parser.add_argument(
        "--reference-resname"
    )
    parser.add_argument(
        "--reference-chain"
    )
    parser.add_argument(
        "--reference-resseq"
    )

    parser.add_argument(
        "--docked-resname"
    )
    parser.add_argument(
        "--docked-chain"
    )
    parser.add_argument(
        "--docked-resseq"
    )

    parser.add_argument(
        "--no-openbabel",
        action="store_true",
        help=(
            "Disable optional Open Babel obrms "
            "for the legacy single-PDB workflow."
        ),
    )

    parser.add_argument(
        "--obrms-bin",
        default="obrms",
    )

    return parser


def main(
    argv: list[str] | None = None,
) -> int:
    args = build_parser().parse_args(
        argv
    )

    if args.poses_sdf is not None:
        summary = evaluate_scored_pose_sdf(
            reference_ligand=(
                args.reference_ligand
            ),
            poses_sdf=args.poses_sdf,
            rmsd_threshold=(
                args.rmsd_threshold
            ),
        )

        outputs = write_scored_pose_outputs(
            summary,
            args.output_dir,
        )

        top_pose = summary[
            "top_cnn_pose"
        ]
        best_pose = summary[
            "best_sampled_pose"
        ]

        print(
            "[POSE_RECOVERY] "
            f"top_cnn_rmsd="
            f"{top_pose['heavy_atom_rmsd']:.3f} Å"
        )
        print(
            "[POSE_RECOVERY] "
            f"best_sampled_rmsd="
            f"{best_pose['heavy_atom_rmsd']:.3f} Å"
        )
        print(
            "[POSE_RECOVERY] "
            f"sampling_pass="
            f"{summary['sampling_pass']}"
        )
        print(
            "[POSE_RECOVERY] "
            f"ranking_pass="
            f"{summary['ranking_pass']}"
        )
        print(
            "[POSE_RECOVERY] "
            f"overall={summary['overall']}"
        )

        for label, output_path in (
            outputs.items()
        ):
            print(
                f"[POSE_RECOVERY] "
                f"{label}: {output_path}"
            )

        return 0

    metrics = compare_pose_recovery(
        reference_ligand=(
            args.reference_ligand
        ),
        docked_pose=args.docked_pose,
        reference_resname=(
            args.reference_resname
        ),
        reference_chain=(
            args.reference_chain
        ),
        reference_resseq=(
            args.reference_resseq
        ),
        docked_resname=(
            args.docked_resname
        ),
        docked_chain=(
            args.docked_chain
        ),
        docked_resseq=(
            args.docked_resseq
        ),
        use_openbabel=(
            not args.no_openbabel
        ),
        obrms_bin=args.obrms_bin,
    )

    outputs = write_outputs(
        metrics,
        args.output_dir,
    )

    print(
        "[POSE_RECOVERY] Metrics:"
    )
    print(
        "[POSE_RECOVERY] "
        f"center_distance="
        f"{metrics.center_distance:.3f} Å"
    )
    print(
        "[POSE_RECOVERY] "
        f"symmetric_nearest_neighbor_rmsd="
        f"{metrics.symmetric_nearest_neighbor_rmsd:.3f} Å"
    )

    if metrics.openbabel_rmsd is not None:
        print(
            "[POSE_RECOVERY] "
            f"openbabel_rmsd="
            f"{metrics.openbabel_rmsd:.3f} Å"
        )

    if (
        metrics.openbabel_minimized_rmsd
        is not None
    ):
        print(
            "[POSE_RECOVERY] "
            f"openbabel_minimized_rmsd="
            f"{metrics.openbabel_minimized_rmsd:.3f} Å"
        )

    print(
        "[POSE_RECOVERY] "
        f"interpretation="
        f"{metrics.interpretation}"
    )

    print(
        "[POSE_RECOVERY] Outputs:"
    )

    for label, output_path in (
        outputs.items()
    ):
        print(
            f"[POSE_RECOVERY] "
            f"{label}: {output_path}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

