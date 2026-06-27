"""Prepare retrieved ligand structures for docking.

Input:
    Stage 4A docking_manifest.csv containing file-based SDF ligands.

Output:
    prepared_ligands/*.sdf
    prepared_docking_manifest.csv
    ligand_preparation_report.json
    ligand_preparation_report.csv

Preparation performed:
    - RDKit parsing and sanitization
    - principal-fragment selection
    - preservation of the existing formal-charge state
    - explicit hydrogen addition
    - ETKDGv3 3D conformer generation
    - MMFF94s optimization, with UFF fallback
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


DEFAULT_RANDOM_SEED = 0xF00D


class LigandPreparationError(RuntimeError):
    """Raised when a ligand cannot be prepared safely."""


def _safe_filename(value: str) -> str:
    text = str(value or "").strip().lower()

    safe = "".join(
        character if character.isalnum() else "_"
        for character in text
    )

    safe = "_".join(
        part
        for part in safe.split("_")
        if part
    )

    return safe or "unnamed_ligand"


def _read_first_sdf_molecule(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
        sanitize=True,
    )

    for molecule in supplier:
        if molecule is not None:
            return molecule

    raise LigandPreparationError(
        f"RDKit could not read a valid molecule from {path}"
    )


def _fragment_rank(molecule: Chem.Mol) -> tuple[int, int, int, float]:
    """Prefer an organic, larger principal component."""
    carbon_count = sum(
        atom.GetAtomicNum() == 6
        for atom in molecule.GetAtoms()
    )
    heavy_atom_count = molecule.GetNumHeavyAtoms()

    return (
        1 if carbon_count > 0 else 0,
        heavy_atom_count,
        carbon_count,
        float(Descriptors.MolWt(molecule)),
    )


def select_principal_fragment(
    molecule: Chem.Mol,
) -> tuple[Chem.Mol, dict[str, Any]]:
    """Select the largest organic fragment without neutralizing it."""
    fragments = list(
        Chem.GetMolFrags(
            molecule,
            asMols=True,
            sanitizeFrags=True,
        )
    )

    if not fragments:
        raise LigandPreparationError(
            "No molecular fragments were available."
        )

    ranked = sorted(
        fragments,
        key=_fragment_rank,
        reverse=True,
    )

    selected = Chem.Mol(ranked[0])
    Chem.SanitizeMol(selected)

    metadata = {
        "input_fragment_count": len(fragments),
        "selected_fragment_heavy_atoms": (
            selected.GetNumHeavyAtoms()
        ),
        "selected_fragment_carbon_atoms": sum(
            atom.GetAtomicNum() == 6
            for atom in selected.GetAtoms()
        ),
        "selected_fragment_formal_charge": (
            Chem.GetFormalCharge(selected)
        ),
        "selected_fragment_molecular_weight": round(
            float(Descriptors.MolWt(selected)),
            6,
        ),
    }

    return selected, metadata


def _embed_3d(
    molecule: Chem.Mol,
    *,
    random_seed: int,
) -> tuple[Chem.Mol, str]:
    """Add hydrogens and generate one deterministic 3D conformer."""
    molecule = Chem.Mol(molecule)

    # Remove the original 2D conformer before embedding.
    molecule.RemoveAllConformers()

    # Existing formal charges are retained. This is not a pH-aware
    # protonation model.
    molecule = Chem.AddHs(
        molecule,
        addCoords=False,
    )

    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = int(random_seed)
    parameters.enforceChirality = True
    parameters.useSmallRingTorsions = True
    parameters.useMacrocycleTorsions = True

    result = AllChem.EmbedMolecule(
        molecule,
        parameters,
    )

    if result < 0:
        retry = AllChem.ETKDGv3()
        retry.randomSeed = int(random_seed)
        retry.enforceChirality = True
        retry.useRandomCoords = True
        retry.useSmallRingTorsions = True
        retry.useMacrocycleTorsions = True

        result = AllChem.EmbedMolecule(
            molecule,
            retry,
        )

        if result < 0:
            raise LigandPreparationError(
                "ETKDGv3 conformer generation failed."
            )

        embedding_route = "ETKDGv3_random_coordinates_retry"
    else:
        embedding_route = "ETKDGv3"

    conformer = molecule.GetConformer()
    conformer.Set3D(True)

    return molecule, embedding_route


def _optimize_geometry(
    molecule: Chem.Mol,
    *,
    max_iterations: int,
) -> dict[str, Any]:
    """Optimize with MMFF94s, falling back to UFF."""
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        result = AllChem.MMFFOptimizeMolecule(
            molecule,
            mmffVariant="MMFF94s",
            maxIters=max_iterations,
        )

        return {
            "force_field": "MMFF94s",
            "optimization_result": int(result),
            "optimization_converged": result == 0,
        }

    if AllChem.UFFHasAllMoleculeParams(molecule):
        result = AllChem.UFFOptimizeMolecule(
            molecule,
            maxIters=max_iterations,
        )

        return {
            "force_field": "UFF",
            "optimization_result": int(result),
            "optimization_converged": result == 0,
        }

    return {
        "force_field": None,
        "optimization_result": None,
        "optimization_converged": False,
        "optimization_warning": (
            "Neither MMFF94s nor UFF had parameters for "
            "the complete molecule."
        ),
    }


def _coordinate_metrics(
    molecule: Chem.Mol,
) -> dict[str, Any]:
    if molecule.GetNumConformers() != 1:
        raise LigandPreparationError(
            "Prepared ligand must contain exactly one conformer."
        )

    conformer = molecule.GetConformer()

    coordinates = [
        conformer.GetAtomPosition(index)
        for index in range(molecule.GetNumAtoms())
    ]

    if not coordinates:
        raise LigandPreparationError(
            "Prepared ligand contains no coordinates."
        )

    x_span = max(point.x for point in coordinates) - min(
        point.x for point in coordinates
    )
    y_span = max(point.y for point in coordinates) - min(
        point.y for point in coordinates
    )
    z_span = max(point.z for point in coordinates) - min(
        point.z for point in coordinates
    )

    if max(x_span, y_span, z_span) <= 0.01:
        raise LigandPreparationError(
            "Prepared conformer has effectively no coordinate span."
        )

    if z_span <= 0.01:
        raise LigandPreparationError(
            "Prepared conformer remains effectively planar in z."
        )

    return {
        "marked_as_3d": bool(conformer.Is3D()),
        "x_span": round(float(x_span), 6),
        "y_span": round(float(y_span), 6),
        "z_span": round(float(z_span), 6),
    }


def prepare_ligand(
    *,
    name: str,
    input_path: Path,
    output_dir: Path,
    random_seed: int = DEFAULT_RANDOM_SEED,
    max_iterations: int = 500,
) -> dict[str, Any]:
    """Prepare one SDF ligand and return its report entry."""
    input_path = input_path.resolve()
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    source = _read_first_sdf_molecule(input_path)

    source_name = (
        source.GetProp("_Name")
        if source.HasProp("_Name")
        else name
    )

    selected, fragment_metadata = (
        select_principal_fragment(source)
    )

    Chem.AssignStereochemistry(
        selected,
        cleanIt=True,
        force=True,
    )

    prepared, embedding_route = _embed_3d(
        selected,
        random_seed=random_seed,
    )

    optimization = _optimize_geometry(
        prepared,
        max_iterations=max_iterations,
    )

    Chem.SanitizeMol(prepared)

    coordinate_metadata = _coordinate_metrics(
        prepared
    )

    explicit_hydrogens = sum(
        atom.GetAtomicNum() == 1
        for atom in prepared.GetAtoms()
    )

    if explicit_hydrogens == 0:
        raise LigandPreparationError(
            "No explicit hydrogens were present after preparation."
        )

    output_path = (
        output_dir
        / f"{_safe_filename(name)}_prepared.sdf"
    )

    prepared.SetProp(
        "_Name",
        str(source_name or name),
    )
    prepared.SetProp(
        "COMPOUNDRANK_INPUT_NAME",
        str(name),
    )
    prepared.SetProp(
        "COMPOUNDRANK_SOURCE_SDF",
        str(input_path),
    )
    prepared.SetProp(
        "COMPOUNDRANK_PREPARATION",
        "fragment_selection;AddHs;ETKDGv3;MMFF94s_or_UFF",
    )
    prepared.SetProp(
        "COMPOUNDRANK_PROTONATION_NOTE",
        (
            "Explicit hydrogens added from stored formal charge; "
            "not pH-aware."
        ),
    )

    writer = Chem.SDWriter(str(output_path))

    try:
        writer.write(prepared)
    finally:
        writer.close()

    if not output_path.exists():
        raise LigandPreparationError(
            f"Prepared SDF was not created: {output_path}"
        )

    return {
        "name": name,
        "status": "prepared",
        "input_path": str(input_path),
        "output_path": str(output_path.resolve()),
        "source_name": source_name,
        "input_atoms": source.GetNumAtoms(),
        "prepared_atoms": prepared.GetNumAtoms(),
        "prepared_heavy_atoms": prepared.GetNumHeavyAtoms(),
        "explicit_hydrogens": explicit_hydrogens,
        "formal_charge": Chem.GetFormalCharge(prepared),
        "embedding_route": embedding_route,
        "random_seed": int(random_seed),
        "protonation_method": (
            "rdkit_explicit_hydrogens_current_formal_charge"
        ),
        "pH_aware_protonation": False,
        **fragment_metadata,
        **optimization,
        **coordinate_metadata,
    }


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def _write_manifest(
    path: Path,
    rows: list[dict[str, str]],
) -> None:
    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "name",
                "source_type",
                "value",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_report_csv(
    path: Path,
    entries: list[dict[str, Any]],
) -> None:
    fieldnames = [
        "name",
        "status",
        "input_path",
        "output_path",
        "error",
        "input_fragment_count",
        "selected_fragment_heavy_atoms",
        "formal_charge",
        "explicit_hydrogens",
        "embedding_route",
        "force_field",
        "optimization_result",
        "optimization_converged",
        "marked_as_3d",
        "x_span",
        "y_span",
        "z_span",
        "pH_aware_protonation",
    ]

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()

        for entry in entries:
            writer.writerow(entry)


def prepare_ligand_manifest(
    *,
    input_manifest: Path,
    output_dir: Path,
    random_seed: int = DEFAULT_RANDOM_SEED,
    max_iterations: int = 500,
) -> dict[str, Path]:
    """Prepare all file-based ligands from a Stage 4A manifest."""
    input_manifest = input_manifest.resolve()
    output_dir = output_dir.resolve()

    prepared_dir = output_dir / "prepared_ligands"
    prepared_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    entries: list[dict[str, Any]] = []
    prepared_manifest_rows: list[dict[str, str]] = []

    for row in _read_manifest(input_manifest):
        name = str(row.get("name") or "").strip()
        source_type = str(
            row.get("source_type") or ""
        ).strip()
        value = str(row.get("value") or "").strip()

        if source_type != "file":
            entries.append(
                {
                    "name": name,
                    "status": "skipped",
                    "input_path": value,
                    "error": (
                        f"Unsupported source_type: {source_type!r}"
                    ),
                }
            )
            continue

        input_path = Path(value)

        if not input_path.exists():
            entries.append(
                {
                    "name": name,
                    "status": "failed",
                    "input_path": value,
                    "error": "Input ligand file does not exist.",
                }
            )
            continue

        try:
            entry = prepare_ligand(
                name=name,
                input_path=input_path,
                output_dir=prepared_dir,
                random_seed=random_seed,
                max_iterations=max_iterations,
            )

            entries.append(entry)

            prepared_manifest_rows.append(
                {
                    "name": name,
                    "source_type": "file",
                    "value": entry["output_path"],
                }
            )

        except Exception as exc:  # noqa: BLE001
            entries.append(
                {
                    "name": name,
                    "status": "failed",
                    "input_path": str(input_path),
                    "error": str(exc),
                }
            )

    prepared_count = sum(
        entry.get("status") == "prepared"
        for entry in entries
    )
    failed_count = sum(
        entry.get("status") == "failed"
        for entry in entries
    )
    skipped_count = sum(
        entry.get("status") == "skipped"
        for entry in entries
    )

    prepared_manifest = (
        output_dir / "prepared_docking_manifest.csv"
    )
    report_json = (
        output_dir / "ligand_preparation_report.json"
    )
    report_csv = (
        output_dir / "ligand_preparation_report.csv"
    )

    _write_manifest(
        prepared_manifest,
        prepared_manifest_rows,
    )

    report_payload = {
        "input_manifest": str(input_manifest),
        "prepared_manifest": str(prepared_manifest),
        "prepared_ligand_directory": str(prepared_dir),
        "preparation_configuration": {
            "random_seed": int(random_seed),
            "max_iterations": int(max_iterations),
            "embedding_method": "ETKDGv3",
            "preferred_force_field": "MMFF94s",
            "fallback_force_field": "UFF",
            "principal_fragment_selection": True,
            "explicit_hydrogens": True,
            "pH_aware_protonation": False,
        },
        "summary": {
            "input_count": len(entries),
            "prepared_count": prepared_count,
            "failed_count": failed_count,
            "skipped_count": skipped_count,
        },
        "ligands": entries,
    }

    report_json.write_text(
        json.dumps(
            report_payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    _write_report_csv(
        report_csv,
        entries,
    )

    if prepared_count == 0:
        raise LigandPreparationError(
            "No ligands were successfully prepared. "
            f"Inspect {report_json}"
        )

    return {
        "prepared_ligands": prepared_dir,
        "prepared_docking_manifest": prepared_manifest,
        "ligand_preparation_report": report_json,
        "ligand_preparation_report_csv": report_csv,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare Stage 4A ligand SDF files for docking."
        )
    )

    parser.add_argument(
        "--input-manifest",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=DEFAULT_RANDOM_SEED,
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=500,
    )

    return parser


def main() -> None:
    args = build_parser().parse_args()

    outputs = prepare_ligand_manifest(
        input_manifest=args.input_manifest,
        output_dir=args.output_dir,
        random_seed=args.random_seed,
        max_iterations=args.max_iterations,
    )

    print("[LIGAND_PREPARATION] Outputs:")

    for label, path in outputs.items():
        print(
            f"[LIGAND_PREPARATION] {label}: {path}"
        )


if __name__ == "__main__":
    main()
