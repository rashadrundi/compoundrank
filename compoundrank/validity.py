from __future__ import annotations

import csv
from pathlib import Path

from .chemistry import write_sdf
from .models import PoseRecord
from .subprocess_utils import resolve_executable, run_command


CRITICAL_BOOLEAN_COLUMNS = {
    "mol_pred_loaded",
    "sanitization",
    "inchi_convertible",
    "all_atoms_connected",
    "no_radicals",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
    "internal_energy",
    "mol_cond_loaded",
    "minimum_distance_to_protein",
    "volume_overlap_with_protein",
    "protein_ligand_maximum_distance",
}


def _as_boolean(value: str) -> bool | None:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    return None


def filter_poses_with_posebusters(
    records: list[PoseRecord],
    receptor_pdb: Path,
    work_dir: Path,
    *,
    posebusters_bin: str,
    skip: bool = False,
) -> tuple[list[PoseRecord], dict[int, list[str]]]:
    if skip:
        return records, {}
    if not records:
        return [], {}

    bust = resolve_executable(posebusters_bin, "PoseBusters")
    work_dir.mkdir(parents=True, exist_ok=True)
    input_sdf = work_dir / "posebusters_input.sdf"
    report_csv = work_dir / "posebusters_report.csv"

    for index, record in enumerate(records):
        record.molecule.SetProp("_Name", f"record_{index:06d}")
    write_sdf(input_sdf, (record.molecule for record in records))

    run_command(
        [
            bust,
            str(input_sdf),
            "-p",
            str(receptor_pdb),
            "--outfmt",
            "csv",
            "--full-report",
            "--output",
            str(report_csv),
            "--max-workers",
            "0",
        ]
    )

    if not report_csv.is_file():
        raise RuntimeError("PoseBusters did not create its CSV report")

    failures: dict[int, list[str]] = {}
    valid_indices: set[int] = set()
    with report_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader):
            molecule_name = row.get("molecule", "")
            if molecule_name.startswith("record_"):
                try:
                    index = int(molecule_name.rsplit("_", 1)[1])
                except ValueError:
                    index = row_number
            else:
                try:
                    index = int(row.get("position", row_number))
                except ValueError:
                    index = row_number

            row_failures: list[str] = []
            for column in CRITICAL_BOOLEAN_COLUMNS:
                if column not in row:
                    continue
                parsed = _as_boolean(row[column])
                if parsed is False:
                    row_failures.append(column)
            if row_failures:
                failures[index] = row_failures
            else:
                valid_indices.add(index)

    valid = [record for index, record in enumerate(records) if index in valid_indices]
    if not valid:
        raise RuntimeError(
            "Every GNINA pose failed the configured PoseBusters checks"
        )
    print(
        f"[VALIDITY] Accepted {len(valid)}/{len(records)} poses; "
        f"rejected {len(records) - len(valid)}"
    )
    return valid, failures
