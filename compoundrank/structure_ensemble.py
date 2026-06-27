from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCHEMA_VERSION = "structure_ensemble.v0.1"
ALIGNMENT_METHOD = "chain_residue_ca_kabsch.v1"


AtomKey = tuple[
    str,
    int,
    str,
    str,
]


def _sha256(
    path: Path,
) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def _read_ca_atoms(
    path: Path,
    *,
    chain_id: str | None = None,
) -> dict[AtomKey, np.ndarray]:
    source = Path(path)

    if (
        not source.is_file()
        or source.stat().st_size == 0
    ):
        raise FileNotFoundError(
            source
        )

    atoms: dict[
        AtomKey,
        np.ndarray,
    ] = {}

    for line in source.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        if not line.startswith(
            "ATOM"
        ):
            continue

        if len(line) < 54:
            continue

        atom_name = line[
            12:16
        ].strip()

        if atom_name != "CA":
            continue

        alternate_location = line[
            16:17
        ].strip()

        if alternate_location not in {
            "",
            "A",
        }:
            continue

        chain = (
            line[21:22].strip()
            or "_"
        )

        if (
            chain_id is not None
            and chain != chain_id
        ):
            continue

        try:
            residue_number = int(
                line[22:26]
            )

            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )

        except ValueError:
            continue

        insertion_code = line[
            26:27
        ].strip()

        residue_name = line[
            17:20
        ].strip().upper()

        key = (
            chain,
            residue_number,
            insertion_code,
            residue_name,
        )

        atoms.setdefault(
            key,
            np.asarray(
                (
                    x,
                    y,
                    z,
                ),
                dtype=float,
            ),
        )

    if not atoms:
        chain_text = (
            f" for chain {chain_id}"
            if chain_id is not None
            else ""
        )

        raise ValueError(
            "No usable alpha-carbon "
            f"coordinates were found{chain_text} "
            f"in {source}"
        )

    return atoms


def calculate_ca_rmsd(
    reference_atoms: dict[
        AtomKey,
        np.ndarray,
    ],
    mobile_atoms: dict[
        AtomKey,
        np.ndarray,
    ],
    *,
    minimum_matched_atoms: int = 3,
) -> dict[str, float | int]:
    common_keys = sorted(
        set(reference_atoms)
        & set(mobile_atoms)
    )

    matched_count = len(
        common_keys
    )

    if (
        matched_count
        < minimum_matched_atoms
    ):
        raise ValueError(
            "Insufficient matching CA atoms: "
            f"{matched_count} found; "
            f"{minimum_matched_atoms} required"
        )

    reference = np.vstack(
        [
            reference_atoms[key]
            for key in common_keys
        ]
    )

    mobile = np.vstack(
        [
            mobile_atoms[key]
            for key in common_keys
        ]
    )

    reference_centered = (
        reference
        - reference.mean(
            axis=0
        )
    )

    mobile_centered = (
        mobile
        - mobile.mean(
            axis=0
        )
    )

    covariance = (
        mobile_centered.T
        @ reference_centered
    )

    left, _, right_transpose = (
        np.linalg.svd(
            covariance
        )
    )

    rotation = (
        left
        @ right_transpose
    )

    if np.linalg.det(
        rotation
    ) < 0:
        left[:, -1] *= -1.0

        rotation = (
            left
            @ right_transpose
        )

    aligned_mobile = (
        mobile_centered
        @ rotation
    )

    differences = (
        aligned_mobile
        - reference_centered
    )

    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(
                    differences
                    * differences,
                    axis=1,
                )
            )
        )
    )

    return {
        "matched_ca_atoms": (
            matched_count
        ),
        "reference_ca_atoms": (
            len(reference_atoms)
        ),
        "mobile_ca_atoms": (
            len(mobile_atoms)
        ),
        "reference_coverage_fraction": (
            matched_count
            / len(reference_atoms)
        ),
        "ca_rmsd_angstrom": rmsd,
    }


def _copy_structure(
    source: Path,
    destination: Path,
) -> Path:
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if (
        source.resolve()
        != destination.resolve()
    ):
        shutil.copy2(
            source,
            destination,
        )

    return destination


def write_structure_ensemble_outputs(
    manifest: dict[str, Any],
    output_dir: Path,
) -> dict[str, Path]:
    destination = Path(
        output_dir
    )

    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    json_path = (
        destination
        / "structure_ensemble.json"
    )

    csv_path = (
        destination
        / "structure_ensemble_snapshots.csv"
    )

    json_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    fieldnames = [
        "snapshot_id",
        "status",
        "source_path",
        "stored_path",
        "checksum_sha256",
        "matched_ca_atoms",
        "reference_ca_atoms",
        "snapshot_ca_atoms",
        "reference_coverage_fraction",
        "ca_rmsd_angstrom",
        "rejection_reason",
    ]

    with csv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for snapshot in manifest[
            "snapshots"
        ]:
            writer.writerow(
                {
                    field: snapshot.get(
                        field,
                        "",
                    )
                    for field in fieldnames
                }
            )

    return {
        "json": json_path,
        "csv": csv_path,
    }


def build_structure_ensemble(
    *,
    reference_pdb: Path,
    snapshot_pdbs: Iterable[Path],
    output_dir: Path,
    source_engine: str = "external",
    chain_id: str | None = None,
    minimum_matched_atoms: int = 3,
    overwrite: bool = False,
) -> dict[str, Any]:
    reference_source = Path(
        reference_pdb
    ).expanduser().resolve()

    snapshot_sources = [
        Path(path)
        .expanduser()
        .resolve()
        for path in snapshot_pdbs
    ]

    destination = Path(
        output_dir
    ).expanduser().resolve()

    if not snapshot_sources:
        raise ValueError(
            "At least one snapshot PDB "
            "is required"
        )

    manifest_path = (
        destination
        / "structure_ensemble.json"
    )

    if (
        manifest_path.exists()
        and not overwrite
    ):
        raise FileExistsError(
            "Structure ensemble output "
            f"already exists: {manifest_path}"
        )

    if overwrite:
        for stale_path in (
            destination
            / "reference",
            destination
            / "snapshots",
        ):
            if stale_path.exists():
                shutil.rmtree(
                    stale_path
                )

        for stale_path in (
            manifest_path,
            destination
            / "structure_ensemble_snapshots.csv",
        ):
            if stale_path.exists():
                stale_path.unlink()

    reference_atoms = (
        _read_ca_atoms(
            reference_source,
            chain_id=chain_id,
        )
    )

    stored_reference = (
        _copy_structure(
            reference_source,
            destination
            / "reference"
            / "reference.pdb",
        )
    )

    snapshot_rows: list[
        dict[str, Any]
    ] = []

    valid_count = 0
    rejected_count = 0

    for index, source in enumerate(
        snapshot_sources,
        start=1,
    ):
        snapshot_id = (
            f"snapshot_{index:04d}"
        )

        base_row: dict[str, Any] = {
            "snapshot_id": snapshot_id,
            "source_path": str(
                source
            ),
            "stored_path": "",
            "checksum_sha256": "",
            "matched_ca_atoms": 0,
            "reference_ca_atoms": (
                len(reference_atoms)
            ),
            "snapshot_ca_atoms": 0,
            "reference_coverage_fraction": 0.0,
            "ca_rmsd_angstrom": "",
            "rejection_reason": "",
        }

        try:
            mobile_atoms = (
                _read_ca_atoms(
                    source,
                    chain_id=chain_id,
                )
            )

            comparison = (
                calculate_ca_rmsd(
                    reference_atoms,
                    mobile_atoms,
                    minimum_matched_atoms=(
                        minimum_matched_atoms
                    ),
                )
            )

            stored_path = (
                _copy_structure(
                    source,
                    destination
                    / "snapshots"
                    / f"{snapshot_id}.pdb",
                )
            )

            base_row.update(
                {
                    "status": "accepted",
                    "stored_path": str(
                        stored_path
                    ),
                    "checksum_sha256": (
                        _sha256(
                            source
                        )
                    ),
                    "matched_ca_atoms": (
                        comparison[
                            "matched_ca_atoms"
                        ]
                    ),
                    "reference_ca_atoms": (
                        comparison[
                            "reference_ca_atoms"
                        ]
                    ),
                    "snapshot_ca_atoms": (
                        comparison[
                            "mobile_ca_atoms"
                        ]
                    ),
                    "reference_coverage_fraction": (
                        comparison[
                            "reference_coverage_fraction"
                        ]
                    ),
                    "ca_rmsd_angstrom": (
                        comparison[
                            "ca_rmsd_angstrom"
                        ]
                    ),
                }
            )

            valid_count += 1

        except Exception as error:
            base_row.update(
                {
                    "status": "rejected",
                    "rejection_reason": (
                        f"{type(error).__name__}: "
                        f"{error}"
                    ),
                }
            )

            rejected_count += 1

        snapshot_rows.append(
            base_row
        )

    if valid_count == 0:
        raise ValueError(
            "No compatible receptor "
            "snapshots were accepted"
        )

    status = (
        "complete"
        if rejected_count == 0
        else "complete_with_rejections"
    )

    manifest: dict[str, Any] = {
        "schema_version": (
            SCHEMA_VERSION
        ),
        "status": status,
        "selection_mode": (
            "report_only"
        ),
        "source_engine": (
            source_engine
        ),
        "alignment_method": (
            ALIGNMENT_METHOD
        ),
        "chain": chain_id,
        "minimum_matched_ca_atoms": (
            minimum_matched_atoms
        ),
        "reference": {
            "source_path": str(
                reference_source
            ),
            "stored_path": str(
                stored_reference
            ),
            "checksum_sha256": (
                _sha256(
                    reference_source
                )
            ),
            "ca_atoms": len(
                reference_atoms
            ),
        },
        "snapshot_count": len(
            snapshot_rows
        ),
        "accepted_snapshot_count": (
            valid_count
        ),
        "rejected_snapshot_count": (
            rejected_count
        ),
        "snapshots": snapshot_rows,
        "limitations": [
            (
                "This manifest records and "
                "geometrically compares receptor "
                "snapshots. It does not establish "
                "that an MD simulation was "
                "equilibrated or physically valid."
            ),
            (
                "CA RMSD is calculated after "
                "least-squares Kabsch alignment "
                "over matching chain, residue "
                "number, insertion code, and "
                "residue name."
            ),
            (
                "Snapshot acceptance does not "
                "yet change docking or pocket "
                "selection."
            ),
        ],
    }

    outputs = (
        write_structure_ensemble_outputs(
            manifest,
            destination,
        )
    )

    manifest["outputs"] = {
        name: str(path)
        for name, path
        in outputs.items()
    }

    # Rewrite once so the JSON contains
    # its own output locations.
    write_structure_ensemble_outputs(
        manifest,
        destination,
    )

    return manifest


def build_cli_parser() -> (
    argparse.ArgumentParser
):
    parser = argparse.ArgumentParser(
        description=(
            "Create a portable receptor-"
            "ensemble manifest from PDB "
            "snapshots generated by OpenMM, "
            "GROMACS, Amber, NAMD, or another "
            "external workflow."
        )
    )

    parser.add_argument(
        "--reference-pdb",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--snapshot",
        action="append",
        type=Path,
        default=[],
    )

    parser.add_argument(
        "--snapshot-dir",
        type=Path,
        default=None,
    )

    parser.add_argument(
        "--snapshot-glob",
        default="*.pdb",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--source-engine",
        default="external",
    )

    parser.add_argument(
        "--chain",
        default=None,
    )

    parser.add_argument(
        "--minimum-matched-atoms",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    return parser


def main() -> int:
    parser = build_cli_parser()
    args = parser.parse_args()

    snapshots = [
        Path(path)
        for path in args.snapshot
    ]

    if args.snapshot_dir is not None:
        snapshots.extend(
            sorted(
                Path(
                    args.snapshot_dir
                ).glob(
                    args.snapshot_glob
                )
            )
        )

    unique_snapshots: list[
        Path
    ] = []

    observed: set[Path] = set()

    for snapshot in snapshots:
        resolved = (
            snapshot
            .expanduser()
            .resolve()
        )

        if resolved in observed:
            continue

        observed.add(
            resolved
        )

        unique_snapshots.append(
            resolved
        )

    manifest = build_structure_ensemble(
        reference_pdb=(
            args.reference_pdb
        ),
        snapshot_pdbs=(
            unique_snapshots
        ),
        output_dir=(
            args.output_dir
        ),
        source_engine=(
            args.source_engine
        ),
        chain_id=args.chain,
        minimum_matched_atoms=(
            args.minimum_matched_atoms
        ),
        overwrite=args.overwrite,
    )

    print(
        json.dumps(
            {
                "status": (
                    manifest["status"]
                ),
                "source_engine": (
                    manifest[
                        "source_engine"
                    ]
                ),
                "accepted_snapshots": (
                    manifest[
                        "accepted_snapshot_count"
                    ]
                ),
                "rejected_snapshots": (
                    manifest[
                        "rejected_snapshot_count"
                    ]
                ),
                "outputs": (
                    manifest["outputs"]
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
