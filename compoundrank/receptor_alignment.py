from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any, Sequence

import gemmi
import numpy as np

from .receptor_ensemble import (
    load_receptor_ensemble_manifest,
)
from .structure_ensemble import (
    _read_ca_atoms,
)


SCHEMA_VERSION = (
    "aligned_receptor_ensemble.v0.1"
)

ALIGNMENT_METHOD = (
    "chain_residue_ca_kabsch.v1"
)


def _sha256(path: Path) -> str:
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


def _matching_coordinate_arrays(
    reference_atoms: dict[Any, np.ndarray],
    mobile_atoms: dict[Any, np.ndarray],
    *,
    minimum_matched_atoms: int,
) -> tuple[
    list[Any],
    np.ndarray,
    np.ndarray,
]:
    common_keys = sorted(
        set(reference_atoms)
        & set(mobile_atoms)
    )

    if (
        len(common_keys)
        < minimum_matched_atoms
    ):
        raise ValueError(
            "Insufficient matching CA atoms: "
            f"{len(common_keys)} found; "
            f"{minimum_matched_atoms} required"
        )

    reference = np.vstack(
        [
            reference_atoms[key]
            for key in common_keys
        ]
    ).astype(float)

    mobile = np.vstack(
        [
            mobile_atoms[key]
            for key in common_keys
        ]
    ).astype(float)

    return (
        common_keys,
        reference,
        mobile,
    )


def _rmsd(
    first: np.ndarray,
    second: np.ndarray,
) -> float:
    difference = (
        first
        - second
    )

    return float(
        math.sqrt(
            float(
                np.mean(
                    np.sum(
                        difference
                        * difference,
                        axis=1,
                    )
                )
            )
        )
    )


def calculate_alignment_transform(
    reference_atoms: dict[Any, np.ndarray],
    mobile_atoms: dict[Any, np.ndarray],
    *,
    minimum_matched_atoms: int = 3,
) -> dict[str, Any]:
    (
        common_keys,
        reference,
        mobile,
    ) = _matching_coordinate_arrays(
        reference_atoms,
        mobile_atoms,
        minimum_matched_atoms=(
            minimum_matched_atoms
        ),
    )

    reference_centroid = (
        reference.mean(
            axis=0
        )
    )

    mobile_centroid = (
        mobile.mean(
            axis=0
        )
    )

    reference_centered = (
        reference
        - reference_centroid
    )

    mobile_centered = (
        mobile
        - mobile_centroid
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
        (
            mobile
            - mobile_centroid
        )
        @ rotation
        + reference_centroid
    )

    return {
        "matched_keys": common_keys,
        "matched_ca_atoms": len(
            common_keys
        ),
        "reference_ca_atoms": len(
            reference_atoms
        ),
        "mobile_ca_atoms": len(
            mobile_atoms
        ),
        "reference_coverage_fraction": (
            len(common_keys)
            / len(reference_atoms)
        ),
        "raw_ca_rmsd_angstrom": (
            _rmsd(
                mobile,
                reference,
            )
        ),
        "aligned_ca_rmsd_angstrom": (
            _rmsd(
                aligned_mobile,
                reference,
            )
        ),
        (
            "centroid_displacement_"
            "angstrom"
        ): float(
            np.linalg.norm(
                mobile_centroid
                - reference_centroid
            )
        ),
        "rotation_matrix": rotation,
        "mobile_centroid": (
            mobile_centroid
        ),
        "reference_centroid": (
            reference_centroid
        ),
    }


def _transform_position(
    coordinate: np.ndarray,
    *,
    rotation: np.ndarray,
    mobile_centroid: np.ndarray,
    reference_centroid: np.ndarray,
) -> np.ndarray:
    return (
        (
            coordinate
            - mobile_centroid
        )
        @ rotation
        + reference_centroid
    )


def _apply_transform_to_pdb(
    source_pdb: Path,
    output_pdb: Path,
    *,
    rotation: np.ndarray,
    mobile_centroid: np.ndarray,
    reference_centroid: np.ndarray,
) -> None:
    structure = gemmi.read_structure(
        str(source_pdb)
    )

    atom_count = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinate = np.asarray(
                        [
                            atom.pos.x,
                            atom.pos.y,
                            atom.pos.z,
                        ],
                        dtype=float,
                    )

                    transformed = (
                        _transform_position(
                            coordinate,
                            rotation=rotation,
                            mobile_centroid=(
                                mobile_centroid
                            ),
                            reference_centroid=(
                                reference_centroid
                            ),
                        )
                    )

                    atom.pos = gemmi.Position(
                        float(
                            transformed[0]
                        ),
                        float(
                            transformed[1]
                        ),
                        float(
                            transformed[2]
                        ),
                    )

                    atom_count += 1

    if atom_count == 0:
        raise ValueError(
            "No atoms were found in the "
            f"mobile structure: {source_pdb}"
        )

    output_pdb.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    structure.write_pdb(
        str(output_pdb)
    )

    if (
        not output_pdb.is_file()
        or output_pdb.stat().st_size == 0
    ):
        raise RuntimeError(
            "Aligned receptor PDB was not "
            f"created: {output_pdb}"
        )


def align_receptor_structure(
    *,
    reference_pdb: Path,
    mobile_pdb: Path,
    output_pdb: Path,
    chain_id: str | None = None,
    minimum_matched_atoms: int = 3,
    overwrite: bool = False,
) -> dict[str, Any]:
    reference = (
        Path(reference_pdb)
        .expanduser()
        .resolve()
    )

    mobile = (
        Path(mobile_pdb)
        .expanduser()
        .resolve()
    )

    output = (
        Path(output_pdb)
        .expanduser()
        .resolve()
    )

    for label, path in (
        (
            "Reference receptor",
            reference,
        ),
        (
            "Mobile receptor",
            mobile,
        ),
    ):
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise FileNotFoundError(
                f"{label} is missing or "
                f"empty: {path}"
            )

    if output.exists() and not overwrite:
        raise FileExistsError(
            "Aligned receptor output already "
            f"exists: {output}"
        )

    reference_atoms = (
        _read_ca_atoms(
            reference,
            chain_id=chain_id,
        )
    )

    mobile_atoms = (
        _read_ca_atoms(
            mobile,
            chain_id=chain_id,
        )
    )

    transform = (
        calculate_alignment_transform(
            reference_atoms,
            mobile_atoms,
            minimum_matched_atoms=(
                minimum_matched_atoms
            ),
        )
    )

    _apply_transform_to_pdb(
        mobile,
        output,
        rotation=transform[
            "rotation_matrix"
        ],
        mobile_centroid=transform[
            "mobile_centroid"
        ],
        reference_centroid=transform[
            "reference_centroid"
        ],
    )

    aligned_atoms = (
        _read_ca_atoms(
            output,
            chain_id=chain_id,
        )
    )

    (
        _,
        aligned_reference,
        aligned_mobile,
    ) = _matching_coordinate_arrays(
        reference_atoms,
        aligned_atoms,
        minimum_matched_atoms=(
            minimum_matched_atoms
        ),
    )

    raw_after_write = _rmsd(
        aligned_mobile,
        aligned_reference,
    )

    return {
        "alignment_method": (
            ALIGNMENT_METHOD
        ),
        "reference_pdb": str(
            reference
        ),
        "mobile_pdb": str(
            mobile
        ),
        "aligned_pdb": str(
            output
        ),
        "aligned_checksum_sha256": (
            _sha256(output)
        ),
        "matched_ca_atoms": (
            transform[
                "matched_ca_atoms"
            ]
        ),
        "reference_ca_atoms": (
            transform[
                "reference_ca_atoms"
            ]
        ),
        "mobile_ca_atoms": (
            transform[
                "mobile_ca_atoms"
            ]
        ),
        (
            "reference_coverage_"
            "fraction"
        ): transform[
            "reference_coverage_fraction"
        ],
        (
            "raw_ca_rmsd_before_"
            "alignment_angstrom"
        ): transform[
            "raw_ca_rmsd_angstrom"
        ],
        (
            "kabsch_ca_rmsd_"
            "angstrom"
        ): transform[
            "aligned_ca_rmsd_angstrom"
        ],
        (
            "raw_ca_rmsd_after_"
            "alignment_angstrom"
        ): raw_after_write,
        (
            "centroid_displacement_"
            "before_alignment_angstrom"
        ): transform[
            (
                "centroid_displacement_"
                "angstrom"
            )
        ],
    }


def build_aligned_receptor_ensemble(
    *,
    ensemble_manifest: Path,
    output_dir: Path,
    verify_checksums: bool = True,
    overwrite: bool = False,
) -> dict[str, Any]:
    source_manifest = (
        Path(ensemble_manifest)
        .expanduser()
        .resolve()
    )

    destination = (
        Path(output_dir)
        .expanduser()
        .resolve()
    )

    if (
        destination.exists()
        and any(destination.iterdir())
        and not overwrite
    ):
        raise FileExistsError(
            "Aligned ensemble directory is "
            f"not empty: {destination}"
        )

    if (
        overwrite
        and destination.exists()
    ):
        shutil.rmtree(
            destination
        )

    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    imported = (
        load_receptor_ensemble_manifest(
            source_manifest,
            verify_checksums=(
                verify_checksums
            ),
        )
    )

    source_copy = (
        destination
        / "source_structure_ensemble.json"
    )

    shutil.copy2(
        source_manifest,
        source_copy,
    )

    reference_source = Path(
        imported[
            "reference"
        ][
            "stored_path"
        ]
    )

    reference_output = (
        destination
        / "reference"
        / "reference.pdb"
    )

    reference_output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    shutil.copy2(
        reference_source,
        reference_output,
    )

    aligned_rows: list[
        dict[str, Any]
    ] = []

    for snapshot in imported[
        "accepted_snapshots"
    ]:
        snapshot_id = str(
            snapshot["snapshot_id"]
        )

        aligned_output = (
            destination
            / "snapshots"
            / f"{snapshot_id}.pdb"
        )

        alignment = (
            align_receptor_structure(
                reference_pdb=(
                    reference_output
                ),
                mobile_pdb=Path(
                    snapshot[
                        "stored_path"
                    ]
                ),
                output_pdb=(
                    aligned_output
                ),
                minimum_matched_atoms=3,
                overwrite=True,
            )
        )

        aligned_rows.append(
            {
                "snapshot_id": (
                    snapshot_id
                ),
                "status": (
                    "aligned"
                ),
                "source_snapshot_path": (
                    snapshot[
                        "stored_path"
                    ]
                ),
                "source_checksum_sha256": (
                    snapshot[
                        "checksum_sha256"
                    ]
                ),
                "aligned_path": (
                    alignment[
                        "aligned_pdb"
                    ]
                ),
                "aligned_checksum_sha256": (
                    alignment[
                        (
                            "aligned_checksum_"
                            "sha256"
                        )
                    ]
                ),
                "matched_ca_atoms": (
                    alignment[
                        "matched_ca_atoms"
                    ]
                ),
                (
                    "reference_coverage_"
                    "fraction"
                ): alignment[
                    (
                        "reference_coverage_"
                        "fraction"
                    )
                ],
                (
                    "raw_ca_rmsd_before_"
                    "alignment_angstrom"
                ): alignment[
                    (
                        "raw_ca_rmsd_before_"
                        "alignment_angstrom"
                    )
                ],
                (
                    "kabsch_ca_rmsd_"
                    "angstrom"
                ): alignment[
                    (
                        "kabsch_ca_rmsd_"
                        "angstrom"
                    )
                ],
                (
                    "raw_ca_rmsd_after_"
                    "alignment_angstrom"
                ): alignment[
                    (
                        "raw_ca_rmsd_after_"
                        "alignment_angstrom"
                    )
                ],
                (
                    "centroid_displacement_"
                    "before_alignment_"
                    "angstrom"
                ): alignment[
                    (
                        "centroid_displacement_"
                        "before_alignment_"
                        "angstrom"
                    )
                ],
            }
        )

    if not aligned_rows:
        raise ValueError(
            "No accepted receptor snapshots "
            "were available for alignment"
        )

    output_manifest = (
        destination
        / "aligned_receptor_ensemble.json"
    )

    manifest: dict[str, Any] = {
        "schema_version": (
            SCHEMA_VERSION
        ),
        "status": "complete",
        "selection_mode": (
            "report_only"
        ),
        "alignment_method": (
            ALIGNMENT_METHOD
        ),
        "source_manifest": str(
            source_manifest
        ),
        (
            "source_manifest_checksum_"
            "sha256"
        ): _sha256(
            source_manifest
        ),
        "copied_source_manifest": str(
            source_copy
        ),
        "reference": {
            "source_path": str(
                reference_source
            ),
            "aligned_reference_path": str(
                reference_output
            ),
            "checksum_sha256": (
                _sha256(
                    reference_output
                )
            ),
            "ca_atoms": imported[
                "reference"
            ][
                "ca_atoms"
            ],
        },
        "snapshot_count": len(
            aligned_rows
        ),
        "snapshots": aligned_rows,
        "docking_behavior": (
            "not_enabled"
        ),
        "limitations": [
            (
                "Snapshots are rigidly aligned "
                "to the reference coordinate "
                "frame using matching CA atoms."
            ),
            (
                "Alignment does not alter the "
                "internal conformation of each "
                "snapshot."
            ),
            (
                "This output does not yet cause "
                "the main pipeline to dock "
                "against multiple receptors."
            ),
        ],
        "outputs": {
            "json": str(
                output_manifest
            ),
        },
    }

    output_manifest.write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return manifest


def build_cli_parser() -> (
    argparse.ArgumentParser
):
    parser = argparse.ArgumentParser(
        description=(
            "Rigidly align accepted receptor "
            "ensemble snapshots to their "
            "reference coordinate frame."
        )
    )

    parser.add_argument(
        "--ensemble-json",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--no-verify-checksums",
        action="store_true",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    return parser


def main(
    argv: Sequence[str] | None = None,
) -> int:
    arguments = (
        build_cli_parser().parse_args(
            argv
        )
    )

    manifest = (
        build_aligned_receptor_ensemble(
            ensemble_manifest=(
                arguments.ensemble_json
            ),
            output_dir=(
                arguments.output_dir
            ),
            verify_checksums=(
                not arguments
                .no_verify_checksums
            ),
            overwrite=(
                arguments.overwrite
            ),
        )
    )

    print(
        json.dumps(
            {
                "status": (
                    manifest["status"]
                ),
                "snapshot_count": (
                    manifest[
                        "snapshot_count"
                    ]
                ),
                "manifest": (
                    manifest[
                        "outputs"
                    ][
                        "json"
                    ]
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
