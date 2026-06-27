from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from compoundrank.receptor_alignment import (
    align_receptor_structure,
    build_aligned_receptor_ensemble,
)


REFERENCE_COORDINATES = np.asarray(
    [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 4.0],
    ],
    dtype=float,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(
        path.read_bytes()
    ).hexdigest()


def _atom_line(
    *,
    serial: int,
    residue_name: str,
    residue_number: int,
    coordinate: np.ndarray,
) -> str:
    return (
        f"ATOM  {serial:5d}  CA  "
        f"{residue_name:>3s} A"
        f"{residue_number:4d}    "
        f"{coordinate[0]:8.3f}"
        f"{coordinate[1]:8.3f}"
        f"{coordinate[2]:8.3f}"
        "  1.00 20.00           C\n"
    )


def _write_pdb(
    path: Path,
    coordinates: np.ndarray,
) -> None:
    residue_names = [
        "ALA",
        "GLY",
        "SER",
        "THR",
    ]

    lines = [
        _atom_line(
            serial=index,
            residue_name=(
                residue_names[
                    index - 1
                ]
            ),
            residue_number=index,
            coordinate=coordinate,
        )
        for index, coordinate
        in enumerate(
            coordinates,
            start=1,
        )
    ]

    path.write_text(
        "".join(
            [
                *lines,
                "TER\n",
                "END\n",
            ]
        ),
        encoding="utf-8",
    )


def _mobile_coordinates() -> np.ndarray:
    rotation = np.asarray(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )

    translation = np.asarray(
        [10.0, -3.0, 5.0],
        dtype=float,
    )

    return (
        REFERENCE_COORDINATES
        @ rotation
        + translation
    )


def _write_source_manifest(
    root: Path,
) -> Path:
    reference = (
        root / "reference.pdb"
    )

    snapshot = (
        root / "snapshot.pdb"
    )

    _write_pdb(
        reference,
        REFERENCE_COORDINATES,
    )

    _write_pdb(
        snapshot,
        _mobile_coordinates(),
    )

    manifest = {
        "schema_version": (
            "structure_ensemble.v0.1"
        ),
        "status": "complete",
        "selection_mode": "report_only",
        "source_engine": "test",
        "snapshot_count": 1,
        "accepted_snapshot_count": 1,
        "rejected_snapshot_count": 0,
        "reference": {
            "stored_path": str(
                reference.resolve()
            ),
            "source_path": str(
                reference.resolve()
            ),
            "checksum_sha256": (
                _sha256(reference)
            ),
            "ca_atoms": 4,
        },
        "snapshots": [
            {
                "snapshot_id": (
                    "snapshot_0001"
                ),
                "status": "accepted",
                "stored_path": str(
                    snapshot.resolve()
                ),
                "source_path": str(
                    snapshot.resolve()
                ),
                "checksum_sha256": (
                    _sha256(snapshot)
                ),
                "matched_ca_atoms": 4,
                "snapshot_ca_atoms": 4,
                "reference_ca_atoms": 4,
                (
                    "reference_coverage_"
                    "fraction"
                ): 1.0,
                "ca_rmsd_angstrom": 0.0,
                "rejection_reason": "",
            }
        ],
    }

    path = (
        root
        / "structure_ensemble.json"
    )

    path.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    return path


class ReceptorAlignmentTests(
    unittest.TestCase
):
    def test_aligns_rigid_transform(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            reference = (
                root / "reference.pdb"
            )

            mobile = (
                root / "mobile.pdb"
            )

            output = (
                root / "aligned.pdb"
            )

            _write_pdb(
                reference,
                REFERENCE_COORDINATES,
            )

            _write_pdb(
                mobile,
                _mobile_coordinates(),
            )

            result = (
                align_receptor_structure(
                    reference_pdb=reference,
                    mobile_pdb=mobile,
                    output_pdb=output,
                )
            )

            self.assertGreater(
                result[
                    (
                        "raw_ca_rmsd_before_"
                        "alignment_angstrom"
                    )
                ],
                1.0,
            )

            self.assertLess(
                result[
                    (
                        "kabsch_ca_rmsd_"
                        "angstrom"
                    )
                ],
                1e-8,
            )

            self.assertLess(
                result[
                    (
                        "raw_ca_rmsd_after_"
                        "alignment_angstrom"
                    )
                ],
                0.002,
            )

    def test_builds_aligned_ensemble(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            source_manifest = (
                _write_source_manifest(
                    root
                )
            )

            output = (
                root / "aligned"
            )

            manifest = (
                build_aligned_receptor_ensemble(
                    ensemble_manifest=(
                        source_manifest
                    ),
                    output_dir=output,
                )
            )

            self.assertEqual(
                manifest[
                    "schema_version"
                ],
                (
                    "aligned_receptor_"
                    "ensemble.v0.1"
                ),
            )

            self.assertEqual(
                manifest[
                    "snapshot_count"
                ],
                1,
            )

            aligned_path = Path(
                manifest[
                    "snapshots"
                ][0][
                    "aligned_path"
                ]
            )

            self.assertTrue(
                aligned_path.is_file()
            )

            self.assertLess(
                manifest[
                    "snapshots"
                ][0][
                    (
                        "raw_ca_rmsd_after_"
                        "alignment_angstrom"
                    )
                ],
                0.002,
            )

    def test_refuses_existing_output(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            source_manifest = (
                _write_source_manifest(
                    root
                )
            )

            output = (
                root / "aligned"
            )

            output.mkdir()
            (
                output / "existing.txt"
            ).write_text(
                "occupied",
                encoding="utf-8",
            )

            with self.assertRaises(
                FileExistsError
            ):
                build_aligned_receptor_ensemble(
                    ensemble_manifest=(
                        source_manifest
                    ),
                    output_dir=output,
                )


if __name__ == "__main__":
    unittest.main()
