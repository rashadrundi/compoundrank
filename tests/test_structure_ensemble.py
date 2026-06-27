from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from compoundrank.structure_ensemble import (
    build_structure_ensemble,
    calculate_ca_rmsd,
    _read_ca_atoms,
)


def write_ca_pdb(
    path: Path,
    coordinates: np.ndarray,
    *,
    chain: str = "A",
) -> None:
    lines: list[str] = []

    for index, point in enumerate(
        coordinates,
        start=1,
    ):
        x, y, z = (
            float(value)
            for value in point
        )

        lines.append(
            f"ATOM  {index:5d}  CA  ALA "
            f"{chain}{index:4d}    "
            f"{x:8.3f}"
            f"{y:8.3f}"
            f"{z:8.3f}"
            "  1.00 20.00           C"
        )

    lines += [
        "TER",
        "END",
    ]

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


class StructureEnsembleTests(
    unittest.TestCase
):
    def setUp(
        self,
    ) -> None:
        self.points = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.5, 0.2, 0.1],
                [0.2, 1.7, 0.3],
                [0.4, 0.6, 2.0],
                [1.7, 1.4, 1.2],
            ],
            dtype=float,
        )

    def test_kabsch_removes_rotation_and_translation(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(
                temporary
            )

            reference = (
                root / "reference.pdb"
            )

            mobile = (
                root / "mobile.pdb"
            )

            angle = np.deg2rad(
                37.0
            )

            rotation = np.asarray(
                [
                    [
                        np.cos(angle),
                        -np.sin(angle),
                        0.0,
                    ],
                    [
                        np.sin(angle),
                        np.cos(angle),
                        0.0,
                    ],
                    [
                        0.0,
                        0.0,
                        1.0,
                    ],
                ]
            )

            transformed = (
                self.points
                @ rotation.T
                + np.asarray(
                    [
                        8.0,
                        -4.0,
                        2.5,
                    ]
                )
            )

            write_ca_pdb(
                reference,
                self.points,
            )

            write_ca_pdb(
                mobile,
                transformed,
            )

            comparison = (
                calculate_ca_rmsd(
                    _read_ca_atoms(
                        reference
                    ),
                    _read_ca_atoms(
                        mobile
                    ),
                )
            )

            self.assertEqual(
                comparison[
                    "matched_ca_atoms"
                ],
                5,
            )

            self.assertLess(
                comparison[
                    "ca_rmsd_angstrom"
                ],
                0.002,
            )

    def test_builds_manifest_and_copies_snapshots(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(
                temporary
            )

            reference = (
                root / "reference.pdb"
            )

            translated = (
                root / "translated.pdb"
            )

            distorted = (
                root / "distorted.pdb"
            )

            output = (
                root / "ensemble"
            )

            write_ca_pdb(
                reference,
                self.points,
            )

            write_ca_pdb(
                translated,
                self.points
                + np.asarray(
                    [
                        4.0,
                        -3.0,
                        1.0,
                    ]
                ),
            )

            distorted_points = (
                self.points.copy()
            )

            distorted_points[
                0
            ] += np.asarray(
                [
                    2.0,
                    0.0,
                    0.0,
                ]
            )

            write_ca_pdb(
                distorted,
                distorted_points,
            )

            manifest = (
                build_structure_ensemble(
                    reference_pdb=reference,
                    snapshot_pdbs=[
                        translated,
                        distorted,
                    ],
                    output_dir=output,
                    source_engine=(
                        "test_engine"
                    ),
                    chain_id="A",
                )
            )

            self.assertEqual(
                manifest[
                    "schema_version"
                ],
                "structure_ensemble.v0.1",
            )

            self.assertEqual(
                manifest[
                    "accepted_snapshot_count"
                ],
                2,
            )

            self.assertEqual(
                manifest[
                    "rejected_snapshot_count"
                ],
                0,
            )

            first = manifest[
                "snapshots"
            ][0]

            second = manifest[
                "snapshots"
            ][1]

            self.assertLess(
                first[
                    "ca_rmsd_angstrom"
                ],
                0.002,
            )

            self.assertGreater(
                second[
                    "ca_rmsd_angstrom"
                ],
                0.2,
            )

            self.assertTrue(
                Path(
                    first[
                        "stored_path"
                    ]
                ).is_file()
            )

            json_path = Path(
                manifest[
                    "outputs"
                ]["json"]
            )

            csv_path = Path(
                manifest[
                    "outputs"
                ]["csv"]
            )

            self.assertTrue(
                json_path.is_file()
            )

            self.assertTrue(
                csv_path.is_file()
            )

            loaded = json.loads(
                json_path.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                loaded[
                    "accepted_snapshot_count"
                ],
                2,
            )

    def test_incompatible_snapshot_is_rejected(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(
                temporary
            )

            reference = (
                root / "reference.pdb"
            )

            valid = root / "valid.pdb"
            invalid = root / "invalid.pdb"

            write_ca_pdb(
                reference,
                self.points,
            )

            write_ca_pdb(
                valid,
                self.points,
            )

            write_ca_pdb(
                invalid,
                self.points[:2],
            )

            manifest = (
                build_structure_ensemble(
                    reference_pdb=reference,
                    snapshot_pdbs=[
                        valid,
                        invalid,
                    ],
                    output_dir=(
                        root / "ensemble"
                    ),
                    source_engine="test",
                    chain_id="A",
                )
            )

            self.assertEqual(
                manifest["status"],
                "complete_with_rejections",
            )

            self.assertEqual(
                manifest[
                    "accepted_snapshot_count"
                ],
                1,
            )

            self.assertEqual(
                manifest[
                    "rejected_snapshot_count"
                ],
                1,
            )

            rejected = manifest[
                "snapshots"
            ][1]

            self.assertEqual(
                rejected["status"],
                "rejected",
            )

            self.assertIn(
                "Insufficient matching CA atoms",
                rejected[
                    "rejection_reason"
                ],
            )


if __name__ == "__main__":
    unittest.main()
