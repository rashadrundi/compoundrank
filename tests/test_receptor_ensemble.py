from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.receptor_ensemble import (
    load_receptor_ensemble_manifest,
    record_receptor_ensemble_input,
)


def checksum(path: Path) -> str:
    return hashlib.sha256(
        path.read_bytes()
    ).hexdigest()


def write_ca_pdb(
    path: Path,
    *,
    x_offset: float = 0.0,
) -> None:
    path.write_text(
        (
            "ATOM      1  CA  ALA A   1    "
            f"{0.0 + x_offset:8.3f}"
            "   0.000   0.000  1.00 20.00"
            "           C\n"
            "ATOM      2  CA  GLY A   2    "
            f"{1.5 + x_offset:8.3f}"
            "   0.000   0.000  1.00 20.00"
            "           C\n"
            "ATOM      3  CA  SER A   3    "
            f"{3.0 + x_offset:8.3f}"
            "   0.000   0.000  1.00 20.00"
            "           C\n"
            "END\n"
        ),
        encoding="utf-8",
    )


def write_manifest(
    root: Path,
) -> Path:
    reference = (
        root / "reference.pdb"
    )

    snapshot = (
        root / "snapshot_0001.pdb"
    )

    write_ca_pdb(
        reference
    )

    write_ca_pdb(
        snapshot,
        x_offset=0.1,
    )

    manifest = {
        "schema_version": (
            "structure_ensemble.v0.1"
        ),
        "status": "complete",
        "selection_mode": "report_only",
        "source_engine": "openmm",
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
                checksum(reference)
            ),
            "ca_atoms": 3,
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
                    checksum(snapshot)
                ),
                "matched_ca_atoms": 3,
                "snapshot_ca_atoms": 3,
                "reference_ca_atoms": 3,
                (
                    "reference_coverage_"
                    "fraction"
                ): 1.0,
                "ca_rmsd_angstrom": 0.1,
                "rejection_reason": "",
            }
        ],
    }

    manifest_path = (
        root
        / "structure_ensemble.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    return manifest_path


class ReceptorEnsembleTests(
    unittest.TestCase
):
    def test_loads_valid_manifest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = (
                write_manifest(root)
            )

            result = (
                load_receptor_ensemble_manifest(
                    manifest_path
                )
            )

            self.assertEqual(
                result["status"],
                "accepted",
            )

            self.assertEqual(
                result[
                    "accepted_snapshot_count"
                ],
                1,
            )

            self.assertEqual(
                result[
                    "docking_behavior"
                ],
                "submitted_receptor_only",
            )

    def test_checksum_mismatch_fails(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = (
                write_manifest(root)
            )

            manifest = json.loads(
                manifest_path.read_text(
                    encoding="utf-8"
                )
            )

            manifest[
                "snapshots"
            ][0][
                "checksum_sha256"
            ] = "0" * 64

            manifest_path.write_text(
                json.dumps(
                    manifest
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "checksum mismatch",
            ):
                load_receptor_ensemble_manifest(
                    manifest_path
                )

    def test_rejects_noncomplete_manifest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = (
                write_manifest(root)
            )

            manifest = json.loads(
                manifest_path.read_text(
                    encoding="utf-8"
                )
            )

            manifest["status"] = "partial"

            manifest_path.write_text(
                json.dumps(
                    manifest
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "status 'complete'",
            ):
                load_receptor_ensemble_manifest(
                    manifest_path
                )

    def test_writes_audit_and_manifest_copy(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = (
                write_manifest(root)
            )

            output_dir = (
                root / "audit"
            )

            audit_path = (
                record_receptor_ensemble_input(
                    manifest_path,
                    output_dir,
                )
            )

            self.assertTrue(
                audit_path.is_file()
            )

            self.assertTrue(
                (
                    output_dir
                    / (
                        "source_structure_"
                        "ensemble.json"
                    )
                ).is_file()
            )

            audit = json.loads(
                audit_path.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                audit[
                    "selection_mode"
                ],
                "report_only",
            )


if __name__ == "__main__":
    unittest.main()
