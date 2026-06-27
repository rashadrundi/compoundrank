from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.aligned_receptor_ensemble import (
    load_aligned_receptor_ensemble_manifest,
    record_aligned_receptor_ensemble_input,
    validate_receptor_ensemble_options,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(
        path.read_bytes()
    ).hexdigest()


def _write_file(
    path: Path,
    content: str,
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        content,
        encoding="utf-8",
    )


def _write_manifest(
    root: Path,
    *,
    snapshot_id: str = "snapshot_0001",
    duplicate_snapshot: bool = False,
    submitted_matches: bool = True,
) -> tuple[Path, Path, Path]:
    submitted = (
        root / "submitted.pdb"
    )

    reference = (
        root / "reference.pdb"
    )

    snapshot = (
        root / "snapshot.pdb"
    )

    submitted_text = (
        "ATOM      1  CA  ALA A   1"
        "       0.000   0.000   0.000\n"
        "END\n"
    )

    reference_text = (
        submitted_text
        if submitted_matches
        else (
            "ATOM      1  CA  GLY A   1"
            "       1.000   0.000   0.000\n"
            "END\n"
        )
    )

    _write_file(
        submitted,
        submitted_text,
    )

    _write_file(
        reference,
        reference_text,
    )

    _write_file(
        snapshot,
        (
            "ATOM      1  CA  ALA A   1"
            "       0.100   0.000   0.000\n"
            "END\n"
        ),
    )

    row = {
        "snapshot_id": snapshot_id,
        "status": "aligned",
        "aligned_path": str(
            snapshot.resolve()
        ),
        "aligned_checksum_sha256": (
            _sha256(snapshot)
        ),
        "matched_ca_atoms": 1,
        (
            "reference_coverage_"
            "fraction"
        ): 1.0,
        (
            "kabsch_ca_rmsd_"
            "angstrom"
        ): 0.1,
        (
            "raw_ca_rmsd_after_"
            "alignment_angstrom"
        ): 0.1,
    }

    rows = [row]

    if duplicate_snapshot:
        rows.append(
            dict(row)
        )

    manifest = {
        "schema_version": (
            "aligned_receptor_ensemble.v0.1"
        ),
        "status": "complete",
        "selection_mode": "report_only",
        "alignment_method": (
            "chain_residue_ca_kabsch.v1"
        ),
        "docking_behavior": "not_enabled",
        "reference": {
            "aligned_reference_path": str(
                reference.resolve()
            ),
            "checksum_sha256": (
                _sha256(reference)
            ),
            "ca_atoms": 1,
        },
        "snapshot_count": len(rows),
        "snapshots": rows,
    }

    manifest_path = (
        root
        / "aligned_receptor_ensemble.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    return (
        manifest_path,
        submitted,
        snapshot,
    )


class AlignedReceptorEnsembleTests(
    unittest.TestCase
):
    def test_loads_valid_manifest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            (
                manifest_path,
                submitted,
                snapshot,
            ) = _write_manifest(
                root
            )

            result = (
                load_aligned_receptor_ensemble_manifest(
                    manifest_path,
                    submitted_receptor_pdb=(
                        submitted
                    ),
                )
            )

            self.assertEqual(
                result[
                    "schema_version"
                ],
                (
                    "aligned_receptor_"
                    "ensemble_input.v0.1"
                ),
            )

            self.assertEqual(
                result["snapshot_count"],
                1,
            )

            self.assertEqual(
                result[
                    "snapshots"
                ][0][
                    "conformer_id"
                ],
                "snapshot_0001",
            )

            self.assertEqual(
                Path(
                    result[
                        "snapshots"
                    ][0][
                        "aligned_path"
                    ]
                ),
                snapshot.resolve(),
            )

    def test_rejects_reference_mismatch(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            (
                manifest_path,
                submitted,
                _,
            ) = _write_manifest(
                root,
                submitted_matches=False,
            )

            with self.assertRaisesRegex(
                ValueError,
                "does not exactly match",
            ):
                load_aligned_receptor_ensemble_manifest(
                    manifest_path,
                    submitted_receptor_pdb=(
                        submitted
                    ),
                )

    def test_rejects_path_separator(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            (
                manifest_path,
                submitted,
                _,
            ) = _write_manifest(
                root,
                snapshot_id=(
                    "../snapshot"
                ),
            )

            with self.assertRaisesRegex(
                ValueError,
                "path separators",
            ):
                load_aligned_receptor_ensemble_manifest(
                    manifest_path,
                    submitted_receptor_pdb=(
                        submitted
                    ),
                )

    def test_rejects_duplicate_ids(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            (
                manifest_path,
                submitted,
                _,
            ) = _write_manifest(
                root,
                duplicate_snapshot=True,
            )

            with self.assertRaisesRegex(
                ValueError,
                "Duplicate",
            ):
                load_aligned_receptor_ensemble_manifest(
                    manifest_path,
                    submitted_receptor_pdb=(
                        submitted
                    ),
                )

    def test_records_audit_and_copy(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            (
                manifest_path,
                submitted,
                _,
            ) = _write_manifest(
                root
            )

            audit, audit_path = (
                record_aligned_receptor_ensemble_input(
                    manifest_path=(
                        manifest_path
                    ),
                    submitted_receptor_pdb=(
                        submitted
                    ),
                    output_dir=(
                        root / "audit"
                    ),
                )
            )

            self.assertTrue(
                audit_path.is_file()
            )

            self.assertTrue(
                Path(
                    audit[
                        "outputs"
                    ][
                        "source_manifest_copy"
                    ]
                ).is_file()
            )

    def test_rejects_both_input_modes(
        self,
    ) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "either",
        ):
            validate_receptor_ensemble_options(
                Path("/tmp/report.json"),
                Path("/tmp/aligned.json"),
            )


if __name__ == "__main__":
    unittest.main()
