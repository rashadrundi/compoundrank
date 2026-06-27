from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from compoundrank.pipeline import (
    _run_receptor_conformer_validations,
)


def _write(
    path: Path,
    text: str = "ATOM\n",
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    path.write_text(
        text,
        encoding="utf-8",
    )
    return path


class PipelineConformerRamachandranTests(
    unittest.TestCase
):
    def test_validates_each_snapshot_once(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)

            submitted_pdb = _write(
                root / "submitted.pdb"
            )
            snapshot_1 = _write(
                root / "snapshot_0001.pdb"
            )
            snapshot_2 = _write(
                root / "snapshot_0002.pdb"
            )

            conformers = [
                (
                    "submitted_receptor",
                    SimpleNamespace(
                        source_pdb=submitted_pdb
                    ),
                ),
                (
                    "snapshot_0001",
                    SimpleNamespace(
                        source_pdb=snapshot_1
                    ),
                ),
                (
                    "snapshot_0002",
                    SimpleNamespace(
                        source_pdb=snapshot_2
                    ),
                ),
            ]

            calls: list[
                dict[str, object]
            ] = []

            def fake_validation(
                *,
                structure_path,
                output_dir,
                chain_id,
                label,
            ):
                output = Path(output_dir)
                output.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                json_path = _write(
                    output
                    / "ramachandran_validation.json",
                    "{}\n",
                )
                csv_path = _write(
                    output
                    / "ramachandran_residues.csv",
                    "residue\n",
                )

                calls.append(
                    {
                        "structure_path": Path(
                            structure_path
                        ),
                        "output_dir": output,
                        "chain_id": chain_id,
                        "label": label,
                    }
                )

                return {
                    "report": {
                        "status": "complete",
                    },
                    "outputs": {
                        "json": str(json_path),
                        "csv": str(csv_path),
                    },
                }

            submitted_validation = {
                "report": {
                    "status": "complete",
                },
                "outputs": {
                    "json": str(
                        root
                        / "submitted_validation.json"
                    ),
                },
            }

            results = (
                _run_receptor_conformer_validations(
                    receptor_conformers=conformers,
                    output_dir=root / "output",
                    chain_id="A",
                    submitted_validation=(
                        submitted_validation
                    ),
                    validation_fn=fake_validation,
                )
            )

            self.assertEqual(
                set(results),
                {
                    "submitted_receptor",
                    "snapshot_0001",
                    "snapshot_0002",
                },
            )

            self.assertIs(
                results["submitted_receptor"],
                submitted_validation,
            )

            self.assertEqual(
                len(calls),
                2,
            )

            self.assertEqual(
                [
                    call["structure_path"]
                    for call in calls
                ],
                [
                    snapshot_1,
                    snapshot_2,
                ],
            )

            self.assertEqual(
                calls[0]["output_dir"],
                (
                    root
                    / "output"
                    / "structure_validation"
                    / "receptor_conformers"
                    / "snapshot_0001"
                ),
            )

            self.assertEqual(
                calls[1]["chain_id"],
                "A",
            )

            self.assertTrue(
                (
                    root
                    / "output"
                    / "structure_validation"
                    / "receptor_conformers"
                    / "snapshot_0002"
                    / "ramachandran_validation.json"
                ).is_file()
            )

    def test_failed_snapshot_is_retained(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)

            submitted = SimpleNamespace(
                source_pdb=_write(
                    root / "submitted.pdb"
                )
            )
            snapshot = SimpleNamespace(
                source_pdb=_write(
                    root / "snapshot.pdb"
                )
            )

            def fake_validation(**kwargs):
                return {
                    "report": {
                        "status": "failed",
                        "error": {
                            "message": "test failure",
                        },
                    },
                    "outputs": {},
                }

            results = (
                _run_receptor_conformer_validations(
                    receptor_conformers=[
                        (
                            "submitted_receptor",
                            submitted,
                        ),
                        (
                            "snapshot_0001",
                            snapshot,
                        ),
                    ],
                    output_dir=root / "output",
                    chain_id=None,
                    submitted_validation={
                        "report": {
                            "status": "complete",
                        },
                        "outputs": {},
                    },
                    validation_fn=fake_validation,
                )
            )

            self.assertEqual(
                results[
                    "snapshot_0001"
                ]["report"]["status"],
                "failed",
            )

    def test_rejects_duplicate_ids(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)

            receptor = SimpleNamespace(
                source_pdb=_write(
                    root / "snapshot.pdb"
                )
            )

            with self.assertRaisesRegex(
                ValueError,
                "Duplicate receptor conformer ID",
            ):
                _run_receptor_conformer_validations(
                    receptor_conformers=[
                        (
                            "snapshot_0001",
                            receptor,
                        ),
                        (
                            "snapshot_0001",
                            receptor,
                        ),
                    ],
                    output_dir=root / "output",
                    chain_id=None,
                    submitted_validation={
                        "report": {
                            "status": "complete",
                        },
                        "outputs": {},
                    },
                    validation_fn=(
                        lambda **kwargs: {}
                    ),
                )


if __name__ == "__main__":
    unittest.main()
