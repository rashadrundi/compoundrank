from __future__ import annotations

import csv
import inspect
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from compoundrank.pipeline import (
    _accepted_records_for_selected_pocket,
    _preserve_posebusters_artifacts,
    _write_docking_attempt_summary,
    run_pipeline,
)
from compoundrank.pocket_selection import (
    rank_pocket_attempts,
    summarize_pocket_attempt,
    write_pocket_selection_summary,
)


def _pose(
    score: float,
    conformer_id: str,
):
    return SimpleNamespace(
        cnn_score=score,
        cnn_affinity=7.0,
        minimized_affinity=-8.0,
        seed=1,
        pose_number=1,
        receptor_conformer_id=(
            conformer_id
        ),
    )


def _pocket():
    return SimpleNamespace(
        pocket_id="pocket_01",
        pocket_rank=1,
        fpocket_score=42.0,
        source="fpocket",
        mode="fpocket",
    )


class MultiConformerExecutionTests(
    unittest.TestCase
):
    def test_attempt_summary_records_conformer(
        self,
    ) -> None:
        row = summarize_pocket_attempt(
            ligand_name="ligand",
            receptor_conformer_id=(
                "snapshot_0001"
            ),
            pocket=_pocket(),
            raw_records=[
                _pose(
                    0.70,
                    "snapshot_0001",
                )
            ],
            accepted_records=[
                _pose(
                    0.70,
                    "snapshot_0001",
                )
            ],
            rejected_pose_count=0,
        )

        self.assertEqual(
            row[
                "receptor_conformer_id"
            ],
            "snapshot_0001",
        )

    def test_ranking_compares_conformer_pocket_pairs(
        self,
    ) -> None:
        submitted = (
            summarize_pocket_attempt(
                ligand_name="ligand",
                receptor_conformer_id=(
                    "submitted_receptor"
                ),
                pocket=_pocket(),
                raw_records=[
                    _pose(
                        0.55,
                        "submitted_receptor",
                    )
                ],
                accepted_records=[
                    _pose(
                        0.55,
                        "submitted_receptor",
                    )
                ],
                rejected_pose_count=0,
            )
        )

        snapshot = (
            summarize_pocket_attempt(
                ligand_name="ligand",
                receptor_conformer_id=(
                    "snapshot_0001"
                ),
                pocket=_pocket(),
                raw_records=[
                    _pose(
                        0.82,
                        "snapshot_0001",
                    )
                ],
                accepted_records=[
                    _pose(
                        0.82,
                        "snapshot_0001",
                    )
                ],
                rejected_pose_count=0,
            )
        )

        ranked = rank_pocket_attempts(
            [
                submitted,
                snapshot,
            ]
        )

        self.assertEqual(
            ranked[0][
                "receptor_conformer_id"
            ],
            "snapshot_0001",
        )

        self.assertTrue(
            ranked[0]["selected"]
        )

    def test_selected_records_use_pair_key(
        self,
    ) -> None:
        submitted_record = object()
        snapshot_record = object()

        records = {
            (
                "submitted_receptor",
                "pocket_01",
            ): [
                submitted_record
            ],
            (
                "snapshot_0001",
                "pocket_01",
            ): [
                snapshot_record
            ],
        }

        selected = (
            _accepted_records_for_selected_pocket(
                valid_records_by_pocket=(
                    records
                ),
                selected_pocket_id=(
                    "pocket_01"
                ),
                selected_receptor_conformer_id=(
                    "snapshot_0001"
                ),
            )
        )

        self.assertEqual(
            selected,
            [
                snapshot_record
            ],
        )

    def test_posebusters_artifacts_are_separated(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            validity = (
                root / "validity"
            )

            validity.mkdir()

            for filename in (
                "posebusters_input.sdf",
                "posebusters_report.csv",
            ):
                (
                    validity
                    / filename
                ).write_text(
                    filename + "\n",
                    encoding="utf-8",
                )

            copied = (
                _preserve_posebusters_artifacts(
                    validity_dir=validity,
                    output_dir=(
                        root / "output"
                    ),
                    ligand_name="ligand",
                    pocket_id="pocket_01",
                    receptor_conformer_id=(
                        "snapshot_0001"
                    ),
                )
            )

            self.assertEqual(
                len(copied),
                2,
            )

            for copied_path in copied:
                self.assertIn(
                    "snapshot_0001",
                    copied_path.parts,
                )

    def test_reports_include_conformer_metadata(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)

            attempt_path = (
                _write_docking_attempt_summary(
                    output,
                    [
                        {
                            "compound": "ligand",
                            (
                                "receptor_"
                                "conformer_id"
                            ): "snapshot_0001",
                            "pocket": "pocket_01",
                            "raw_poses": 1,
                            "accepted_poses": 1,
                            "rejected_poses": 0,
                            "status": "accepted",
                            (
                                "best_raw_"
                                "cnn_score"
                            ): "0.8",
                            (
                                "best_accepted_"
                                "cnn_score"
                            ): "0.8",
                        }
                    ],
                )
            )

            self.assertIsNotNone(
                attempt_path
            )

            with Path(
                attempt_path
            ).open(
                newline="",
                encoding="utf-8",
            ) as handle:
                rows = list(
                    csv.DictReader(
                        handle
                    )
                )

            self.assertEqual(
                rows[0][
                    "receptor_conformer_id"
                ],
                "snapshot_0001",
            )

            selection_paths = (
                write_pocket_selection_summary(
                    output,
                    [
                        {
                            "compound": "ligand",
                            (
                                "receptor_"
                                "conformer_id"
                            ): "snapshot_0001",
                            "selection_rank": 1,
                            "selected": True,
                            "pocket_id": "pocket_01",
                        }
                    ],
                )
            )

            self.assertIsNotNone(
                selection_paths
            )

            selection_csv, selection_json = (
                selection_paths
            )

            with selection_csv.open(
                newline="",
                encoding="utf-8",
            ) as handle:
                selection_rows = list(
                    csv.DictReader(
                        handle
                    )
                )

            self.assertEqual(
                selection_rows[0][
                    "receptor_conformer_id"
                ],
                "snapshot_0001",
            )

            payload = json.loads(
                selection_json.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                payload[
                    "receptor_conformer_count"
                ],
                1,
            )

            self.assertEqual(
                payload["selection_unit"],
                (
                    "receptor_conformer_"
                    "pocket_pair"
                ),
            )

    def test_pipeline_source_executes_each_conformer(
        self,
    ) -> None:
        source = inspect.getsource(
            run_pipeline
        )

        self.assertIn(
            "in receptor_conformers:",
            source,
        )

        self.assertIn(
            "receptor_conformer_id=(",
            source,
        )

        self.assertIn(
            "conformer_receptor",
            source,
        )

        self.assertIn(
            ".display_pdb",
            source,
        )

        self.assertIn(
            (
                "selected_receptor_"
                "conformer_id"
            ),
            source,
        )


if __name__ == "__main__":
    unittest.main()
