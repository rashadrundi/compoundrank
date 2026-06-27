from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from compoundrank.pipeline import (
    _run_selected_conformer_structure_pocket_quality,
)


class PipelineConformerStructurePocketQualityTests(
    unittest.TestCase
):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

        self.pockets = self.root / "pockets.json"
        self.selection = self.root / "selection.json"

        self.pockets.write_text(
            "{}\n",
            encoding="utf-8",
        )
        self.selection.write_text(
            "{}\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _receptor(
        self,
        conformer_id: str,
    ):
        path = (
            self.root
            / f"{conformer_id}.pdb"
        )
        path.write_text(
            "ATOM\n",
            encoding="utf-8",
        )
        return SimpleNamespace(
            source_pdb=path
        )

    def _validation(
        self,
        conformer_id: str,
        *,
        status: str = "complete",
    ) -> dict:
        path = (
            self.root
            / f"{conformer_id}_rama.json"
        )
        path.write_text(
            "{}\n",
            encoding="utf-8",
        )

        return {
            "report": {
                "status": status,
            },
            "outputs": {
                "json": str(path),
            },
        }

    def _fake_runner(
        self,
        calls: list[dict],
        verdicts: dict[str, str],
    ):
        def run(**kwargs):
            conformer_id = kwargs[
                "receptor_conformer_id"
            ]

            calls.append(kwargs)

            report = {
                "status": "complete",
                "receptor_conformer_id": (
                    conformer_id
                ),
                "selected_pocket_ids": [
                    f"{conformer_id}_pocket"
                ],
                "outlier_count": 2,
                "selected_box_local_outliers": [
                    "ALA:A:1"
                ],
                "selected_pose_local_outliers": [],
                "box_edge_only_outliers": [
                    "ALA:A:1"
                ],
                "selected_pose_available": True,
                "verdict": verdicts[
                    conformer_id
                ],
            }

            output_path = (
                Path(kwargs["output_dir"])
                / "structure_pocket_quality"
                / conformer_id
                / "structure_pocket_quality.json"
            )
            output_path.parent.mkdir(
                parents=True,
                exist_ok=True,
            )
            output_path.write_text(
                json.dumps(report) + "\n",
                encoding="utf-8",
            )

            return {
                "report": report,
                "output_path": output_path,
            }

        return run

    def test_evaluates_only_selected_conformer(
        self,
    ) -> None:
        calls: list[dict] = []

        pose_summary = (
            self.root
            / "pose_set_recovery_summary.json"
        )
        pose_summary.write_text(
            "{}\n",
            encoding="utf-8",
        )

        result = (
            _run_selected_conformer_structure_pocket_quality(
                selected_conformer_ids={
                    "snapshot_0002"
                },
                receptor_conformers=[
                    (
                        "submitted_receptor",
                        self._receptor(
                            "submitted_receptor"
                        ),
                    ),
                    (
                        "snapshot_0001",
                        self._receptor(
                            "snapshot_0001"
                        ),
                    ),
                    (
                        "snapshot_0002",
                        self._receptor(
                            "snapshot_0002"
                        ),
                    ),
                ],
                receptor_structure_validations={
                    "submitted_receptor": (
                        self._validation(
                            "submitted_receptor"
                        )
                    ),
                    "snapshot_0001": (
                        self._validation(
                            "snapshot_0001"
                        )
                    ),
                    "snapshot_0002": (
                        self._validation(
                            "snapshot_0002"
                        )
                    ),
                },
                pocket_definitions_path=(
                    self.pockets
                ),
                pocket_selection_summary_path=(
                    self.selection
                ),
                pose_recovery_summary_path=(
                    pose_summary
                ),
                output_dir=self.root / "output",
                quality_runner=self._fake_runner(
                    calls,
                    {
                        "snapshot_0002": (
                            "usable_with_global_"
                            "geometry_caution"
                        )
                    },
                ),
            )
        )

        self.assertEqual(
            len(calls),
            1,
        )
        self.assertEqual(
            calls[0][
                "receptor_conformer_id"
            ],
            "snapshot_0002",
        )
        self.assertEqual(
            Path(
                calls[0]["structure_path"]
            ).name,
            "snapshot_0002.pdb",
        )
        self.assertEqual(
            calls[0][
                "pose_recovery_summary_path"
            ],
            pose_summary,
        )
        self.assertEqual(
            result["report"][
                "selected_conformer_count"
            ],
            1,
        )

        conformer_record = (
            result["report"]["conformers"][0]
        )

        self.assertEqual(
            conformer_record[
                "selected_box_local_outlier_count"
            ],
            1,
        )
        self.assertEqual(
            conformer_record[
                "selected_pose_local_outlier_count"
            ],
            0,
        )
        self.assertEqual(
            conformer_record[
                "box_edge_only_outlier_count"
            ],
            1,
        )
        self.assertEqual(
            result["report"][
                "completed_conformer_count"
            ],
            1,
        )

    def test_aggregate_uses_worst_verdict(
        self,
    ) -> None:
        calls: list[dict] = []

        result = (
            _run_selected_conformer_structure_pocket_quality(
                selected_conformer_ids={
                    "snapshot_0001",
                    "snapshot_0002",
                },
                receptor_conformers=[
                    (
                        "snapshot_0001",
                        self._receptor(
                            "snapshot_0001"
                        ),
                    ),
                    (
                        "snapshot_0002",
                        self._receptor(
                            "snapshot_0002"
                        ),
                    ),
                ],
                receptor_structure_validations={
                    "snapshot_0001": (
                        self._validation(
                            "snapshot_0001"
                        )
                    ),
                    "snapshot_0002": (
                        self._validation(
                            "snapshot_0002"
                        )
                    ),
                },
                pocket_definitions_path=(
                    self.pockets
                ),
                pocket_selection_summary_path=(
                    self.selection
                ),
                output_dir=self.root / "output",
                quality_runner=self._fake_runner(
                    calls,
                    {
                        "snapshot_0001": "strong",
                        "snapshot_0002": (
                            "selected_pocket_"
                            "geometry_concern"
                        ),
                    },
                ),
            )
        )

        self.assertEqual(
            result["report"]["status"],
            "complete",
        )
        self.assertEqual(
            result["report"][
                "overall_verdict"
            ],
            (
                "selected_pocket_"
                "geometry_concern"
            ),
        )
        self.assertEqual(
            len(calls),
            2,
        )

    def test_incomplete_validation_is_recorded(
        self,
    ) -> None:
        calls: list[dict] = []

        result = (
            _run_selected_conformer_structure_pocket_quality(
                selected_conformer_ids={
                    "snapshot_0001"
                },
                receptor_conformers=[
                    (
                        "snapshot_0001",
                        self._receptor(
                            "snapshot_0001"
                        ),
                    ),
                ],
                receptor_structure_validations={
                    "snapshot_0001": (
                        self._validation(
                            "snapshot_0001",
                            status="failed",
                        )
                    ),
                },
                pocket_definitions_path=(
                    self.pockets
                ),
                pocket_selection_summary_path=(
                    self.selection
                ),
                output_dir=self.root / "output",
                quality_runner=self._fake_runner(
                    calls,
                    {
                        "snapshot_0001": "strong"
                    },
                ),
            )
        )

        self.assertEqual(calls, [])
        self.assertEqual(
            result["report"]["status"],
            "failed",
        )
        self.assertEqual(
            result["report"][
                "overall_verdict"
            ],
            (
                "manual_review_incomplete_"
                "conformer_quality"
            ),
        )
        self.assertEqual(
            result["report"]["conformers"][0][
                "reason_code"
            ],
            (
                "ramachandran_validation_"
                "incomplete"
            ),
        )

    def test_single_receptor_writes_legacy_alias(
        self,
    ) -> None:
        calls: list[dict] = []
        output_dir = self.root / "output"

        _run_selected_conformer_structure_pocket_quality(
            selected_conformer_ids={
                "submitted_receptor"
            },
            receptor_conformers=[
                (
                    "submitted_receptor",
                    self._receptor(
                        "submitted_receptor"
                    ),
                ),
            ],
            receptor_structure_validations={
                "submitted_receptor": (
                    self._validation(
                        "submitted_receptor"
                    )
                ),
            },
            pocket_definitions_path=self.pockets,
            pocket_selection_summary_path=(
                self.selection
            ),
            output_dir=output_dir,
            quality_runner=self._fake_runner(
                calls,
                {
                    "submitted_receptor": "strong"
                },
            ),
        )

        legacy_path = (
            output_dir
            / "structure_pocket_quality.json"
        )

        self.assertTrue(
            legacy_path.is_file()
        )

        saved = json.loads(
            legacy_path.read_text(
                encoding="utf-8"
            )
        )

        self.assertEqual(
            saved["receptor_conformer_id"],
            "submitted_receptor",
        )


if __name__ == "__main__":
    unittest.main()
