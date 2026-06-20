from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.models import (
    PocketDefinition,
    PoseRecord,
)
from compoundrank.pocket_selection import (
    rank_pocket_attempts,
    summarize_pocket_attempt,
    write_pocket_selection_summary,
)


def pose(
    *,
    pocket_id: str,
    score: float,
    affinity: float,
    energy: float,
    pose_number: int = 1,
) -> PoseRecord:
    return PoseRecord(
        ligand_name="darunavir",
        seed=61453,
        pose_number=pose_number,
        molecule=None,
        cnn_score=score,
        cnn_affinity=affinity,
        minimized_affinity=energy,
        source_sdf=Path("poses.sdf"),
        pocket_id=pocket_id,
    )


class PocketSelectionTests(unittest.TestCase):
    def test_accepted_poses_are_preferred(
        self,
    ) -> None:
        pocket = PocketDefinition(
            mode="explicit",
            pocket_id="pocket_1",
            pocket_rank=1,
        )

        row = summarize_pocket_attempt(
            ligand_name="darunavir",
            pocket=pocket,
            raw_records=[
                pose(
                    pocket_id="pocket_1",
                    score=0.99,
                    affinity=10.0,
                    energy=-10.0,
                )
            ],
            accepted_records=[
                pose(
                    pocket_id="pocket_1",
                    score=0.80,
                    affinity=9.0,
                    energy=-9.0,
                )
            ],
            rejected_pose_count=1,
        )

        self.assertEqual(
            row["score_source"],
            "accepted_poses",
        )
        self.assertAlmostEqual(
            row["top_cnn_score"],
            0.80,
        )

    def test_raw_fallback_is_recorded(
        self,
    ) -> None:
        pocket = PocketDefinition(
            mode="explicit",
            pocket_id="pocket_2",
            pocket_rank=2,
        )

        row = summarize_pocket_attempt(
            ligand_name="darunavir",
            pocket=pocket,
            raw_records=[
                pose(
                    pocket_id="pocket_2",
                    score=0.60,
                    affinity=7.0,
                    energy=-2.0,
                )
            ],
            accepted_records=[],
            rejected_pose_count=1,
        )

        self.assertEqual(
            row["score_source"],
            "raw_pose_fallback",
        )
        self.assertAlmostEqual(
            row["top_cnn_score"],
            0.60,
        )

    def test_highest_cnn_score_selects_pocket(
        self,
    ) -> None:
        rows = [
            {
                "compound": "darunavir",
                "pocket_id": "pocket_1",
                "pocket_rank": 1,
                "top_cnn_score": 0.61,
                "top_cnn_affinity": 4.9,
                "top_minimized_affinity": 0.6,
            },
            {
                "compound": "darunavir",
                "pocket_id": "pocket_3",
                "pocket_rank": 3,
                "top_cnn_score": 0.98,
                "top_cnn_affinity": 10.3,
                "top_minimized_affinity": -9.4,
            },
        ]

        ranked = rank_pocket_attempts(rows)

        self.assertEqual(
            ranked[0]["pocket_id"],
            "pocket_3",
        )
        self.assertTrue(ranked[0]["selected"])
        self.assertEqual(
            ranked[0]["selection_rank"],
            1,
        )
        self.assertFalse(ranked[1]["selected"])

    def test_summary_writer_creates_csv_and_json(
        self,
    ) -> None:
        rows = rank_pocket_attempts(
            [
                {
                    "compound": "darunavir",
                    "pocket_id": "pocket_3",
                    "pocket_rank": 3,
                    "fpocket_score": 0.126,
                    "pocket_source": "fpocket",
                    "raw_poses": 20,
                    "accepted_poses": 20,
                    "rejected_poses": 0,
                    "score_source": "accepted_poses",
                    "top_cnn_score": 0.9818,
                    "top_cnn_affinity": 10.363,
                    "top_minimized_affinity": -9.428,
                    "top_seed": 61453,
                    "top_pose_number": 1,
                }
            ]
        )

        with tempfile.TemporaryDirectory() as directory:
            result = write_pocket_selection_summary(
                Path(directory),
                rows,
            )

            self.assertIsNotNone(result)
            csv_path, json_path = result

            self.assertTrue(csv_path.is_file())
            self.assertTrue(json_path.is_file())

            payload = json.loads(
                json_path.read_text(
                    encoding="utf-8"
                )
            )

            self.assertFalse(
                payload[
                    "reference_ligand_used_for_selection"
                ]
            )
            self.assertEqual(
                payload["selected_pockets"][0][
                    "pocket_id"
                ],
                "pocket_3",
            )


if __name__ == "__main__":
    unittest.main()
