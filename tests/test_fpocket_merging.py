from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.cli import build_parser
from compoundrank.models import PocketDefinition
from compoundrank.pocket import (
    _merge_nearby_fpocket_definitions,
    write_pocket_definitions,
)


def atom_line(
    serial: int,
    x: float,
    y: float,
    z: float,
) -> str:
    return (
        f"ATOM  {serial:5d}  C   UNK A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        "  1.00  0.00           C\n"
    )


class FpocketMergingTests(unittest.TestCase):
    def test_nearby_pair_is_appended_and_far_pair_is_not(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            out_dir = Path(directory)
            pockets_dir = out_dir / "pockets"
            pockets_dir.mkdir()

            (
                pockets_dir
                / "pocket4_vert.pqr"
            ).write_text(
                atom_line(
                    1,
                    0.0,
                    0.0,
                    0.0,
                )
                + atom_line(
                    2,
                    1.0,
                    0.0,
                    0.0,
                ),
                encoding="utf-8",
            )

            (
                pockets_dir
                / "pocket7_vert.pqr"
            ).write_text(
                atom_line(
                    1,
                    4.0,
                    0.0,
                    0.0,
                )
                + atom_line(
                    2,
                    5.0,
                    0.0,
                    0.0,
                ),
                encoding="utf-8",
            )

            (
                pockets_dir
                / "pocket9_vert.pqr"
            ).write_text(
                atom_line(
                    1,
                    30.0,
                    0.0,
                    0.0,
                ),
                encoding="utf-8",
            )

            selected = [
                {
                    "number": 4,
                    "metrics": {
                        "Score": 0.8,
                    },
                },
                {
                    "number": 7,
                    "metrics": {
                        "Score": 0.7,
                    },
                },
                {
                    "number": 9,
                    "metrics": {
                        "Score": 0.6,
                    },
                },
            ]

            independent = [
                PocketDefinition(
                    mode="explicit",
                    pocket_id=(
                        "fpocket_01_pocket_4"
                    ),
                    pocket_rank=1,
                ),
                PocketDefinition(
                    mode="explicit",
                    pocket_id=(
                        "fpocket_02_pocket_7"
                    ),
                    pocket_rank=2,
                ),
                PocketDefinition(
                    mode="explicit",
                    pocket_id=(
                        "fpocket_03_pocket_9"
                    ),
                    pocket_rank=3,
                ),
            ]

            merged = (
                _merge_nearby_fpocket_definitions(
                    selected_pockets=selected,
                    independent_definitions=(
                        independent
                    ),
                    out_dir=out_dir,
                    padding=4.0,
                    distance_threshold=4.0,
                    starting_rank=4,
                )
            )

            self.assertEqual(
                len(merged),
                1,
            )

            candidate = merged[0]

            self.assertEqual(
                candidate.pocket_id,
                (
                    "fpocket_merge_01_02_"
                    "pockets_4_7"
                ),
            )

            self.assertEqual(
                candidate.pocket_rank,
                4,
            )

            self.assertEqual(
                candidate.merged_from,
                (
                    "fpocket_01_pocket_4",
                    "fpocket_02_pocket_7",
                ),
            )

            self.assertAlmostEqual(
                candidate.merge_distance,
                3.0,
            )

            self.assertAlmostEqual(
                candidate.center_x,
                2.5,
            )

            self.assertAlmostEqual(
                candidate.size_x,
                20.0,
            )

            self.assertIsNone(
                candidate.fpocket_score
            )

            output = (
                out_dir
                / "pocket_definitions.json"
            )

            write_pocket_definitions(
                output,
                independent + merged,
            )

            payload = json.loads(
                output.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                payload[
                    "independent_pocket_count"
                ],
                3,
            )

            self.assertEqual(
                payload[
                    "merged_pocket_count"
                ],
                1,
            )

            self.assertIn(
                "merged candidates appended",
                payload["ranking_method"],
            )

            self.assertEqual(
                payload["pockets"][3][
                    "merged_from"
                ],
                [
                    "fpocket_01_pocket_4",
                    "fpocket_02_pocket_7",
                ],
            )

    def test_cli_accepts_merge_options(
        self,
    ) -> None:
        arguments = (
            build_parser().parse_args(
                [
                    "--receptor",
                    "/tmp/receptor.pdb",
                    "--data-root",
                    "/tmp/data",
                    "--fpocket-top-n",
                    "7",
                    "--fpocket-merge-nearby",
                    "--fpocket-merge-distance",
                    "4.5",
                ]
            )
        )

        self.assertEqual(
            arguments.fpocket_top_n,
            7,
        )

        self.assertTrue(
            arguments.fpocket_merge_nearby
        )

        self.assertAlmostEqual(
            arguments.fpocket_merge_distance,
            4.5,
        )


if __name__ == "__main__":
    unittest.main()
