from __future__ import annotations

import unittest

from compoundrank.cli import build_parser


class CliPocketEvidenceTests(
    unittest.TestCase
):
    def test_parser_accepts_pocket_evidence_json(
        self,
    ) -> None:
        arguments = (
            build_parser().parse_args(
                [
                    "--receptor",
                    "/tmp/receptor.pdb",
                    "--data-root",
                    "/tmp/data",
                    "--pocket-evidence-json",
                    "/tmp/evidence.json",
                ]
            )
        )

        self.assertEqual(
            arguments.pocket_evidence_json,
            "/tmp/evidence.json",
        )


if __name__ == "__main__":
    unittest.main()
