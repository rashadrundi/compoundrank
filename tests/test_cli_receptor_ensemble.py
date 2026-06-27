from __future__ import annotations

import inspect
import unittest

from compoundrank.cli import (
    build_parser,
)
from compoundrank.pipeline import (
    run_pipeline,
)


class CliReceptorEnsembleTests(
    unittest.TestCase
):
    def test_parser_accepts_manifest(
        self,
    ) -> None:
        arguments = (
            build_parser().parse_args(
                [
                    "--receptor",
                    "/tmp/receptor.pdb",
                    "--data-root",
                    "/tmp/data",
                    (
                        "--receptor-"
                        "ensemble-json"
                    ),
                    "/tmp/ensemble.json",
                ]
            )
        )

        self.assertEqual(
            arguments
            .receptor_ensemble_json,
            "/tmp/ensemble.json",
        )

    def test_pipeline_signature_accepts_manifest(
        self,
    ) -> None:
        parameter = (
            inspect.signature(
                run_pipeline
            )
            .parameters[
                "receptor_ensemble_json"
            ]
        )

        self.assertIsNone(
            parameter.default
        )


if __name__ == "__main__":
    unittest.main()
