from __future__ import annotations

import unittest
from pathlib import Path

from compoundrank.openmm_md import (
    OpenMMMDConfig,
    build_pdb2pqr_command,
    choose_platform_name,
    validate_config,
)


class OpenMMMDTests(
    unittest.TestCase
):
    def test_auto_platform_prefers_cuda(
        self,
    ) -> None:
        result = choose_platform_name(
            [
                "Reference",
                "CPU",
                "CUDA",
            ],
            "auto",
        )

        self.assertEqual(
            result,
            "CUDA",
        )

    def test_auto_platform_falls_back_to_cpu(
        self,
    ) -> None:
        result = choose_platform_name(
            [
                "Reference",
                "CPU",
            ],
            "auto",
        )

        self.assertEqual(
            result,
            "CPU",
        )

    def test_missing_requested_platform_fails(
        self,
    ) -> None:
        with self.assertRaises(
            RuntimeError
        ):
            choose_platform_name(
                [
                    "Reference",
                    "CPU",
                ],
                "CUDA",
            )

    def test_snapshot_interval_must_divide_production(
        self,
    ) -> None:
        config = OpenMMMDConfig(
            production_steps=1001,
            snapshot_interval=500,
        )

        with self.assertRaises(
            ValueError
        ):
            validate_config(
                config
            )

    def test_builds_pdb2pqr_command(
        self,
    ) -> None:
        command = build_pdb2pqr_command(
            pdb2pqr_bin=(
                "/tools/pdb2pqr"
            ),
            receptor_pdb=Path(
                "/data/receptor.pdb"
            ),
            output_pdb=Path(
                "/output/prepared.pdb"
            ),
            output_pqr=Path(
                "/output/prepared.pqr"
            ),
            ph=7.4,
        )

        self.assertEqual(
            command[0],
            "/tools/pdb2pqr",
        )

        self.assertIn(
            "--ff=AMBER",
            command,
        )

        self.assertIn(
            "--keep-chain",
            command,
        )

        self.assertIn(
            "--with-ph=7.4",
            command,
        )

        self.assertEqual(
            command[-2:],
            [
                "/data/receptor.pdb",
                "/output/prepared.pqr",
            ],
        )


if __name__ == "__main__":
    unittest.main()
