from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.pipeline import (
    _run_structure_validation,
)


TEST_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.460   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.450   1.100   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       2.200   2.300   0.000  1.00 20.00           O
ATOM      5  N   ALA A   2       3.700   0.800   0.200  1.00 20.00           N
ATOM      6  CA  ALA A   2       4.600   1.900   0.400  1.00 20.00           C
ATOM      7  C   ALA A   2       6.000   1.500   0.600  1.00 20.00           C
ATOM      8  O   ALA A   2       6.300   0.300   0.600  1.00 20.00           O
ATOM      9  N   GLY A   3       6.900   2.400   0.800  1.00 20.00           N
ATOM     10  CA  GLY A   3       8.300   2.000   1.000  1.00 20.00           C
ATOM     11  C   GLY A   3       9.100   3.200   1.500  1.00 20.00           C
ATOM     12  O   GLY A   3       8.700   4.300   1.700  1.00 20.00           O
TER
END
"""


class PipelineRamachandranIntegrationTests(
    unittest.TestCase
):
    def test_writes_receptor_report(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            receptor = root / "receptor.pdb"
            output = root / "validation"

            receptor.write_text(
                TEST_PDB,
                encoding="utf-8",
            )

            result = _run_structure_validation(
                structure_path=receptor,
                output_dir=output,
                chain_id="A",
                label="submitted receptor",
            )

            report_path = Path(
                result["outputs"]["json"]
            )

            report = json.loads(
                report_path.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                report["status"],
                "complete",
            )

            self.assertEqual(
                report["selection_mode"],
                "report_only",
            )

    def test_failure_is_nonblocking(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            receptor = root / "receptor.pdb"

            receptor.write_text(
                TEST_PDB,
                encoding="utf-8",
            )

            result = _run_structure_validation(
                structure_path=receptor,
                output_dir=root / "validation",
                chain_id="Z",
                label="submitted receptor",
            )

            self.assertEqual(
                result["report"]["status"],
                "failed",
            )

            self.assertTrue(
                Path(
                    result["outputs"]["json"]
                ).is_file()
            )


if __name__ == "__main__":
    unittest.main()
