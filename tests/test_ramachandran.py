from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.ramachandran import (
    GRID_FILES,
    classify_ramachandran,
    classify_score,
    load_top8000_grids,
    lookup_top8000_score,
    validate_ramachandran,
    write_ramachandran_outputs,
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


class RamachandranTests(
    unittest.TestCase
):
    def test_loads_all_top8000_grids(
        self,
    ) -> None:
        grids = load_top8000_grids()

        self.assertEqual(
            set(grids),
            set(GRID_FILES),
        )

        for category in GRID_FILES:
            self.assertGreater(
                len(grids[category]),
                6000,
            )

    def test_official_score_cutoffs(
        self,
    ) -> None:
        self.assertEqual(
            classify_score(
                "general",
                0.02,
            ),
            "favored",
        )

        self.assertEqual(
            classify_score(
                "general",
                0.0005,
            ),
            "allowed",
        )

        self.assertEqual(
            classify_score(
                "general",
                0.00049,
            ),
            "outlier",
        )

        self.assertEqual(
            classify_score(
                "cis_proline",
                0.0020,
            ),
            "allowed",
        )

        self.assertEqual(
            classify_score(
                "trans_proline",
                0.0010,
            ),
            "allowed",
        )

    def test_classifies_top8000_peaks(
        self,
    ) -> None:
        favored_points = (
            (
                "general",
                -63.0,
                -43.0,
            ),
            (
                "glycine",
                63.0,
                41.0,
            ),
            (
                "cis_proline",
                -75.0,
                155.0,
            ),
            (
                "trans_proline",
                -59.0,
                143.0,
            ),
            (
                "pre_proline",
                -57.0,
                -45.0,
            ),
            (
                "isoleucine_valine",
                -121.0,
                129.0,
            ),
        )

        for (
            category,
            phi,
            psi,
        ) in favored_points:
            with self.subTest(
                category=category
            ):
                result = (
                    classify_ramachandran(
                        phi,
                        psi,
                        category,
                    )
                )

                self.assertEqual(
                    result[
                        "classification"
                    ],
                    "favored",
                )

                self.assertGreaterEqual(
                    result["score"],
                    0.02,
                )

    def test_sparse_missing_bin_is_zero(
        self,
    ) -> None:
        result = (
            lookup_top8000_score(
                179.0,
                -1.0,
                "cis_proline",
            )
        )

        self.assertFalse(
            result[
                "grid_point_present"
            ]
        )

        self.assertEqual(
            result["score"],
            0.0,
        )

    def test_validates_connected_backbone(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(
                temporary
            )

            structure = (
                root / "peptide.pdb"
            )

            structure.write_text(
                TEST_PDB,
                encoding="utf-8",
            )

            report = (
                validate_ramachandran(
                    structure
                )
            )

            self.assertEqual(
                report[
                    "schema_version"
                ],
                (
                    "ramachandran_"
                    "validation.v0.2"
                ),
            )

            self.assertEqual(
                report["status"],
                "complete",
            )

            self.assertEqual(
                report[
                    "selection_mode"
                ],
                "report_only",
            )

            self.assertEqual(
                report[
                    "total_polymer_residues"
                ],
                3,
            )

            self.assertEqual(
                report[
                    "evaluable_residues"
                ],
                1,
            )

            self.assertEqual(
                report["residues"][0][
                    "residue"
                ],
                "ALA:A:2",
            )

            self.assertIn(
                report["residues"][0][
                    "classification"
                ],
                {
                    "favored",
                    "allowed",
                    "outlier",
                },
            )

    def test_writes_json_and_csv(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(
                temporary
            )

            structure = (
                root / "peptide.pdb"
            )

            output = (
                root / "output"
            )

            structure.write_text(
                TEST_PDB,
                encoding="utf-8",
            )

            report = (
                validate_ramachandran(
                    structure
                )
            )

            paths = (
                write_ramachandran_outputs(
                    report,
                    output,
                )
            )

            json_path = Path(
                paths["json"]
            )

            csv_path = Path(
                paths["csv"]
            )

            self.assertTrue(
                json_path.is_file()
            )

            self.assertTrue(
                csv_path.is_file()
            )

            loaded = json.loads(
                json_path.read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                loaded[
                    "evaluable_residues"
                ],
                1,
            )

            self.assertIn(
                "score_percent",
                csv_path.read_text(
                    encoding="utf-8"
                ),
            )


if __name__ == "__main__":
    unittest.main()
