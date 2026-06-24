from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.models import (
    PocketDefinition,
)
from compoundrank.pocket_evidence import (
    load_pocket_evidence,
    score_pocket_biological_evidence,
)


def atom_line(
    *,
    serial: int,
    residue: str,
    chain: str,
    number: int,
) -> str:
    return (
        f"ATOM  {serial:5d}  CA  "
        f"{residue:>3s} {chain:1s}"
        f"{number:4d}    "
        f"{0.0:8.3f}"
        f"{0.0:8.3f}"
        f"{0.0:8.3f}"
        "  1.00  0.00           C\n"
    )


class PocketEvidenceTests(
    unittest.TestCase
):
    def test_merged_pocket_unions_component_residues(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pockets_dir = (
                root / "pockets"
            )
            pockets_dir.mkdir()

            (
                pockets_dir
                / "pocket4_atm.pdb"
            ).write_text(
                atom_line(
                    serial=1,
                    residue="ARG",
                    chain="B",
                    number=118,
                )
                + atom_line(
                    serial=2,
                    residue="GLU",
                    chain="B",
                    number=277,
                ),
                encoding="utf-8",
            )

            (
                pockets_dir
                / "pocket7_atm.pdb"
            ).write_text(
                atom_line(
                    serial=1,
                    residue="GLU",
                    chain="B",
                    number=277,
                )
                + atom_line(
                    serial=2,
                    residue="TYR",
                    chain="B",
                    number=406,
                ),
                encoding="utf-8",
            )

            pockets = [
                PocketDefinition(
                    mode="explicit",
                    pocket_id=(
                        "fpocket_04_pocket_4"
                    ),
                    pocket_rank=4,
                ),
                PocketDefinition(
                    mode="explicit",
                    pocket_id=(
                        "fpocket_07_pocket_7"
                    ),
                    pocket_rank=7,
                ),
                PocketDefinition(
                    mode="explicit",
                    pocket_id=(
                        "fpocket_merge_04_07_"
                        "pockets_4_7"
                    ),
                    pocket_rank=9,
                    merged_from=(
                        "fpocket_04_pocket_4",
                        "fpocket_07_pocket_7",
                    ),
                ),
            ]

            evidence = {
                "schema_version": (
                    "pocket_evidence.v0.1"
                ),
                "evidence_id": "test",
                "evidence_origin": (
                    "curated_annotation"
                ),
                "selection_mode": (
                    "prioritize_supported"
                ),
                "confidence": "high",
                "residues": [
                    "ARG:B:118",
                    "GLU:B:277",
                    "TYR:B:406",
                ],
            }

            scores, report = (
                score_pocket_biological_evidence(
                    pockets,
                    fpocket_output_dir=root,
                    evidence=evidence,
                )
            )

            merged = scores[
                "fpocket_merge_04_07_"
                "pockets_4_7"
            ]

            self.assertEqual(
                merged[
                    "biological_evidence_overlap_count"
                ],
                3,
            )

            self.assertAlmostEqual(
                merged[
                    "biological_evidence_recall"
                ],
                1.0,
            )

            self.assertEqual(
                merged[
                    "biological_evidence_rank"
                ],
                1,
            )

            self.assertTrue(
                report[
                    "used_for_selection"
                ]
            )

    def test_reference_ligand_posthoc_cannot_select(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "evidence.json"
            )

            path.write_text(
                json.dumps(
                    {
                        "schema_version": (
                            "pocket_evidence.v0.1"
                        ),
                        "evidence_origin": (
                            "reference_ligand_posthoc"
                        ),
                        "selection_mode": (
                            "prioritize_supported"
                        ),
                        "residues": [
                            "ARG:B:118",
                        ],
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "cannot be used",
            ):
                load_pocket_evidence(path)

    def test_nested_target_evidence_extension_is_supported(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "target_evidence.json"
            )

            path.write_text(
                json.dumps(
                    {
                        "target_interpretation": {
                            "target_class": (
                                "viral enzyme"
                            )
                        },
                        "pocket_residue_evidence": {
                            "schema_version": (
                                "pocket_evidence.v0.1"
                            ),
                            "evidence_id": (
                                "nested-test"
                            ),
                            "evidence_origin": (
                                "homolog_transfer"
                            ),
                            "selection_mode": (
                                "report_only"
                            ),
                            "residues": [
                                "asp:b:151",
                            ],
                        },
                    }
                ),
                encoding="utf-8",
            )

            evidence = (
                load_pocket_evidence(path)
            )

            self.assertEqual(
                evidence["residues"],
                ["ASP:B:151"],
            )

            self.assertEqual(
                evidence["selection_mode"],
                "report_only",
            )


if __name__ == "__main__":
    unittest.main()
