from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from compoundrank.functional_site_transfer import (
    build_homolog_pocket_evidence,
    write_homolog_pocket_evidence,
)
from compoundrank.pocket_evidence import (
    load_pocket_evidence,
)


AA1_TO_3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "N": "ASN",
    "R": "ARG",
    "Y": "TYR",
}


def atom_line(
    *,
    serial: int,
    amino_acid: str,
    chain: str,
    residue_number: int,
) -> str:
    return (
        f"ATOM  {serial:5d}  CA  "
        f"{AA1_TO_3[amino_acid]:>3s} "
        f"{chain:1s}"
        f"{residue_number:4d}    "
        f"{0.0:8.3f}"
        f"{0.0:8.3f}"
        f"{0.0:8.3f}"
        "  1.00  0.00           C\n"
    )


def write_chain(
    path: Path,
    *,
    sequence: str,
    chain: str,
    start_number: int,
) -> None:
    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for index, amino_acid in enumerate(
            sequence,
            start=0,
        ):
            handle.write(
                atom_line(
                    serial=index + 1,
                    amino_acid=(
                        amino_acid
                    ),
                    chain=chain,
                    residue_number=(
                        start_number
                        + index
                    ),
                )
            )


def reference_record(
    *,
    selection_mode: str = (
        "prioritize_supported"
    ),
) -> dict:
    return {
        "schema_version": (
            "functional_site_reference.v0.1"
        ),
        "evidence_id": (
            "example_functional_site"
        ),
        "selection_mode": (
            selection_mode
        ),
        "confidence": "high",
        "reference_chain": "R",
        "residues": [
            {
                "residue_number": 85,
                "amino_acid": "D",
                "label": "catalytic",
            },
            {
                "residue_number": 87,
                "amino_acid": "F",
                "label": "binding",
            },
        ],
        "source": {
            "database": "test",
        },
    }


class FunctionalSiteTransferTests(
    unittest.TestCase
):
    def test_strong_transfer_can_prioritize(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            reference_pdb = (
                root / "reference.pdb"
            )

            receptor_pdb = (
                root / "receptor.pdb"
            )

            write_chain(
                reference_pdb,
                sequence="ACDEFG",
                chain="R",
                start_number=83,
            )

            write_chain(
                receptor_pdb,
                sequence="ACDEFG",
                chain="B",
                start_number=101,
            )

            evidence = (
                build_homolog_pocket_evidence(
                    reference_record=(
                        reference_record()
                    ),
                    reference_sequence=(
                        "ACDEFG"
                    ),
                    submitted_sequence=(
                        "ACDEFG"
                    ),
                    reference_pdb=(
                        reference_pdb
                    ),
                    receptor_pdb=(
                        receptor_pdb
                    ),
                    receptor_chain_id="B",
                )
            )

            self.assertEqual(
                evidence[
                    "selection_mode"
                ],
                "prioritize_supported",
            )

            self.assertEqual(
                evidence[
                    "evidence_origin"
                ],
                "homolog_transfer",
            )

            self.assertEqual(
                evidence["residues"],
                [
                    "ASP:B:103",
                    "PHE:B:105",
                ],
            )

            self.assertEqual(
                evidence[
                    "transfer_summary"
                ][
                    "mapped_conserved_count"
                ],
                2,
            )

    def test_weak_homolog_is_report_only(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            reference_pdb = (
                root / "reference.pdb"
            )

            receptor_pdb = (
                root / "receptor.pdb"
            )

            write_chain(
                reference_pdb,
                sequence="ACDEFG",
                chain="R",
                start_number=83,
            )

            write_chain(
                receptor_pdb,
                sequence="ACNEYG",
                chain="B",
                start_number=101,
            )

            evidence = (
                build_homolog_pocket_evidence(
                    reference_record=(
                        reference_record()
                    ),
                    reference_sequence=(
                        "ACDEFG"
                    ),
                    submitted_sequence=(
                        "ACNEYG"
                    ),
                    reference_pdb=(
                        reference_pdb
                    ),
                    receptor_pdb=(
                        receptor_pdb
                    ),
                    receptor_chain_id="B",
                )
            )

            self.assertEqual(
                evidence[
                    "selection_mode"
                ],
                "report_only",
            )

            self.assertFalse(
                evidence[
                    "transfer_summary"
                ][
                    "selection_checks"
                ][
                    "homolog_alignment_supported"
                ]
            )

    def test_written_output_loads_as_pocket_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            reference_pdb = (
                root / "reference.pdb"
            )

            receptor_pdb = (
                root / "receptor.pdb"
            )

            output_path = (
                root
                / "pocket_evidence.json"
            )

            write_chain(
                reference_pdb,
                sequence="ACDEFG",
                chain="R",
                start_number=83,
            )

            write_chain(
                receptor_pdb,
                sequence="ACDEFG",
                chain="B",
                start_number=101,
            )

            write_homolog_pocket_evidence(
                output_path,
                reference_record=(
                    reference_record()
                ),
                reference_sequence=(
                    "ACDEFG"
                ),
                submitted_sequence=(
                    "ACDEFG"
                ),
                reference_pdb=(
                    reference_pdb
                ),
                receptor_pdb=(
                    receptor_pdb
                ),
                receptor_chain_id="B",
            )

            loaded = (
                load_pocket_evidence(
                    output_path
                )
            )

            self.assertEqual(
                loaded[
                    "evidence_origin"
                ],
                "homolog_transfer",
            )

            self.assertEqual(
                loaded[
                    "selection_mode"
                ],
                "prioritize_supported",
            )

            self.assertEqual(
                loaded["residues"],
                [
                    "ASP:B:103",
                    "PHE:B:105",
                ],
            )


if __name__ == "__main__":
    unittest.main()
