from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.functional_site_transfer import (
    main,
    read_single_fasta,
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
}


def write_fasta(
    path: Path,
    *,
    header: str,
    sequence: str,
) -> None:
    path.write_text(
        f">{header}\n{sequence}\n",
        encoding="utf-8",
    )


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


class FunctionalSiteTransferCliTests(
    unittest.TestCase
):
    def test_read_single_fasta_rejects_multiple_records(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "multiple.faa"
            )

            path.write_text(
                ">one\nACDE\n>two\nFG\n",
                encoding="utf-8",
            )

            with self.assertRaises(
                ValueError
            ):
                read_single_fasta(path)

    def test_cli_writes_loadable_pocket_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            record_path = (
                root / "reference.json"
            )

            reference_fasta = (
                root / "reference.faa"
            )

            submitted_fasta = (
                root / "submitted.faa"
            )

            reference_pdb = (
                root / "reference.pdb"
            )

            receptor_pdb = (
                root / "receptor.pdb"
            )

            output_path = (
                root / "pocket_evidence.json"
            )

            record_path.write_text(
                json.dumps(
                    {
                        "schema_version": (
                            "functional_site_reference.v0.1"
                        ),
                        "evidence_id": (
                            "cli_test_site"
                        ),
                        "selection_mode": (
                            "prioritize_supported"
                        ),
                        "confidence": "high",
                        "reference_chain": "R",
                        "residues": [
                            {
                                "residue_number": 85,
                                "amino_acid": "D",
                            },
                            {
                                "residue_number": 87,
                                "amino_acid": "F",
                            },
                        ],
                    }
                ),
                encoding="utf-8",
            )

            write_fasta(
                reference_fasta,
                header="reference",
                sequence="ACDEFG",
            )

            write_fasta(
                submitted_fasta,
                header="submitted",
                sequence="ACDEFG",
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

            exit_code = main(
                [
                    "--reference-record",
                    str(record_path),
                    "--reference-fasta",
                    str(reference_fasta),
                    "--submitted-fasta",
                    str(submitted_fasta),
                    "--reference-pdb",
                    str(reference_pdb),
                    "--receptor-pdb",
                    str(receptor_pdb),
                    "--reference-chain",
                    "R",
                    "--receptor-chain",
                    "B",
                    "--output",
                    str(output_path),
                ]
            )

            self.assertEqual(
                exit_code,
                0,
            )

            loaded = (
                load_pocket_evidence(
                    output_path
                )
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
