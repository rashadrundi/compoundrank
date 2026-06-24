from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from compoundrank.uniprot_accession import (
    resolve_uniprot_accession_from_fasta,
)


class UniProtAccessionResolutionTests(
    unittest.TestCase
):
    def test_resolves_canonical_swissprot_header(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "protein.faa"
            )

            path.write_text(
                (
                    ">sp|Q6DPL2|NA_INFA5 "
                    "Neuraminidase\n"
                    "ACDEFG\n"
                ),
                encoding="utf-8",
            )

            result = (
                resolve_uniprot_accession_from_fasta(
                    path
                )
            )

            self.assertEqual(
                result[
                    "selected_accession"
                ],
                "Q6DPL2",
            )

            self.assertEqual(
                result[
                    "resolution_method"
                ],
                "canonical_uniprot_header",
            )

    def test_resolves_strict_plain_header_token(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "protein.faa"
            )

            path.write_text(
                (
                    ">influenza_neuraminidase "
                    "Q6DPL2 chain B\n"
                    "ACDEFG\n"
                ),
                encoding="utf-8",
            )

            result = (
                resolve_uniprot_accession_from_fasta(
                    path
                )
            )

            self.assertEqual(
                result[
                    "selected_accession"
                ],
                "Q6DPL2",
            )

            self.assertEqual(
                result[
                    "resolution_method"
                ],
                "header_token",
            )

    def test_filename_is_conservative_fallback(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / (
                    "2HU0_neuraminidase_"
                    "Q6DPL2.faa"
                )
            )

            path.write_text(
                ">chain_B\nACDEFG\n",
                encoding="utf-8",
            )

            result = (
                resolve_uniprot_accession_from_fasta(
                    path
                )
            )

            self.assertEqual(
                result[
                    "selected_accession"
                ],
                "Q6DPL2",
            )

            self.assertEqual(
                result[
                    "resolution_method"
                ],
                "filename_token",
            )

    def test_ambiguous_header_is_rejected(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "protein.faa"
            )

            path.write_text(
                (
                    ">possible Q6DPL2 "
                    "or P12345\n"
                    "ACDEFG\n"
                ),
                encoding="utf-8",
            )

            with self.assertRaises(
                ValueError
            ):
                resolve_uniprot_accession_from_fasta(
                    path
                )

    def test_multiple_fasta_records_are_rejected(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "protein_Q6DPL2.faa"
            )

            path.write_text(
                (
                    ">one\nACD\n"
                    ">two\nEFG\n"
                ),
                encoding="utf-8",
            )

            with self.assertRaises(
                ValueError
            ):
                resolve_uniprot_accession_from_fasta(
                    path
                )


if __name__ == "__main__":
    unittest.main()
