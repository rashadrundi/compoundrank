from __future__ import annotations

import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from compoundrank.uniprot_accession import (
    parse_ebi_blast_tsv,
    parse_ebi_blast_xml,
    resolve_uniprot_accession_by_sequence,
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

    def test_parses_ebi_blast_xml(
        self,
    ) -> None:
        xml_text = """
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>sp|Q6DPL2|NA_INFA5</Hit_id>
          <Hit_def>Neuraminidase</Hit_def>
          <Hit_accession>Q6DPL2</Hit_accession>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>500</Hsp_bit-score>
              <Hsp_evalue>1e-100</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>20</Hsp_hit-to>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_align-len>20</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

        candidates = parse_ebi_blast_xml(
            xml_text
        )

        self.assertEqual(
            len(candidates),
            1,
        )

        self.assertEqual(
            candidates[0]["accession"],
            "Q6DPL2",
        )

        self.assertEqual(
            candidates[0][
                "best_hsp"
            ]["identity"],
            1.0,
        )

    def test_sequence_search_accepts_strong_match(
        self,
    ) -> None:
        sequence = (
            "ACDEFGHIKLMNPQRSTVWY"
        )

        xml_text = """
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>sp|Q6DPL2|NA_INFA5</Hit_id>
          <Hit_def>Neuraminidase</Hit_def>
          <Hit_accession>Q6DPL2</Hit_accession>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>500</Hsp_bit-score>
              <Hsp_evalue>1e-100</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>20</Hsp_hit-to>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_align-len>20</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

        payload = {
            "primaryAccession": "Q6DPL2",
            "entryType": (
                "UniProtKB reviewed "
                "(Swiss-Prot)"
            ),
            "annotationScore": 5.0,
            "sequence": {
                "value": sequence,
            },
            "features": [
                {
                    "type": "Active site",
                }
            ],
            "uniProtKBCrossReferences": [
                {
                    "database": "PDB",
                    "id": "3CKZ",
                }
            ],
        }

        with patch(
            "compoundrank.uniprot_accession."
            "_submit_ebi_blast",
            return_value=(
                "job-test",
                {"status": 200},
            ),
        ), patch(
            "compoundrank.uniprot_accession."
            "_poll_ebi_blast",
            return_value={
                "job_id": "job-test",
                "status": "FINISHED",
            },
        ), patch(
            "compoundrank.uniprot_accession."
            "_fetch_ebi_blast_xml",
            return_value=(
                xml_text,
                {"status": 200},
            ),
        ), patch(
            "compoundrank.uniprot_accession."
            "fetch_uniprot_entry",
            return_value=(
                payload,
                {"source": "test"},
            ),
        ):
            result = (
                resolve_uniprot_accession_by_sequence(
                    sequence,
                    email="test@example.org",
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
            (
                "ebi_blast_"
                "sequence_alignment"
            ),
        )

    def test_sequence_search_rejects_weak_match(
        self,
    ) -> None:
        query = (
            "ACDEFGHIKLMNPQRSTVWY"
        )

        xml_text = """
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>tr|P12345|WEAK</Hit_id>
          <Hit_def>Weak match</Hit_def>
          <Hit_accession>P12345</Hit_accession>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>20</Hsp_bit-score>
              <Hsp_evalue>0.5</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>5</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>5</Hsp_hit-to>
              <Hsp_identity>5</Hsp_identity>
              <Hsp_align-len>5</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

        weak_payload = {
            "primaryAccession": "P12345",
            "entryType": (
                "UniProtKB unreviewed "
                "(TrEMBL)"
            ),
            "sequence": {
                "value": (
                    "AAAAAAAAAAAAAAAAAAAA"
                ),
            },
        }

        with patch(
            "compoundrank.uniprot_accession."
            "_submit_ebi_blast",
            return_value=(
                "job-weak",
                {},
            ),
        ), patch(
            "compoundrank.uniprot_accession."
            "_poll_ebi_blast",
            return_value={
                "status": "FINISHED",
            },
        ), patch(
            "compoundrank.uniprot_accession."
            "_fetch_ebi_blast_xml",
            return_value=(
                xml_text,
                {},
            ),
        ), patch(
            "compoundrank.uniprot_accession."
            "fetch_uniprot_entry",
            return_value=(
                weak_payload,
                {},
            ),
        ):
            with self.assertRaises(
                ValueError
            ):
                resolve_uniprot_accession_by_sequence(
                    query,
                    email="test@example.org",
                )

    def test_parses_ebi_blast_tsv(
        self,
    ) -> None:
        tsv_text = (
            "Hit\tDB\tAccession\tDescription\t"
            "Organism\tLength\tScore(Bits)\t"
            "Identities(%)\tPositives(%)\tE()\n"
            "1\tSP\tQ6DPL2\tNeuraminidase\t"
            "Influenza A virus\t449\t2140\t"
            "99.7\t100.0\t0.0\n"
            "2\tTR\tQ5SDA6\tNeuraminidase\t"
            "Influenza A virus\t449\t2135\t"
            "99.5\t100.0\t0.0\n"
        )

        candidates = (
            parse_ebi_blast_tsv(
                tsv_text
            )
        )

        self.assertEqual(
            len(candidates),
            2,
        )

        self.assertEqual(
            candidates[0][
                "accession"
            ],
            "Q6DPL2",
        )

        self.assertAlmostEqual(
            candidates[0][
                "best_hsp"
            ]["identity"],
            0.997,
        )

        self.assertEqual(
            candidates[0][
                "best_hsp"
            ]["bit_score"],
            2140.0,
        )


if __name__ == "__main__":
    unittest.main()
