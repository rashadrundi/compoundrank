from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from compoundrank.pipeline import (
    _resolve_pocket_evidence_json,
)


class PipelineAutoReferenceEvidenceTests(
    unittest.TestCase
):
    def test_manual_evidence_is_returned_unchanged(
        self,
    ) -> None:
        evidence = Path(
            "/tmp/manual_evidence.json"
        )

        resolved = (
            _resolve_pocket_evidence_json(
                pocket_evidence_json=(
                    evidence
                ),
                auto_reference_evidence=False,
                reference_uniprot_accession=None,
                reference_uniprot_json=None,
                reference_pdb_id=None,
                reference_chain_id=None,
                receptor_chain_id=None,
                reference_pdb=None,
                reference_evidence_timeout_seconds=(
                    60.0
                ),
                fasta_path=None,
                receptor_pdb=Path(
                    "/tmp/receptor.pdb"
                ),
                output_dir=Path(
                    "/tmp/output"
                ),
            )
        )

        self.assertEqual(
            resolved,
            evidence,
        )

    def test_automatic_workflow_result_becomes_effective_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            uniprot_json = (
                root / "uniprot.json"
            )

            fasta = root / "protein.faa"
            receptor = root / "receptor.pdb"
            generated = (
                root
                / "output"
                / "automatic_reference_evidence"
                / "pocket_evidence.json"
            )

            uniprot_json.write_text(
                json.dumps(
                    {
                        "primaryAccession": (
                            "P12345"
                        )
                    }
                ),
                encoding="utf-8",
            )

            fasta.write_text(
                ">protein\nAAAA\n",
                encoding="utf-8",
            )

            receptor.write_text(
                "END\n",
                encoding="utf-8",
            )

            with patch(
                "compoundrank.pipeline."
                "run_reference_evidence_workflow",
                return_value={
                    "pocket_evidence_path": (
                        str(generated)
                    )
                },
            ) as mocked:
                resolved = (
                    _resolve_pocket_evidence_json(
                        pocket_evidence_json=None,
                        auto_reference_evidence=True,
                        reference_uniprot_accession=None,
                        reference_uniprot_json=(
                            uniprot_json
                        ),
                        reference_pdb_id="1ABC",
                        reference_chain_id="A",
                        receptor_chain_id="B",
                        reference_pdb=None,
                        reference_evidence_timeout_seconds=(
                            30.0
                        ),
                        fasta_path=fasta,
                        receptor_pdb=receptor,
                        output_dir=(
                            root / "output"
                        ),
                    )
                )

            self.assertEqual(
                resolved,
                generated,
            )

            kwargs = mocked.call_args.kwargs

            self.assertEqual(
                kwargs["requested_pdb_id"],
                "1ABC",
            )

            self.assertEqual(
                kwargs["reference_chain_id"],
                "A",
            )

            self.assertEqual(
                kwargs["receptor_chain_id"],
                "B",
            )

            self.assertEqual(
                kwargs["submitted_fasta"],
                fasta,
            )

    def test_fasta_accession_is_resolved_when_source_omitted(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            fasta = (
                root / "protein_Q6DPL2.faa"
            )
            receptor = root / "receptor.pdb"
            output = root / "output"
            generated = (
                output
                / "automatic_reference_evidence"
                / "pocket_evidence.json"
            )

            fasta.write_text(
                ">protein\nAAAA\n",
                encoding="utf-8",
            )

            receptor.write_text(
                "END\n",
                encoding="utf-8",
            )

            resolution = {
                "schema_version": (
                    "uniprot_accession_"
                    "resolution.v0.1"
                ),
                "status": "resolved",
                "selected_accession": (
                    "Q6DPL2"
                ),
                "resolution_method": (
                    "filename_token"
                ),
            }

            payload = {
                "primaryAccession": (
                    "Q6DPL2"
                )
            }

            with patch(
                "compoundrank.pipeline."
                "resolve_uniprot_accession_from_fasta",
                return_value=resolution,
            ) as resolver, patch(
                "compoundrank.pipeline."
                "fetch_uniprot_entry",
                return_value=(
                    payload,
                    {
                        "source": (
                            "uniprot_rest"
                        )
                    },
                ),
            ) as fetcher, patch(
                "compoundrank.pipeline."
                "run_reference_evidence_workflow",
                return_value={
                    "pocket_evidence_path": (
                        str(generated)
                    )
                },
            ) as workflow:
                resolved = (
                    _resolve_pocket_evidence_json(
                        pocket_evidence_json=None,
                        auto_reference_evidence=True,
                        reference_uniprot_accession=None,
                        reference_uniprot_json=None,
                        reference_pdb_id=None,
                        reference_chain_id=None,
                        receptor_chain_id="B",
                        reference_pdb=None,
                        reference_evidence_timeout_seconds=(
                            60.0
                        ),
                        fasta_path=fasta,
                        receptor_pdb=receptor,
                        output_dir=output,
                    )
                )

            self.assertEqual(
                resolved,
                generated,
            )

            resolver.assert_called_once_with(
                fasta,
                sequence_search_email=None,
                sequence_search_timeout_seconds=(
                    600.0
                ),
            )

            fetcher.assert_called_once_with(
                "Q6DPL2",
                timeout_seconds=60.0,
            )

            self.assertEqual(
                workflow.call_args.kwargs[
                    "payload"
                ],
                payload,
            )

            resolution_path = (
                output
                / "automatic_reference_evidence"
                / (
                    "uniprot_accession_"
                    "resolution.json"
                )
            )

            self.assertTrue(
                resolution_path.is_file()
            )


if __name__ == "__main__":
    unittest.main()
