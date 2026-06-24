from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from compoundrank.cli import (
    build_parser,
    main,
)


class CliAutoReferenceEvidenceTests(
    unittest.TestCase
):
    def test_parser_accepts_automatic_reference_flags(
        self,
    ) -> None:
        arguments = (
            build_parser().parse_args(
                [
                    "--receptor",
                    "/tmp/receptor.pdb",
                    "--data-root",
                    "/tmp/data",
                    "--fasta",
                    "/tmp/protein.faa",
                    "--auto-reference-evidence",
                    "--reference-uniprot-accession",
                    "Q6DPL2",
                    "--reference-pdb-id",
                    "3CKZ",
                    "--reference-chain",
                    "A",
                    "--receptor-chain",
                    "B",
                ]
            )
        )

        self.assertTrue(
            arguments.auto_reference_evidence
        )

        self.assertEqual(
            arguments.reference_uniprot_accession,
            "Q6DPL2",
        )

        self.assertEqual(
            arguments.reference_pdb_id,
            "3CKZ",
        )

    def test_manual_and_automatic_evidence_conflict(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            receptor = root / "receptor.pdb"
            fasta = root / "protein.faa"

            receptor.write_text(
                "END\n",
                encoding="utf-8",
            )

            fasta.write_text(
                ">protein\nAAAA\n",
                encoding="utf-8",
            )

            with self.assertRaises(
                ValueError
            ):
                main(
                    [
                        "--receptor",
                        str(receptor),
                        "--data-root",
                        str(root / "data"),
                        "--fasta",
                        str(fasta),
                        "--auto-reference-evidence",
                        "--reference-uniprot-accession",
                        "Q6DPL2",
                        "--pocket-evidence-json",
                        str(
                            root / "manual.json"
                        ),
                    ]
                )

    def test_automatic_reference_arguments_are_forwarded(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            receptor = root / "receptor.pdb"
            fasta = root / "protein.faa"
            ligand = root / "ligand.sdf"
            uniprot_json = (
                root / "uniprot.json"
            )
            reference_pdb = (
                root / "reference.pdb"
            )

            receptor.write_text(
                "END\n",
                encoding="utf-8",
            )

            fasta.write_text(
                ">protein\nAAAA\n",
                encoding="utf-8",
            )

            ligand.write_text(
                "test ligand\n",
                encoding="utf-8",
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

            reference_pdb.write_text(
                "END\n",
                encoding="utf-8",
            )

            with patch(
                "compoundrank.cli.run_pipeline"
            ) as mocked:
                main(
                    [
                        "--receptor",
                        str(receptor),
                        "--data-root",
                        str(root / "data"),
                        "--fasta",
                        str(fasta),
                        "--ligand-file",
                        str(ligand),
                        "--auto-reference-evidence",
                        "--reference-uniprot-json",
                        str(uniprot_json),
                        "--reference-pdb-id",
                        "1ABC",
                        "--reference-chain",
                        "A",
                        "--receptor-chain",
                        "B",
                        "--reference-pdb",
                        str(reference_pdb),
                        "--reference-evidence-timeout-seconds",
                        "45",
                    ]
                )

            kwargs = mocked.call_args.kwargs

            self.assertTrue(
                kwargs[
                    "auto_reference_evidence"
                ]
            )

            self.assertEqual(
                kwargs[
                    "reference_uniprot_json"
                ],
                uniprot_json.resolve(),
            )

            self.assertEqual(
                kwargs["reference_pdb_id"],
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
                kwargs["reference_pdb"],
                reference_pdb.resolve(),
            )

            self.assertEqual(
                kwargs[
                    "reference_evidence_timeout_seconds"
                ],
                45.0,
            )


if __name__ == "__main__":
    unittest.main()
