from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.functional_site_transfer import (
    load_functional_site_reference,
)
from compoundrank.uniprot_acquisition import (
    build_uniprot_acquisition,
    main,
)


def example_payload(
    *,
    entry_type: str = (
        "UniProtKB reviewed (Swiss-Prot)"
    ),
) -> dict:
    return {
        "primaryAccession": "P12345",
        "uniProtkbId": "TEST_VIRUS",
        "entryType": entry_type,
        "annotationScore": 5.0,
        "organism": {
            "scientificName": (
                "Example virus"
            ),
            "taxonId": 123,
        },
        "sequence": {
            "value": "ACDEFGHIK",
            "length": 9,
        },
        "features": [
            {
                "type": "Active site",
                "description": (
                    "Proton donor"
                ),
                "location": {
                    "start": {
                        "value": 3
                    },
                    "end": {
                        "value": 3
                    },
                },
                "evidences": [
                    {
                        "evidenceCode": (
                            "ECO:0000255"
                        ),
                        "source": (
                            "HAMAP-Rule"
                        ),
                    }
                ],
            },
            {
                "type": "Binding site",
                "location": {
                    "start": {
                        "value": 5
                    },
                    "end": {
                        "value": 5
                    },
                },
                "ligand": {
                    "name": "substrate"
                },
                "evidences": [
                    {
                        "evidenceCode": (
                            "ECO:0000269"
                        ),
                        "source": "PubMed",
                    }
                ],
            },
            {
                "type": "Binding site",
                "location": {
                    "start": {
                        "value": 6
                    },
                    "end": {
                        "value": 6
                    },
                },
                "ligand": {
                    "name": "Ca(2+)"
                },
                "evidences": [
                    {
                        "evidenceCode": (
                            "ECO:0007744"
                        ),
                        "source": "PDB",
                    }
                ],
            },
        ],
        "uniProtKBCrossReferences": [
            {
                "database": "PDB",
                "id": "1AAA",
                "properties": [
                    {
                        "key": "Method",
                        "value": "X-ray",
                    },
                    {
                        "key": "Resolution",
                        "value": "2.50 A",
                    },
                    {
                        "key": "Chains",
                        "value": "A=1-9",
                    },
                ],
            },
            {
                "database": "PDB",
                "id": "2BBB",
                "properties": [
                    {
                        "key": "Method",
                        "value": "X-ray",
                    },
                    {
                        "key": "Resolution",
                        "value": "1.80 A",
                    },
                    {
                        "key": "Chains",
                        "value": "B=1-9",
                    },
                ],
            },
        ],
    }


class UniProtAcquisitionTests(
    unittest.TestCase
):
    def test_builds_reference_and_excludes_metal_site(
        self,
    ) -> None:
        result = (
            build_uniprot_acquisition(
                example_payload()
            )
        )

        reference = result[
            "reference_record"
        ]

        self.assertEqual(
            [
                residue[
                    "sequence_position"
                ]
                for residue in reference[
                    "residues"
                ]
            ],
            [3, 5],
        )

        self.assertEqual(
            [
                residue[
                    "amino_acid"
                ]
                for residue in reference[
                    "residues"
                ]
            ],
            ["D", "F"],
        )

        self.assertEqual(
            reference[
                "selection_mode"
            ],
            "prioritize_supported",
        )

        self.assertEqual(
            result[
                "pdb_candidates"
            ][0]["pdb_id"],
            "2BBB",
        )

        self.assertEqual(
            result["summary"][
                "excluded_feature_count"
            ],
            1,
        )

    def test_unreviewed_entry_is_report_only(
        self,
    ) -> None:
        result = (
            build_uniprot_acquisition(
                example_payload(
                    entry_type=(
                        "UniProtKB unreviewed "
                        "(TrEMBL)"
                    )
                )
            )
        )

        self.assertEqual(
            result[
                "reference_record"
            ]["selection_mode"],
            "report_only",
        )

        self.assertEqual(
            result[
                "reference_record"
            ]["confidence"],
            "low",
        )

    def test_cli_writes_loadable_reference_packet(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            input_path = (
                root / "entry.json"
            )

            output_dir = (
                root / "output"
            )

            input_path.write_text(
                json.dumps(
                    example_payload()
                ),
                encoding="utf-8",
            )

            exit_code = main(
                [
                    "--input-json",
                    str(input_path),
                    "--output-dir",
                    str(output_dir),
                ]
            )

            self.assertEqual(
                exit_code,
                0,
            )

            reference = (
                load_functional_site_reference(
                    output_dir
                    / "functional_site_reference.json"
                )
            )

            self.assertEqual(
                reference[
                    "selection_mode"
                ],
                "prioritize_supported",
            )

            self.assertEqual(
                len(
                    reference[
                        "residues"
                    ]
                ),
                2,
            )

            self.assertTrue(
                (
                    output_dir
                    / "P12345.fasta"
                ).is_file()
            )


if __name__ == "__main__":
    unittest.main()
