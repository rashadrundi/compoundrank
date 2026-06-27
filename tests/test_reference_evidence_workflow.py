from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from compoundrank.pocket_evidence import (
    load_pocket_evidence,
)
from compoundrank.reference_evidence_workflow import (
    choose_pdb_candidate,
    choose_reference_chain,
    run_reference_evidence_workflow,
)


AA1_TO_3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
}


def write_fasta(
    path: Path,
    sequence: str,
) -> None:
    path.write_text(
        f">submitted\n{sequence}\n",
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
            sequence
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


def example_payload() -> dict:
    return {
        "primaryAccession": "P12345",
        "uniProtkbId": "TEST_VIRUS",
        "entryType": (
            "UniProtKB reviewed "
            "(Swiss-Prot)"
        ),
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
                        "value": "2.00 A",
                    },
                    {
                        "key": "Chains",
                        "value": "A=1-9",
                    },
                ],
            }
        ],
    }


class ReferenceEvidenceWorkflowTests(
    unittest.TestCase
):
    def test_candidate_override_is_respected(
        self,
    ) -> None:
        candidates = [
            {
                "pdb_id": "1AAA",
            },
            {
                "pdb_id": "2BBB",
            },
        ]

        selected = choose_pdb_candidate(
            candidates,
            requested_pdb_id="2bbb",
        )

        self.assertEqual(
            selected["pdb_id"],
            "2BBB",
        )

    def test_chain_with_best_site_coverage_is_selected(
        self,
    ) -> None:
        candidate = {
            "chain_ranges": [
                {
                    "chain": "A",
                    "uniprot_start": 1,
                    "uniprot_end": 4,
                    "mapped_length": 4,
                },
                {
                    "chain": "B",
                    "uniprot_start": 1,
                    "uniprot_end": 9,
                    "mapped_length": 9,
                },
            ]
        }

        selected = (
            choose_reference_chain(
                candidate,
                functional_positions=[
                    3,
                    5,
                ],
            )
        )

        self.assertEqual(
            selected,
            "B",
        )

    def test_workflow_writes_loadable_pocket_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            submitted_fasta = (
                root / "submitted.faa"
            )

            reference_pdb = (
                root / "reference.pdb"
            )

            receptor_pdb = (
                root / "receptor.pdb"
            )

            output_dir = root / "output"

            write_fasta(
                submitted_fasta,
                "ACDEFGHIK",
            )

            write_chain(
                reference_pdb,
                sequence="ACDEFGHIK",
                chain="A",
                start_number=11,
            )

            write_chain(
                receptor_pdb,
                sequence="ACDEFGHIK",
                chain="B",
                start_number=101,
            )

            result = (
                run_reference_evidence_workflow(
                    payload=(
                        example_payload()
                    ),
                    submitted_fasta=(
                        submitted_fasta
                    ),
                    receptor_pdb=(
                        receptor_pdb
                    ),
                    output_dir=(
                        output_dir
                    ),
                    response_metadata={
                        "source": "test"
                    },
                    reference_pdb=(
                        reference_pdb
                    ),
                    receptor_chain_id="B",
                )
            )

            loaded = (
                load_pocket_evidence(
                    output_dir
                    / "pocket_evidence.json"
                )
            )

            self.assertEqual(
                result[
                    "selected_reference"
                ]["pdb_id"],
                "1AAA",
            )

            self.assertEqual(
                result[
                    "selected_reference"
                ]["chain"],
                "A",
            )

            self.assertEqual(
                loaded["residues"],
                [
                    "ASP:B:103",
                    "PHE:B:105",
                ],
            )

            self.assertEqual(
                loaded[
                    "selection_mode"
                ],
                "prioritize_supported",
            )

            summary = json.loads(
                (
                    output_dir
                    / "workflow_summary.json"
                ).read_text(
                    encoding="utf-8"
                )
            )

            self.assertEqual(
                summary["status"],
                "complete",
            )


if __name__ == "__main__":
    unittest.main()
