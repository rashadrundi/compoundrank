from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from compoundrank.residue_transfer import (
    build_reference_numbering_map,
    extract_pdb_chain_sequences,
    map_submitted_sequence_to_structure,
    smith_waterman_residue_map,
    transfer_reference_residues_to_structure,
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
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def atom_line(
    *,
    serial: int,
    amino_acid: str,
    chain: str,
    residue_number: int,
) -> str:
    residue_name = AA1_TO_3[
        amino_acid
    ]

    return (
        f"ATOM  {serial:5d}  CA  "
        f"{residue_name:>3s} "
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
        "a",
        encoding="utf-8",
    ) as handle:
        for offset, amino_acid in enumerate(
            sequence
        ):
            handle.write(
                atom_line(
                    serial=(
                        offset + 1
                    ),
                    amino_acid=(
                        amino_acid
                    ),
                    chain=chain,
                    residue_number=(
                        start_number
                        + offset
                    ),
                )
            )


class ResidueTransferTests(
    unittest.TestCase
):
    def test_exact_alignment_maps_positions(
        self,
    ) -> None:
        alignment = (
            smith_waterman_residue_map(
                query_sequence="ACDEFG",
                reference_sequence=(
                    "ACDEFG"
                ),
            )
        )

        self.assertEqual(
            alignment[
                "reference_to_query"
            ],
            {
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                6: 6,
            },
        )

        self.assertAlmostEqual(
            alignment["identity"],
            1.0,
        )

        self.assertAlmostEqual(
            alignment[
                "reference_coverage"
            ],
            1.0,
        )

        self.assertEqual(
            alignment[
                "alignment_grade"
            ],
            "strong",
        )

    def test_terminal_insertion_shifts_mapping(
        self,
    ) -> None:
        alignment = (
            smith_waterman_residue_map(
                query_sequence=(
                    "TTACDEFGAA"
                ),
                reference_sequence=(
                    "ACDEFG"
                ),
            )
        )

        self.assertEqual(
            alignment[
                "reference_to_query"
            ][1],
            3,
        )

        self.assertEqual(
            alignment[
                "reference_to_query"
            ][6],
            8,
        )

    def test_extracts_ordered_pdb_chain_sequence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "receptor.pdb"
            )

            write_chain(
                path,
                sequence="ACDE",
                chain="B",
                start_number=101,
            )

            chains = (
                extract_pdb_chain_sequences(
                    path
                )
            )

            self.assertEqual(
                chains["B"][
                    "sequence"
                ],
                "ACDE",
            )

            self.assertEqual(
                chains["B"][
                    "residue_ids"
                ],
                [
                    "ALA:B:101",
                    "CYS:B:102",
                    "ASP:B:103",
                    "GLU:B:104",
                ],
            )

    def test_best_matching_structure_chain_is_selected(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "receptor.pdb"
            )

            write_chain(
                path,
                sequence="AAAA",
                chain="A",
                start_number=1,
            )

            write_chain(
                path,
                sequence="ACDE",
                chain="B",
                start_number=101,
            )

            mapping = (
                map_submitted_sequence_to_structure(
                    submitted_sequence=(
                        "ACDE"
                    ),
                    receptor_pdb=path,
                )
            )

            self.assertEqual(
                mapping[
                    "selected_chain"
                ],
                "B",
            )

            self.assertEqual(
                mapping[
                    "submitted_to_structure"
                ][3],
                "ASP:B:103",
            )

    def test_reference_residue_transfers_through_sequence_offset(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "receptor.pdb"
            )

            write_chain(
                path,
                sequence="ACDEFG",
                chain="B",
                start_number=101,
            )

            result = (
                transfer_reference_residues_to_structure(
                    reference_sequence=(
                        "ACDEFG"
                    ),
                    submitted_sequence=(
                        "TTACDEFGAA"
                    ),
                    receptor_pdb=path,
                    reference_residues=[
                        {
                            "position": 3,
                            "amino_acid": "D",
                            "label": (
                                "catalytic"
                            ),
                        }
                    ],
                )
            )

            self.assertEqual(
                result["residues"],
                ["ASP:B:103"],
            )

            transfer = result[
                "transfers"
            ][0]

            self.assertEqual(
                transfer[
                    "submitted_position"
                ],
                5,
            )

            self.assertEqual(
                transfer["status"],
                "mapped_conserved",
            )

            self.assertTrue(
                transfer[
                    "selection_eligible"
                ]
            )

    def test_substituted_residue_is_not_selection_eligible(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "receptor.pdb"
            )

            write_chain(
                path,
                sequence="ACNEFG",
                chain="B",
                start_number=101,
            )

            result = (
                transfer_reference_residues_to_structure(
                    reference_sequence=(
                        "ACDEFG"
                    ),
                    submitted_sequence=(
                        "ACNEFG"
                    ),
                    receptor_pdb=path,
                    reference_residues=[
                        {
                            "position": 3,
                            "amino_acid": "D",
                        }
                    ],
                )
            )

            self.assertEqual(
                result["residues"],
                [],
            )

            self.assertEqual(
                result["transfers"][0][
                    "status"
                ],
                "mapped_with_substitution",
            )

            self.assertFalse(
                result["transfers"][0][
                    "selection_eligible"
                ]
            )

    def test_source_amino_acid_mismatch_is_rejected(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "receptor.pdb"
            )

            write_chain(
                path,
                sequence="ACDEFG",
                chain="B",
                start_number=101,
            )

            result = (
                transfer_reference_residues_to_structure(
                    reference_sequence=(
                        "ACDEFG"
                    ),
                    submitted_sequence=(
                        "ACDEFG"
                    ),
                    receptor_pdb=path,
                    reference_residues=[
                        {
                            "position": 3,
                            "amino_acid": "E",
                        }
                    ],
                )
            )

            self.assertEqual(
                result["residues"],
                [],
            )

            self.assertEqual(
                result["transfers"][0][
                    "status"
                ],
                "reference_amino_acid_mismatch",
            )

    def test_reference_pdb_number_maps_through_sequence(
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

            result = (
                transfer_reference_residues_to_structure(
                    reference_sequence="ACDEFG",
                    submitted_sequence=(
                        "TTACDEFGAA"
                    ),
                    reference_pdb=(
                        reference_pdb
                    ),
                    reference_chain_id="R",
                    receptor_pdb=(
                        receptor_pdb
                    ),
                    chain_id="B",
                    reference_residues=[
                        {
                            "residue_number": 85,
                            "amino_acid": "D",
                            "label": "catalytic",
                        }
                    ],
                )
            )

            self.assertEqual(
                result["residues"],
                ["ASP:B:103"],
            )

            transfer = result[
                "transfers"
            ][0]

            self.assertEqual(
                transfer[
                    "reference_sequence_position"
                ],
                3,
            )

            self.assertEqual(
                transfer[
                    "submitted_position"
                ],
                5,
            )

            self.assertEqual(
                transfer[
                    "reference_residue_id"
                ],
                "ASP:R:85",
            )

            self.assertEqual(
                transfer[
                    "reference_numbering_mode"
                ],
                "pdb_residue_number",
            )

            self.assertEqual(
                transfer["status"],
                "mapped_conserved",
            )

    def test_numbering_gap_does_not_use_constant_offset(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = (
                Path(directory)
                / "reference.pdb"
            )

            with path.open(
                "w",
                encoding="utf-8",
            ) as handle:
                handle.write(
                    atom_line(
                        serial=1,
                        amino_acid="A",
                        chain="R",
                        residue_number=83,
                    )
                )

                handle.write(
                    atom_line(
                        serial=2,
                        amino_acid="C",
                        chain="R",
                        residue_number=84,
                    )
                )

                handle.write(
                    atom_line(
                        serial=3,
                        amino_acid="D",
                        chain="R",
                        residue_number=86,
                    )
                )

                handle.write(
                    atom_line(
                        serial=4,
                        amino_acid="E",
                        chain="R",
                        residue_number=87,
                    )
                )

            numbering = (
                build_reference_numbering_map(
                    reference_sequence=(
                        "ACDE"
                    ),
                    reference_pdb=path,
                    chain_id="R",
                )
            )

            self.assertEqual(
                numbering[
                    "residues_by_number"
                ]["86"][
                    "reference_sequence_position"
                ],
                3,
            )

            self.assertEqual(
                numbering[
                    "residues_by_number"
                ]["87"][
                    "reference_sequence_position"
                ],
                4,
            )

    def test_biological_number_requires_reference_map(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            receptor = (
                Path(directory)
                / "receptor.pdb"
            )

            write_chain(
                receptor,
                sequence="ACDE",
                chain="B",
                start_number=101,
            )

            result = (
                transfer_reference_residues_to_structure(
                    reference_sequence="ACDE",
                    submitted_sequence="ACDE",
                    receptor_pdb=receptor,
                    reference_residues=[
                        {
                            "residue_number": 85,
                            "amino_acid": "D",
                        }
                    ],
                )
            )

            self.assertEqual(
                result["residues"],
                [],
            )

            self.assertEqual(
                result[
                    "transfers"
                ][0]["status"],
                "reference_numbering_map_required",
            )



if __name__ == "__main__":
    unittest.main()
