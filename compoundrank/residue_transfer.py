from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Iterable


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
}

SCHEMA_VERSION = "residue_transfer.v0.1"


def normalize_protein_sequence(
    sequence: object,
) -> str:
    normalized = "".join(
        character
        for character in str(
            sequence or ""
        ).upper()
        if character.isalpha()
    )

    if not normalized:
        raise ValueError(
            "Protein sequence is empty."
        )

    return normalized


def _alignment_grade(
    *,
    identity: float,
    reference_coverage: float,
) -> str:
    if (
        identity >= 0.90
        and reference_coverage >= 0.80
    ):
        return "strong"

    if (
        identity >= 0.70
        and reference_coverage >= 0.70
    ):
        return "moderate"

    if (
        identity >= 0.40
        and reference_coverage >= 0.50
    ):
        return "weak"

    return "exploratory"


def smith_waterman_residue_map(
    *,
    query_sequence: str,
    reference_sequence: str,
) -> dict[str, Any]:
    """Map 1-based reference positions onto query positions.

    The scoring scheme intentionally matches the existing
    CompoundRank local-alignment implementation:

    - match: +2
    - mismatch: -1
    - gap: -2
    """

    query = normalize_protein_sequence(
        query_sequence
    )
    reference = normalize_protein_sequence(
        reference_sequence
    )

    rows = len(query) + 1
    columns = len(reference) + 1

    scores = [
        [0] * columns
        for _ in range(rows)
    ]

    directions = [
        [0] * columns
        for _ in range(rows)
    ]

    best_score = 0
    best_position = (0, 0)

    for query_index in range(
        1,
        rows,
    ):
        for reference_index in range(
            1,
            columns,
        ):
            diagonal = (
                scores[
                    query_index - 1
                ][
                    reference_index - 1
                ]
                + (
                    2
                    if (
                        query[
                            query_index - 1
                        ]
                        == reference[
                            reference_index - 1
                        ]
                    )
                    else -1
                )
            )

            upward = (
                scores[
                    query_index - 1
                ][reference_index]
                - 2
            )

            leftward = (
                scores[query_index][
                    reference_index - 1
                ]
                - 2
            )

            score = max(
                0,
                diagonal,
                upward,
                leftward,
            )

            scores[
                query_index
            ][reference_index] = score

            if score == 0:
                direction = 0
            elif score == diagonal:
                direction = 1
            elif score == upward:
                direction = 2
            else:
                direction = 3

            directions[
                query_index
            ][reference_index] = (
                direction
            )

            if score > best_score:
                best_score = score
                best_position = (
                    query_index,
                    reference_index,
                )

    query_index, reference_index = (
        best_position
    )

    query_end = query_index
    reference_end = reference_index

    matches = 0
    aligned_columns = 0
    query_residues = 0
    reference_residues = 0

    mapping: dict[int, int] = {}

    aligned_pairs: list[
        dict[str, Any]
    ] = []

    while (
        query_index > 0
        and reference_index > 0
    ):
        direction = directions[
            query_index
        ][reference_index]

        if direction == 0:
            break

        aligned_columns += 1

        if direction == 1:
            query_amino_acid = query[
                query_index - 1
            ]

            reference_amino_acid = (
                reference[
                    reference_index - 1
                ]
            )

            mapping[
                reference_index
            ] = query_index

            conserved = (
                query_amino_acid
                == reference_amino_acid
            )

            if conserved:
                matches += 1

            aligned_pairs.append(
                {
                    "reference_position": (
                        reference_index
                    ),
                    "reference_amino_acid": (
                        reference_amino_acid
                    ),
                    "query_position": (
                        query_index
                    ),
                    "query_amino_acid": (
                        query_amino_acid
                    ),
                    "conserved": conserved,
                }
            )

            query_index -= 1
            reference_index -= 1
            query_residues += 1
            reference_residues += 1

        elif direction == 2:
            query_index -= 1
            query_residues += 1

        else:
            reference_index -= 1
            reference_residues += 1

    aligned_pairs.reverse()

    query_start = (
        query_index + 1
        if query_residues
        else None
    )

    reference_start = (
        reference_index + 1
        if reference_residues
        else None
    )

    identity = (
        matches / aligned_columns
        if aligned_columns
        else 0.0
    )

    query_coverage = (
        query_residues / len(query)
    )

    reference_coverage = (
        reference_residues
        / len(reference)
    )

    return {
        "score": best_score,
        "identity": identity,
        "query_coverage": (
            query_coverage
        ),
        "reference_coverage": (
            reference_coverage
        ),
        "alignment_grade": (
            _alignment_grade(
                identity=identity,
                reference_coverage=(
                    reference_coverage
                ),
            )
        ),
        "matches": matches,
        "aligned_columns": (
            aligned_columns
        ),
        "query_residues": query_residues,
        "reference_residues": (
            reference_residues
        ),
        "query_start": query_start,
        "query_end": (
            query_end
            if query_residues
            else None
        ),
        "reference_start": (
            reference_start
        ),
        "reference_end": (
            reference_end
            if reference_residues
            else None
        ),
        "reference_to_query": mapping,
        "aligned_pairs": aligned_pairs,
    }


def extract_pdb_chain_sequences(
    pdb_path: Path,
) -> dict[str, dict[str, Any]]:
    """Extract ordered residue sequences from receptor chains."""

    path = Path(pdb_path)

    if not path.is_file():
        raise FileNotFoundError(path)

    chain_records: dict[
        str,
        list[dict[str, str]],
    ] = {}

    seen: set[
        tuple[str, str, str]
    ] = set()

    for line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        if not line.startswith(
            ("ATOM", "HETATM")
        ):
            continue

        residue_name = (
            line[17:20].strip().upper()
        )

        if residue_name not in AA3_TO_1:
            continue

        chain = (
            line[21:22].strip().upper()
            or "_"
        )

        residue_number = (
            line[22:26].strip()
        )

        insertion_code = (
            line[26:27].strip().upper()
        )

        if not residue_number:
            continue

        residue_key = (
            chain,
            residue_number,
            insertion_code,
        )

        if residue_key in seen:
            continue

        seen.add(residue_key)

        residue_id = (
            f"{residue_name}:"
            f"{chain}:"
            f"{residue_number}"
            f"{insertion_code}"
        )

        chain_records.setdefault(
            chain,
            [],
        ).append(
            {
                "amino_acid": (
                    AA3_TO_1[
                        residue_name
                    ]
                ),
                "residue_id": (
                    residue_id
                ),
            }
        )

    if not chain_records:
        raise ValueError(
            "No protein residues were found "
            f"in {path}."
        )

    return {
        chain: {
            "sequence": "".join(
                record["amino_acid"]
                for record in records
            ),
            "residue_ids": [
                record["residue_id"]
                for record in records
            ],
        }
        for chain, records
        in chain_records.items()
    }




def _normalize_residue_number(
    value: object,
    *,
    insertion_code: object = "",
) -> str:
    text = str(value or "").strip().upper()

    match = re.fullmatch(
        r"(-?\d+)([A-Z]?)",
        text,
    )

    if match is None:
        raise ValueError(
            "Residue numbers must contain an integer "
            "with an optional insertion code; "
            f"received {value!r}"
        )

    number = match.group(1)
    embedded_insertion = match.group(2)

    supplied_insertion = str(
        insertion_code or ""
    ).strip().upper()

    if (
        embedded_insertion
        and supplied_insertion
        and embedded_insertion
        != supplied_insertion
    ):
        raise ValueError(
            "Conflicting residue insertion codes: "
            f"{embedded_insertion!r} and "
            f"{supplied_insertion!r}"
        )

    return (
        number
        + (
            supplied_insertion
            or embedded_insertion
        )
    )


def build_reference_numbering_map(
    *,
    reference_sequence: str,
    reference_pdb: Path,
    chain_id: str | None = None,
) -> dict[str, Any]:
    """Map biological/PDB residue numbers to FASTA positions.

    The reference PDB chain is aligned to the complete reference
    sequence. This handles truncated structures, unresolved residues,
    numbering offsets, numbering gaps, and insertion codes without
    assuming that a PDB residue number equals a FASTA index.
    """

    reference = normalize_protein_sequence(
        reference_sequence
    )

    chains = extract_pdb_chain_sequences(
        reference_pdb
    )

    normalized_chain = (
        chain_id.strip().upper()
        if chain_id
        else None
    )

    if normalized_chain is not None:
        if normalized_chain not in chains:
            raise ValueError(
                "Requested reference chain was not "
                f"found: {normalized_chain}"
            )

        candidate_chains = {
            normalized_chain: chains[
                normalized_chain
            ]
        }
    else:
        candidate_chains = chains

    candidates = []

    for chain, chain_data in (
        candidate_chains.items()
    ):
        alignment = (
            smith_waterman_residue_map(
                query_sequence=(
                    chain_data["sequence"]
                ),
                reference_sequence=(
                    reference
                ),
            )
        )

        candidates.append(
            {
                "chain": chain,
                "chain_data": chain_data,
                "alignment": alignment,
            }
        )

    selected = max(
        candidates,
        key=lambda candidate: (
            candidate["alignment"][
                "reference_coverage"
            ],
            candidate["alignment"][
                "identity"
            ],
            candidate["alignment"][
                "score"
            ],
            candidate["alignment"][
                "query_coverage"
            ],
        ),
    )

    query_to_reference = {
        int(pair["query_position"]): int(
            pair["reference_position"]
        )
        for pair in selected[
            "alignment"
        ]["aligned_pairs"]
    }

    chain_data = selected[
        "chain_data"
    ]

    residues_by_number = {}
    unresolved_structure_residues = []

    for chain_position, residue_id in enumerate(
        chain_data["residue_ids"],
        start=1,
    ):
        residue_number = (
            residue_id.rsplit(
                ":",
                1,
            )[1].upper()
        )

        reference_position = (
            query_to_reference.get(
                chain_position
            )
        )

        record = {
            "reference_chain": (
                selected["chain"]
            ),
            "reference_residue_number": (
                residue_number
            ),
            "reference_residue_id": (
                residue_id
            ),
            "reference_chain_position": (
                chain_position
            ),
            "reference_sequence_position": (
                reference_position
            ),
            "amino_acid": (
                chain_data["sequence"][
                    chain_position - 1
                ]
            ),
        }

        if reference_position is None:
            unresolved_structure_residues.append(
                record
            )
            continue

        residues_by_number[
            residue_number
        ] = record

    return {
        "schema_version": (
            "reference_numbering_map.v0.1"
        ),
        "reference_pdb": str(
            Path(reference_pdb)
        ),
        "selected_chain": (
            selected["chain"]
        ),
        "alignment": (
            selected["alignment"]
        ),
        "residue_count": len(
            residues_by_number
        ),
        "residues_by_number": (
            residues_by_number
        ),
        "unresolved_structure_residues": (
            unresolved_structure_residues
        ),
        "candidate_chains": [
            {
                "chain": (
                    candidate["chain"]
                ),
                "identity": (
                    candidate["alignment"][
                        "identity"
                    ]
                ),
                "reference_coverage": (
                    candidate["alignment"][
                        "reference_coverage"
                    ]
                ),
                "structure_coverage": (
                    candidate["alignment"][
                        "query_coverage"
                    ]
                ),
                "score": (
                    candidate["alignment"][
                        "score"
                    ]
                ),
            }
            for candidate in candidates
        ],
    }

def map_submitted_sequence_to_structure(
    *,
    submitted_sequence: str,
    receptor_pdb: Path,
    chain_id: str | None = None,
) -> dict[str, Any]:
    """Map submitted-sequence positions to receptor PDB residues."""

    submitted = (
        normalize_protein_sequence(
            submitted_sequence
        )
    )

    chains = extract_pdb_chain_sequences(
        receptor_pdb
    )

    normalized_chain = (
        chain_id.strip().upper()
        if chain_id
        else None
    )

    if normalized_chain is not None:
        if normalized_chain not in chains:
            raise ValueError(
                "Requested receptor chain was "
                f"not found: {normalized_chain}"
            )

        candidate_chains = {
            normalized_chain: chains[
                normalized_chain
            ]
        }
    else:
        candidate_chains = chains

    candidates = []

    for chain, chain_data in (
        candidate_chains.items()
    ):
        alignment = (
            smith_waterman_residue_map(
                query_sequence=(
                    chain_data["sequence"]
                ),
                reference_sequence=(
                    submitted
                ),
            )
        )

        candidates.append(
            {
                "chain": chain,
                "chain_data": chain_data,
                "alignment": alignment,
            }
        )

    selected = max(
        candidates,
        key=lambda candidate: (
            candidate["alignment"][
                "reference_coverage"
            ],
            candidate["alignment"][
                "identity"
            ],
            candidate["alignment"][
                "score"
            ],
            candidate["alignment"][
                "query_coverage"
            ],
        ),
    )

    selected_chain_data = selected[
        "chain_data"
    ]

    position_mapping = {}

    for (
        submitted_position,
        chain_sequence_position,
    ) in selected["alignment"][
        "reference_to_query"
    ].items():
        residue_index = (
            chain_sequence_position - 1
        )

        position_mapping[
            submitted_position
        ] = selected_chain_data[
            "residue_ids"
        ][residue_index]

    return {
        "selected_chain": (
            selected["chain"]
        ),
        "alignment": (
            selected["alignment"]
        ),
        "submitted_to_structure": (
            position_mapping
        ),
        "candidate_chains": [
            {
                "chain": (
                    candidate["chain"]
                ),
                "identity": (
                    candidate[
                        "alignment"
                    ]["identity"]
                ),
                "submitted_coverage": (
                    candidate[
                        "alignment"
                    ][
                        "reference_coverage"
                    ]
                ),
                "structure_coverage": (
                    candidate[
                        "alignment"
                    ]["query_coverage"]
                ),
                "score": (
                    candidate[
                        "alignment"
                    ]["score"]
                ),
            }
            for candidate in candidates
        ],
    }


def transfer_reference_residues_to_structure(
    *,
    reference_sequence: str,
    submitted_sequence: str,
    receptor_pdb: Path,
    reference_residues: Iterable[
        dict[str, Any]
    ],
    chain_id: str | None = None,
    reference_pdb: Path | None = None,
    reference_chain_id: str | None = None,
) -> dict[str, Any]:
    """Transfer homolog residue positions onto receptor PDB IDs.

    Residues may be specified in either of two explicit ways:

    - sequence_position: 1-based index in reference_sequence
    - residue_number: biological/PDB numbering, requiring reference_pdb

    The legacy field ``position`` remains an alias for
    ``sequence_position`` for backward compatibility.
    """

    reference = (
        normalize_protein_sequence(
            reference_sequence
        )
    )

    submitted = (
        normalize_protein_sequence(
            submitted_sequence
        )
    )

    residue_requests = list(
        reference_residues
    )

    homolog_alignment = (
        smith_waterman_residue_map(
            query_sequence=submitted,
            reference_sequence=reference,
        )
    )

    structure_mapping = (
        map_submitted_sequence_to_structure(
            submitted_sequence=submitted,
            receptor_pdb=receptor_pdb,
            chain_id=chain_id,
        )
    )

    needs_numbering_map = any(
        isinstance(raw_residue, dict)
        and "residue_number" in raw_residue
        for raw_residue in residue_requests
    )

    reference_numbering_map = None

    if (
        needs_numbering_map
        and reference_pdb is not None
    ):
        reference_numbering_map = (
            build_reference_numbering_map(
                reference_sequence=reference,
                reference_pdb=reference_pdb,
                chain_id=(
                    reference_chain_id
                ),
            )
        )

    reference_to_submitted = (
        homolog_alignment[
            "reference_to_query"
        ]
    )

    submitted_to_structure = (
        structure_mapping[
            "submitted_to_structure"
        ]
    )

    transfers = []
    mapped_structure_residues = []

    for raw_residue in residue_requests:
        if not isinstance(
            raw_residue,
            dict,
        ):
            raise ValueError(
                "Each reference residue must "
                "be a JSON-like object."
            )

        expected_amino_acid = str(
            raw_residue.get(
                "amino_acid"
            )
            or ""
        ).strip().upper()

        label = raw_residue.get(
            "label"
        )

        transfer = {
            "reference_numbering_mode": None,
            "reference_sequence_position": None,
            "reference_position": None,
            "reference_residue_number": (
                raw_residue.get(
                    "residue_number"
                )
            ),
            "reference_residue_id": None,
            "reference_chain": None,
            "expected_amino_acid": (
                expected_amino_acid
                or None
            ),
            "label": label,
            "submitted_position": None,
            "submitted_amino_acid": None,
            "structure_residue_id": None,
            "status": None,
            "selection_eligible": False,
        }

        reference_position = None

        if "residue_number" in raw_residue:
            transfer[
                "reference_numbering_mode"
            ] = "pdb_residue_number"

            if reference_numbering_map is None:
                transfer["status"] = (
                    "reference_numbering_map_required"
                )
                transfers.append(transfer)
                continue

            requested_chain = str(
                raw_residue.get("chain")
                or reference_chain_id
                or ""
            ).strip().upper()

            selected_reference_chain = str(
                reference_numbering_map[
                    "selected_chain"
                ]
            )

            if (
                requested_chain
                and requested_chain
                != selected_reference_chain
            ):
                transfer["status"] = (
                    "reference_chain_mismatch"
                )
                transfer[
                    "reference_chain"
                ] = selected_reference_chain
                transfers.append(transfer)
                continue

            try:
                normalized_number = (
                    _normalize_residue_number(
                        raw_residue[
                            "residue_number"
                        ],
                        insertion_code=(
                            raw_residue.get(
                                "insertion_code"
                            )
                        ),
                    )
                )
            except ValueError:
                transfer["status"] = (
                    "invalid_reference_residue_number"
                )
                transfers.append(transfer)
                continue

            transfer[
                "reference_residue_number"
            ] = normalized_number

            numbering_record = (
                reference_numbering_map[
                    "residues_by_number"
                ].get(normalized_number)
            )

            if numbering_record is None:
                transfer["status"] = (
                    "reference_residue_number_not_found"
                )
                transfers.append(transfer)
                continue

            reference_position = (
                numbering_record[
                    "reference_sequence_position"
                ]
            )

            transfer[
                "reference_residue_id"
            ] = numbering_record[
                "reference_residue_id"
            ]

            transfer[
                "reference_chain"
            ] = numbering_record[
                "reference_chain"
            ]

        elif "sequence_position" in raw_residue:
            transfer[
                "reference_numbering_mode"
            ] = "sequence_position"

            reference_position = int(
                raw_residue[
                    "sequence_position"
                ]
            )

        elif "position" in raw_residue:
            transfer[
                "reference_numbering_mode"
            ] = (
                "legacy_sequence_position"
            )

            reference_position = int(
                raw_residue["position"]
            )

        else:
            transfer["status"] = (
                "reference_position_missing"
            )
            transfers.append(transfer)
            continue

        transfer[
            "reference_sequence_position"
        ] = reference_position

        transfer[
            "reference_position"
        ] = reference_position

        if (
            reference_position < 1
            or reference_position
            > len(reference)
        ):
            transfer["status"] = (
                "reference_position_out_of_range"
            )
            transfers.append(transfer)
            continue

        reference_amino_acid = (
            reference[
                reference_position - 1
            ]
        )

        transfer[
            "reference_amino_acid"
        ] = reference_amino_acid

        if (
            expected_amino_acid
            and expected_amino_acid
            != reference_amino_acid
        ):
            transfer["status"] = (
                "reference_amino_acid_mismatch"
            )
            transfers.append(transfer)
            continue

        submitted_position = (
            reference_to_submitted.get(
                reference_position
            )
        )

        if submitted_position is None:
            transfer["status"] = (
                "not_aligned_to_submitted_sequence"
            )
            transfers.append(transfer)
            continue

        submitted_amino_acid = (
            submitted[
                submitted_position - 1
            ]
        )

        transfer[
            "submitted_position"
        ] = submitted_position

        transfer[
            "submitted_amino_acid"
        ] = submitted_amino_acid

        if (
            submitted_amino_acid
            != reference_amino_acid
        ):
            transfer["status"] = (
                "mapped_with_substitution"
            )
            transfers.append(transfer)
            continue

        structure_residue_id = (
            submitted_to_structure.get(
                submitted_position
            )
        )

        if structure_residue_id is None:
            transfer["status"] = (
                "not_mapped_to_receptor_structure"
            )
            transfers.append(transfer)
            continue

        structure_amino_acid = (
            AA3_TO_1.get(
                structure_residue_id.split(
                    ":",
                    1,
                )[0]
            )
        )

        if (
            structure_amino_acid
            != submitted_amino_acid
        ):
            transfer["status"] = (
                "structure_amino_acid_mismatch"
            )
            transfers.append(transfer)
            continue

        transfer[
            "structure_residue_id"
        ] = structure_residue_id

        transfer["status"] = (
            "mapped_conserved"
        )

        transfer[
            "selection_eligible"
        ] = True

        mapped_structure_residues.append(
            structure_residue_id
        )

        transfers.append(transfer)

    eligible_transfers = [
        transfer
        for transfer in transfers
        if transfer[
            "selection_eligible"
        ]
    ]

    status_counts = {}

    for transfer in transfers:
        status = str(
            transfer["status"]
        )

        status_counts[status] = (
            status_counts.get(
                status,
                0,
            )
            + 1
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "reference_to_submitted_alignment": (
            homolog_alignment
        ),
        "reference_numbering_map": (
            reference_numbering_map
        ),
        "submitted_to_structure_alignment": (
            structure_mapping[
                "alignment"
            ]
        ),
        "selected_structure_chain": (
            structure_mapping[
                "selected_chain"
            ]
        ),
        "structure_chain_candidates": (
            structure_mapping[
                "candidate_chains"
            ]
        ),
        "requested_residue_count": (
            len(transfers)
        ),
        "mapped_conserved_count": (
            len(eligible_transfers)
        ),
        "mapping_fraction": (
            len(eligible_transfers)
            / len(transfers)
            if transfers
            else 0.0
        ),
        "status_counts": status_counts,
        "residues": sorted(
            set(
                mapped_structure_residues
            )
        ),
        "transfers": transfers,
    }
