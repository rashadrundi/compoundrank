from __future__ import annotations

import argparse
import json
import re
import urllib.request
from pathlib import Path
from typing import Any, Sequence


SCHEMA_VERSION = (
    "uniprot_acquisition.v0.1"
)

REFERENCE_SCHEMA_VERSION = (
    "functional_site_reference.v0.1"
)

EXPERIMENTAL_EVIDENCE_CODES = {
    "ECO:0000269",
    "ECO:0007744",
}


def _location_value(
    feature: dict[str, Any],
    side: str,
) -> int | None:
    location = feature.get(
        "location"
    )

    if not isinstance(location, dict):
        return None

    boundary = location.get(side)

    if not isinstance(boundary, dict):
        return None

    value = boundary.get("value")

    if isinstance(value, bool):
        return None

    try:
        return int(value)
    except (
        TypeError,
        ValueError,
    ):
        return None


def _feature_evidence(
    feature: dict[str, Any],
) -> tuple[list[str], list[str]]:
    raw_evidences = feature.get(
        "evidences"
    )

    if not isinstance(
        raw_evidences,
        list,
    ):
        raw_evidences = []

    codes = sorted(
        {
            str(
                evidence.get(
                    "evidenceCode"
                )
                or ""
            )
            for evidence in raw_evidences
            if isinstance(
                evidence,
                dict,
            )
            and evidence.get(
                "evidenceCode"
            )
        }
    )

    sources = sorted(
        {
            str(
                evidence.get(
                    "source"
                )
                or ""
            )
            for evidence in raw_evidences
            if isinstance(
                evidence,
                dict,
            )
            and evidence.get("source")
        }
    )

    return codes, sources


def _ligand_block(
    feature: dict[str, Any],
) -> dict[str, Any] | None:
    ligand = feature.get("ligand")

    if not isinstance(ligand, dict):
        return None

    return {
        key: ligand.get(key)
        for key in (
            "name",
            "id",
            "label",
        )
        if ligand.get(key) is not None
    }


def _selected_feature_label(
    feature: dict[str, Any],
) -> tuple[str | None, str | None]:
    feature_type = str(
        feature.get("type")
        or ""
    )

    if feature_type == "Active site":
        return (
            "active_site",
            None,
        )

    if feature_type != "Binding site":
        return (
            None,
            "feature_type_not_selected",
        )

    ligand = _ligand_block(feature)

    ligand_name = str(
        (
            ligand or {}
        ).get("name")
        or ""
    ).strip().lower()

    if ligand_name == "substrate":
        return (
            "substrate_binding",
            None,
        )

    if ligand_name:
        return (
            None,
            "non_substrate_binding_site:"
            + ligand_name,
        )

    return (
        None,
        "binding_site_without_substrate_label",
    )


def _expand_feature_positions(
    feature: dict[str, Any],
    *,
    maximum_span: int = 3,
) -> tuple[list[int], str | None]:
    start = _location_value(
        feature,
        "start",
    )

    end = _location_value(
        feature,
        "end",
    )

    if start is None or end is None:
        return (
            [],
            "feature_location_missing",
        )

    if start < 1 or end < start:
        return (
            [],
            "feature_location_invalid",
        )

    span = end - start + 1

    if span > maximum_span:
        return (
            [],
            "feature_span_exceeds_limit",
        )

    return (
        list(
            range(
                start,
                end + 1,
            )
        ),
        None,
    )


def _parse_resolution(
    value: object,
) -> float | None:
    text = str(
        value or ""
    ).strip()

    match = re.search(
        r"(\d+(?:\.\d+)?)",
        text,
    )

    if match is None:
        return None

    return float(
        match.group(1)
    )


def _parse_chain_ranges(
    value: object,
) -> list[dict[str, Any]]:
    text = str(
        value or ""
    ).strip()

    if not text:
        return []

    mappings = []

    for group in text.split(","):
        group = group.strip()

        match = re.fullmatch(
            r"([^=]+)="
            r"(\d+)-(\d+)",
            group,
        )

        if match is None:
            continue

        chains = [
            chain.strip().upper()
            for chain in (
                match.group(1)
                .replace(";", "/")
                .split("/")
            )
            if chain.strip()
        ]

        start = int(match.group(2))
        end = int(match.group(3))

        if end < start:
            continue

        for chain in chains:
            mappings.append(
                {
                    "chain": chain,
                    "uniprot_start": (
                        start
                    ),
                    "uniprot_end": end,
                    "mapped_length": (
                        end - start + 1
                    ),
                }
            )

    return mappings


def parse_pdb_candidates(
    payload: dict[str, Any],
    *,
    selected_positions: Sequence[int],
) -> list[dict[str, Any]]:
    cross_references = payload.get(
        "uniProtKBCrossReferences"
    )

    if not isinstance(
        cross_references,
        list,
    ):
        cross_references = []

    candidates = []

    for cross_reference in (
        cross_references
    ):
        if not isinstance(
            cross_reference,
            dict,
        ):
            continue

        if (
            cross_reference.get(
                "database"
            )
            != "PDB"
        ):
            continue

        properties = {}

        for property_record in (
            cross_reference.get(
                "properties"
            )
            or []
        ):
            if not isinstance(
                property_record,
                dict,
            ):
                continue

            key = str(
                property_record.get("key")
                or ""
            )

            if key:
                properties[key] = (
                    property_record.get(
                        "value"
                    )
                )

        chain_ranges = (
            _parse_chain_ranges(
                properties.get("Chains")
            )
        )

        covered_positions = sorted(
            {
                position
                for position
                in selected_positions
                if any(
                    mapping[
                        "uniprot_start"
                    ]
                    <= position
                    <= mapping[
                        "uniprot_end"
                    ]
                    for mapping
                    in chain_ranges
                )
            }
        )

        method = str(
            properties.get("Method")
            or ""
        )

        resolution = (
            _parse_resolution(
                properties.get(
                    "Resolution"
                )
            )
        )

        candidates.append(
            {
                "pdb_id": str(
                    cross_reference.get(
                        "id"
                    )
                    or ""
                ).upper(),
                "method": method,
                "resolution_angstrom": (
                    resolution
                ),
                "chain_ranges": (
                    chain_ranges
                ),
                "functional_site_positions_covered": (
                    covered_positions
                ),
                "functional_site_coverage_count": (
                    len(covered_positions)
                ),
                "functional_site_coverage_fraction": (
                    len(covered_positions)
                    / len(selected_positions)
                    if selected_positions
                    else 0.0
                ),
                "properties": properties,
            }
        )

    def method_rank(
        candidate: dict[str, Any],
    ) -> int:
        method = str(
            candidate.get("method")
            or ""
        ).lower()

        if "x-ray" in method:
            return 0

        if "electron microscopy" in method:
            return 1

        if "nmr" in method:
            return 2

        return 3

    candidates.sort(
        key=lambda candidate: (
            -int(
                candidate[
                    "functional_site_coverage_count"
                ]
            ),
            method_rank(candidate),
            (
                float(
                    candidate[
                        "resolution_angstrom"
                    ]
                )
                if candidate[
                    "resolution_angstrom"
                ]
                is not None
                else float("inf")
            ),
            -max(
                (
                    int(
                        mapping[
                            "mapped_length"
                        ]
                    )
                    for mapping
                    in candidate[
                        "chain_ranges"
                    ]
                ),
                default=0,
            ),
            candidate["pdb_id"],
        )
    )

    for rank, candidate in enumerate(
        candidates,
        start=1,
    ):
        candidate["rank"] = rank

    return candidates


def build_uniprot_acquisition(
    payload: dict[str, Any],
    *,
    requested_selection_mode: str = (
        "prioritize_supported"
    ),
) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise ValueError(
            "UniProt entry must be a JSON object."
        )

    accession = str(
        payload.get("primaryAccession")
        or ""
    ).strip()

    if not accession:
        raise ValueError(
            "UniProt entry has no primary accession."
        )

    sequence_block = payload.get(
        "sequence"
    )

    if not isinstance(
        sequence_block,
        dict,
    ):
        raise ValueError(
            "UniProt entry has no sequence block."
        )

    sequence = "".join(
        character
        for character in str(
            sequence_block.get("value")
            or ""
        ).upper()
        if character.isalpha()
    )

    if not sequence:
        raise ValueError(
            "UniProt entry has no sequence."
        )

    entry_type = str(
        payload.get("entryType")
        or ""
    )

    reviewed = (
        "reviewed" in entry_type.lower()
        and "unreviewed"
        not in entry_type.lower()
    )

    features = payload.get("features")

    if not isinstance(features, list):
        features = []

    residue_records: dict[
        int,
        dict[str, Any],
    ] = {}

    excluded_features = []

    for feature_index, feature in enumerate(
        features,
        start=1,
    ):
        if not isinstance(feature, dict):
            continue

        feature_type = str(
            feature.get("type")
            or ""
        )

        if feature_type not in {
            "Active site",
            "Binding site",
        }:
            continue

        label, exclusion_reason = (
            _selected_feature_label(
                feature
            )
        )

        positions, location_error = (
            _expand_feature_positions(
                feature
            )
        )

        ligand = _ligand_block(feature)
        codes, sources = (
            _feature_evidence(feature)
        )

        feature_summary = {
            "feature_index": (
                feature_index
            ),
            "feature_type": (
                feature_type
            ),
            "start": _location_value(
                feature,
                "start",
            ),
            "end": _location_value(
                feature,
                "end",
            ),
            "description": str(
                feature.get(
                    "description"
                )
                or ""
            ),
            "ligand": ligand,
            "evidence_codes": codes,
            "evidence_sources": (
                sources
            ),
        }

        if exclusion_reason is not None:
            feature_summary["reason"] = (
                exclusion_reason
            )

            excluded_features.append(
                feature_summary
            )
            continue

        if location_error is not None:
            feature_summary["reason"] = (
                location_error
            )

            excluded_features.append(
                feature_summary
            )
            continue

        for position in positions:
            if position > len(sequence):
                excluded_features.append(
                    {
                        **feature_summary,
                        "reason": (
                            "position_outside_sequence"
                        ),
                        "position": (
                            position
                        ),
                    }
                )
                continue

            amino_acid = sequence[
                position - 1
            ]

            existing = (
                residue_records.get(
                    position
                )
            )

            if existing is None:
                residue_records[
                    position
                ] = {
                    "sequence_position": (
                        position
                    ),
                    "amino_acid": (
                        amino_acid
                    ),
                    "label": label,
                    "feature_types": [
                        feature_type
                    ],
                    "descriptions": (
                        [
                            feature_summary[
                                "description"
                            ]
                        ]
                        if feature_summary[
                            "description"
                        ]
                        else []
                    ),
                    "ligands": (
                        [ligand]
                        if ligand
                        else []
                    ),
                    "evidence_codes": (
                        list(codes)
                    ),
                    "evidence_sources": (
                        list(sources)
                    ),
                }
                continue

            existing[
                "feature_types"
            ] = sorted(
                {
                    *existing[
                        "feature_types"
                    ],
                    feature_type,
                }
            )

            existing[
                "evidence_codes"
            ] = sorted(
                {
                    *existing[
                        "evidence_codes"
                    ],
                    *codes,
                }
            )

            existing[
                "evidence_sources"
            ] = sorted(
                {
                    *existing[
                        "evidence_sources"
                    ],
                    *sources,
                }
            )

    residues = [
        residue_records[position]
        for position in sorted(
            residue_records
        )
    ]

    selected_positions = [
        int(
            residue[
                "sequence_position"
            ]
        )
        for residue in residues
    ]

    pdb_candidates = (
        parse_pdb_candidates(
            payload,
            selected_positions=(
                selected_positions
            ),
        )
    )

    has_experimental_evidence = any(
        EXPERIMENTAL_EVIDENCE_CODES
        & set(
            residue[
                "evidence_codes"
            ]
        )
        for residue in residues
    )

    if (
        reviewed
        and len(residues) >= 3
        and has_experimental_evidence
    ):
        confidence = "high"
    elif reviewed and len(residues) >= 2:
        confidence = "moderate"
    else:
        confidence = "low"

    selection_allowed = (
        requested_selection_mode
        == "prioritize_supported"
        and reviewed
        and len(residues) >= 2
    )

    effective_selection_mode = (
        "prioritize_supported"
        if selection_allowed
        else "report_only"
    )

    organism = payload.get(
        "organism"
    )

    if not isinstance(organism, dict):
        organism = {}

    source = {
        "database": "UniProtKB",
        "primary_accession": (
            accession
        ),
        "uniprot_id": payload.get(
            "uniProtkbId"
        ),
        "entry_type": entry_type,
        "reviewed": reviewed,
        "annotation_score": (
            payload.get(
                "annotationScore"
            )
        ),
        "organism": {
            "scientific_name": (
                organism.get(
                    "scientificName"
                )
            ),
            "taxon_id": organism.get(
                "taxonId"
            ),
        },
        "coordinate_system": (
            "UniProt full-sequence positions"
        ),
        "feature_policy": {
            "included": [
                "Active site",
                (
                    "Binding site with "
                    "ligand name substrate"
                ),
            ],
            "excluded_from_selection": [
                (
                    "Non-substrate ligand "
                    "binding sites"
                ),
                (
                    "Long or unresolved "
                    "feature ranges"
                ),
            ],
        },
    }

    notes = [
        (
            "Residue coordinates are UniProt "
            "full-sequence positions, not PDB "
            "author residue numbers."
        ),
        (
            "Non-substrate ligand-binding sites "
            "were preserved in acquisition "
            "metadata but excluded from the "
            "selection residue list."
        ),
    ]

    if (
        requested_selection_mode
        == "prioritize_supported"
        and not selection_allowed
    ):
        notes.append(
            "Requested prioritization was "
            "downgraded to report_only because "
            "the entry was unreviewed or fewer "
            "than two usable functional residues "
            "were found."
        )

    reference_record = {
        "schema_version": (
            REFERENCE_SCHEMA_VERSION
        ),
        "evidence_id": (
            f"uniprot_{accession}_"
            "functional_sites"
        ),
        "selection_mode": (
            effective_selection_mode
        ),
        "confidence": confidence,
        "residues": residues,
        "source": source,
        "notes": notes,
    }

    summary = {
        "schema_version": (
            SCHEMA_VERSION
        ),
        "primary_accession": (
            accession
        ),
        "entry_type": entry_type,
        "reviewed": reviewed,
        "sequence_length": (
            len(sequence)
        ),
        "requested_selection_mode": (
            requested_selection_mode
        ),
        "effective_selection_mode": (
            effective_selection_mode
        ),
        "confidence": confidence,
        "selected_residue_count": (
            len(residues)
        ),
        "selected_sequence_positions": (
            selected_positions
        ),
        "excluded_feature_count": (
            len(excluded_features)
        ),
        "excluded_features": (
            excluded_features
        ),
        "pdb_candidate_count": (
            len(pdb_candidates)
        ),
        "selected_pdb_candidate": (
            pdb_candidates[0]
            if pdb_candidates
            else None
        ),
    }

    return {
        "entry": payload,
        "sequence": sequence,
        "reference_record": (
            reference_record
        ),
        "pdb_candidates": (
            pdb_candidates
        ),
        "summary": summary,
    }


def fetch_uniprot_entry(
    accession: str,
    *,
    timeout_seconds: float = 60.0,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
]:
    normalized_accession = str(
        accession or ""
    ).strip().upper()

    if not normalized_accession:
        raise ValueError(
            "UniProt accession is empty."
        )

    request = urllib.request.Request(
        (
            "https://rest.uniprot.org/"
            "uniprotkb/"
            f"{normalized_accession}.json"
        ),
        headers={
            "Accept": "application/json",
            "User-Agent": (
                "CompoundRank-EXORCIST/0.1"
            ),
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=timeout_seconds,
    ) as response:
        raw = response.read()

        payload = json.loads(
            raw.decode("utf-8")
        )

        metadata = {
            "url": response.geturl(),
            "status": getattr(
                response,
                "status",
                None,
            ),
            "content_type": (
                response.headers.get(
                    "Content-Type"
                )
            ),
            "uniprot_release": (
                response.headers.get(
                    "X-UniProt-Release"
                )
            ),
            "uniprot_release_date": (
                response.headers.get(
                    "X-UniProt-Release-Date"
                )
            ),
            "etag": response.headers.get(
                "ETag"
            ),
        }

    return payload, metadata


def write_acquisition_outputs(
    output_dir: Path,
    *,
    payload: dict[str, Any],
    response_metadata: (
        dict[str, Any] | None
    ) = None,
    requested_selection_mode: str = (
        "prioritize_supported"
    ),
) -> dict[str, str]:
    result = build_uniprot_acquisition(
        payload,
        requested_selection_mode=(
            requested_selection_mode
        ),
    )

    output = Path(output_dir)

    output.mkdir(
        parents=True,
        exist_ok=True,
    )

    accession = result[
        "summary"
    ]["primary_accession"]

    paths = {
        "entry": (
            output / "uniprot_entry.json"
        ),
        "fasta": (
            output / f"{accession}.fasta"
        ),
        "reference_record": (
            output
            / "functional_site_reference.json"
        ),
        "pdb_candidates": (
            output / "pdb_candidates.json"
        ),
        "summary": (
            output / "acquisition_summary.json"
        ),
    }

    paths["entry"].write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    sequence = result["sequence"]

    paths["fasta"].write_text(
        (
            f">{accession}|"
            f"{payload.get('uniProtkbId') or ''}\n"
        )
        + "\n".join(
            sequence[
                index:index + 60
            ]
            for index in range(
                0,
                len(sequence),
                60,
            )
        )
        + "\n",
        encoding="utf-8",
    )

    paths[
        "reference_record"
    ].write_text(
        json.dumps(
            result[
                "reference_record"
            ],
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    paths[
        "pdb_candidates"
    ].write_text(
        json.dumps(
            {
                "schema_version": (
                    SCHEMA_VERSION
                ),
                "primary_accession": (
                    accession
                ),
                "candidates": result[
                    "pdb_candidates"
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    summary = dict(
        result["summary"]
    )

    summary[
        "response_metadata"
    ] = response_metadata or {}

    paths["summary"].write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        key: str(path)
        for key, path in paths.items()
    }


def build_cli_parser() -> (
    argparse.ArgumentParser
):
    parser = argparse.ArgumentParser(
        description=(
            "Acquire a UniProt entry and build "
            "functional_site_reference.v0.1."
        )
    )

    parser.add_argument(
        "--accession",
        default=None,
    )

    parser.add_argument(
        "--input-json",
        type=Path,
        default=None,
        help=(
            "Use a previously downloaded "
            "UniProt JSON record."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--selection-mode",
        choices=(
            "report_only",
            "prioritize_supported",
        ),
        default="prioritize_supported",
    )

    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=60.0,
    )

    return parser


def main(
    argv: Sequence[str] | None = None,
) -> int:
    parser = build_cli_parser()

    arguments = parser.parse_args(
        argv
    )

    if (
        arguments.input_json is None
        and not arguments.accession
    ):
        parser.error(
            "Provide --accession or "
            "--input-json."
        )

    if (
        arguments.input_json is not None
        and arguments.accession
    ):
        parser.error(
            "Use only one of --accession "
            "or --input-json."
        )

    response_metadata = {}

    if arguments.input_json is not None:
        payload = json.loads(
            arguments.input_json.read_text(
                encoding="utf-8"
            )
        )

        response_metadata = {
            "source": "input_json",
            "input_path": str(
                arguments.input_json
            ),
        }
    else:
        payload, response_metadata = (
            fetch_uniprot_entry(
                arguments.accession,
                timeout_seconds=(
                    arguments.timeout_seconds
                ),
            )
        )

    paths = write_acquisition_outputs(
        arguments.output_dir,
        payload=payload,
        response_metadata=(
            response_metadata
        ),
        requested_selection_mode=(
            arguments.selection_mode
        ),
    )

    summary = json.loads(
        Path(
            paths["summary"]
        ).read_text(
            encoding="utf-8"
        )
    )

    print(
        json.dumps(
            {
                "primary_accession": (
                    summary[
                        "primary_accession"
                    ]
                ),
                "reviewed": (
                    summary["reviewed"]
                ),
                "sequence_length": (
                    summary[
                        "sequence_length"
                    ]
                ),
                "selection_mode": (
                    summary[
                        "effective_selection_mode"
                    ]
                ),
                "confidence": (
                    summary["confidence"]
                ),
                "selected_residue_count": (
                    summary[
                        "selected_residue_count"
                    ]
                ),
                "selected_sequence_positions": (
                    summary[
                        "selected_sequence_positions"
                    ]
                ),
                "selected_pdb_candidate": (
                    summary[
                        "selected_pdb_candidate"
                    ]
                ),
                "outputs": paths,
            },
            indent=2,
            sort_keys=True,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
