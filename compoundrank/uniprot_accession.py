from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


SCHEMA_VERSION = (
    "uniprot_accession_resolution.v0.1"
)

_ACCESSION_BODY = (
    r"(?:"
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|"
    r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]"
    r"|"
    r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]"
    r"[A-Z0-9]{3}[0-9]"
    r")"
)

_ACCESSION_PATTERN = re.compile(
    rf"(?<![A-Z0-9])"
    rf"(?P<accession>{_ACCESSION_BODY})"
    rf"(?:-\d+)?"
    rf"(?![A-Z0-9])",
    re.IGNORECASE,
)

_PIPE_PATTERN = re.compile(
    rf"(?:^|\s)(?:sp|tr)\|"
    rf"(?P<accession>{_ACCESSION_BODY})"
    rf"(?:-\d+)?\|",
    re.IGNORECASE,
)

_LABELED_PATTERN = re.compile(
    rf"(?:"
    rf"uniprot(?:kb)?"
    rf"|accession"
    rf")"
    rf"\s*[:=|]\s*"
    rf"(?P<accession>{_ACCESSION_BODY})"
    rf"(?:-\d+)?",
    re.IGNORECASE,
)


def normalize_uniprot_accession(
    value: str,
) -> str:
    normalized = str(
        value or ""
    ).strip().upper()

    if "-" in normalized:
        base, suffix = (
            normalized.split(
                "-",
                1,
            )
        )

        if suffix.isdigit():
            normalized = base

    if not re.fullmatch(
        _ACCESSION_BODY,
        normalized,
        re.IGNORECASE,
    ):
        raise ValueError(
            "Value is not a recognized "
            f"UniProt accession: {value!r}"
        )

    return normalized


def _read_single_fasta_header(
    path: Path,
) -> str:
    source = Path(path)

    if not source.is_file():
        raise FileNotFoundError(
            source
        )

    headers: list[str] = []

    with source.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                headers.append(
                    line[1:].strip()
                )

    if not headers:
        raise ValueError(
            f"FASTA has no header: {source}"
        )

    if len(headers) != 1:
        raise ValueError(
            "Automatic UniProt accession "
            "resolution requires a "
            "single-record FASTA."
        )

    return headers[0]


def _add_matches(
    candidates: dict[
        str,
        dict[str, Any],
    ],
    *,
    text: str,
    pattern: re.Pattern[str],
    source: str,
    priority: int,
) -> None:
    for match in pattern.finditer(
        text
    ):
        accession = (
            normalize_uniprot_accession(
                match.group(
                    "accession"
                )
            )
        )

        candidate = candidates.setdefault(
            accession,
            {
                "accession": accession,
                "sources": [],
                "highest_priority": (
                    priority
                ),
            },
        )

        candidate[
            "highest_priority"
        ] = max(
            int(
                candidate[
                    "highest_priority"
                ]
            ),
            priority,
        )

        source_record = {
            "source": source,
            "matched_text": (
                match.group(0)
            ),
            "priority": priority,
        }

        if (
            source_record
            not in candidate["sources"]
        ):
            candidate[
                "sources"
            ].append(
                source_record
            )


def resolve_uniprot_accession_from_fasta(
    fasta_path: Path,
) -> dict[str, Any]:
    source = Path(fasta_path)
    header = (
        _read_single_fasta_header(
            source
        )
    )

    candidates: dict[
        str,
        dict[str, Any],
    ] = {}

    _add_matches(
        candidates,
        text=header,
        pattern=_PIPE_PATTERN,
        source=(
            "canonical_uniprot_header"
        ),
        priority=400,
    )

    _add_matches(
        candidates,
        text=header,
        pattern=_LABELED_PATTERN,
        source="labeled_header",
        priority=300,
    )

    _add_matches(
        candidates,
        text=header,
        pattern=_ACCESSION_PATTERN,
        source="header_token",
        priority=200,
    )

    _add_matches(
        candidates,
        text=source.stem,
        pattern=_ACCESSION_PATTERN,
        source="filename_token",
        priority=100,
    )

    if not candidates:
        raise ValueError(
            "No strict UniProt accession "
            "was found in the FASTA header "
            "or filename. Supply "
            "--reference-uniprot-accession."
        )

    highest_priority = max(
        int(
            candidate[
                "highest_priority"
            ]
        )
        for candidate
        in candidates.values()
    )

    strongest = [
        candidate
        for candidate
        in candidates.values()
        if int(
            candidate[
                "highest_priority"
            ]
        )
        == highest_priority
    ]

    if len(strongest) != 1:
        raise ValueError(
            "Automatic UniProt accession "
            "resolution was ambiguous: "
            + ", ".join(
                sorted(
                    candidate[
                        "accession"
                    ]
                    for candidate
                    in strongest
                )
            )
        )

    selected = strongest[0]

    methods = sorted(
        {
            str(record["source"])
            for record
            in selected["sources"]
            if int(
                record["priority"]
            )
            == highest_priority
        }
    )

    return {
        "schema_version": (
            SCHEMA_VERSION
        ),
        "status": "resolved",
        "selected_accession": (
            selected["accession"]
        ),
        "resolution_method": (
            methods[0]
        ),
        "source_fasta": str(
            source
        ),
        "fasta_header": header,
        "candidate_count": len(
            candidates
        ),
        "candidates": sorted(
            candidates.values(),
            key=lambda candidate: (
                -int(
                    candidate[
                        "highest_priority"
                    ]
                ),
                str(
                    candidate[
                        "accession"
                    ]
                ),
            ),
        ),
        "safety_checks": {
            "strict_accession_format": (
                True
            ),
            "single_highest_priority_candidate": (
                True
            ),
            "single_fasta_record": True,
            "downstream_sequence_validation_required": (
                True
            ),
        },
    }


def write_uniprot_accession_resolution(
    path: Path,
    resolution: dict[str, Any],
) -> Path:
    destination = Path(path)

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    destination.write_text(
        json.dumps(
            resolution,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    return destination
