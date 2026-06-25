from __future__ import annotations

import csv
import io
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from .chembl_search import (
    local_alignment_metrics,
    read_fasta_sequence,
)
from .uniprot_acquisition import (
    fetch_uniprot_entry,
)


SCHEMA_VERSION = (
    "uniprot_accession_resolution.v0.2"
)

EBI_BLAST_BASE_URL = (
    "https://www.ebi.ac.uk/Tools/"
    "services/rest/ncbiblast"
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

    normalized = re.sub(
        r"\.\d+$",
        "",
        normalized,
    )

    if "-" in normalized:
        base, suffix = normalized.split(
            "-",
            1,
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


def _request_text(
    request: urllib.request.Request,
    *,
    timeout_seconds: float,
) -> tuple[
    str,
    dict[str, Any],
]:
    with urllib.request.urlopen(
        request,
        timeout=min(
            float(timeout_seconds),
            60.0,
        ),
    ) as response:
        raw = response.read()

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
        }

    return (
        raw.decode(
            "utf-8",
            errors="replace",
        ),
        metadata,
    )


def _submit_ebi_blast(
    sequence: str,
    *,
    email: str,
    timeout_seconds: float,
    maximum_hits: int,
) -> tuple[
    str,
    dict[str, Any],
]:
    normalized_email = str(
        email or ""
    ).strip()

    if (
        not normalized_email
        or "@" not in normalized_email
    ):
        raise ValueError(
            "Sequence-based UniProt discovery "
            "requires a valid Job Dispatcher "
            "contact email."
        )

    normalized_sequence = re.sub(
        r"[^A-Za-z]",
        "",
        str(sequence or ""),
    ).upper()

    if not normalized_sequence:
        raise ValueError(
            "Protein sequence is empty."
        )

    payload = urllib.parse.urlencode(
        {
            "email": normalized_email,
            "title": (
                "EXORCIST UniProt "
                "accession discovery"
            ),
            "program": "blastp",
            "database": "uniprotkb",
            "stype": "protein",
            "sequence": (
                ">exorcist_query\n"
                f"{normalized_sequence}\n"
            ),
            "matrix": "BLOSUM62",
            "exp": "1e-10",
            "filter": "F",
            "gapalign": "true",
            "scores": str(
                maximum_hits
            ),
            "alignments": str(
                maximum_hits
            ),
        }
    ).encode(
        "utf-8"
    )

    request = urllib.request.Request(
        f"{EBI_BLAST_BASE_URL}/run",
        data=payload,
        headers={
            "Content-Type": (
                "application/"
                "x-www-form-urlencoded"
            ),
            "Accept": "text/plain",
            "User-Agent": (
                "CompoundRank-EXORCIST/0.1"
            ),
        },
        method="POST",
    )

    job_id, metadata = _request_text(
        request,
        timeout_seconds=timeout_seconds,
    )

    normalized_job_id = job_id.strip()

    if not normalized_job_id:
        raise RuntimeError(
            "EMBL-EBI BLAST submission "
            "returned no job identifier."
        )

    return normalized_job_id, metadata


def _poll_ebi_blast(
    job_id: str,
    *,
    timeout_seconds: float,
    poll_interval_seconds: float,
) -> dict[str, Any]:
    if timeout_seconds <= 0:
        raise ValueError(
            "Sequence-search timeout must "
            "be greater than zero."
        )

    if poll_interval_seconds <= 0:
        raise ValueError(
            "Sequence-search polling interval "
            "must be greater than zero."
        )

    started = time.monotonic()
    history: list[str] = []

    while True:
        request = urllib.request.Request(
            (
                f"{EBI_BLAST_BASE_URL}/"
                f"status/{job_id}"
            ),
            headers={
                "Accept": "text/plain",
                "User-Agent": (
                    "CompoundRank-"
                    "EXORCIST/0.1"
                ),
            },
        )

        status_text, _ = _request_text(
            request,
            timeout_seconds=(
                timeout_seconds
            ),
        )

        status = status_text.strip().upper()

        if (
            not history
            or history[-1] != status
        ):
            history.append(status)

        if status == "FINISHED":
            return {
                "job_id": job_id,
                "status": status,
                "status_history": history,
                "elapsed_seconds": (
                    time.monotonic()
                    - started
                ),
            }

        if status in {
            "ERROR",
            "FAILURE",
            "NOT_FOUND",
        }:
            raise RuntimeError(
                "EMBL-EBI BLAST job "
                f"{job_id} ended with "
                f"status {status}."
            )

        if status not in {
            "QUEUED",
            "RUNNING",
            "PENDING",
        }:
            raise RuntimeError(
                "EMBL-EBI BLAST returned "
                f"unknown status {status!r}."
            )

        elapsed = (
            time.monotonic()
            - started
        )

        if elapsed >= timeout_seconds:
            raise TimeoutError(
                "EMBL-EBI BLAST did not "
                "finish before the "
                f"{timeout_seconds:g}-second "
                "timeout."
            )

        time.sleep(
            min(
                poll_interval_seconds,
                max(
                    0.0,
                    timeout_seconds
                    - elapsed,
                ),
            )
        )


def _fetch_ebi_blast_xml(
    job_id: str,
    *,
    timeout_seconds: float,
) -> tuple[
    str,
    dict[str, Any],
]:
    request = urllib.request.Request(
        (
            f"{EBI_BLAST_BASE_URL}/"
            f"result/{job_id}/xml"
        ),
        headers={
            "Accept": "application/xml",
            "User-Agent": (
                "CompoundRank-EXORCIST/0.1"
            ),
        },
    )

    return _request_text(
        request,
        timeout_seconds=timeout_seconds,
    )


def _fetch_ebi_blast_tsv(
    job_id: str,
    *,
    timeout_seconds: float,
) -> tuple[
    str,
    dict[str, Any],
]:
    request = urllib.request.Request(
        (
            f"{EBI_BLAST_BASE_URL}/"
            f"result/{job_id}/tsv"
        ),
        headers={
            "Accept": (
                "text/tab-separated-values"
            ),
            "User-Agent": (
                "CompoundRank-EXORCIST/0.1"
            ),
        },
    )

    return _request_text(
        request,
        timeout_seconds=timeout_seconds,
    )


def _safe_float(
    value: Any,
    *,
    default: float = 0.0,
) -> float:
    try:
        return float(
            str(value or "").strip()
        )
    except (
        TypeError,
        ValueError,
    ):
        return default


def _safe_int(
    value: Any,
    *,
    default: int = 0,
) -> int:
    try:
        return int(
            float(
                str(value or "").strip()
            )
        )
    except (
        TypeError,
        ValueError,
    ):
        return default


def parse_ebi_blast_tsv(
    tsv_text: str,
) -> list[dict[str, Any]]:
    reader = csv.DictReader(
        io.StringIO(
            str(tsv_text or "")
        ),
        delimiter="\t",
    )

    candidates: dict[
        str,
        dict[str, Any],
    ] = {}

    for raw_row in reader:
        row = {
            str(key or "").strip(): (
                str(value or "").strip()
            )
            for key, value
            in raw_row.items()
        }

        raw_accession = row.get(
            "Accession",
            "",
        )

        if not raw_accession:
            continue

        try:
            accession = (
                normalize_uniprot_accession(
                    raw_accession
                )
            )
        except ValueError:
            continue

        bit_score = _safe_float(
            row.get(
                "Score(Bits)"
            )
        )

        identity_percent = _safe_float(
            row.get(
                "Identities(%)"
            )
        )

        positives_percent = _safe_float(
            row.get(
                "Positives(%)"
            )
        )

        evalue = _safe_float(
            row.get(
                "E()"
            ),
            default=float(
                "inf"
            ),
        )

        target_length = _safe_int(
            row.get(
                "Length"
            )
        )

        rank = _safe_int(
            row.get(
                "Hit"
            )
        )

        database = row.get(
            "DB",
            "",
        )

        candidate = {
            "accession": accession,
            "hit_id": (
                f"{database}:{accession}"
                if database
                else accession
            ),
            "hit_definition": row.get(
                "Description",
                "",
            ),
            "reported_accession": (
                raw_accession
            ),
            "best_hsp": {
                "aligned_length": 0,
                "identities": 0,
                "identity": (
                    identity_percent
                    / 100.0
                ),
                "bit_score": bit_score,
                "evalue": evalue,
                "query_start": 0,
                "query_end": 0,
                "target_start": 0,
                "target_end": 0,
            },
            "hsp_count": 1,
            "remote_search": {
                "rank": rank,
                "database": database,
                "organism": row.get(
                    "Organism",
                    "",
                ),
                "target_length": (
                    target_length
                ),
                "identity_percent": (
                    identity_percent
                ),
                "positives_percent": (
                    positives_percent
                ),
            },
        }

        existing = candidates.get(
            accession
        )

        if (
            existing is None
            or bit_score
            > float(
                existing[
                    "best_hsp"
                ]["bit_score"]
            )
        ):
            candidates[
                accession
            ] = candidate

    return sorted(
        candidates.values(),
        key=lambda candidate: (
            -float(
                candidate[
                    "best_hsp"
                ]["bit_score"]
            ),
            float(
                candidate[
                    "best_hsp"
                ]["evalue"]
            ),
            int(
                candidate.get(
                    "remote_search",
                    {},
                ).get(
                    "rank",
                    0,
                )
                or 0
            ),
            str(
                candidate[
                    "accession"
                ]
            ),
        ),
    )


def _extract_accession_from_hit(
    *values: str,
) -> str | None:
    for value in values:
        text = str(
            value or ""
        )

        for pattern in (
            _PIPE_PATTERN,
            _ACCESSION_PATTERN,
        ):
            match = pattern.search(text)

            if match is None:
                continue

            try:
                return (
                    normalize_uniprot_accession(
                        match.group(
                            "accession"
                        )
                    )
                )
            except ValueError:
                continue

    return None


def parse_ebi_blast_xml(
    xml_text: str,
) -> list[dict[str, Any]]:
    try:
        root = ET.fromstring(
            xml_text
        )
    except ET.ParseError as error:
        raise ValueError(
            "Could not parse EMBL-EBI "
            "BLAST XML."
        ) from error

    candidates: dict[
        str,
        dict[str, Any],
    ] = {}

    for hit in root.findall(
        ".//Hit"
    ):
        hit_id = (
            hit.findtext(
                "Hit_id"
            )
            or ""
        )

        hit_definition = (
            hit.findtext(
                "Hit_def"
            )
            or ""
        )

        hit_accession = (
            hit.findtext(
                "Hit_accession"
            )
            or ""
        )

        accession = (
            _extract_accession_from_hit(
                hit_accession,
                hit_id,
                hit_definition,
            )
        )

        if accession is None:
            continue

        hsps: list[dict[str, Any]] = []

        for hsp in hit.findall(
            "./Hit_hsps/Hsp"
        ):
            aligned_length = int(
                hsp.findtext(
                    "Hsp_align-len"
                )
                or 0
            )

            identities = int(
                hsp.findtext(
                    "Hsp_identity"
                )
                or 0
            )

            bit_score = float(
                hsp.findtext(
                    "Hsp_bit-score"
                )
                or 0.0
            )

            evalue = float(
                hsp.findtext(
                    "Hsp_evalue"
                )
                or "inf"
            )

            hsps.append(
                {
                    "aligned_length": (
                        aligned_length
                    ),
                    "identities": (
                        identities
                    ),
                    "identity": (
                        identities
                        / aligned_length
                        if aligned_length
                        else 0.0
                    ),
                    "bit_score": (
                        bit_score
                    ),
                    "evalue": evalue,
                    "query_start": int(
                        hsp.findtext(
                            "Hsp_query-from"
                        )
                        or 0
                    ),
                    "query_end": int(
                        hsp.findtext(
                            "Hsp_query-to"
                        )
                        or 0
                    ),
                    "target_start": int(
                        hsp.findtext(
                            "Hsp_hit-from"
                        )
                        or 0
                    ),
                    "target_end": int(
                        hsp.findtext(
                            "Hsp_hit-to"
                        )
                        or 0
                    ),
                }
            )

        best_hsp = max(
            hsps,
            key=lambda record: (
                float(
                    record[
                        "bit_score"
                    ]
                ),
                -float(
                    record[
                        "evalue"
                    ]
                ),
                int(
                    record[
                        "aligned_length"
                    ]
                ),
            ),
            default={
                "aligned_length": 0,
                "identities": 0,
                "identity": 0.0,
                "bit_score": 0.0,
                "evalue": float(
                    "inf"
                ),
                "query_start": 0,
                "query_end": 0,
                "target_start": 0,
                "target_end": 0,
            },
        )

        candidate = {
            "accession": accession,
            "hit_id": hit_id,
            "hit_definition": (
                hit_definition
            ),
            "reported_accession": (
                hit_accession
            ),
            "best_hsp": best_hsp,
            "hsp_count": len(
                hsps
            ),
        }

        existing = candidates.get(
            accession
        )

        if (
            existing is None
            or float(
                best_hsp["bit_score"]
            )
            > float(
                existing[
                    "best_hsp"
                ]["bit_score"]
            )
        ):
            candidates[
                accession
            ] = candidate

    return sorted(
        candidates.values(),
        key=lambda candidate: (
            -float(
                candidate[
                    "best_hsp"
                ]["bit_score"]
            ),
            float(
                candidate[
                    "best_hsp"
                ]["evalue"]
            ),
            str(
                candidate[
                    "accession"
                ]
            ),
        ),
    )


def _entry_sequence(
    payload: dict[str, Any],
) -> str:
    sequence_record = payload.get(
        "sequence"
    )

    if not isinstance(
        sequence_record,
        dict,
    ):
        return ""

    return re.sub(
        r"[^A-Za-z]",
        "",
        str(
            sequence_record.get(
                "value"
            )
            or ""
        ),
    ).upper()


def _entry_reviewed(
    payload: dict[str, Any],
) -> bool:
    entry_type = str(
        payload.get(
            "entryType"
        )
        or ""
    ).casefold()

    return (
        "reviewed" in entry_type
        and "unreviewed"
        not in entry_type
    )


def _pdb_reference_count(
    payload: dict[str, Any],
) -> int:
    references = payload.get(
        "uniProtKBCrossReferences"
    )

    if not isinstance(
        references,
        list,
    ):
        return 0

    return sum(
        1
        for reference
        in references
        if (
            isinstance(
                reference,
                dict,
            )
            and str(
                reference.get(
                    "database"
                )
                or ""
            ).upper()
            == "PDB"
        )
    )


def _functional_site_count(
    payload: dict[str, Any],
) -> int:
    features = payload.get(
        "features"
    )

    if not isinstance(
        features,
        list,
    ):
        return 0

    supported_types = {
        "active site",
        "binding site",
    }

    return sum(
        1
        for feature
        in features
        if (
            isinstance(
                feature,
                dict,
            )
            and str(
                feature.get(
                    "type"
                )
                or ""
            ).casefold()
            in supported_types
        )
    )


def _candidate_rank_key(
    candidate: dict[str, Any],
) -> tuple[Any, ...]:
    alignment = candidate[
        "alignment"
    ]

    return (
        not bool(
            candidate[
                "eligible"
            ]
        ),
        -float(
            alignment[
                "combined_score"
            ]
        ),
        -float(
            alignment[
                "identity"
            ]
        ),
        -float(
            alignment[
                "query_coverage"
            ]
        ),
        -float(
            alignment[
                "target_coverage"
            ]
        ),
        -int(
            bool(
                candidate[
                    "reviewed"
                ]
            )
        ),
        -int(
            candidate[
                "pdb_cross_reference_count"
            ]
        ),
        -int(
            candidate[
                "functional_site_feature_count"
            ]
        ),
        -float(
            candidate[
                "annotation_score"
            ]
        ),
        -float(
            candidate[
                "blast"
            ][
                "best_hsp"
            ][
                "bit_score"
            ]
        ),
        str(
            candidate[
                "accession"
            ]
        ),
    )


def _candidates_are_ambiguous(
    first: dict[str, Any],
    second: dict[str, Any],
) -> bool:
    if not (
        first["eligible"]
        and second["eligible"]
    ):
        return False

    first_alignment = first[
        "alignment"
    ]

    second_alignment = second[
        "alignment"
    ]

    sequence_scores_equivalent = (
        abs(
            float(
                first_alignment[
                    "identity"
                ]
            )
            - float(
                second_alignment[
                    "identity"
                ]
            )
        )
        <= 0.005
        and abs(
            float(
                first_alignment[
                    "query_coverage"
                ]
            )
            - float(
                second_alignment[
                    "query_coverage"
                ]
            )
        )
        <= 0.01
        and abs(
            float(
                first_alignment[
                    "target_coverage"
                ]
            )
            - float(
                second_alignment[
                    "target_coverage"
                ]
            )
        )
        <= 0.01
    )

    evidence_equivalent = (
        bool(
            first["reviewed"]
        )
        == bool(
            second["reviewed"]
        )
        and int(
            first[
                "pdb_cross_reference_count"
            ]
        )
        == int(
            second[
                "pdb_cross_reference_count"
            ]
        )
        and int(
            first[
                "functional_site_feature_count"
            ]
        )
        == int(
            second[
                "functional_site_feature_count"
            ]
        )
        and abs(
            float(
                first[
                    "annotation_score"
                ]
            )
            - float(
                second[
                    "annotation_score"
                ]
            )
        )
        <= 0.01
    )

    return (
        sequence_scores_equivalent
        and evidence_equivalent
    )


def resolve_uniprot_accession_by_sequence(
    sequence: str,
    *,
    email: str,
    timeout_seconds: float = 600.0,
    poll_interval_seconds: float = 3.0,
    maximum_hits: int = 20,
    minimum_identity: float = 0.90,
    minimum_query_coverage: float = 0.90,
    minimum_target_coverage: float = 0.80,
) -> dict[str, Any]:
    if maximum_hits < 1:
        raise ValueError(
            "maximum_hits must be "
            "at least one."
        )

    for label, threshold in (
        (
            "minimum_identity",
            minimum_identity,
        ),
        (
            "minimum_query_coverage",
            minimum_query_coverage,
        ),
        (
            "minimum_target_coverage",
            minimum_target_coverage,
        ),
    ):
        if not 0 < threshold <= 1:
            raise ValueError(
                f"{label} must be "
                "greater than zero and "
                "no greater than one."
            )

    query_sequence = re.sub(
        r"[^A-Za-z]",
        "",
        str(sequence or ""),
    ).upper()

    if not query_sequence:
        raise ValueError(
            "Protein sequence is empty."
        )

    job_id, submission_metadata = (
        _submit_ebi_blast(
            query_sequence,
            email=email,
            timeout_seconds=(
                timeout_seconds
            ),
            maximum_hits=maximum_hits,
        )
    )

    poll_metadata = _poll_ebi_blast(
        job_id,
        timeout_seconds=timeout_seconds,
        poll_interval_seconds=(
            poll_interval_seconds
        ),
    )

    xml_text, xml_metadata = (
        _fetch_ebi_blast_xml(
            job_id,
            timeout_seconds=(
                timeout_seconds
            ),
        )
    )

    blast_candidates = (
        parse_ebi_blast_xml(
            xml_text
        )
    )

    if blast_candidates:
        result_metadata = {
            "selected_format": "xml",
            "xml": xml_metadata,
        }
    else:
        tsv_text, tsv_metadata = (
            _fetch_ebi_blast_tsv(
                job_id,
                timeout_seconds=(
                    timeout_seconds
                ),
            )
        )

        blast_candidates = (
            parse_ebi_blast_tsv(
                tsv_text
            )
        )

        result_metadata = {
            "selected_format": "tsv",
            "xml": xml_metadata,
            "tsv": tsv_metadata,
            "xml_candidate_count": 0,
        }

    if not blast_candidates:
        raise ValueError(
            "Sequence search returned no "
            "parseable UniProtKB candidates "
            "from either XML or TSV results."
        )

    evaluated: list[
        dict[str, Any]
    ] = []

    for blast_candidate in (
        blast_candidates[
            :maximum_hits
        ]
    ):
        accession = str(
            blast_candidate[
                "accession"
            ]
        )

        try:
            payload, entry_metadata = (
                fetch_uniprot_entry(
                    accession,
                    timeout_seconds=min(
                        timeout_seconds,
                        60.0,
                    ),
                )
            )

            candidate_sequence = (
                _entry_sequence(
                    payload
                )
            )

            if not candidate_sequence:
                raise ValueError(
                    "UniProt entry contains "
                    "no sequence."
                )

            alignment = (
                local_alignment_metrics(
                    query_sequence,
                    candidate_sequence,
                )
            )

            identity = float(
                alignment.get(
                    "identity"
                )
                or 0.0
            )

            query_coverage = float(
                alignment.get(
                    "query_coverage"
                )
                or 0.0
            )

            target_coverage = float(
                alignment.get(
                    "target_coverage"
                )
                or 0.0
            )

            checks = {
                "identity_supported": (
                    identity
                    >= minimum_identity
                ),
                "query_coverage_supported": (
                    query_coverage
                    >= minimum_query_coverage
                ),
                "target_coverage_supported": (
                    target_coverage
                    >= minimum_target_coverage
                ),
            }

            primary_accession = (
                normalize_uniprot_accession(
                    str(
                        payload.get(
                            "primaryAccession"
                        )
                        or accession
                    )
                )
            )

            evaluated.append(
                {
                    "accession": (
                        primary_accession
                    ),
                    "reviewed": (
                        _entry_reviewed(
                            payload
                        )
                    ),
                    "entry_type": (
                        payload.get(
                            "entryType"
                        )
                    ),
                    "annotation_score": float(
                        payload.get(
                            "annotationScore"
                        )
                        or 0.0
                    ),
                    "pdb_cross_reference_count": (
                        _pdb_reference_count(
                            payload
                        )
                    ),
                    "functional_site_feature_count": (
                        _functional_site_count(
                            payload
                        )
                    ),
                    "sequence_length": len(
                        candidate_sequence
                    ),
                    "alignment": alignment,
                    "eligibility_checks": (
                        checks
                    ),
                    "eligible": all(
                        checks.values()
                    ),
                    "blast": (
                        blast_candidate
                    ),
                    "entry_metadata": (
                        entry_metadata
                    ),
                }
            )
        except Exception as error:
            evaluated.append(
                {
                    "accession": accession,
                    "eligible": False,
                    "retrieval_error": str(
                        error
                    ),
                    "blast": (
                        blast_candidate
                    ),
                }
            )

    eligible_candidates = [
        candidate
        for candidate in evaluated
        if candidate.get(
            "eligible"
        )
    ]

    if not eligible_candidates:
        best_available = next(
            (
                candidate
                for candidate
                in sorted(
                    evaluated,
                    key=lambda value: (
                        -float(
                            (
                                value.get(
                                    "alignment"
                                )
                                or {}
                            ).get(
                                "combined_score"
                            )
                            or 0.0
                        ),
                        str(
                            value.get(
                                "accession"
                            )
                            or ""
                        ),
                    ),
                )
                if candidate.get(
                    "alignment"
                )
            ),
            None,
        )

        raise ValueError(
            "No sequence-search candidate "
            "met the strict identity and "
            "coverage gates. Best candidate: "
            + json.dumps(
                best_available,
                sort_keys=True,
                default=str,
            )
        )

    ranked = sorted(
        eligible_candidates,
        key=_candidate_rank_key,
    )

    selected = ranked[0]

    if (
        len(ranked) > 1
        and _candidates_are_ambiguous(
            ranked[0],
            ranked[1],
        )
    ):
        raise ValueError(
            "Sequence-based UniProt "
            "resolution was ambiguous "
            "between equally supported "
            "candidates: "
            f"{ranked[0]['accession']}, "
            f"{ranked[1]['accession']}."
        )

    selected_alignment = (
        selected["alignment"]
    )

    confidence = (
        "high"
        if (
            float(
                selected_alignment[
                    "identity"
                ]
            )
            >= 0.98
            and float(
                selected_alignment[
                    "query_coverage"
                ]
            )
            >= 0.95
            and float(
                selected_alignment[
                    "target_coverage"
                ]
            )
            >= 0.90
        )
        else "moderate"
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
            "ebi_blast_sequence_alignment"
        ),
        "confidence": confidence,
        "candidate_count": len(
            evaluated
        ),
        "eligible_candidate_count": len(
            ranked
        ),
        "selected_candidate": (
            selected
        ),
        "candidates": sorted(
            evaluated,
            key=lambda candidate: (
                _candidate_rank_key(
                    candidate
                )
                if candidate.get(
                    "alignment"
                )
                else (
                    True,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0,
                    0,
                    0,
                    0.0,
                    0.0,
                    str(
                        candidate.get(
                            "accession"
                        )
                        or ""
                    ),
                )
            ),
        ),
        "thresholds": {
            "minimum_identity": (
                minimum_identity
            ),
            "minimum_query_coverage": (
                minimum_query_coverage
            ),
            "minimum_target_coverage": (
                minimum_target_coverage
            ),
        },
        "sequence_search": {
            "provider": (
                "EMBL-EBI Job Dispatcher"
            ),
            "tool": "NCBI BLAST+",
            "program": "blastp",
            "database": "uniprotkb",
            "job_id": job_id,
            "submission": (
                submission_metadata
            ),
            "polling": poll_metadata,
            "result": result_metadata,
        },
        "safety_checks": {
            "sequence_alignment_supported": (
                True
            ),
            "single_best_candidate": True,
            "strict_identity_gate": True,
            "strict_query_coverage_gate": (
                True
            ),
            "strict_target_coverage_gate": (
                True
            ),
            "downstream_sequence_validation_required": (
                True
            ),
        },
    }


def resolve_uniprot_accession_from_fasta(
    fasta_path: Path,
    *,
    sequence_search_email: (
        str | None
    ) = None,
    sequence_search_timeout_seconds: (
        float
    ) = 600.0,
    sequence_search_poll_interval_seconds: (
        float
    ) = 3.0,
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
        if not sequence_search_email:
            raise ValueError(
                "No strict UniProt accession "
                "was found in the FASTA header "
                "or filename. Supply "
                "--reference-uniprot-accession "
                "or configure "
                "--reference-sequence-search-email "
                "for sequence-based discovery."
            )

        sequence_resolution = (
            resolve_uniprot_accession_by_sequence(
                read_fasta_sequence(
                    source
                ),
                email=(
                    sequence_search_email
                ),
                timeout_seconds=(
                    sequence_search_timeout_seconds
                ),
                poll_interval_seconds=(
                    sequence_search_poll_interval_seconds
                ),
            )
        )

        return {
            **sequence_resolution,
            "source_fasta": str(
                source
            ),
            "fasta_header": header,
        }

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
        "confidence": "high",
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
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )

    return destination
