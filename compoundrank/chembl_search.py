"""Generic ChEMBL target and bioactivity retrieval.

This module must not contain target-specific compound names.

Its responsibility is:

generic target queries
    -> ChEMBL target records
    -> measured activity records
    -> normalized ligand candidates

The module does not use the local ligand rule registry.
"""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from collections.abc import Callable
from typing import Any


CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

TARGET_SCORE_WINDOW = 45.0


SUPPORTED_ACTIVITY_TYPES = {
    "IC50",
    "EC50",
    "Ki",
    "Kd",
    "Potency",
    "Inhibition",
}

RequestJson = Callable[
    [str, dict[str, Any], int],
    dict[str, Any],
]


def chembl_request_json(
    resource: str,
    params: dict[str, Any],
    timeout_seconds: int,
) -> dict[str, Any]:
    """Request one JSON page from a ChEMBL web-service resource."""
    query = urllib.parse.urlencode(params, doseq=True)
    url = f"{CHEMBL_BASE}/{resource}.json"

    if query:
        url = f"{url}?{query}"

    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": (
                "CompoundRank-stage4A/0.2 "
                "educational-research-use"
            ),
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=timeout_seconds,
    ) as response:
        raw = response.read().decode(
            "utf-8",
            errors="replace",
        )

    payload = json.loads(raw)

    if not isinstance(payload, dict):
        raise ValueError(
            f"Expected ChEMBL JSON object from {resource}"
        )

    return payload


def _flatten_text(value: Any) -> str:
    parts: list[str] = []

    def visit(item: Any) -> None:
        if item is None:
            return

        if isinstance(item, dict):
            for key, nested in item.items():
                parts.append(str(key))
                visit(nested)
            return

        if isinstance(item, (list, tuple, set)):
            for nested in item:
                visit(nested)
            return

        parts.append(str(item))

    visit(value)
    return " ".join(parts)


def _normalized_text(value: Any) -> str:
    text = str(value or "").lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def _target_noun(target_context: dict[str, Any]) -> str:
    """Extract the principal protein-function noun from target context."""
    values = [
        target_context.get("target_name"),
        target_context.get("enzyme_class"),
        target_context.get("target_class"),
        target_context.get("special_domain_label"),
    ]
    blob = _normalized_text(" ".join(str(value or "") for value in values))

    supported_nouns = (
        "neuraminidase",
        "methyltransferase",
        "endonuclease",
        "exonuclease",
        "polymerase",
        "integrase",
        "protease",
        "peptidase",
        "helicase",
        "hydrolase",
        "kinase",
    )

    for noun in supported_nouns:
        if noun in blob:
            return noun

    return ""



def _context_acronyms(
    target_context: dict[str, Any] | None,
) -> set[str]:
    """Extract explicit organism acronyms such as HIV from evidence."""
    if not target_context:
        return set()

    values = [
        target_context.get("target_name"),
        target_context.get("viral_family"),
    ]

    acronyms: set[str] = set()

    for value in values:
        text = str(value or "")
        for match in re.findall(
            r"\b[A-Z][A-Z0-9]{1,9}\b",
            text,
        ):
            # Ignore accessions and terms containing no letters.
            if any(character.isalpha() for character in match):
                acronyms.add(match.lower())

    return acronyms


def _organism_identity_aliases(
    target: dict[str, Any],
) -> set[str]:
    """Build normalized aliases and initialisms for a target organism."""
    organism = _normalized_text(
        target.get("organism")
    )
    pref_name = _normalized_text(
        target.get("pref_name")
    )

    aliases: set[str] = set()

    for text in (organism, pref_name):
        if not text:
            continue

        aliases.add(text)

        for token in text.split():
            if len(token) >= 2:
                aliases.add(token)

    words = [
        word
        for word in organism.split()
        if word.isalpha()
        and word not in {
            "type",
            "subtype",
            "strain",
            "isolate",
            "group",
        }
        and len(word) > 1
    ]

    if words:
        aliases.add(
            "".join(word[0] for word in words)
        )

    # Also calculate an initialism before trailing subtype numbers.
    leading_words: list[str] = []
    for word in organism.split():
        if word.isdigit():
            break
        if word in {
            "type",
            "subtype",
            "strain",
            "isolate",
            "group",
        }:
            continue
        if word.isalpha() and len(word) > 1:
            leading_words.append(word)

    if leading_words:
        aliases.add(
            "".join(word[0] for word in leading_words)
        )

    return aliases


def _target_identity_matches(
    target: dict[str, Any],
    target_context: dict[str, Any] | None,
) -> bool:
    expected = _context_acronyms(target_context)

    if not expected:
        return True

    aliases = _organism_identity_aliases(target)
    target_blob = _normalized_text(_flatten_text(target))

    return any(
        acronym in aliases
        or acronym in target_blob.split()
        for acronym in expected
    )


def build_context_target_terms(
    target_context: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Generate prioritized ChEMBL target searches from evidence."""
    if not target_context:
        return []

    target_name = str(
        target_context.get("target_name") or ""
    ).strip()
    domain_label = str(
        target_context.get("special_domain_label") or ""
    ).strip()
    viral_family = str(
        target_context.get("viral_family") or ""
    ).strip()
    target_noun = _target_noun(target_context)
    acronyms = _context_acronyms(target_context)

    raw_terms: list[tuple[str, int, str]] = []

    for acronym in sorted(acronyms):
        raw_terms.append(
            (
                acronym.upper(),
                140,
                "context_organism_acronym_search",
            )
        )

        if target_noun:
            raw_terms.append(
                (
                    f"{acronym.upper()} {target_noun}",
                    145,
                    "context_acronym_target_search",
                )
            )

    if target_name and target_name.lower() != "unknown":
        cleaned_name = re.sub(
            r"-like\b",
            "",
            target_name,
            flags=re.IGNORECASE,
        )
        cleaned_name = " ".join(cleaned_name.split())

        raw_terms.append(
            (
                cleaned_name,
                135,
                "context_target_name_search",
            )
        )

    if domain_label:
        cleaned_domain = re.sub(
            r"\bdomain\b",
            "",
            domain_label,
            flags=re.IGNORECASE,
        )
        cleaned_domain = " ".join(cleaned_domain.split())

        if cleaned_domain:
            raw_terms.append(
                (
                    cleaned_domain,
                    125,
                    "context_domain_label_search",
                )
            )

    if viral_family and target_noun:
        for token in re.findall(
            r"[A-Za-z]+",
            viral_family,
        ):
            normalized = token.lower()

            if normalized.endswith("viral"):
                raw_terms.append(
                    (
                        f"{token} {target_noun}",
                        120,
                        "context_viral_family_search",
                    )
                )

    output: list[dict[str, Any]] = []
    seen: set[str] = set()

    for term, specificity, route in raw_terms:
        key = _normalized_text(term)

        if not key or key in seen:
            continue

        seen.add(key)

        output.append(
            {
                "term": term,
                "specificity": specificity,
                "retrieval_route": route,
                "original_query": term,
            }
        )

    return output

def _expects_viral_target(
    target_context: dict[str, Any] | None,
) -> bool:
    if not target_context:
        return False

    context_blob = _normalized_text(target_context)

    viral_markers = (
        "virus",
        "viral",
        "viridae",
        "retrovir",
        "hiv",
        "influenza",
        "coronavirus",
        "bacteriophage",
        "phage",
    )

    return any(
        marker in context_blob
        for marker in viral_markers
    )



def target_is_compatible(
    target: dict[str, Any],
    target_context: dict[str, Any] | None,
) -> bool:
    """Reject targets that contradict organism or enzyme context."""
    if not target_context:
        return True

    target_type = _normalized_text(
        target.get("target_type")
    )

    if target_type and "protein" not in target_type:
        return False

    target_blob = _normalized_text(
        _flatten_text(target)
    )
    expected_noun = _target_noun(target_context)

    if expected_noun and expected_noun not in target_blob:
        return False

    if _expects_viral_target(target_context):
        organism = _normalized_text(
            target.get("organism")
        )

        viral_markers = (
            "virus",
            "viridae",
            "retrovir",
            "hiv",
            "influenza",
            "coronavirus",
            "bacteriophage",
            "phage",
        )

        if not any(
            marker in organism
            for marker in viral_markers
        ):
            return False

    if not _target_identity_matches(
        target,
        target_context,
    ):
        return False

    expected_enzyme = _normalized_text(
        target_context.get("enzyme_class")
    )

    if (
        "aspartyl protease" in expected_enzyme
        or "aspartic protease" in expected_enzyme
    ):
        expected_markers = (
            "aspartic",
            "aspartyl",
            "aspartate",
            "retropepsin",
        )
        conflicting_markers = (
            "cysteine protease",
            "serine protease",
            "metallo protease",
            "metalloprotease",
            "threonine protease",
        )

        has_expected_marker = any(
            marker in target_blob
            for marker in expected_markers
        )
        has_conflicting_marker = any(
            marker in target_blob
            for marker in conflicting_markers
        )

        if has_conflicting_marker and not has_expected_marker:
            return False

    return True

def _target_search_term(query: str) -> str:
    """Convert a ligand-style query into a ChEMBL target term."""
    text = str(query).strip()

    text = re.sub(
        r"\s+(inhibitor|ligand)\s*$",
        "",
        text,
        flags=re.IGNORECASE,
    )

    return " ".join(text.split())


def build_chembl_target_terms(
    generic_queries: list[dict[str, Any]],
    *,
    max_terms: int = 6,
) -> list[dict[str, Any]]:
    """Create prioritized ChEMBL target-search terms."""
    ordered = sorted(
        generic_queries,
        key=lambda item: int(item.get("specificity") or 0),
        reverse=True,
    )

    output: list[dict[str, Any]] = []
    seen: set[str] = set()

    for item in ordered:
        route = str(item.get("retrieval_route") or "")

        # Pfam/domain accessions are useful provenance, but ChEMBL's
        # full-text target search is primarily name-based.
        if "accession" in route:
            continue

        term = _target_search_term(
            str(item.get("query") or "")
        )
        key = _normalized_text(term)

        if not key or key == "unknown":
            continue
        if key in seen:
            continue

        seen.add(key)

        output.append(
            {
                "term": term,
                "specificity": int(
                    item.get("specificity") or 0
                ),
                "retrieval_route": route,
                "original_query": item.get("query"),
            }
        )

        if max_terms > 0 and len(output) >= max_terms:
            break

    return output



def _score_target(
    target: dict[str, Any],
    search_item: dict[str, Any],
    target_context: dict[str, Any] | None = None,
) -> float:
    """Score a biologically compatible raw ChEMBL target."""
    term = _normalized_text(search_item["term"])
    pref_name = _normalized_text(
        target.get("pref_name")
    )
    target_blob = _normalized_text(
        _flatten_text(target)
    )

    term_tokens = set(term.split())
    target_tokens = set(target_blob.split())

    overlap = 0.0
    if term_tokens:
        overlap = len(
            term_tokens & target_tokens
        ) / len(term_tokens)

    score = float(search_item["specificity"])

    if term and term == pref_name:
        score += 100.0
    elif term and term in target_blob:
        score += 50.0

    score += overlap * 40.0

    if "protein" in _normalized_text(
        target.get("target_type")
    ):
        score += 10.0

    expected_acronyms = _context_acronyms(
        target_context
    )
    target_aliases = _organism_identity_aliases(
        target
    )

    if expected_acronyms & target_aliases:
        score += 200.0

    generic_names = {
        "protease",
        "polymerase",
        "integrase",
        "helicase",
        "kinase",
        "hydrolase",
        "enzyme",
        "protein",
    }

    if pref_name in generic_names:
        score -= 60.0
    else:
        descriptive_tokens = {
            token
            for token in pref_name.split()
            if token not in generic_names
            and token not in {
                "type",
                "subtype",
                "protein",
                "enzyme",
            }
        }
        score += min(
            60.0,
            len(descriptive_tokens) * 10.0,
        )

    if target_context:
        expected_noun = _target_noun(
            target_context
        )
        if expected_noun and expected_noun in pref_name:
            score += 30.0

        expected_enzyme = _normalized_text(
            target_context.get("enzyme_class")
        )

        if (
            "aspartyl protease" in expected_enzyme
            or "aspartic protease" in expected_enzyme
        ):
            if any(
                marker in target_blob
                for marker in (
                    "aspartic",
                    "aspartyl",
                    "aspartate",
                    "retropepsin",
                )
            ):
                score += 50.0

    return round(score, 4)


def search_chembl_targets(
    generic_queries: list[dict[str, Any]],
    *,
    target_context: dict[str, Any] | None = None,
    max_terms: int = 8,
    per_query_limit: int = 50,
    max_targets: int = 5,
    timeout_seconds: int = 60,
    request_json: RequestJson = chembl_request_json,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Search broadly, then apply strict biological target ranking."""
    combined_terms = (
        build_context_target_terms(target_context)
        + build_chembl_target_terms(
            generic_queries,
            max_terms=0,
        )
    )

    combined_terms.sort(
        key=lambda item: int(
            item.get("specificity") or 0
        ),
        reverse=True,
    )

    search_terms: list[dict[str, Any]] = []
    seen_terms: set[str] = set()

    for item in combined_terms:
        key = _normalized_text(item.get("term"))

        if not key or key in seen_terms:
            continue

        seen_terms.add(key)
        search_terms.append(item)

        if max_terms > 0 and len(search_terms) >= max_terms:
            break

    targets_by_id: dict[str, dict[str, Any]] = {}

    for search_item in search_terms:
        payload = request_json(
            "target/search",
            {
                "q": search_item["term"],
                "limit": per_query_limit,
            },
            timeout_seconds,
        )

        records = payload.get("targets") or []

        if not isinstance(records, list):
            continue

        for record in records:
            if not isinstance(record, dict):
                continue

            target_id = str(
                record.get("target_chembl_id") or ""
            ).strip()

            if not target_id:
                continue

            if not target_is_compatible(
                record,
                target_context,
            ):
                continue

            score = _score_target(
                record,
                search_item,
                target_context,
            )

            candidate = dict(record)
            candidate["discovery_query"] = search_item[
                "original_query"
            ]
            candidate["target_search_term"] = search_item[
                "term"
            ]
            candidate["query_specificity"] = search_item[
                "specificity"
            ]
            candidate["query_retrieval_route"] = search_item[
                "retrieval_route"
            ]
            candidate["target_match_score"] = score

            existing = targets_by_id.get(target_id)

            if (
                existing is None
                or score
                > float(
                    existing.get(
                        "target_match_score"
                    )
                    or 0
                )
            ):
                targets_by_id[target_id] = candidate

    targets = list(targets_by_id.values())
    targets.sort(
        key=lambda item: (
            -float(
                item.get("target_match_score") or 0
            ),
            str(
                item.get("target_chembl_id") or ""
            ),
        )
    )

    # Keep only targets close to the best biologically specific
    # target instead of filling the quota with distant viral proteins.
    if targets:
        top_score = float(
            targets[0].get("target_match_score") or 0
        )
        targets = [
            target
            for target in targets
            if float(
                target.get("target_match_score") or 0
            )
            >= top_score - TARGET_SCORE_WINDOW
        ]

    if max_targets > 0:
        targets = targets[:max_targets]

    return targets, search_terms


SEQUENCE_SCORE_WINDOW = 0.03


def normalize_protein_sequence(value: str) -> str:
    """Normalize and validate a canonical amino-acid sequence."""
    sequence = re.sub(
        r"\s+",
        "",
        str(value or ""),
    ).upper()

    allowed = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = sorted(set(sequence) - allowed)

    if invalid:
        raise ValueError(
            f"Invalid amino-acid characters: {invalid}"
        )

    return sequence


def read_fasta_sequence(path: Any) -> str:
    """Read all sequence lines from a single-protein FASTA."""
    from pathlib import Path

    fasta_path = Path(path)
    parts: list[str] = []

    for line in fasta_path.read_text(
        encoding="utf-8"
    ).splitlines():
        line = line.strip()

        if not line or line.startswith(">"):
            continue

        parts.append(line)

    sequence = normalize_protein_sequence(
        "".join(parts)
    )

    if not sequence:
        raise ValueError(
            f"No protein sequence found in {fasta_path}"
        )

    return sequence


def local_alignment_metrics(
    query: str,
    target: str,
) -> dict[str, float | int]:
    """Calculate Smith-Waterman local-alignment metrics."""
    query = normalize_protein_sequence(query)
    target = normalize_protein_sequence(target)

    rows = len(query) + 1
    columns = len(target) + 1

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

    for i in range(1, rows):
        for j in range(1, columns):
            diagonal = (
                scores[i - 1][j - 1]
                + (
                    2
                    if query[i - 1] == target[j - 1]
                    else -1
                )
            )
            upward = scores[i - 1][j] - 2
            leftward = scores[i][j - 1] - 2

            score = max(
                0,
                diagonal,
                upward,
                leftward,
            )
            scores[i][j] = score

            if score == 0:
                direction = 0
            elif score == diagonal:
                direction = 1
            elif score == upward:
                direction = 2
            else:
                direction = 3

            directions[i][j] = direction

            if score > best_score:
                best_score = score
                best_position = (i, j)

    i, j = best_position

    matches = 0
    aligned_columns = 0
    query_residues = 0
    target_residues = 0

    while i > 0 and j > 0:
        direction = directions[i][j]

        if direction == 0:
            break

        aligned_columns += 1

        if direction == 1:
            query_residues += 1
            target_residues += 1

            if query[i - 1] == target[j - 1]:
                matches += 1

            i -= 1
            j -= 1

        elif direction == 2:
            query_residues += 1
            i -= 1

        else:
            target_residues += 1
            j -= 1

    identity = (
        matches / aligned_columns
        if aligned_columns
        else 0.0
    )
    query_coverage = (
        query_residues / len(query)
        if query
        else 0.0
    )
    target_coverage = (
        target_residues / len(target)
        if target
        else 0.0
    )

    return {
        "alignment_score": best_score,
        "matches": matches,
        "aligned_columns": aligned_columns,
        "identity": identity,
        "query_coverage": query_coverage,
        "target_coverage": target_coverage,
        "combined_score": identity * query_coverage,
    }


def _component_records(
    payload: dict[str, Any],
) -> list[dict[str, Any]]:
    records = payload.get("target_components")

    if isinstance(records, list):
        return [
            item
            for item in records
            if isinstance(item, dict)
        ]

    if (
        "component_id" in payload
        or "accession" in payload
        or "sequence" in payload
    ):
        return [payload]

    return []


def _fetch_component_records(
    component: dict[str, Any],
    *,
    timeout_seconds: int,
    request_json: RequestJson,
) -> list[dict[str, Any]]:
    """Hydrate a target-component stub."""
    if str(component.get("sequence") or "").strip():
        return [component]

    component_id = component.get("component_id")
    accession = str(
        component.get("accession") or ""
    ).strip()

    if component_id not in (None, ""):
        try:
            payload = request_json(
                f"target_component/{component_id}",
                {},
                timeout_seconds,
            )
            records = _component_records(payload)

            if records:
                return records
        except Exception:
            pass

    if accession:
        try:
            payload = request_json(
                "target_component",
                {
                    "accession": accession,
                    "limit": 100,
                },
                timeout_seconds,
            )
            records = _component_records(payload)

            if records:
                return records
        except Exception:
            pass

    return []


def _combined_target_search_terms(
    generic_queries: list[dict[str, Any]],
    target_context: dict[str, Any] | None,
    *,
    max_terms: int,
) -> list[dict[str, Any]]:
    combined = (
        build_context_target_terms(target_context)
        + build_chembl_target_terms(
            generic_queries,
            max_terms=0,
        )
    )

    combined.sort(
        key=lambda item: int(
            item.get("specificity") or 0
        ),
        reverse=True,
    )

    output: list[dict[str, Any]] = []
    seen: set[str] = set()

    for item in combined:
        key = _normalized_text(item.get("term"))

        if not key or key in seen:
            continue

        seen.add(key)
        output.append(item)

        if max_terms > 0 and len(output) >= max_terms:
            break

    return output


def search_chembl_targets_sequence_first(
    generic_queries: list[dict[str, Any]],
    *,
    target_context: dict[str, Any] | None,
    query_sequence: str,
    max_terms: int = 8,
    per_query_limit: int = 50,
    max_targets: int = 5,
    minimum_identity: float = 0.70,
    minimum_query_coverage: float = 0.70,
    timeout_seconds: int = 60,
    request_json: RequestJson = chembl_request_json,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Resolve ChEMBL targets using sequence evidence before text rank."""
    query_sequence = normalize_protein_sequence(
        query_sequence
    )

    search_terms = _combined_target_search_terms(
        generic_queries,
        target_context,
        max_terms=max_terms,
    )

    targets_by_id: dict[str, dict[str, Any]] = {}

    for search_item in search_terms:
        payload = request_json(
            "target/search",
            {
                "q": search_item["term"],
                "limit": per_query_limit,
            },
            timeout_seconds,
        )

        records = payload.get("targets") or []

        if not isinstance(records, list):
            continue

        for record in records:
            if not isinstance(record, dict):
                continue

            target_id = str(
                record.get("target_chembl_id") or ""
            ).strip()

            if not target_id:
                continue

            if not target_is_compatible(
                record,
                target_context,
            ):
                continue

            text_score = _score_target(
                record,
                search_item,
                target_context,
            )

            candidate = dict(record)
            candidate["discovery_query"] = search_item[
                "original_query"
            ]
            candidate["target_search_term"] = search_item[
                "term"
            ]
            candidate["query_specificity"] = search_item[
                "specificity"
            ]
            candidate["query_retrieval_route"] = search_item[
                "retrieval_route"
            ]
            candidate["target_match_score"] = text_score

            existing = targets_by_id.get(target_id)

            if (
                existing is None
                or text_score
                > float(
                    existing.get("target_match_score")
                    or 0
                )
            ):
                targets_by_id[target_id] = candidate

    sequence_ranked: list[dict[str, Any]] = []

    for target_id, search_record in targets_by_id.items():
        try:
            full_target = request_json(
                f"target/{target_id}",
                {},
                timeout_seconds,
            )
        except Exception:
            full_target = search_record

        if not isinstance(full_target, dict):
            full_target = search_record

        components = (
            full_target.get("target_components")
            or search_record.get("target_components")
            or []
        )

        component_matches: list[dict[str, Any]] = []

        for component_stub in components:
            if not isinstance(component_stub, dict):
                continue

            records = _fetch_component_records(
                component_stub,
                timeout_seconds=timeout_seconds,
                request_json=request_json,
            )

            for component in records:
                component_sequence = str(
                    component.get("sequence") or ""
                )

                try:
                    component_sequence = (
                        normalize_protein_sequence(
                            component_sequence
                        )
                    )
                except ValueError:
                    continue

                if not component_sequence:
                    continue

                metrics = local_alignment_metrics(
                    query_sequence,
                    component_sequence,
                )

                component_matches.append(
                    {
                        "component_id": component.get(
                            "component_id"
                        ),
                        "accession": component.get(
                            "accession"
                        ),
                        "description": component.get(
                            "component_description"
                        ),
                        "component_sequence_length": len(
                            component_sequence
                        ),
                        **metrics,
                    }
                )

        component_matches.sort(
            key=lambda item: (
                -float(item["combined_score"]),
                -float(item["identity"]),
                -float(item["query_coverage"]),
                str(item.get("accession") or ""),
            )
        )

        if not component_matches:
            continue

        best = component_matches[0]

        enriched = dict(full_target)

        for key, value in search_record.items():
            enriched.setdefault(key, value)

        enriched["sequence_resolution"] = {
            "route": "local_component_alignment",
            "matched_component_id": best.get(
                "component_id"
            ),
            "matched_component_accession": best.get(
                "accession"
            ),
            "matched_component_description": best.get(
                "description"
            ),
            "component_sequence_length": best.get(
                "component_sequence_length"
            ),
            "identity": best["identity"],
            "query_coverage": best[
                "query_coverage"
            ],
            "target_coverage": best[
                "target_coverage"
            ],
            "alignment_score": best[
                "alignment_score"
            ],
            "combined_score": best[
                "combined_score"
            ],
        }
        enriched["sequence_identity"] = best["identity"]
        enriched["sequence_query_coverage"] = best[
            "query_coverage"
        ]
        enriched["sequence_match_score"] = best[
            "combined_score"
        ]
        enriched["matched_component_accession"] = (
            best.get("accession")
        )
        enriched["target_resolution_route"] = (
            "sequence_alignment"
        )

        sequence_ranked.append(enriched)

    sequence_ranked.sort(
        key=lambda item: (
            -float(
                item.get("sequence_match_score")
                or 0
            ),
            -float(
                item.get("sequence_identity")
                or 0
            ),
            -float(
                item.get("sequence_query_coverage")
                or 0
            ),
            -float(
                item.get("target_match_score")
                or 0
            ),
            str(
                item.get("target_chembl_id")
                or ""
            ),
        )
    )

    strong = [
        target
        for target in sequence_ranked
        if float(
            target.get("sequence_identity") or 0
        )
        >= minimum_identity
        and float(
            target.get(
                "sequence_query_coverage"
            )
            or 0
        )
        >= minimum_query_coverage
    ]

    if strong:
        best_score = float(
            strong[0].get("sequence_match_score")
            or 0
        )

        selected = [
            target
            for target in strong
            if float(
                target.get("sequence_match_score")
                or 0
            )
            >= best_score - SEQUENCE_SCORE_WINDOW
        ]

        if max_targets > 0:
            selected = selected[:max_targets]

        return selected, search_terms

    fallback_targets, _ = search_chembl_targets(
        generic_queries,
        target_context=target_context,
        max_terms=max_terms,
        per_query_limit=per_query_limit,
        max_targets=max_targets,
        timeout_seconds=timeout_seconds,
        request_json=request_json,
    )

    for target in fallback_targets:
        target["target_resolution_route"] = (
            "text_search_fallback"
        )
        target["sequence_resolution"] = {
            "route": "text_search_fallback",
            "reason": (
                "No candidate component met the minimum "
                "sequence identity and query-coverage thresholds."
            ),
        }

    return fallback_targets, search_terms


def fetch_chembl_activities(
    target_chembl_id: str,
    *,
    limit: int = 200,
    timeout_seconds: int = 60,
    request_json: RequestJson = chembl_request_json,
) -> list[dict[str, Any]]:
    """Retrieve measured activities for one ChEMBL target."""
    payload = request_json(
        "activity",
        {
            "target_chembl_id": target_chembl_id,
            "standard_type__in": ",".join(
                sorted(SUPPORTED_ACTIVITY_TYPES)
            ),
            "pchembl_value__isnull": "false",
            "order_by": "-pchembl_value",
            "limit": limit,
        },
        timeout_seconds,
    )

    records = payload.get("activities") or []

    if not isinstance(records, list):
        return []

    return [
        record
        for record in records
        if isinstance(record, dict)
    ]


def _as_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def activity_is_usable(
    activity: dict[str, Any],
) -> bool:
    """Reject incomplete or clearly questionable activity records."""
    molecule_id = str(
        activity.get("molecule_chembl_id") or ""
    ).strip()
    if not molecule_id:
        return False

    activity_type = str(
        activity.get("standard_type") or ""
    ).strip()
    if activity_type not in SUPPORTED_ACTIVITY_TYPES:
        return False

    if _as_float(activity.get("pchembl_value")) is None:
        return False

    if activity.get("potential_duplicate") in {
        1,
        "1",
        True,
    }:
        return False

    validity_comment = str(
        activity.get("data_validity_comment") or ""
    ).strip()
    if validity_comment:
        return False

    relation = str(
        activity.get("standard_relation") or ""
    ).strip()
    if relation not in {"", "=", "<", "<=", "~"}:
        return False

    return True


def _evidence_level(pchembl_value: float) -> str:
    if pchembl_value >= 7.0:
        return "strong"
    if pchembl_value >= 6.0:
        return "moderate"
    if pchembl_value >= 5.0:
        return "weak"
    return "exploratory"


def make_chembl_candidate(
    activity: dict[str, Any],
    target: dict[str, Any],
) -> dict[str, Any]:
    """Normalize one ChEMBL activity into the Stage 4A schema."""
    molecule_id = str(
        activity["molecule_chembl_id"]
    ).strip()
    preferred_name = str(
        activity.get("molecule_pref_name") or ""
    ).strip()
    pchembl_value = float(activity["pchembl_value"])

    compound_name = preferred_name or molecule_id

    return {
        "compound_name": compound_name,
        "retrieval_rank": None,
        "design_status": "database_observed_ligand",
        "discovery_source": "chembl_activity",
        "hardcoded_seed": False,
        "retrieval_route": (
            "generic_chembl_target_activity_search"
        ),
        "retrieval_rule_id": None,
        "retrieval_rule_label": None,
        "retrieval_reason": (
            "Retrieved from a measured ChEMBL activity "
            f"associated with target "
            f"{target.get('target_chembl_id')}."
        ),
        "discovery_query": target.get(
            "discovery_query"
        ),
        "target_search_term": target.get(
            "target_search_term"
        ),
        "target_match_score": target.get(
            "target_match_score"
        ),
        "target_name": target.get("pref_name"),
        "target_class": target.get("target_type"),
        "enzyme_class": None,
        "viral_family": target.get("organism"),
        "target_evidence_confidence": None,
        "special_domain_label": None,
        "special_domain_accession": None,
        "target_family_basis": target.get(
            "pref_name"
        ),
        "domain_basis": [],
        "evidence_level": _evidence_level(
            pchembl_value
        ),
        "retrieval_terms": [
            target.get("target_search_term")
        ],
        "source_databases": ["ChEMBL"],
        "source_notes": "",
        "chembl_target_id": target.get(
            "target_chembl_id"
        ),
        "chembl_molecule_id": molecule_id,
        "chembl_activity_id": activity.get(
            "activity_id"
        ),
        "chembl_assay_id": activity.get(
            "assay_chembl_id"
        ),
        "chembl_document_id": activity.get(
            "document_chembl_id"
        ),
        "activity_type": activity.get(
            "standard_type"
        ),
        "activity_relation": activity.get(
            "standard_relation"
        ),
        "activity_value": activity.get(
            "standard_value"
        ),
        "activity_units": activity.get(
            "standard_units"
        ),
        "pchembl_value": pchembl_value,
        "assay_type": activity.get("assay_type"),
        "canonical_smiles": activity.get(
            "canonical_smiles"
        ),
        "pubchem_cid": None,
        "smiles": activity.get(
            "canonical_smiles"
        ),
        "inchi_key": None,
        "structure_source": None,
        "structure_fetch_status": "not_attempted",
        "structure_fetch_error": None,
        "local_sdf_path": None,
        "selected_for_docking": True,
        "mutation_relevance": {
            "status": "not_evaluated",
            "pocket_mutation_overlap": None,
            "expected_interaction_change": None,
            "adaptation_rationale": None,
        },
        "novel_design_relevance": {
            "status": "not_generated",
            "pocket_feature_basis": None,
            "pharmacophore_basis": None,
            "hypothetical_modification": None,
        },
        "limitations": [
            "Database activity does not prove inhibition of "
            "the submitted sequence.",
            "Assay conditions and target transferability require "
            "later evidence review.",
        ],
        "supporting_activities": [
            {
                "activity_id": activity.get(
                    "activity_id"
                ),
                "assay_chembl_id": activity.get(
                    "assay_chembl_id"
                ),
                "standard_type": activity.get(
                    "standard_type"
                ),
                "standard_value": activity.get(
                    "standard_value"
                ),
                "standard_units": activity.get(
                    "standard_units"
                ),
                "pchembl_value": pchembl_value,
            }
        ],
    }


def retrieve_chembl_candidates(
    generic_queries: list[dict[str, Any]],
    *,
    target_context: dict[str, Any] | None = None,
    query_sequence: str | None = None,
    max_targets: int = 5,
    activities_per_target: int = 200,
    max_candidates: int = 20,
    timeout_seconds: int = 60,
    request_json: RequestJson = chembl_request_json,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run generic target search and measured-activity retrieval."""
    if query_sequence:
        targets, search_terms = (
            search_chembl_targets_sequence_first(
                generic_queries,
                target_context=target_context,
                query_sequence=query_sequence,
                max_targets=max_targets,
                timeout_seconds=timeout_seconds,
                request_json=request_json,
            )
        )
    else:
        targets, search_terms = search_chembl_targets(
            generic_queries,
            target_context=target_context,
            max_targets=max_targets,
            timeout_seconds=timeout_seconds,
            request_json=request_json,
        )

    candidates_by_molecule: dict[
        str,
        dict[str, Any],
    ] = {}

    activity_count = 0
    accepted_activity_count = 0

    for target in targets:
        target_id = str(
            target.get("target_chembl_id") or ""
        ).strip()
        if not target_id:
            continue

        activities = fetch_chembl_activities(
            target_id,
            limit=activities_per_target,
            timeout_seconds=timeout_seconds,
            request_json=request_json,
        )
        activity_count += len(activities)

        for activity in activities:
            if not activity_is_usable(activity):
                continue

            accepted_activity_count += 1

            candidate = make_chembl_candidate(
                activity,
                target,
            )
            molecule_id = str(
                candidate["chembl_molecule_id"]
            )

            existing = candidates_by_molecule.get(
                molecule_id
            )
            if existing is None:
                candidates_by_molecule[
                    molecule_id
                ] = candidate
                continue

            existing["supporting_activities"].extend(
                candidate["supporting_activities"]
            )

            if float(
                candidate["pchembl_value"]
            ) > float(existing["pchembl_value"]):
                candidate["supporting_activities"] = (
                    existing["supporting_activities"]
                )
                candidates_by_molecule[
                    molecule_id
                ] = candidate

    candidates = list(candidates_by_molecule.values())
    candidates.sort(
        key=lambda item: (
            -float(item.get("pchembl_value") or 0),
            -float(
                item.get("target_match_score") or 0
            ),
            str(item.get("chembl_molecule_id") or ""),
        )
    )

    if max_candidates > 0:
        candidates = candidates[:max_candidates]

    for rank, candidate in enumerate(
        candidates,
        start=1,
    ):
        candidate["retrieval_rank"] = rank

    trace = {
        "backend": "ChEMBL",
        "search_terms": search_terms,
        "target_count": len(targets),
        "targets": targets,
        "activity_count": activity_count,
        "accepted_activity_count": (
            accepted_activity_count
        ),
        "candidate_count": len(candidates),
    }

    return candidates, trace
