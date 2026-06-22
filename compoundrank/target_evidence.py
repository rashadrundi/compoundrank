from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


TARGET_RULES = [
    {
        "target_class": "viral protease",
        "enzyme_class": None,
        "keywords": [
            "protease",
            "peptidase",
            "proteinase",
        ],
        "recommended_compound_classes": [
            "viral protease inhibitors",
        ],
        "ligand_database_query_terms": [
            "viral protease inhibitor",
            "antiviral protease inhibitor",
        ],
        "docking_priority": "high",
        "notes": [
            "Proteases are common antiviral drug targets.",
            "Generic protease wording does not establish an aspartyl, cysteine, serine, or metalloprotease subclass.",
            "Catalytic subclass and family-specific ligand transfer require specific annotation or sequence evidence.",
        ],
    },
    {
        "target_class": "viral polymerase",
        "enzyme_class": "polymerase",
        "keywords": [
            "polymerase",
            "rna-dependent",
            "rdrp",
            "replicase",
            "reverse transcriptase",
            "transcriptase",
        ],
        "recommended_compound_classes": [
            "polymerase inhibitors",
            "nucleoside analog inhibitors",
            "non-nucleoside polymerase inhibitors",
        ],
        "ligand_database_query_terms": [
            "viral polymerase inhibitor",
            "RNA dependent RNA polymerase inhibitor",
            "nucleoside analog antiviral",
        ],
        "docking_priority": "high",
        "notes": [
            "Viral polymerases are major antiviral targets.",
            "Docking interpretation should distinguish active-site nucleotide analogs from allosteric inhibitors.",
        ],
    },
    {
        "target_class": "viral neuraminidase",
        "enzyme_class": "sialidase",
        "keywords": [
            "neuraminidase",
            "sialidase",
        ],
        "recommended_compound_classes": [
            "neuraminidase inhibitors",
        ],
        "ligand_database_query_terms": [
            "neuraminidase inhibitor",
            "oseltamivir",
            "zanamivir",
            "peramivir",
        ],
        "docking_priority": "high",
        "notes": [
            "Neuraminidase is a known antiviral target in influenza.",
        ],
    },
    {
        "target_class": "viral helicase",
        "enzyme_class": "helicase",
        "keywords": [
            "helicase",
            "nucleoside triphosphatase",
            "ntpase",
            "atpase",
        ],
        "recommended_compound_classes": [
            "viral helicase inhibitors",
            "ATPase-site inhibitors",
        ],
        "ligand_database_query_terms": [
            "viral helicase inhibitor",
            "helicase ATPase inhibitor",
        ],
        "docking_priority": "medium",
        "notes": [
            "Helicases can be druggable, but binding-site interpretation is usually more complex than known protease/polymerase systems.",
        ],
    },
    {
        "target_class": "viral integrase",
        "enzyme_class": "integrase",
        "keywords": [
            "integrase",
            "strand transfer",
            "retroviral integrase",
        ],
        "recommended_compound_classes": [
            "integrase inhibitors",
            "strand-transfer inhibitors",
        ],
        "ligand_database_query_terms": [
            "viral integrase inhibitor",
            "HIV integrase strand transfer inhibitor",
        ],
        "docking_priority": "high",
        "notes": [
            "Retroviral integrases are established antiviral targets.",
        ],
    },
    {
        "target_class": "viral methyltransferase",
        "enzyme_class": "methyltransferase",
        "keywords": [
            "methyltransferase",
            "mtase",
            "s-adenosylmethionine",
            "sam-dependent",
        ],
        "recommended_compound_classes": [
            "viral methyltransferase inhibitors",
            "SAM-competitive inhibitors",
        ],
        "ligand_database_query_terms": [
            "viral methyltransferase inhibitor",
            "SAM competitive methyltransferase inhibitor",
        ],
        "docking_priority": "medium",
        "notes": [
            "Methyltransferases may be targetable, but ligand selection needs careful cofactor-site context.",
        ],
    },
    {
        "target_class": "viral structural protein",
        "enzyme_class": None,
        "keywords": [
            "capsid",
            "coat protein",
            "nucleocapsid",
            "matrix protein",
            "envelope glycoprotein",
            "spike",
        ],
        "recommended_compound_classes": [
            "entry inhibitors",
            "capsid assembly modulators",
            "protein-protein interaction inhibitors",
        ],
        "ligand_database_query_terms": [
            "viral capsid inhibitor",
            "viral entry inhibitor",
            "viral protein protein interaction inhibitor",
        ],
        "docking_priority": "low",
        "notes": [
            "Structural proteins may be druggable, but simple small-molecule docking is often less reliable without a known pocket or mechanism.",
        ],
    },
]


SPECIAL_DOMAIN_RULES = [
    {
        "match_names": [
            "RVP",
        ],
        "match_accessions": [
            "pfam00077",
            "PF00077",
        ],
        "match_terms": [
            "retroviral aspartyl protease",
            "retroviral protease domain",
            "retropepsin",
        ],
        "label": "Retroviral aspartyl protease domain",
        "target_name": "HIV-like retroviral aspartyl protease",
        "target_class": "viral protease",
        "enzyme_class": "aspartyl protease",
        "viral_family": "Retroviridae-like / retroviral",
        "matched_terms": [
            "RVP",
            "pfam00077",
            "retroviral protease domain",
        ],
        "recommended_compound_classes": [
            "viral protease inhibitors",
            "HIV protease inhibitors",
            "aspartyl protease inhibitors",
        ],
        "ligand_database_query_terms": [
            "HIV protease inhibitor",
            "viral aspartyl protease inhibitor",
            "retroviral protease inhibitor",
        ],
        "docking_priority": "high",
        "notes": [
            "Catalytic aspartate and flap-region contacts are important for retroviral aspartyl proteases.",
        ],
        "confidence": "high",
        "confidence_reason": (
            "Specific CDD/Pfam retroviral protease domain evidence "
            "was detected. This is stronger than generic protease "
            "keyword matching."
        ),
    },
    {
        "match_names": [
            "3C-like protease",
            "3CLpro",
            "Mpro",
            "Peptidase C30",
            "Peptidase_C30",
            "main protease",
            "main proteinase",
        ],
        "match_accessions": [],
        "match_terms": [
            "3c-like protease",
            "3cl protease",
            "3clpro",
            "coronavirus main protease",
            "coronavirus main proteinase",
            "peptidase c30",
            "peptidase_c30",
            "coronavirus endopeptidase",
            "coronavirus polyprotein-processing protease",
        ],
        "label": "Coronavirus 3C-like main protease domain",
        "target_name": "coronavirus 3C-like main protease",
        "target_class": "viral protease",
        "enzyme_class": "cysteine protease",
        "viral_family": "Coronaviridae-like / coronavirus",
        "matched_terms": [
            "3C-like protease",
            "3CLpro",
            "main protease",
            "peptidase C30",
        ],
        "recommended_compound_classes": [
            "coronavirus main protease inhibitors",
            "viral cysteine protease inhibitors",
            "3CLpro inhibitors",
        ],
        "ligand_database_query_terms": [
            "coronavirus main protease inhibitor",
            "SARS-CoV-2 Mpro inhibitor",
            "3CLpro inhibitor",
            "peptidase C30 inhibitor",
        ],
        "docking_priority": "high",
        "notes": [
            "Coronavirus 3C-like proteases use a catalytic cysteine-histidine dyad.",
            "Covalent and noncovalent inhibitors should be distinguished during docking interpretation.",
        ],
        "confidence": "high",
        "confidence_reason": (
            "Specific coronavirus 3C-like protease or peptidase "
            "C30 annotation evidence was detected."
        ),
    },
]


def _flatten_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (int, float, bool)):
        return str(value)
    if isinstance(value, dict):
        return " ".join(_flatten_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(_flatten_text(item) for item in value)
    return str(value)


def _shorten(text: str, limit: int = 220) -> str:
    cleaned = " ".join(text.split())
    if len(cleaned) <= limit:
        return cleaned
    return cleaned[: limit - 3] + "..."


def _row_label(row: dict[str, Any]) -> str:
    for key in (
        "name",
        "description",
        "signature_desc",
        "interpro_description",
        "title",
        "accession",
        "hit_id",
        "query",
        "domain",
        "family",
        "member_database",
    ):
        value = row.get(key)
        if value:
            return _shorten(str(value), 140)
    return _shorten(_flatten_text(row), 140)


def _collect_hits(summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = summary.get("rows", {})
    hits: list[dict[str, Any]] = []

    if not isinstance(rows, dict):
        return hits

    for tool_name in ("cdd", "interpro", "vogdb"):
        tool_rows = rows.get(tool_name, [])
        if not isinstance(tool_rows, list):
            continue

        for index, row in enumerate(tool_rows[:20], start=1):
            if not isinstance(row, dict):
                continue

            hits.append(
                {
                    "tool": tool_name,
                    "rank": index,
                    "label": _row_label(row),
                    "text": _shorten(_flatten_text(row)),
                }
            )

    return hits


def _score_rules(search_text: str) -> list[dict[str, Any]]:
    text = search_text.lower()
    scored: list[dict[str, Any]] = []

    for rule in TARGET_RULES:
        matched_keywords = [
            keyword
            for keyword in rule["keywords"]
            if keyword.lower() in text
        ]
        scored.append(
            {
                "rule": rule,
                "score": len(matched_keywords),
                "matched_keywords": matched_keywords,
            }
        )

    scored.sort(
        key=lambda item: item["score"],
        reverse=True,
    )
    return scored


def _detect_special_domain(
    summary: dict[str, Any],
) -> dict[str, Any] | None:
    rows = summary.get("rows", {})

    if not isinstance(rows, dict):
        return None

    for tool_name, tool_rows in rows.items():
        if not isinstance(tool_rows, list):
            continue

        for row in tool_rows:
            if not isinstance(row, dict):
                continue

            row_name = str(
                row.get("name")
                or row.get("signature_desc")
                or row.get("description")
                or ""
            ).strip()

            row_accession = str(
                row.get("accession")
                or row.get("signature_accession")
                or row.get("interpro_accession")
                or ""
            ).strip()

            row_text = _flatten_text(row).casefold()

            for rule in SPECIAL_DOMAIN_RULES:
                match_names = [
                    str(value)
                    for value in rule.get(
                        "match_names",
                        [],
                    )
                ]
                match_accessions = [
                    str(value)
                    for value in rule.get(
                        "match_accessions",
                        [],
                    )
                ]
                match_terms = [
                    str(value)
                    for value in rule.get(
                        "match_terms",
                        [],
                    )
                ]

                name_match = any(
                    row_name.casefold()
                    == value.casefold()
                    for value in match_names
                )

                accession_match = any(
                    row_accession.casefold()
                    == value.casefold()
                    for value in match_accessions
                )

                text_match = any(
                    value.casefold() in row_text
                    for value in (
                        match_names
                        + match_accessions
                        + match_terms
                    )
                    if value
                )

                if not (
                    name_match
                    or accession_match
                    or text_match
                ):
                    continue

                return {
                    "label": rule["label"],
                    "tool": tool_name,
                    "hit_name": (
                        row_name
                        or _row_label(row)
                    ),
                    "accession": row_accession,
                    "start": row.get("start"),
                    "end": row.get("end"),
                    "evalue": row.get("evalue"),
                    "score": row.get("score"),
                    "target_name": rule[
                        "target_name"
                    ],
                    "target_class": rule[
                        "target_class"
                    ],
                    "enzyme_class": rule[
                        "enzyme_class"
                    ],
                    "viral_family": rule[
                        "viral_family"
                    ],
                    "matched_terms": rule[
                        "matched_terms"
                    ],
                    "recommended_compound_classes": (
                        rule[
                            "recommended_compound_classes"
                        ]
                    ),
                    "ligand_database_query_terms": (
                        rule[
                            "ligand_database_query_terms"
                        ]
                    ),
                    "docking_priority": rule[
                        "docking_priority"
                    ],
                    "notes": rule["notes"],
                    "confidence": rule[
                        "confidence"
                    ],
                    "confidence_reason": rule[
                        "confidence_reason"
                    ],
                }

    return None


def _confidence(
    best_score: int,
    total_hits: int,
    *,
    special_domain: dict[str, Any] | None = None,
) -> str:
    if special_domain is not None:
        return str(special_domain.get("confidence", "high"))
    if best_score >= 3 and total_hits >= 2:
        return "high"
    if best_score >= 2 or total_hits >= 3:
        return "medium"
    if best_score >= 1:
        return "low"
    return "unknown"


def _annotation_tool_statuses(
    summary: dict[str, Any],
) -> dict[str, str]:
    tool_names = (
        "cdd",
        "interpro",
        "vogdb",
    )

    statuses: dict[str, str] = {}

    explicit_statuses = summary.get(
        "tool_statuses",
        {},
    )

    tools = summary.get(
        "tools",
        {},
    )

    for tool_name in tool_names:
        status = None

        if isinstance(explicit_statuses, dict):
            explicit_status = explicit_statuses.get(
                tool_name
            )

            if explicit_status is not None:
                status = str(
                    explicit_status
                ).strip().lower()

        if status is None and isinstance(
            tools,
            dict,
        ):
            tool_result = tools.get(
                tool_name,
                {},
            )

            if isinstance(tool_result, dict):
                tool_status = tool_result.get(
                    "status"
                )

                if tool_status is not None:
                    status = str(
                        tool_status
                    ).strip().lower()

        statuses[tool_name] = (
            status or "unknown"
        )

    return statuses


def _annotation_tool_errors(
    summary: dict[str, Any],
) -> dict[str, str]:
    errors: dict[str, str] = {}

    explicit_errors = summary.get(
        "tool_errors",
        {},
    )

    if isinstance(explicit_errors, dict):
        for tool_name, error in (
            explicit_errors.items()
        ):
            if error:
                errors[str(tool_name)] = str(
                    error
                )

    tools = summary.get(
        "tools",
        {},
    )

    if isinstance(tools, dict):
        for tool_name, tool_result in (
            tools.items()
        ):
            if not isinstance(
                tool_result,
                dict,
            ):
                continue

            error = tool_result.get("error")

            if error and tool_name not in errors:
                errors[str(tool_name)] = str(
                    error
                )

    return errors


def _summarize_tool_error(
    error: Any,
    *,
    maximum_length: int = 300,
) -> str | None:
    if not error:
        return None

    lines = [
        line.strip()
        for line in str(error).splitlines()
        if line.strip()
    ]

    if not lines:
        return None

    preferred_fragments = (
        "bad file format",
        "no such file",
        "permission denied",
        "out of memory",
        "timed out",
        "timeout",
        "command not found",
        "exit code:",
        "runtimeerror:",
        "error:",
    )

    selected = None

    for line in reversed(lines):
        lowered = line.lower()

        if any(
            fragment in lowered
            for fragment in preferred_fragments
        ):
            selected = line
            break

    if selected is None:
        selected = lines[-1]

    if len(selected) > maximum_length:
        return (
            selected[: maximum_length - 3]
            + "..."
        )

    return selected


def _annotation_source_block(
    summary: dict[str, Any],
    *,
    source_fasta: str | None = None,
) -> dict[str, Any]:
    return {
        "job_id": summary.get("job_id"),
        "status": summary.get("status"),
        "source_fasta": source_fasta,
        "result_counts": summary.get(
            "result_counts",
            {},
        ),
        "tool_statuses": (
            _annotation_tool_statuses(
                summary
            )
        ),
        "tool_errors": (
            _annotation_tool_errors(
                summary
            )
        ),
    }


def _annotation_limitations(
    summary: dict[str, Any],
) -> list[str]:
    statuses = _annotation_tool_statuses(
        summary
    )
    errors = _annotation_tool_errors(
        summary
    )

    display_names = {
        "cdd": "CDD",
        "interpro": "InterPro",
        "vogdb": "VOGDB",
    }

    limitations: list[str] = []

    for tool_name, status in statuses.items():
        display_name = display_names.get(
            tool_name,
            tool_name,
        )

        if status == "failed":
            error_summary = (
                _summarize_tool_error(
                    errors.get(tool_name)
                )
            )

            message = (
                f"{display_name} failed; its "
                "zero usable results must not "
                "be interpreted as a successful "
                "no-hit result."
            )

            if error_summary:
                message += (
                    f" Error summary: "
                    f"{error_summary}"
                )

            limitations.append(message)

        elif status == "partial":
            limitations.append(
                f"{display_name} completed only "
                "partially; evidence from this "
                "tool may be incomplete."
            )

        elif status == "unknown":
            limitations.append(
                f"{display_name} did not provide "
                "an explicit execution status; "
                "its result count alone cannot "
                "distinguish no hits from failure."
            )

    return limitations


def _annotation_status_table_lines(
    source: dict[str, Any],
) -> list[str]:
    statuses = source.get(
        "tool_statuses",
        {},
    )
    counts = source.get(
        "result_counts",
        {},
    )
    errors = source.get(
        "tool_errors",
        {},
    )

    if not isinstance(statuses, dict):
        statuses = {}

    if not isinstance(counts, dict):
        counts = {}

    if not isinstance(errors, dict):
        errors = {}

    display_names = {
        "cdd": "CDD",
        "interpro": "InterPro",
        "vogdb": "VOGDB",
    }

    lines = [
        "| Tool | Status | Usable results | Error summary |",
        "|---|---|---:|---|",
    ]

    for tool_name in (
        "cdd",
        "interpro",
        "vogdb",
    ):
        status = statuses.get(
            tool_name,
            "unknown",
        )
        count = counts.get(
            tool_name,
            0,
        )
        error_summary = (
            _summarize_tool_error(
                errors.get(tool_name)
            )
            or ""
        )

        error_summary = (
            error_summary
            .replace("|", "\\|")
            .replace("\n", " ")
        )

        lines.append(
            f"| {display_names[tool_name]} "
            f"| {status} "
            f"| {count} "
            f"| {error_summary} |"
        )

    return lines


def _confidence_reasoning(
    *,
    confidence: str,
    best_score: int,
    total_hits: int,
    special_domain: dict[str, Any] | None,
    summary: dict[str, Any],
) -> list[str]:
    reasons: list[str] = []

    if special_domain is not None:
        reasons.append(str(special_domain["confidence_reason"]))
        reasons.append(
            "Special domain evidence: "
            f"{special_domain.get('label')} "
            f"({special_domain.get('hit_name')}/{special_domain.get('accession')}), "
            f"e-value={special_domain.get('evalue')}, "
            f"score={special_domain.get('score')}."
        )

    result_counts = summary.get(
        "result_counts",
        {},
    )

    if isinstance(result_counts, dict):
        reasons.append(
            "Usable annotation rows retained: "
            f"CDD={result_counts.get('cdd', 0)}, "
            f"InterPro={result_counts.get('interpro', 0)}, "
            f"VOGDB={result_counts.get('vogdb', 0)}."
        )

    tool_statuses = (
        _annotation_tool_statuses(
            summary
        )
    )
    tool_errors = (
        _annotation_tool_errors(
            summary
        )
    )

    display_names = {
        "cdd": "CDD",
        "interpro": "InterPro",
        "vogdb": "VOGDB",
    }

    if all(
        status == "unknown"
        for status in tool_statuses.values()
    ):
        reasons.append(
            "Per-tool execution statuses were "
            "not recorded; result counts alone "
            "cannot distinguish successful no-hit "
            "results from tool failures."
        )
    else:
        for tool_name in (
            "cdd",
            "interpro",
            "vogdb",
        ):
            status = tool_statuses.get(
                tool_name,
                "unknown",
            )
            count = (
                result_counts.get(
                    tool_name,
                    0,
                )
                if isinstance(
                    result_counts,
                    dict,
                )
                else 0
            )
            display_name = display_names[
                tool_name
            ]

            if status == "complete":
                reasons.append(
                    f"{display_name} completed "
                    f"successfully with {count} "
                    "usable result(s)."
                )

            elif status == "complete_no_hits":
                reasons.append(
                    f"{display_name} completed "
                    "successfully and returned no "
                    "hits."
                )

            elif status == "failed":
                error_summary = (
                    _summarize_tool_error(
                        tool_errors.get(
                            tool_name
                        )
                    )
                )

                reason = (
                    f"{display_name} failed; "
                    f"{count} usable result(s) "
                    "were retained."
                )

                if error_summary:
                    reason += (
                        f" Error summary: "
                        f"{error_summary}"
                    )

                reasons.append(reason)

            elif status == "partial":
                reasons.append(
                    f"{display_name} completed "
                    f"partially with {count} "
                    "usable result(s)."
                )

            elif status == "skipped":
                reasons.append(
                    f"{display_name} was skipped."
                )

            else:
                reasons.append(
                    f"{display_name} execution "
                    "status is unknown."
                )

    status = summary.get("status")

    if status and status != "complete":
        reasons.append(
            f"Overall CPU annotation status was "
            f"{status}; confidence is based only "
            "on available successful tool outputs."
        )

    if not reasons:
        reasons.append(
            f"Confidence={confidence} based on {best_score} matched keyword(s) across {total_hits} supporting hit(s)."
        )

    return reasons


def _default_unknown_evidence(summary: dict[str, Any], hits: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "schema_version": "target_evidence.v0.2",
        "source": _annotation_source_block(
            summary,
        ),
        "target_interpretation": {
            "target_name": "unknown",
            "target_class": "unknown",
            "enzyme_class": None,
            "viral_family": "unknown",
            "predicted_function": "No confident functional class inferred from current annotation summary.",
            "docking_priority": "low",
            "evidence_confidence": "unknown",
        },
        "evidence": {
            "matched_keywords": [],
            "supporting_hits": hits[:10],
        },
        "future_ligand_database_query": {
            "recommended_compound_classes": [],
            "query_terms": [],
            "status": "not_queried",
            "notes": [
                "No ligand database query was performed.",
                "Target evidence did not identify a confident known antiviral target class.",
            ],
        },
        "active_site_or_motif_notes": [],
        "limitations": [
            "This target evidence packet is generated from computational annotation only.",
            "It does not prove target identity, binding, inhibition, or antiviral efficacy.",
            "Experimental validation and literature review remain required.",
        ] + _annotation_limitations(
            summary
        ),
        "recommended_next_action": "Review annotation hits manually before selecting ligands for docking.",
    }


def build_target_evidence(
    summary: dict[str, Any],
    *,
    source_fasta: str | None = None,
) -> dict[str, Any]:
    hits = _collect_hits(summary)
    search_text = _flatten_text(summary)
    scored_rules = _score_rules(search_text)

    best = (
        scored_rules[0]
        if scored_rules
        else {
            "score": 0,
            "rule": None,
            "matched_keywords": [],
        }
    )

    best_score = int(
        best.get("score") or 0
    )
    special_domain = _detect_special_domain(
        summary
    )

    if (
        best_score <= 0
        and special_domain is None
    ):
        return _default_unknown_evidence(
            summary,
            hits,
        )

    if special_domain is not None:
        interpretation_source = special_domain

        target_name = str(
            special_domain["target_name"]
        )
        target_class = str(
            special_domain["target_class"]
        )
        enzyme_class = special_domain.get(
            "enzyme_class"
        )
        viral_family = str(
            special_domain["viral_family"]
        )
        docking_priority = str(
            special_domain["docking_priority"]
        )

        predicted_function = (
            f"Likely {target_name} based on "
            "specific annotation-domain evidence."
        )
    else:
        rule = best.get("rule")

        if not isinstance(rule, dict):
            return _default_unknown_evidence(
                summary,
                hits,
            )

        interpretation_source = rule
        target_name = str(
            rule["target_class"]
        )
        target_class = str(
            rule["target_class"]
        )
        enzyme_class = rule.get(
            "enzyme_class"
        )
        viral_family = "unknown"
        docking_priority = str(
            rule["docking_priority"]
        )

        predicted_function = (
            f"Likely {target_class} based on "
            "annotation and homology keyword "
            "evidence. Catalytic subclass remains "
            "unresolved unless specific evidence "
            "is present."
        )

    confidence = _confidence(
        best_score,
        len(hits),
        special_domain=special_domain,
    )

    matched_keywords = list(
        dict.fromkeys(
            list(
                best.get(
                    "matched_keywords",
                    [],
                )
            )
            + (
                list(
                    special_domain.get(
                        "matched_terms",
                        [],
                    )
                )
                if special_domain is not None
                else []
            )
        )
    )

    return {
        "schema_version": (
            "target_evidence.v0.2"
        ),
        "source": _annotation_source_block(
            summary,
            source_fasta=source_fasta,
        ),
        "target_interpretation": {
            "target_name": target_name,
            "target_class": target_class,
            "enzyme_class": enzyme_class,
            "viral_family": viral_family,
            "predicted_function": (
                predicted_function
            ),
            "docking_priority": (
                docking_priority
            ),
            "evidence_confidence": (
                confidence
            ),
        },
        "evidence": {
            "matched_keywords": (
                matched_keywords
            ),
            "special_domain_evidence": (
                special_domain
            ),
            "confidence_reasoning": (
                _confidence_reasoning(
                    confidence=confidence,
                    best_score=best_score,
                    total_hits=len(hits),
                    special_domain=(
                        special_domain
                    ),
                    summary=summary,
                )
            ),
            "supporting_hits": hits[:10],
        },
        "future_ligand_database_query": {
            "recommended_compound_classes": (
                interpretation_source[
                    "recommended_compound_classes"
                ]
            ),
            "query_terms": (
                interpretation_source[
                    "ligand_database_query_terms"
                ]
            ),
            "status": "not_queried",
            "notes": [
                "No ligand database query was performed in this stage.",
                "These terms describe reference ligand-evidence searches; the submitted FASTA remains the docking target.",
            ],
        },
        "active_site_or_motif_notes": (
            interpretation_source["notes"]
        ),
        "limitations": [
            "This target evidence packet is generated from computational annotation only.",
            "It does not replace or redefine the submitted FASTA target.",
            "It does not prove binding, inhibition, or antiviral efficacy.",
            "Docking should be interpreted as hypothesis generation, not validation.",
            "Experimental validation and literature review remain required.",
        ] + _annotation_limitations(
            summary
        ),
        "recommended_next_action": (
            "Use annotation and sequence similarity "
            "to retrieve transferable ligand evidence, "
            "then dock candidates against the structure "
            "generated from the submitted FASTA."
        ),
    }


def write_target_evidence_outputs(
    evidence: dict[str, Any],
    output_dir: Path,
) -> dict[str, str]:
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    json_path = output_dir / "target_evidence.json"
    report_path = output_dir / "target_evidence_report.md"

    json_path.write_text(
        json.dumps(
            evidence,
            indent=2,
            sort_keys=True,
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )

    report_path.write_text(
        render_target_evidence_report(evidence),
        encoding="utf-8",
    )

    return {
        "target_evidence": str(json_path),
        "target_evidence_report": str(report_path),
    }


def render_target_evidence_report(evidence: dict[str, Any]) -> str:
    interpretation = evidence.get("target_interpretation", {})
    future_query = evidence.get("future_ligand_database_query", {})
    evidence_block = evidence.get("evidence", {})
    source = evidence.get("source", {})

    lines = [
        "# Target Evidence Report",
        "",
        "## Interpretation",
        "",
        f"- Target name: {interpretation.get('target_name')}",
        f"- Target class: {interpretation.get('target_class')}",
        f"- Enzyme class: {interpretation.get('enzyme_class')}",
        f"- Viral family evidence: {interpretation.get('viral_family')}",
        f"- Docking priority: {interpretation.get('docking_priority')}",
        f"- Evidence confidence: {interpretation.get('evidence_confidence')}",
        f"- Predicted function: {interpretation.get('predicted_function')}",
        "",
        "## Annotation source status",
        "",
        f"- Overall CPU annotation status: {source.get('status')}",
        "",
        *_annotation_status_table_lines(
            source
        ),
        "",
        "## Evidence assessment",
        "",
    ]

    confidence_reasoning = evidence_block.get("confidence_reasoning", [])
    if confidence_reasoning:
        for reason in confidence_reasoning:
            lines.append(f"- {reason}")
    else:
        lines.append("- No confidence reasoning was recorded.")

    special_domain = evidence_block.get("special_domain_evidence")
    if special_domain:
        lines += [
            "",
            "## Specific domain evidence",
            "",
            f"- Label: {special_domain.get('label')}",
            f"- Tool: {special_domain.get('tool')}",
            f"- Hit: {special_domain.get('hit_name')}",
            f"- Accession: {special_domain.get('accession')}",
            f"- Coordinates: {special_domain.get('start')}–{special_domain.get('end')}",
            f"- E-value: {special_domain.get('evalue')}",
            f"- Score: {special_domain.get('score')}",
        ]

    lines += [
        "",
        "## Matched keywords and terms",
        "",
    ]

    matched_keywords = evidence_block.get("matched_keywords", [])
    if matched_keywords:
        for keyword in matched_keywords:
            lines.append(f"- {keyword}")
    else:
        lines.append("- None")

    lines += [
        "",
        "## Supporting annotation hits",
        "",
    ]

    supporting_hits = evidence_block.get("supporting_hits", [])
    if supporting_hits:
        for hit in supporting_hits:
            lines.append(
                f"- {hit.get('tool')} hit {hit.get('rank')}: {hit.get('label')}"
            )
    else:
        lines.append("- No supporting hits were available.")

    lines += [
        "",
        "## Future ligand database query terms",
        "",
    ]

    query_terms = future_query.get("query_terms", [])
    if query_terms:
        for term in query_terms:
            lines.append(f"- {term}")
    else:
        lines.append("- No ligand query terms recommended.")

    lines += [
        "",
        "## Limitations",
        "",
    ]

    for limitation in evidence.get("limitations", []):
        lines.append(f"- {limitation}")

    lines.append("")
    return "\n".join(lines)


def build_and_write_target_evidence(
    *,
    summary_path: Path,
    output_dir: Path,
    source_fasta: Path | None = None,
) -> dict[str, str]:
    summary = json.loads(
        Path(summary_path).read_text(
            encoding="utf-8",
        )
    )

    if not isinstance(summary, dict):
        raise RuntimeError("Target evidence summary input must be a JSON object.")

    evidence = build_target_evidence(
        summary,
        source_fasta=str(source_fasta) if source_fasta is not None else None,
    )

    return write_target_evidence_outputs(
        evidence,
        Path(output_dir),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build target_evidence.json from homolog_search_summary.json."
    )
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, default=None)
    args = parser.parse_args()

    outputs = build_and_write_target_evidence(
        summary_path=args.summary,
        output_dir=args.output_dir,
        source_fasta=args.fasta,
    )

    for label, path in outputs.items():
        print(f"{label}: {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
