from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


TARGET_RULES = [
    {
        "target_class": "viral protease",
        "enzyme_class": "aspartyl protease",
        "keywords": [
            "protease",
            "peptidase",
            "aspartyl",
            "aspartic",
            "retropepsin",
            "hiv protease",
            "proteinase",
        ],
        "recommended_compound_classes": [
            "viral protease inhibitors",
            "HIV protease inhibitors",
            "aspartyl protease inhibitors",
        ],
        "ligand_database_query_terms": [
            "HIV protease inhibitor",
            "viral aspartyl protease inhibitor",
            "protease inhibitor antiviral",
        ],
        "docking_priority": "high",
        "notes": [
            "Proteases are common antiviral drug targets.",
            "For HIV-like aspartyl proteases, catalytic aspartate and flap-region contacts are biologically important.",
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


def _confidence(best_score: int, total_hits: int) -> str:
    if best_score >= 3 and total_hits >= 2:
        return "high"
    if best_score >= 2 or total_hits >= 3:
        return "medium"
    if best_score >= 1:
        return "low"
    return "unknown"


def _default_unknown_evidence(summary: dict[str, Any], hits: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "schema_version": "target_evidence.v0.1",
        "source": {
            "job_id": summary.get("job_id"),
            "status": summary.get("status"),
            "result_counts": summary.get("result_counts", {}),
        },
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
        ],
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
    best = scored_rules[0] if scored_rules else {"score": 0, "rule": None, "matched_keywords": []}

    if best["score"] <= 0:
        evidence = _default_unknown_evidence(summary, hits)
    else:
        rule = best["rule"]
        confidence = _confidence(best["score"], len(hits))

        target_name = rule["target_class"]
        if rule["target_class"] == "viral protease" and any(
            keyword in search_text.lower()
            for keyword in ("hiv", "retropepsin", "aspartyl", "aspartic")
        ):
            target_name = "HIV-like viral aspartyl protease"

        evidence = {
            "schema_version": "target_evidence.v0.1",
            "source": {
                "job_id": summary.get("job_id"),
                "status": summary.get("status"),
                "source_fasta": source_fasta,
                "result_counts": summary.get("result_counts", {}),
            },
            "target_interpretation": {
                "target_name": target_name,
                "target_class": rule["target_class"],
                "enzyme_class": rule["enzyme_class"],
                "viral_family": "unknown",
                "predicted_function": (
                    f"Likely {rule['target_class']} based on annotation and homology keyword evidence."
                ),
                "docking_priority": rule["docking_priority"],
                "evidence_confidence": confidence,
            },
            "evidence": {
                "matched_keywords": best["matched_keywords"],
                "supporting_hits": hits[:10],
            },
            "future_ligand_database_query": {
                "recommended_compound_classes": rule["recommended_compound_classes"],
                "query_terms": rule["ligand_database_query_terms"],
                "status": "not_queried",
                "notes": [
                    "No ligand database query was performed in this stage.",
                    "These query terms are intended for a future external ligand database service.",
                ],
            },
            "active_site_or_motif_notes": rule["notes"],
            "limitations": [
                "This target evidence packet is generated from computational annotation only.",
                "It does not prove target identity, binding, inhibition, or antiviral efficacy.",
                "Docking should be interpreted as hypothesis generation, not validation.",
                "Experimental validation and literature review remain required.",
            ],
            "recommended_next_action": (
                "Use this target evidence to guide ligand database queries and docking setup, "
                "then manually review whether docking contacts match biologically important regions."
            ),
        }

    return evidence


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

    lines = [
        "# Target Evidence Report",
        "",
        "## Interpretation",
        "",
        f"- Target name: {interpretation.get('target_name')}",
        f"- Target class: {interpretation.get('target_class')}",
        f"- Enzyme class: {interpretation.get('enzyme_class')}",
        f"- Docking priority: {interpretation.get('docking_priority')}",
        f"- Evidence confidence: {interpretation.get('evidence_confidence')}",
        "",
        "## Matched keywords",
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
