from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None

    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return None

    return data


def _format_value(value: Any) -> str:
    if value is None:
        return "unknown"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def _extract_float(text: str) -> float | None:
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", text)
    if not match:
        return None

    try:
        return float(match.group(0))
    except ValueError:
        return None


def _extract_float_after_keyword(text: str, keyword: str) -> float | None:
    upper = text.upper()
    keyword_upper = keyword.upper()

    index = upper.find(keyword_upper)
    if index == -1:
        return None

    value_text = text[index + len(keyword):].strip()
    return _extract_float(value_text)


def _parse_pdb_remark_metadata(path: Path) -> dict[str, Any]:
    metadata: dict[str, Any] = {}

    try:
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()
    except OSError:
        return metadata

    for line in lines[:250]:
        if not line.startswith("REMARK"):
            continue

        clean = " ".join(line.split())
        upper = clean.upper()

        if "GNINA CNN SCORE" in upper:
            value = _extract_float_after_keyword(clean, "GNINA CNN SCORE")
            if value is not None:
                metadata["gnina_cnn_score"] = value

        elif "POSE CONFIDENCE" in upper:
            metadata["pose_confidence"] = clean.split("POSE CONFIDENCE", 1)[-1].strip()

        elif "POCKET ID" in upper:
            metadata["pocket_id"] = clean.split("POCKET ID", 1)[-1].strip()

        elif "POCKET SOURCE" in upper:
            metadata["pocket_source"] = clean.split("POCKET SOURCE", 1)[-1].strip()

        elif "FPOCKET SCORE" in upper:
            value = _extract_float_after_keyword(clean, "FPOCKET SCORE")
            if value is not None:
                metadata["fpocket_score"] = value

    return metadata


def _parse_hypothesis_filename(path: Path) -> dict[str, Any]:
    stem = path.stem
    pieces = stem.split("__")

    parsed: dict[str, Any] = {
        "file": path.name,
        "compound": "unknown",
        "hypothesis": "unknown",
        "pocket": "unknown",
    }

    if len(pieces) >= 2:
        compound_piece = pieces[1]
        parsed["compound"] = compound_piece.replace("_", " ")

    for piece in pieces:
        if piece.startswith("hypothesis_"):
            parsed["hypothesis"] = piece.replace("hypothesis_", "")
        elif piece.startswith("fpocket_"):
            parsed["pocket"] = piece

    rank_match = re.match(r"^(\d+)$", pieces[0]) if pieces else None
    if rank_match:
        parsed["compound_rank"] = int(rank_match.group(1))

    return parsed


def collect_pdb_hypotheses(output_dir: Path) -> list[dict[str, Any]]:
    hypotheses: list[dict[str, Any]] = []

    for pdb_path in sorted(output_dir.glob("*.pdb")):
        item = _parse_hypothesis_filename(pdb_path)
        item.update(_parse_pdb_remark_metadata(pdb_path))
        item["path"] = str(pdb_path)
        hypotheses.append(item)

    hypotheses.sort(
        key=lambda item: (
            item.get("compound_rank", 9999),
            item.get("compound", ""),
            item.get("pocket", ""),
            item.get("hypothesis", ""),
            item.get("file", ""),
        )
    )

    return hypotheses


def _summarize_annotation_error(
    error: Any,
    *,
    maximum_length: int = 220,
) -> str:
    if not error:
        return ""

    lines = [
        line.strip()
        for line in str(error).splitlines()
        if line.strip()
    ]

    if not lines:
        return ""

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

    selected = (
        selected
        .replace("|", "\\|")
        .replace("\n", " ")
    )

    if len(selected) > maximum_length:
        selected = (
            selected[: maximum_length - 3]
            + "..."
        )

    return selected


def _render_annotation_status_table(
    source: dict[str, Any],
) -> list[str]:
    overall_status = source.get("status")
    result_counts = source.get(
        "result_counts",
        {},
    )
    tool_statuses = source.get(
        "tool_statuses",
        {},
    )
    tool_errors = source.get(
        "tool_errors",
        {},
    )

    if not isinstance(result_counts, dict):
        result_counts = {}

    if not isinstance(tool_statuses, dict):
        tool_statuses = {}

    if not isinstance(tool_errors, dict):
        tool_errors = {}

    lines = [
        "### Annotation Execution Status",
        "",
        (
            "- Overall CPU annotation status: "
            f"{_format_value(overall_status)}"
        ),
        "",
    ]

    if not tool_statuses:
        lines += [
            (
                "- Per-tool execution statuses were "
                "not available in this legacy evidence file."
            ),
            (
                "- Result counts: "
                f"{_format_value(result_counts)}"
            ),
            "",
        ]
        return lines

    display_names = {
        "cdd": "CDD",
        "interpro": "InterPro",
        "vogdb": "VOGDB",
    }

    lines += [
        "| Tool | Status | Usable results | Error summary |",
        "|---|---|---:|---|",
    ]

    failed_tools: list[str] = []
    unknown_tools: list[str] = []

    for tool_name in (
        "cdd",
        "interpro",
        "vogdb",
    ):
        display_name = display_names[tool_name]
        status = str(
            tool_statuses.get(
                tool_name,
                "unknown",
            )
        )
        count = result_counts.get(
            tool_name,
            0,
        )
        error_summary = (
            _summarize_annotation_error(
                tool_errors.get(tool_name)
            )
        )

        lines.append(
            f"| {display_name} "
            f"| {status} "
            f"| {count} "
            f"| {error_summary} |"
        )

        if status == "failed":
            failed_tools.append(display_name)

        if status == "unknown":
            unknown_tools.append(display_name)

    lines.append("")

    if failed_tools:
        lines.append(
            "- Failed annotation tools: "
            + ", ".join(failed_tools)
            + ". Their zero usable-result counts "
            "must not be interpreted as successful "
            "no-hit results."
        )
        lines.append("")

    if unknown_tools:
        lines.append(
            "- Unknown annotation status: "
            + ", ".join(unknown_tools)
            + ". Result counts alone cannot confirm "
            "successful execution."
        )
        lines.append("")

    return lines


def _render_target_section(
    target_evidence: dict[str, Any] | None,
) -> list[str]:
    if target_evidence is None:
        return [
            "## Target Evidence",
            "",
            (
                "No target evidence file was "
                "available for this run."
            ),
            "",
        ]

    interpretation = target_evidence.get(
        "target_interpretation",
        {},
    )
    evidence = target_evidence.get(
        "evidence",
        {},
    )
    future_query = target_evidence.get(
        "future_ligand_database_query",
        {},
    )
    source = target_evidence.get(
        "source",
        {},
    )
    limitations = target_evidence.get(
        "limitations",
        [],
    )

    if not isinstance(interpretation, dict):
        interpretation = {}

    if not isinstance(evidence, dict):
        evidence = {}

    if not isinstance(future_query, dict):
        future_query = {}

    if not isinstance(source, dict):
        source = {}

    if not isinstance(limitations, list):
        limitations = []

    lines = [
        "## Target Evidence",
        "",
        (
            "- Target name: "
            f"{_format_value(interpretation.get('target_name'))}"
        ),
        (
            "- Target class: "
            f"{_format_value(interpretation.get('target_class'))}"
        ),
        (
            "- Enzyme class: "
            f"{_format_value(interpretation.get('enzyme_class'))}"
        ),
        (
            "- Viral family evidence: "
            f"{_format_value(interpretation.get('viral_family'))}"
        ),
        (
            "- Evidence confidence: "
            f"{_format_value(interpretation.get('evidence_confidence'))}"
        ),
        (
            "- Docking priority: "
            f"{_format_value(interpretation.get('docking_priority'))}"
        ),
        "",
    ]

    lines.extend(
        _render_annotation_status_table(
            source
        )
    )

    confidence_reasoning = evidence.get(
        "confidence_reasoning",
        [],
    )

    if isinstance(
        confidence_reasoning,
        list,
    ) and confidence_reasoning:
        lines += [
            "### Evidence Reasoning",
            "",
        ]

        for reason in confidence_reasoning:
            lines.append(f"- {reason}")

        lines.append("")

    annotation_limitations = [
        str(limitation)
        for limitation in limitations
        if any(
            term in str(limitation).lower()
            for term in (
                "cdd",
                "interpro",
                "vogdb",
                "annotation",
                "no-hit",
                "execution status",
            )
        )
    ]

    if annotation_limitations:
        lines += [
            "### Annotation Limitations",
            "",
        ]

        for limitation in annotation_limitations:
            lines.append(f"- {limitation}")

        lines.append("")

    special_domain = evidence.get(
        "special_domain_evidence"
    )

    if isinstance(
        special_domain,
        dict,
    ) and special_domain:
        lines += [
            "### Specific Domain Evidence",
            "",
            (
                "- Label: "
                f"{_format_value(special_domain.get('label'))}"
            ),
            (
                "- Tool: "
                f"{_format_value(special_domain.get('tool'))}"
            ),
            (
                "- Hit: "
                f"{_format_value(special_domain.get('hit_name'))}"
            ),
            (
                "- Accession: "
                f"{_format_value(special_domain.get('accession'))}"
            ),
            (
                "- Coordinates: "
                f"{_format_value(special_domain.get('start'))}"
                "–"
                f"{_format_value(special_domain.get('end'))}"
            ),
            (
                "- E-value: "
                f"{_format_value(special_domain.get('evalue'))}"
            ),
            (
                "- Score: "
                f"{_format_value(special_domain.get('score'))}"
            ),
            "",
        ]

    query_terms = future_query.get(
        "query_terms",
        [],
    )

    lines += [
        "### Future Ligand Database Query Terms",
        "",
    ]

    if isinstance(query_terms, list) and query_terms:
        for term in query_terms:
            lines.append(f"- {term}")
    else:
        lines.append(
            "- No ligand database query terms "
            "were recommended."
        )

    lines.append("")

    return lines


def _relative_or_original(output_dir: Path, value: Any) -> str:
    if value is None:
        return ""

    text = str(value).strip()
    if not text:
        return ""

    path = Path(text)
    try:
        return str(path.relative_to(output_dir))
    except ValueError:
        return text


def _read_candidate_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def _render_ligand_retrieval_section(output_dir: Path) -> list[str]:
    stage4a_dir = output_dir / "stage4a_compound_retrieval"
    candidate_csv = stage4a_dir / "candidate_ligands.csv"
    docking_manifest = stage4a_dir / "docking_manifest.csv"
    ligand_report = stage4a_dir / "ligand_search_report.md"

    candidates = _read_candidate_csv(candidate_csv)

    if not candidates and not stage4a_dir.exists():
        return []

    lines = [
        "## Stage 4A Ligand Retrieval",
        "",
    ]

    if not candidates:
        lines += [
            "Stage 4A output directory was found, but no candidate ligand CSV was available.",
            "",
        ]
        return lines

    lines += [
        f"- Candidate table: `{_relative_or_original(output_dir, candidate_csv)}`",
        f"- Docking manifest: `{_relative_or_original(output_dir, docking_manifest)}`" if docking_manifest.exists() else "- Docking manifest: unavailable",
        f"- Ligand search report: `{_relative_or_original(output_dir, ligand_report)}`" if ligand_report.exists() else "- Ligand search report: unavailable",
        f"- Retrieved candidate count: {len(candidates)}",
        "",
        "| Retrieval rank | Compound | Design status | Evidence | Rule | PubChem CID | Structure status |",
        "|---:|---|---|---|---|---|---|",
    ]

    def _rank_value(row: dict[str, str]) -> int:
        try:
            return int(str(row.get("retrieval_rank", "9999")))
        except ValueError:
            return 9999

    candidates = sorted(candidates, key=_rank_value)

    for row in candidates:
        lines.append(
            "| "
            f"{_format_value(row.get('retrieval_rank'))} | "
            f"{_format_value(row.get('compound_name'))} | "
            f"{_format_value(row.get('design_status'))} | "
            f"{_format_value(row.get('evidence_level'))} | "
            f"{_format_value(row.get('retrieval_rule_id'))} | "
            f"{_format_value(row.get('pubchem_cid'))} | "
            f"{_format_value(row.get('structure_fetch_status'))} |"
        )

    lines += [
        "",
        "### Retrieval Basis",
        "",
    ]

    seen_basis: set[tuple[str, str, str, str]] = set()
    for row in candidates:
        basis = (
            str(row.get("retrieval_rule_id", "")),
            str(row.get("target_family_basis", "")),
            str(row.get("special_domain_label", "")),
            str(row.get("special_domain_accession", "")),
        )
        if basis in seen_basis:
            continue
        seen_basis.add(basis)

        rule, family, label, accession = basis
        lines += [
            f"- Rule: `{_format_value(rule)}`",
            f"  - Target family basis: {_format_value(family)}",
            f"  - Special domain: {_format_value(label)} ({_format_value(accession)})",
        ]

    lines += [
        "",
        "### Candidate Reasoning",
        "",
    ]

    for row in candidates[:10]:
        lines += [
            f"#### {_format_value(row.get('compound_name'))}",
            "",
            f"- Retrieval reason: {_format_value(row.get('retrieval_reason'))}",
            f"- Local SDF: `{_relative_or_original(output_dir, row.get('local_sdf_path'))}`",
            "",
        ]

    if len(candidates) > 10:
        lines += [
            f"_Only the first 10 candidates are shown here. See `{_relative_or_original(output_dir, candidate_csv)}` for the full table._",
            "",
        ]

    return lines

def _render_docking_section(hypotheses: list[dict[str, Any]]) -> list[str]:
    lines = [
        "## Docking Hypotheses",
        "",
    ]

    if not hypotheses:
        lines += [
            "No final PDB hypothesis files were found.",
            "",
        ]
        return lines

    compound_best: dict[str, dict[str, Any]] = {}
    for item in hypotheses:
        compound = str(item.get("compound", "unknown"))
        current = compound_best.get(compound)
        score = item.get("gnina_cnn_score")

        if current is None:
            compound_best[compound] = item
        elif isinstance(score, float) and score > float(current.get("gnina_cnn_score", -9999)):
            compound_best[compound] = item

    ranked_compounds = sorted(
        compound_best.values(),
        key=lambda item: float(item.get("gnina_cnn_score", -9999)),
        reverse=True,
    )

    lines += [
        "### Compound Priority",
        "",
        "| Rank | Compound | Best GNINA CNN score | Confidence | Best hypothesis file |",
        "|---:|---|---:|---|---|",
    ]

    for index, item in enumerate(ranked_compounds, start=1):
        lines.append(
            "| "
            f"{index} | "
            f"{_format_value(item.get('compound'))} | "
            f"{_format_value(item.get('gnina_cnn_score'))} | "
            f"{_format_value(item.get('pose_confidence'))} | "
            f"`{_format_value(item.get('file'))}` |"
        )

    lines += [
        "",
        "### Final Hypothesis Files",
        "",
        "| File | Compound | Hypothesis | Pocket | GNINA CNN score | Pose confidence |",
        "|---|---|---:|---|---:|---|",
    ]

    for item in hypotheses:
        lines.append(
            "| "
            f"`{_format_value(item.get('file'))}` | "
            f"{_format_value(item.get('compound'))} | "
            f"{_format_value(item.get('hypothesis'))} | "
            f"{_format_value(item.get('pocket_id') or item.get('pocket'))} | "
            f"{_format_value(item.get('gnina_cnn_score'))} | "
            f"{_format_value(item.get('pose_confidence'))} |"
        )

    lines.append("")
    return lines



def _render_docking_attempt_summary_section(output_dir: Path) -> list[str]:
    summary_csv = output_dir / "docking_attempt_summary.csv"
    rows = _read_candidate_csv(summary_csv)

    if not rows:
        return []

    lines = [
        "## Docking Attempt Summary",
        "",
        f"- Attempt summary table: `{_relative_or_original(output_dir, summary_csv)}`",
        "",
        "| Compound | Pocket | Raw poses | Accepted | Rejected | Status | Best raw CNN | Best accepted CNN |",
        "|---|---|---:|---:|---:|---|---:|---:|",
    ]

    def _as_int(row: dict[str, str], key: str) -> int:
        try:
            return int(str(row.get(key, "0")))
        except ValueError:
            return 0

    failed_rows: list[dict[str, str]] = []

    for row in rows:
        accepted = _as_int(row, "accepted_poses")
        if accepted == 0 or str(row.get("status", "")) != "accepted":
            failed_rows.append(row)

        lines.append(
            "| "
            f"{_format_value(row.get('compound'))} | "
            f"{_format_value(row.get('pocket'))} | "
            f"{_format_value(row.get('raw_poses'))} | "
            f"{_format_value(row.get('accepted_poses'))} | "
            f"{_format_value(row.get('rejected_poses'))} | "
            f"{_format_value(row.get('status'))} | "
            f"{_format_value(row.get('best_raw_cnn_score'))} | "
            f"{_format_value(row.get('best_accepted_cnn_score'))} |"
        )

    lines.append("")

    if failed_rows:
        lines += [
            "### Attempted Ligands With No Accepted Final Pose",
            "",
            "| Compound | Pocket | Raw poses | Accepted | Rejected | Status | Best raw CNN |",
            "|---|---|---:|---:|---:|---|---:|",
        ]

        for row in failed_rows:
            lines.append(
                "| "
                f"{_format_value(row.get('compound'))} | "
                f"{_format_value(row.get('pocket'))} | "
                f"{_format_value(row.get('raw_poses'))} | "
                f"{_format_value(row.get('accepted_poses'))} | "
                f"{_format_value(row.get('rejected_poses'))} | "
                f"{_format_value(row.get('status'))} | "
                f"{_format_value(row.get('best_raw_cnn_score'))} |"
            )

        lines += [
            "",
            "These ligands were attempted but did not produce a final PoseBusters-accepted hypothesis file.",
            "",
        ]

    return lines


def _coerce_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None

    if isinstance(value, (int, float)):
        return float(value)

    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _coerce_int(value: Any) -> int | None:
    number = _coerce_float(value)

    if number is None:
        return None

    return int(number)


def _selection_rank_value(row: dict[str, Any]) -> int:
    rank = _coerce_int(row.get("selection_rank"))

    if rank is None:
        return 9999

    return rank


def _selected_attempt(
    rows: list[dict[str, Any]],
) -> dict[str, Any] | None:
    for row in rows:
        if row.get("selected") is True:
            return row

    for row in rows:
        if _selection_rank_value(row) == 1:
            return row

    return rows[0] if rows else None


def _selection_confidence(
    rows: list[dict[str, Any]],
) -> tuple[float | None, str]:
    selected = _selected_attempt(rows)

    if selected is None:
        return None, "unavailable"

    selected_score = _coerce_float(
        selected.get("top_cnn_score")
    )

    alternative_scores = [
        score
        for row in rows
        if row is not selected
        for score in [
            _coerce_float(
                row.get("top_cnn_score")
            )
        ]
        if score is not None
    ]

    if selected_score is None:
        return None, "unavailable"

    if not alternative_scores:
        return None, "single pocket"

    margin = selected_score - max(
        alternative_scores
    )

    if margin >= 0.20:
        label = "high separation"
    elif margin >= 0.10:
        label = "moderate separation"
    else:
        label = "low separation"

    return margin, label


def _format_box_triplet(
    pocket: dict[str, Any],
    prefix: str,
) -> str:
    values = [
        _coerce_float(
            pocket.get(f"{prefix}_{axis}")
        )
        for axis in ("x", "y", "z")
    ]

    if any(value is None for value in values):
        return "unknown"

    return ", ".join(
        f"{value:.3f}"
        for value in values
        if value is not None
    )


def _render_pocket_selection_section(
    output_dir: Path,
) -> list[str]:
    summary_path = (
        output_dir
        / "pocket_selection_summary.json"
    )
    definitions_path = (
        output_dir
        / "pocket_definitions.json"
    )
    summary_csv = (
        output_dir
        / "pocket_selection_summary.csv"
    )

    summary = _load_json(summary_path)
    definitions = _load_json(definitions_path)

    if summary is None and definitions is None:
        return []

    lines = [
        "## Pocket Detection and Selection",
        "",
    ]

    if summary is not None:
        attempts_value = summary.get(
            "attempts",
            [],
        )
        attempts = [
            row
            for row in attempts_value
            if isinstance(row, dict)
        ]

        grouped: dict[
            str,
            list[dict[str, Any]],
        ] = {}

        for row in attempts:
            compound = str(
                row.get(
                    "compound",
                    "unknown",
                )
            )
            grouped.setdefault(
                compound,
                [],
            ).append(row)

        for rows in grouped.values():
            rows.sort(
                key=_selection_rank_value
            )

        reference_used = bool(
            summary.get(
                "reference_ligand_used_for_selection",
                False,
            )
        )

        lines += [
            f"- Selection summary: `{_relative_or_original(output_dir, summary_path)}`",
            (
                f"- Selection table: "
                f"`{_relative_or_original(output_dir, summary_csv)}`"
                if summary_csv.exists()
                else "- Selection table: unavailable"
            ),
            f"- Selection method: {_format_value(summary.get('selection_method'))}",
            (
                "- Reference ligand used for pocket selection: "
                + ("yes" if reference_used else "no")
            ),
            f"- Compounds evaluated: {_format_value(summary.get('compound_count'))}",
            f"- Total pocket attempts: {_format_value(summary.get('attempt_count'))}",
            "",
        ]

        if grouped:
            lines += [
                "### Selected Pocket per Compound",
                "",
                "| Compound | Selected pocket | fpocket rank | fpocket score | Top CNN score | CNN affinity | Score source | CNN-score margin | Selection confidence |",
                "|---|---|---:|---:|---:|---:|---|---:|---|",
            ]

            for compound in sorted(grouped):
                rows = grouped[compound]
                selected = _selected_attempt(
                    rows
                )

                if selected is None:
                    continue

                margin, confidence = (
                    _selection_confidence(rows)
                )

                lines.append(
                    "| "
                    f"{_format_value(compound)} | "
                    f"{_format_value(selected.get('pocket_id'))} | "
                    f"{_format_value(selected.get('pocket_rank'))} | "
                    f"{_format_value(selected.get('fpocket_score'))} | "
                    f"{_format_value(selected.get('top_cnn_score'))} | "
                    f"{_format_value(selected.get('top_cnn_affinity'))} | "
                    f"{_format_value(selected.get('score_source'))} | "
                    f"{_format_value(margin)} | "
                    f"{confidence} |"
                )

            lines += [
                "",
                "### All Pocket Attempts",
                "",
                "| Compound | Selection rank | Pocket | fpocket rank | fpocket score | Raw poses | Accepted | Rejected | Top CNN score | CNN affinity | Affinity | Score source |",
                "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
            ]

            for compound in sorted(grouped):
                for row in grouped[compound]:
                    lines.append(
                        "| "
                        f"{_format_value(compound)} | "
                        f"{_format_value(row.get('selection_rank'))} | "
                        f"{_format_value(row.get('pocket_id'))} | "
                        f"{_format_value(row.get('pocket_rank'))} | "
                        f"{_format_value(row.get('fpocket_score'))} | "
                        f"{_format_value(row.get('raw_poses'))} | "
                        f"{_format_value(row.get('accepted_poses'))} | "
                        f"{_format_value(row.get('rejected_poses'))} | "
                        f"{_format_value(row.get('top_cnn_score'))} | "
                        f"{_format_value(row.get('top_cnn_affinity'))} | "
                        f"{_format_value(row.get('top_minimized_affinity'))} | "
                        f"{_format_value(row.get('score_source'))} |"
                    )

            lines.append("")

            warnings: list[str] = []

            for compound in sorted(grouped):
                rows = grouped[compound]
                selected = _selected_attempt(
                    rows
                )

                if selected is None:
                    continue

                margin, confidence = (
                    _selection_confidence(rows)
                )

                if (
                    selected.get("score_source")
                    == "raw_pose_fallback"
                ):
                    warnings.append(
                        f"{compound}: selected pocket "
                        "was ranked using raw GNINA poses "
                        "because no PoseBusters-accepted "
                        "pose was available."
                    )

                if confidence == "low separation":
                    warnings.append(
                        f"{compound}: the selected pocket "
                        f"led the next alternative by only "
                        f"{_format_value(margin)} CNNscore; "
                        "the site assignment should be "
                        "treated as uncertain."
                    )

                fpocket_rank = _coerce_int(
                    selected.get("pocket_rank")
                )

                if (
                    fpocket_rank is not None
                    and fpocket_rank > 1
                ):
                    warnings.append(
                        f"{compound}: GNINA selected "
                        f"fpocket rank {fpocket_rank}, "
                        "overriding fpocket's geometry-only "
                        "rank 1 result."
                    )

            if warnings:
                lines += [
                    "### Pocket-Selection Notes",
                    "",
                ]

                for warning in warnings:
                    lines.append(
                        f"- {warning}"
                    )

                lines.append("")

            lines += [
                "_Selection confidence is a descriptive "
                "heuristic based on the difference between "
                "the selected pocket's top CNNscore and the "
                "next-best tested pocket. It is not a "
                "calibrated probability of binding-site "
                "correctness._",
                "",
            ]

    if definitions is not None:
        pockets_value = definitions.get(
            "pockets",
            [],
        )
        pockets = [
            pocket
            for pocket in pockets_value
            if isinstance(pocket, dict)
        ]

        lines += [
            "### Docking Pocket Definitions",
            "",
            f"- Pocket-definition file: `{_relative_or_original(output_dir, definitions_path)}`",
            f"- Pocket count: {_format_value(definitions.get('pocket_count'))}",
            f"- Ranking method: {_format_value(definitions.get('ranking_method'))}",
            "",
        ]

        if pockets:
            lines += [
                "| Pocket | fpocket rank | fpocket score | Center x, y, z | Size x, y, z | Source |",
                "|---|---:|---:|---|---|---|",
            ]

            pockets.sort(
                key=lambda pocket: (
                    _coerce_int(
                        pocket.get("pocket_rank")
                    )
                    or 9999
                )
            )

            for pocket in pockets:
                lines.append(
                    "| "
                    f"{_format_value(pocket.get('pocket_id'))} | "
                    f"{_format_value(pocket.get('pocket_rank'))} | "
                    f"{_format_value(pocket.get('fpocket_score'))} | "
                    f"{_format_box_triplet(pocket, 'center')} | "
                    f"{_format_box_triplet(pocket, 'size')} | "
                    f"{_format_value(pocket.get('source'))} |"
                )

            lines.append("")

    return lines



def _render_pose_recovery_section(
    output_dir: Path,
) -> list[str]:
    """Render an optional cognate-redocking benchmark section."""
    summary_path = (
        output_dir
        / "pose_set_recovery_summary.json"
    )
    metrics_path = (
        output_dir
        / "pose_set_recovery_metrics.csv"
    )
    standalone_report_path = (
        output_dir
        / "pose_set_recovery_report.md"
    )

    summary = _load_json(summary_path)

    if summary is None:
        return []

    top_pose_value = summary.get(
        "top_cnn_pose",
        {},
    )
    best_pose_value = summary.get(
        "best_sampled_pose",
        {},
    )

    top_pose = (
        top_pose_value
        if isinstance(top_pose_value, dict)
        else {}
    )
    best_pose = (
        best_pose_value
        if isinstance(best_pose_value, dict)
        else {}
    )

    def format_float(
        value: Any,
        *,
        decimals: int = 3,
    ) -> str:
        number = _coerce_float(value)

        if number is None:
            return "unknown"

        return f"{number:.{decimals}f}"

    def yes_no_unknown(value: Any) -> str:
        if value is True:
            return "yes"
        if value is False:
            return "no"
        return "unknown"

    def filename_or_unknown(value: Any) -> str:
        if value is None:
            return "unknown"

        text_value = str(value).strip()

        if not text_value:
            return "unknown"

        return Path(text_value).name

    def format_string_list(value: Any) -> str:
        if not isinstance(value, list):
            return "unknown"

        items = [
            str(item).strip()
            for item in value
            if str(item).strip()
        ]

        return ", ".join(items) if items else "unknown"

    threshold = _coerce_float(
        summary.get(
            "rmsd_threshold_angstrom"
        )
    )

    threshold_text = (
        f"{threshold:.3f} Å"
        if threshold is not None
        else "unknown"
    )

    mapping_failures = _coerce_int(
        summary.get(
            "mapping_failure_count"
        )
    )

    lines = [
        "## Cognate Pose-Recovery Benchmark",
        "",
        (
            f"- Summary: "
            f"`{_relative_or_original(output_dir, summary_path)}`"
        ),
        (
            f"- Pose metrics: "
            f"`{_relative_or_original(output_dir, metrics_path)}`"
            if metrics_path.exists()
            else "- Pose metrics: unavailable"
        ),
        (
            f"- Standalone report: "
            f"`{_relative_or_original(output_dir, standalone_report_path)}`"
            if standalone_report_path.exists()
            else "- Standalone report: unavailable"
        ),
        (
            "- Reference ligand: "
            f"`{filename_or_unknown(summary.get('reference_ligand'))}`"
        ),
        (
            "- GNINA pose set: "
            f"`{filename_or_unknown(summary.get('poses_sdf'))}`"
        ),
        (
            "- RMSD method: "
            f"{_format_value(summary.get('rmsd_method'))}"
        ),
        (
            "- Complete atom mapping policy: "
            f"{_format_value(summary.get('bond_order_mapping'))}"
        ),
        (
            "- Evaluated compound: "
            f"{_format_value(summary.get('evaluated_compound'))}"
        ),
        (
            "- Evaluated pocket: "
            f"{_format_value(summary.get('evaluated_pocket_id'))}"
        ),
        (
            "- Normally selected receptor conformer: "
            f"{_format_value(summary.get('normally_selected_receptor_conformer_id'))}"
        ),
        (
            "- Evaluated receptor conformers: "
            f"{format_string_list(summary.get('evaluated_receptor_conformer_ids'))}"
        ),
        (
            "- Evaluated receptor conformer count: "
            f"{_format_value(summary.get('evaluated_receptor_conformer_count'))}"
        ),
        (
            "- Evaluation stage: "
            f"{_format_value(summary.get('evaluation_stage'))}"
        ),
        (
            "- Reference ligand used for post hoc RMSD evaluation: "
            f"{yes_no_unknown(summary.get('reference_ligand_used_for_posthoc_evaluation'))}"
        ),
        (
            "- Same reference file supplied as the GNINA autobox ligand: "
            f"{yes_no_unknown(summary.get('reference_ligand_also_supplied_as_autobox_ligand'))}"
        ),
        (
            "- Reference ligand used to define the docking box: "
            f"{yes_no_unknown(summary.get('reference_ligand_used_for_box_definition'))}"
        ),
        (
            "- Reference ligand used to choose among detected pockets: "
            f"{yes_no_unknown(summary.get('reference_ligand_used_for_pocket_selection'))}"
        ),
        "",
        "| Metric | Value |",
        "|---|---:|",
        (
            "| Cognate RMSD threshold | "
            f"{threshold_text} |"
        ),
        (
            "| Chemically mapped poses | "
            f"{_format_value(summary.get('mapped_pose_count'))} |"
        ),
        (
            "| Mapping failures | "
            f"{_format_value(summary.get('mapping_failure_count'))} |"
        ),
        (
            "| Top CNN pose index | "
            f"{_format_value(top_pose.get('pose_index'))} |"
        ),
        (
            "| Top CNN receptor conformer | "
            f"{_format_value(top_pose.get('receptor_conformer_id'))} |"
        ),
        (
            "| Top CNN seed | "
            f"{_format_value(top_pose.get('seed'))} |"
        ),
        (
            "| Top CNN source pose | "
            f"{_format_value(top_pose.get('source_pose_number'))} |"
        ),
        (
            "| Top CNN pocket | "
            f"{_format_value(top_pose.get('pocket_id'))} |"
        ),
        (
            "| Top CNN score | "
            f"{format_float(top_pose.get('cnnscore'), decimals=6)} |"
        ),
        (
            "| Top CNN pose RMSD | "
            f"{format_float(top_pose.get('heavy_atom_rmsd'))} Å |"
        ),
        (
            "| Best sampled pose index | "
            f"{_format_value(best_pose.get('pose_index'))} |"
        ),
        (
            "| Best sampled receptor conformer | "
            f"{_format_value(best_pose.get('receptor_conformer_id'))} |"
        ),
        (
            "| Best sampled seed | "
            f"{_format_value(best_pose.get('seed'))} |"
        ),
        (
            "| Best sampled source pose | "
            f"{_format_value(best_pose.get('source_pose_number'))} |"
        ),
        (
            "| Best sampled pocket | "
            f"{_format_value(best_pose.get('pocket_id'))} |"
        ),
        (
            "| Best sampled pose RMSD | "
            f"{format_float(best_pose.get('heavy_atom_rmsd'))} Å |"
        ),
        (
            "| Sampling pass | "
            f"{yes_no_unknown(summary.get('sampling_pass'))} |"
        ),
        (
            "| Ranking pass | "
            f"{yes_no_unknown(summary.get('ranking_pass'))} |"
        ),
        "",
        "### Benchmark Interpretation",
        "",
        (
            f"- Overall result: "
            f"`{_format_value(summary.get('overall'))}`"
        ),
    ]

    if (
        mapping_failures is not None
        and mapping_failures > 0
    ):
        lines.append(
            "- Warning: one or more poses could not be "
            "chemically mapped to the reference ligand."
        )

    if summary.get("sampling_pass") is False:
        lines.append(
            "- Sampling failure: no evaluated pose met "
            "the configured cognate RMSD threshold."
        )

    if (
        summary.get("sampling_pass") is True
        and summary.get("ranking_pass") is False
    ):
        lines.append(
            "- Ranking failure: GNINA sampled a qualifying "
            "pose but did not rank one first."
        )

    if (
        summary.get("sampling_pass") is True
        and summary.get("ranking_pass") is True
    ):
        lines.append(
            "- The highest-CNN-scoring pose also met the "
            "configured cognate RMSD threshold."
        )

    lines += [
        "",
        "_This cognate-redocking benchmark evaluates pose "
        "recovery under a known reference-complex condition. "
        "It does not establish biological activity or clinical "
        "efficacy._",
        "",
    ]

    return lines

def _render_structure_validation_section(
    output_dir: Path,
) -> list[str]:
    candidates = (
        (
            "Submitted receptor",
            (
                output_dir
                / "structure_validation"
                / "receptor"
                / "ramachandran_validation.json"
            ),
        ),
        (
            "Automatic reference structure",
            (
                output_dir
                / "automatic_reference_evidence"
                / "structure_validation"
                / "reference"
                / "ramachandran_validation.json"
            ),
        ),
        (
            "Reference structure",
            (
                output_dir
                / "structure_validation"
                / "reference"
                / "ramachandran_validation.json"
            ),
        ),
    )

    reports: list[
        tuple[
            str,
            Path,
            dict[str, Any],
        ]
    ] = []

    observed_paths: set[Path] = set()

    for label, path in candidates:
        resolved = path.resolve()

        if resolved in observed_paths:
            continue

        observed_paths.add(
            resolved
        )

        report = _load_json(
            path
        )

        if report is not None:
            reports.append(
                (
                    label,
                    path,
                    report,
                )
            )

    if not reports:
        return []

    lines = [
        "## Structure Geometry Validation",
        "",
        (
            "Ramachandran validation is "
            "reported for structural review "
            "and does not currently reject "
            "or rerank structures."
        ),
        "",
        (
            "| Structure | Status | Chain | "
            "Evaluable | Favored | Allowed | "
            "Outliers | Screening flag |"
        ),
        (
            "|---|---|---|---:|---:|---:|"
            "---:|---|"
        ),
    ]

    artifact_lines: list[str] = []
    warnings: list[str] = []

    def fraction_text(
        value: Any,
    ) -> str:
        try:
            return f"{float(value):.2%}"
        except (
            TypeError,
            ValueError,
        ):
            return "unknown"

    for label, path, report in reports:
        summary = report.get(
            "summary",
            {},
        )

        if not isinstance(
            summary,
            dict,
        ):
            summary = {}

        status = str(
            report.get(
                "status",
                "unknown",
            )
        )

        flag = str(
            summary.get(
                "screening_flag",
                "unknown",
            )
        )

        favored_text = (
            f"{_format_value(summary.get('favored'))} "
            f"({fraction_text(summary.get('favored_fraction'))})"
        )

        allowed_text = (
            f"{_format_value(summary.get('allowed'))} "
            f"({fraction_text(summary.get('allowed_fraction'))})"
        )

        outlier_text = (
            f"{_format_value(summary.get('outliers'))} "
            f"({fraction_text(summary.get('outlier_fraction'))})"
        )

        lines.append(
            "| "
            f"{label} | "
            f"{status} | "
            f"{_format_value(report.get('requested_chain'))} | "
            f"{_format_value(report.get('evaluable_residues'))} | "
            f"{favored_text} | "
            f"{allowed_text} | "
            f"{outlier_text} | "
            f"{flag} |"
        )

        artifact_lines.append(
            f"- {label} validation file: "
            f"`{_relative_or_original(output_dir, path)}`"
        )

        csv_path = (
            path.parent
            / "ramachandran_residues.csv"
        )

        if csv_path.is_file():
            artifact_lines.append(
                f"- {label} residue table: "
                f"`{_relative_or_original(output_dir, csv_path)}`"
            )

        if status == "failed":
            error = report.get(
                "error",
                {},
            )

            message = (
                error.get("message")
                if isinstance(
                    error,
                    dict,
                )
                else error
            )

            warnings.append(
                f"{label}: Ramachandran "
                "validation failed: "
                f"{_format_value(message)}"
            )

        elif flag not in {
            "meets_ramalyze_goals",
            "unknown",
        }:
            warnings.append(
                f"{label}: structure geometry "
                f"was flagged as `{flag}`. "
                "Review residue-level outliers "
                "before relying on this model."
            )

    lines.append("")

    if artifact_lines:
        lines += [
            "### Validation Artifacts",
            "",
            *artifact_lines,
            "",
        ]

    if warnings:
        lines += [
            "### Geometry Review Notes",
            "",
        ]

        for warning in warnings:
            lines.append(
                f"- {warning}"
            )

        lines.append("")

    lines += [
        (
            "_Ramachandran results describe "
            "backbone conformational geometry "
            "only. They do not replace clash, "
            "rotamer, bond-length, bond-angle, "
            "electron-density, or dynamics "
            "validation._"
        ),
        "",
    ]

    return lines


def render_run_report(
    *,
    output_dir: Path,
    target_evidence: dict[str, Any] | None,
    hypotheses: list[dict[str, Any]],
) -> str:
    lines = [
        "# CompoundRank Run Report",
        "",
        f"- Output directory: `{output_dir}`",
        "",
    ]

    lines.extend(_render_target_section(target_evidence))
    lines.extend(
        _render_structure_validation_section(
            output_dir
        )
    )
    lines.extend(_render_ligand_retrieval_section(output_dir))
    lines.extend(_render_docking_section(hypotheses))
    lines.extend(_render_pocket_selection_section(output_dir))
    lines.extend(_render_docking_attempt_summary_section(output_dir))

    lines += _render_pose_recovery_section(
        output_dir
    )

    lines += [
        "## Interpretation Limits",
        "",
        "- This report summarizes computational target annotation and docking hypotheses.",
        "- GNINA/PoseBusters outputs are not proof of binding, inhibition, safety, or antiviral activity.",
        "- Target evidence is based on computational annotation and should be reviewed against literature.",
        "- Experimental validation remains required before biological or translational claims.",
        "",
    ]

    return "\n".join(lines)


def write_run_report(
    *,
    output_dir: Path,
    target_evidence_path: Path | None = None,
    report_name: str = "compoundrank_run_report.md",
) -> Path:
    output_dir = Path(output_dir)

    if target_evidence_path is None:
        target_evidence_path = output_dir / "target_evidence.json"

    target_evidence = _load_json(Path(target_evidence_path))
    hypotheses = collect_pdb_hypotheses(output_dir)

    report_text = render_run_report(
        output_dir=output_dir,
        target_evidence=target_evidence,
        hypotheses=hypotheses,
    )

    report_path = output_dir / report_name
    report_path.write_text(
        report_text,
        encoding="utf-8",
    )

    return report_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a combined CompoundRank run report."
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target-evidence", type=Path, default=None)
    parser.add_argument("--report-name", default="compoundrank_run_report.md")
    args = parser.parse_args()

    report_path = write_run_report(
        output_dir=args.output_dir,
        target_evidence_path=args.target_evidence,
        report_name=args.report_name,
    )

    print(f"Run report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
