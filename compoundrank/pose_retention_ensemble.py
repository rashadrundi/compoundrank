"""Tier-aware ensemble summary for production-style pose retention.

Reads production_pose_retention_candidates.csv files across retrieved receptor
snapshot output directories and summarizes compound-level support.

This does not use RMSD. It summarizes:
- primary_supported poses
- strong_alternative poses
- diagnostic_alternative poses
- physically_warned poses
- alternatives that exceed primary evidence
- recurrent contacts across receptor snapshots
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


def parse_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null", "unknown", "unavailable"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def parse_bool(value) -> bool:
    return str(value).strip().lower() in {"true", "yes", "1", "retain", "retained"}


def split_semicolon(value) -> set[str]:
    if not value:
        return set()
    return {part.strip() for part in str(value).split(";") if part.strip()}


def mean(values):
    values = [value for value in values if value is not None]
    if not values:
        return ""
    return sum(values) / len(values)


def infer_confidence_tier(row: dict) -> str:
    """Fallback for older CSVs that lack confidence_tier."""
    existing = str(row.get("confidence_tier") or "").strip()
    if existing:
        return existing

    retained = parse_bool(row.get("retain_recommended"))
    if not retained:
        return "not_retained"

    reasons = split_semicolon(row.get("retain_reasons"))
    tags = split_semicolon(row.get("evidence_tags"))

    if {"possible_clash", "extreme_minimized_affinity", "weak_contact_distance"} & tags:
        return "physically_warned"

    if row.get("candidate_type") == "primary_top_cnn":
        return "primary_supported"

    if "evidence_score_exceeds_primary" in reasons or "known_contact_supported" in reasons:
        return "strong_alternative"

    if "ensemble_recurrent_contacts" in reasons or "near_top_cnn" in reasons:
        return "diagnostic_alternative"

    return "diagnostic_alternative"


def find_candidate_csvs(runroot: Path) -> list[Path]:
    runroot = Path(runroot)

    direct = runroot / "production_pose_retention_candidates.csv"
    if direct.exists():
        return [direct]

    snapshots_root = runroot / "retrieved_on_top_g39_snapshots"
    if snapshots_root.exists():
        return sorted(snapshots_root.glob("*/production_pose_retention_candidates.csv"))

    return sorted(runroot.glob("*/production_pose_retention_candidates.csv"))


def read_candidate_rows(candidate_csvs: list[Path]) -> list[dict]:
    rows = []

    for csv_path in candidate_csvs:
        snapshot_run = csv_path.parent.name

        with csv_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)

            for row in reader:
                row = dict(row)
                row["snapshot_run"] = snapshot_run
                row["candidate_csv"] = str(csv_path)

                row["_retain"] = parse_bool(row.get("retain_recommended"))
                row["_evidence_score"] = parse_float(row.get("evidence_score"))
                row["_primary_evidence_score"] = parse_float(row.get("primary_evidence_score"))
                row["_cnn_score"] = parse_float(row.get("cnn_score"))
                row["_cnn_rank"] = parse_float(row.get("cnn_rank_within_compound"))
                row["_known_overlap_count"] = parse_float(row.get("known_contact_overlap_count")) or 0
                row["_recurrent_overlap_count"] = parse_float(row.get("recurrent_contact_overlap_count")) or 0

                row["_retain_reasons"] = split_semicolon(row.get("retain_reasons"))
                row["_evidence_tags"] = split_semicolon(row.get("evidence_tags"))
                row["_all_contacts"] = split_semicolon(row.get("all_contacts"))
                row["_recurrent_contacts"] = split_semicolon(row.get("recurrent_contact_overlap"))
                row["_known_contacts"] = split_semicolon(row.get("known_contact_overlap"))

                tier = infer_confidence_tier(row)
                row["_confidence_tier"] = tier
                row["confidence_tier"] = tier

                rows.append(row)

    return rows


def summarize_compound(compound: str, rows: list[dict]) -> dict:
    snapshots = sorted({row["snapshot_run"] for row in rows})

    retained = [row for row in rows if row["_retain"]]
    primary = [row for row in rows if row.get("candidate_type") == "primary_top_cnn"]
    retained_alt = [
        row for row in retained if row.get("candidate_type") != "primary_top_cnn"
    ]

    tier_counter = Counter(row["_confidence_tier"] for row in retained)

    evidence_exceeds_primary = [
        row for row in retained_alt
        if "evidence_score_exceeds_primary" in row["_retain_reasons"]
    ]

    near_top = [
        row for row in retained_alt
        if "near_top_cnn" in row["_retain_reasons"]
    ]

    recurrent_supported = [
        row for row in retained
        if "ensemble_recurrent_contacts" in row["_retain_reasons"]
        or row["_recurrent_overlap_count"] > 0
    ]

    known_supported = [
        row for row in retained
        if "known_contact_supported" in row["_retain_reasons"]
        or row["_known_overlap_count"] > 0
    ]

    contact_snapshot_counter: dict[str, set[str]] = defaultdict(set)
    for row in retained:
        for residue in row["_all_contacts"]:
            contact_snapshot_counter[residue].add(row["snapshot_run"])

    contact_counts = {
        residue: len(snapshot_set)
        for residue, snapshot_set in contact_snapshot_counter.items()
    }

    top_contacts = sorted(
        contact_counts.items(),
        key=lambda item: (-item[1], item[0]),
    )[:12]

    top_contacts_text = ";".join(
        f"{residue}:{count}" for residue, count in top_contacts
    )

    best = None
    for row in retained:
        if row["_evidence_score"] is None:
            continue
        if best is None or row["_evidence_score"] > best["_evidence_score"]:
            best = row

    max_evidence_score = best["_evidence_score"] if best else ""

    primary_supported_count = tier_counter.get("primary_supported", 0)
    strong_alternative_count = tier_counter.get("strong_alternative", 0)
    diagnostic_alternative_count = tier_counter.get("diagnostic_alternative", 0)
    physically_warned_count = tier_counter.get("physically_warned", 0)
    not_retained_count = tier_counter.get("not_retained", 0)
    unknown_tier_count = tier_counter.get("unknown", 0)

    clean_support_count = primary_supported_count + strong_alternative_count
    warning_ratio = physically_warned_count / len(retained) if retained else 0.0

    if clean_support_count and warning_ratio < 0.34:
        tier_assessment = "cleaner_support"
        interpretation = (
            "Compound has cleaner support: primary/strong retained poses outweigh "
            "physical warnings."
        )
    elif clean_support_count and warning_ratio < 0.67:
        tier_assessment = "mixed_support"
        interpretation = (
            "Compound has mixed support: useful retained poses exist, but physical "
            "warnings require caution."
        )
    elif physically_warned_count:
        tier_assessment = "warning_dominated"
        interpretation = (
            "Compound is warning-dominated: recurrence exists, but most retained "
            "poses carry physical-sanity warnings."
        )
    elif diagnostic_alternative_count:
        tier_assessment = "diagnostic_only"
        interpretation = (
            "Compound is mostly diagnostic: retained alternatives need manual review "
            "before confidence increases."
        )
    else:
        tier_assessment = "low_support"
        interpretation = "Compound has low retained-pose support under the current policy."

    return {
        "compound": compound,
        "snapshots_with_compound": len(snapshots),
        "snapshots": ";".join(snapshots),
        "total_rows": len(rows),
        "retained_rows": len(retained),
        "primary_rows": len(primary),
        "retained_alternative_rows": len(retained_alt),
        "primary_supported_count": primary_supported_count,
        "strong_alternative_count": strong_alternative_count,
        "diagnostic_alternative_count": diagnostic_alternative_count,
        "physically_warned_count": physically_warned_count,
        "not_retained_count": not_retained_count,
        "unknown_tier_count": unknown_tier_count,
        "evidence_exceeds_primary_count": len(evidence_exceeds_primary),
        "near_top_cnn_alternative_count": len(near_top),
        "recurrent_supported_retained_count": len(recurrent_supported),
        "known_supported_retained_count": len(known_supported),
        "avg_primary_evidence_score": mean(row["_evidence_score"] for row in primary),
        "avg_retained_alternative_evidence_score": mean(row["_evidence_score"] for row in retained_alt),
        "max_evidence_score": max_evidence_score,
        "best_snapshot": best["snapshot_run"] if best else "",
        "best_hypothesis_file": best.get("hypothesis_file", "") if best else "",
        "best_candidate_type": best.get("candidate_type", "") if best else "",
        "best_confidence_tier": best.get("_confidence_tier", "") if best else "",
        "best_cnn_score": best.get("cnn_score", "") if best else "",
        "top_cross_snapshot_contacts": top_contacts_text,
        "tier_assessment": tier_assessment,
        "ensemble_interpretation": interpretation,
    }


def tier_priority(row: dict) -> tuple:
    priority = {
        "strong_alternative": 0,
        "diagnostic_alternative": 1,
        "physically_warned": 2,
        "primary_supported": 3,
        "not_retained": 4,
        "unknown": 5,
    }
    tier = row.get("_confidence_tier", "unknown")
    evidence_score = row.get("_evidence_score")
    score_sort = -(evidence_score if evidence_score is not None else -1e9)
    return (priority.get(tier, 5), score_sort, row.get("compound", ""), row.get("snapshot_run", ""))


def write_ensemble_artifacts(runroot: Path) -> tuple[Path, Path]:
    runroot = Path(runroot)
    candidate_csvs = find_candidate_csvs(runroot)

    if not candidate_csvs:
        raise FileNotFoundError(
            f"No production_pose_retention_candidates.csv files found under {runroot}"
        )

    rows = read_candidate_rows(candidate_csvs)

    by_compound: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_compound[row.get("compound", "unknown")].append(row)

    summaries = [
        summarize_compound(compound, compound_rows)
        for compound, compound_rows in sorted(by_compound.items())
    ]

    summaries = sorted(
        summaries,
        key=lambda row: (
            -row["snapshots_with_compound"],
            -row["primary_supported_count"],
            -row["strong_alternative_count"],
            row["physically_warned_count"],
            row["compound"],
        ),
    )

    out_csv = runroot / "production_pose_retention_ensemble_summary.csv"
    out_md = runroot / "production_pose_retention_ensemble_report.md"

    fieldnames = [
        "compound",
        "snapshots_with_compound",
        "snapshots",
        "total_rows",
        "retained_rows",
        "primary_rows",
        "retained_alternative_rows",
        "primary_supported_count",
        "strong_alternative_count",
        "diagnostic_alternative_count",
        "physically_warned_count",
        "not_retained_count",
        "unknown_tier_count",
        "evidence_exceeds_primary_count",
        "near_top_cnn_alternative_count",
        "recurrent_supported_retained_count",
        "known_supported_retained_count",
        "avg_primary_evidence_score",
        "avg_retained_alternative_evidence_score",
        "max_evidence_score",
        "best_snapshot",
        "best_hypothesis_file",
        "best_candidate_type",
        "best_confidence_tier",
        "best_cnn_score",
        "top_cross_snapshot_contacts",
        "tier_assessment",
        "ensemble_interpretation",
    ]

    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)

    retained_alts = [
        row for row in rows
        if row["_retain"] and row.get("candidate_type") != "primary_top_cnn"
    ]

    retained_alts = sorted(retained_alts, key=tier_priority)

    md = [
        "# Production Pose Retention Ensemble Report",
        "",
        f"- Run root: `{runroot}`",
        f"- Candidate CSV files summarized: {len(candidate_csvs)}",
        f"- Candidate rows summarized: {len(rows)}",
        f"- Compounds summarized: {len(summaries)}",
        "",
        "## Purpose",
        "",
        "This report summarizes production-style pose-retention evidence across receptor snapshots. It does not use RMSD. It asks which compounds repeatedly produce retained primary or alternative hypotheses across conformational variants, then separates clean support from physically warned support.",
        "",
        "## Compound-Level Confidence Tier Summary",
        "",
        "| Compound | Snapshots | Primary supported | Strong alternatives | Diagnostic alternatives | Physically warned | Alternatives > primary | Best evidence | Assessment | Interpretation |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]

    for row in summaries:
        best_score = row["max_evidence_score"]
        best_score_text = f"{best_score:.3f}" if isinstance(best_score, float) else str(best_score)

        md.append(
            "| {compound} | {snapshots_with_compound} | {primary_supported_count} | {strong_alternative_count} | {diagnostic_alternative_count} | {physically_warned_count} | {evidence_exceeds_primary_count} | ".format(
                **row
            )
            + f"{best_score_text} | "
            + f"{row['tier_assessment']} | "
            + f"{row['ensemble_interpretation']} |"
        )

    md.extend(
        [
            "",
            "## Top Cross-Snapshot Contacts",
            "",
            "| Compound | Top recurring contacts |",
            "|---|---|",
        ]
    )

    for row in summaries:
        md.append(f"| {row['compound']} | {row['top_cross_snapshot_contacts']} |")

    md.extend(
        [
            "",
            "## Top Retained Alternatives by Confidence Tier",
            "",
            "| Compound | Snapshot | Tier | File | CNN rank | CNN score | Evidence score | Reasons | Interpretation |",
            "|---|---|---|---|---:|---:|---:|---|---|",
        ]
    )

    for row in retained_alts[:25]:
        evidence_score = row["_evidence_score"]
        evidence_text = f"{evidence_score:.3f}" if evidence_score is not None else ""

        md.append(
            "| {compound} | {snapshot_run} | {confidence_tier} | `{hypothesis_file}` | {cnn_rank_within_compound} | {cnn_score} | ".format(
                **row
            )
            + f"{evidence_text} | "
            + f"{row.get('retain_reasons', '')} | "
            + f"{row.get('interpretation', '')} |"
        )

    md.extend(
        [
            "",
            "## Interpretation Notes",
            "",
            "- This is a production-style summary, so it does not use ligand-reference RMSD.",
            "- `primary_supported` and `strong_alternative` are cleaner support tiers.",
            "- `diagnostic_alternative` means the pose is worth reviewing but should not be elevated yet.",
            "- `physically_warned` means recurrence or CNN support exists, but clash/extreme-affinity warnings reduce confidence.",
            "- Alternatives that exceed the primary evidence score deserve manual review because they may represent cases where top-CNN ranking is not the best structural hypothesis.",
            "",
            f"CSV output: `{out_csv.name}`",
            "",
        ]
    )

    out_md.write_text("\n".join(md), encoding="utf-8")

    return out_csv, out_md


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runroot", required=True, type=Path)
    args = parser.parse_args(argv)

    out_csv, out_md = write_ensemble_artifacts(args.runroot)

    print("Wrote:", out_csv)
    print("Wrote:", out_md)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
