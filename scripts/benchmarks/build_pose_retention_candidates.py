#!/usr/bin/env python3
"""
Build alternative-pose retention candidates from pose-recovery benchmark metrics.

Purpose:
- Keep GNINA top-CNN pose as the primary hypothesis.
- Retain alternative poses when they have stronger benchmark geometry,
  near-threshold RMSD, or are the best sampled pose in a receptor snapshot.

This is a benchmark/report tool. It does not change production ranking yet.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null", "unknown", "unavailable"}:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def pick_column(fieldnames, exact=(), contains=()):
    fields = list(fieldnames or [])
    lower = {field.lower(): field for field in fields}

    for candidate in exact:
        if candidate.lower() in lower:
            return lower[candidate.lower()]

    for field in fields:
        fl = field.lower()
        if all(piece.lower() in fl for piece in contains):
            return field

    return None


def read_metrics(metrics_path: Path):
    with metrics_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    rmsd_col = pick_column(
        fieldnames,
        exact=("heavy_atom_rmsd", "rmsd", "rmsd_angstrom", "pose_rmsd"),
        contains=("rmsd",),
    )

    cnn_col = pick_column(
        fieldnames,
        exact=("cnnscore", "cnn_score", "CNNscore", "cnn_pose_score"),
        contains=("cnn",),
    )

    if cnn_col and "affinity" in cnn_col.lower():
        non_affinity = [
            field for field in fieldnames
            if "cnn" in field.lower() and "affinity" not in field.lower()
        ]
        if non_affinity:
            cnn_col = non_affinity[0]

    pose_col = pick_column(
        fieldnames,
        exact=("pose_index", "pose", "pose_number", "mode", "mode_index"),
    )

    seed_col = pick_column(
        fieldnames,
        exact=("seed", "random_seed"),
    )

    if rmsd_col is None:
        raise ValueError(f"No RMSD column found in {metrics_path}. Columns: {fieldnames}")

    parsed = []
    for file_row_index, row in enumerate(rows, start=1):
        rmsd = parse_float(row.get(rmsd_col))
        if rmsd is None:
            continue

        cnnscore = parse_float(row.get(cnn_col)) if cnn_col else None

        parsed.append(
            {
                "file_row_index": file_row_index,
                "raw": row,
                "rmsd": rmsd,
                "cnnscore": cnnscore,
                "pose_index": row.get(pose_col, "") if pose_col else "",
                "seed": row.get(seed_col, "") if seed_col else "",
            }
        )

    if not parsed:
        raise ValueError(f"No numeric RMSD rows found in {metrics_path}")

    # GNINA CNN pose score is higher-is-better.
    ranked = sorted(
        parsed,
        key=lambda item: (
            item["cnnscore"] is None,
            -(item["cnnscore"] if item["cnnscore"] is not None else -1e9),
            item["file_row_index"],
        ),
    )

    for rank, item in enumerate(ranked, start=1):
        item["cnn_rank"] = rank

    return parsed, ranked


def add_candidate(
    candidates,
    *,
    snapshot_run,
    metrics_path,
    candidate_type,
    item,
    top_item,
    best_item,
    threshold,
    near_threshold,
    reason,
):
    cnn_gap = ""
    if top_item.get("cnnscore") is not None and item.get("cnnscore") is not None:
        cnn_gap = top_item["cnnscore"] - item["cnnscore"]

    rmsd_delta_from_top = item["rmsd"] - top_item["rmsd"]

    tags = []
    if item is top_item:
        tags.append("primary_top_cnn")
    if item is best_item:
        tags.append("best_geometry")
    if item["rmsd"] <= threshold:
        tags.append("sub_threshold_rmsd")
    elif item["rmsd"] <= near_threshold:
        tags.append("near_threshold_rmsd")
    if item.get("cnn_rank", 9999) <= 5:
        tags.append("top5_cnn")
    elif item.get("cnn_rank", 9999) <= 10:
        tags.append("top10_cnn")
    elif item.get("cnn_rank", 9999) <= 20:
        tags.append("top20_cnn")

    retain = (
        item is top_item
        or item is best_item
        or item["rmsd"] <= threshold
        or item["rmsd"] <= near_threshold
    )

    if item is best_item and item is not top_item and item["rmsd"] <= threshold:
        interpretation = (
            "Benchmark-supported alternative: this pose recovered sub-threshold "
            "reference geometry but was not selected as top CNN."
        )
    elif item is best_item and item is not top_item:
        interpretation = (
            "Best sampled geometry for this snapshot, but above strict threshold. "
            "Retain for diagnostic review, not as a benchmark pass."
        )
    elif item is top_item and item["rmsd"] > threshold:
        interpretation = (
            "Primary GNINA/CNN hypothesis, but benchmark geometry does not pass "
            "the cognate RMSD threshold."
        )
    elif item is top_item:
        interpretation = "Primary GNINA/CNN hypothesis and benchmark geometry passes."
    else:
        interpretation = "Retained as a supporting alternative."

    candidates.append(
        {
            "snapshot_run": snapshot_run,
            "candidate_type": candidate_type,
            "retain_recommended": retain,
            "pose_index": item.get("pose_index", ""),
            "seed": item.get("seed", ""),
            "cnn_rank": item.get("cnn_rank", ""),
            "cnnscore": item.get("cnnscore", ""),
            "heavy_atom_rmsd": item.get("rmsd", ""),
            "rmsd_threshold": threshold,
            "near_threshold": near_threshold,
            "cnn_score_gap_from_top": cnn_gap,
            "rmsd_delta_from_top_cnn_pose": rmsd_delta_from_top,
            "evidence_tags": ";".join(tags),
            "reason": reason,
            "interpretation": interpretation,
            "metrics_path": str(metrics_path),
        }
    )


def unique_candidates(candidates):
    seen = set()
    unique = []
    for row in candidates:
        key = (
            row["snapshot_run"],
            row["pose_index"],
            row["seed"],
            row["candidate_type"],
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(row)
    return unique


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runroot", required=True, type=Path)
    parser.add_argument("--threshold", type=float, default=2.0)
    parser.add_argument("--near-threshold", type=float, default=3.0)
    parser.add_argument("--retain-top-n", type=int, default=20)
    args = parser.parse_args()

    runroot = args.runroot
    metrics_files = sorted(runroot.glob("g39_md_snapshots/*/pose_set_recovery_metrics.csv"))

    if not metrics_files:
        raise SystemExit(f"No pose_set_recovery_metrics.csv files found under {runroot}")

    candidates = []
    snapshot_summaries = []

    for metrics_path in metrics_files:
        snapshot_run = metrics_path.parent.name
        parsed, ranked = read_metrics(metrics_path)

        top_item = ranked[0]
        best_item = min(parsed, key=lambda item: item["rmsd"])

        add_candidate(
            candidates,
            snapshot_run=snapshot_run,
            metrics_path=metrics_path,
            candidate_type="primary_top_cnn",
            item=top_item,
            top_item=top_item,
            best_item=best_item,
            threshold=args.threshold,
            near_threshold=args.near_threshold,
            reason="Highest GNINA CNN pose score.",
        )

        if best_item is not top_item:
            add_candidate(
                candidates,
                snapshot_run=snapshot_run,
                metrics_path=metrics_path,
                candidate_type="best_geometry",
                item=best_item,
                top_item=top_item,
                best_item=best_item,
                threshold=args.threshold,
                near_threshold=args.near_threshold,
                reason="Lowest heavy-atom RMSD to reference ligand in benchmark mode.",
            )

        # Retain every sub-threshold pose, even if it is not the single best.
        for item in parsed:
            if item["rmsd"] <= args.threshold:
                add_candidate(
                    candidates,
                    snapshot_run=snapshot_run,
                    metrics_path=metrics_path,
                    candidate_type="sub_threshold_pose",
                    item=item,
                    top_item=top_item,
                    best_item=best_item,
                    threshold=args.threshold,
                    near_threshold=args.near_threshold,
                    reason="Pose passes the configured cognate RMSD threshold.",
                )

        # Retain near-threshold top-N diagnostic alternatives.
        for item in ranked[: args.retain_top_n]:
            if item["rmsd"] <= args.near_threshold and item["rmsd"] > args.threshold:
                add_candidate(
                    candidates,
                    snapshot_run=snapshot_run,
                    metrics_path=metrics_path,
                    candidate_type="near_threshold_topn",
                    item=item,
                    top_item=top_item,
                    best_item=best_item,
                    threshold=args.threshold,
                    near_threshold=args.near_threshold,
                    reason=(
                        f"Pose is within top {args.retain_top_n} by CNN and "
                        f"within {args.near_threshold:.2f} Å RMSD."
                    ),
                )

        snapshot_summaries.append(
            {
                "snapshot_run": snapshot_run,
                "top_cnn_rmsd": top_item["rmsd"],
                "top_cnn_score": top_item.get("cnnscore", ""),
                "best_rmsd": best_item["rmsd"],
                "best_rmsd_cnn_rank": best_item.get("cnn_rank", ""),
                "best_rmsd_cnn_score": best_item.get("cnnscore", ""),
                "sampling_pass": best_item["rmsd"] <= args.threshold,
                "ranking_pass": top_item["rmsd"] <= args.threshold,
                "cnn_score_gap_top_minus_best": (
                    top_item["cnnscore"] - best_item["cnnscore"]
                    if top_item.get("cnnscore") is not None and best_item.get("cnnscore") is not None
                    else ""
                ),
            }
        )

    candidates = unique_candidates(candidates)

    out_csv = runroot / "alternative_pose_retention_candidates.csv"
    out_md = runroot / "alternative_pose_retention_report.md"

    fieldnames = [
        "snapshot_run",
        "candidate_type",
        "retain_recommended",
        "pose_index",
        "seed",
        "cnn_rank",
        "cnnscore",
        "heavy_atom_rmsd",
        "rmsd_threshold",
        "near_threshold",
        "cnn_score_gap_from_top",
        "rmsd_delta_from_top_cnn_pose",
        "evidence_tags",
        "reason",
        "interpretation",
        "metrics_path",
    ]

    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(candidates)

    sampling_passes = [row for row in snapshot_summaries if row["sampling_pass"]]
    ranking_passes = [row for row in snapshot_summaries if row["ranking_pass"]]
    retained_subthreshold = [
        row for row in candidates
        if "sub_threshold_rmsd" in row["evidence_tags"]
    ]

    best_overall = min(snapshot_summaries, key=lambda row: row["best_rmsd"])

    md = [
        "# Alternative Pose Retention Report",
        "",
        f"- Run root: `{runroot}`",
        f"- RMSD threshold: {args.threshold:.3f} Å",
        f"- Near-threshold diagnostic cutoff: {args.near_threshold:.3f} Å",
        f"- Retained top-N search window: {args.retain_top_n}",
        "",
        "## Benchmark Summary",
        "",
        f"- Snapshot analyses: {len(snapshot_summaries)}",
        f"- Sampling passes: {len(sampling_passes)} / {len(snapshot_summaries)}",
        f"- Ranking passes: {len(ranking_passes)} / {len(snapshot_summaries)}",
        f"- Retained sub-threshold candidate poses: {len(retained_subthreshold)}",
        "",
        "## Best Overall Geometry",
        "",
        f"- Snapshot: `{best_overall['snapshot_run']}`",
        f"- Best RMSD: {best_overall['best_rmsd']:.3f} Å",
        f"- CNN rank of best-RMSD pose: {best_overall['best_rmsd_cnn_rank']}",
        f"- Top-CNN RMSD in same snapshot: {best_overall['top_cnn_rmsd']:.3f} Å",
        f"- CNN score gap, top minus best-RMSD pose: {best_overall['cnn_score_gap_top_minus_best']}",
        "",
        "## Retention Policy Suggested by This Run",
        "",
        "Use the top GNINA CNN pose as the primary hypothesis, but retain alternative poses when they show stronger benchmark geometry or other structural evidence.",
        "",
        "This run specifically supports an evidence-aware retention policy because at least one sub-threshold pose was sampled but not ranked first by CNN.",
        "",
        "## Retained Candidate Poses",
        "",
        "| Snapshot | Type | Pose | Seed | CNN rank | CNN score | RMSD Å | Evidence tags | Interpretation |",
        "|---|---|---:|---:|---:|---:|---:|---|---|",
    ]

    for row in sorted(
        candidates,
        key=lambda r: (
            r["snapshot_run"],
            float(r["heavy_atom_rmsd"]) if str(r["heavy_atom_rmsd"]) else 999.0,
            int(r["cnn_rank"]) if str(r["cnn_rank"]).isdigit() else 999,
        ),
    ):
        md.append(
            "| {snapshot_run} | {candidate_type} | {pose_index} | {seed} | {cnn_rank} | {cnnscore} | {heavy_atom_rmsd} | {evidence_tags} | {interpretation} |".format(
                **row
            )
        )

    md.extend(
        [
            "",
            "## Production Translation",
            "",
            "In production runs without a reference ligand, EXORCIST should not use RMSD. Instead, analogous evidence channels should include:",
            "",
            "- active-site/contact recovery",
            "- known catalytic residue overlap",
            "- pose clustering stability",
            "- recurrence across MD/receptor snapshots",
            "- physical sanity filters",
            "- CNN score as one signal, not the only signal",
            "",
            f"CSV output: `{out_csv.name}`",
            "",
        ]
    )

    out_md.write_text("\n".join(md), encoding="utf-8")

    print("Wrote:", out_csv)
    print("Wrote:", out_md)
    print()
    print("=== SUMMARY ===")
    print(f"Snapshots: {len(snapshot_summaries)}")
    print(f"Sampling passes: {len(sampling_passes)} / {len(snapshot_summaries)}")
    print(f"Ranking passes: {len(ranking_passes)} / {len(snapshot_summaries)}")
    print(f"Retained candidates: {len(candidates)}")
    print(
        "Best overall:",
        f"{best_overall['best_rmsd']:.3f} Å",
        "in",
        best_overall["snapshot_run"],
        "| CNN rank",
        best_overall["best_rmsd_cnn_rank"],
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
