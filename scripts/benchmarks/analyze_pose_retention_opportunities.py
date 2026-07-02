#!/usr/bin/env python3
"""
Analyze whether alternative-pose retention would rescue benchmark-quality poses.

This script reads a completed 2HU0 MD-ensemble benchmark run and asks:

1. What was the top-CNN pose RMSD?
2. What was the best-RMSD pose?
3. Where did the best-RMSD pose rank by CNN?
4. Would top-3, top-5, top-10, or top-20 retention have captured a sub-2 Å pose?
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


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


def choose_column(fieldnames: Iterable[str], candidates: list[str], contains_all: list[str] | None = None):
    fields = list(fieldnames or [])
    lower_map = {field.lower(): field for field in fields}

    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]

    if contains_all:
        needles = [item.lower() for item in contains_all]
        for field in fields:
            fl = field.lower()
            if all(needle in fl for needle in needles):
                return field

    return None


def truthy(value):
    return str(value).strip().lower() in {"true", "yes", "1", "pass", "passed"}


def summarize_metrics(metrics_path: Path, threshold: float, top_ns: list[int]) -> dict:
    with metrics_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    if not rows:
        return {
            "metrics_path": str(metrics_path),
            "status": "empty_metrics",
        }

    rmsd_col = choose_column(
        fieldnames,
        [
            "heavy_atom_rmsd",
            "rmsd",
            "rmsd_angstrom",
            "heavy_atom_rmsd_angstrom",
            "pose_rmsd",
        ],
        contains_all=["rmsd"],
    )

    cnn_col = choose_column(
        fieldnames,
        [
            "cnn_score",
            "cnnscore",
            "CNNscore",
            "cnn_pose_score",
            "pose_cnn_score",
            "CNN",
        ],
        contains_all=["cnn"],
    )

    # Avoid choosing CNN affinity if a pose-score column exists.
    if cnn_col and "affinity" in cnn_col.lower():
        non_affinity = [
            f for f in fieldnames
            if "cnn" in f.lower() and "affinity" not in f.lower()
        ]
        if non_affinity:
            cnn_col = non_affinity[0]

    rank_col = choose_column(
        fieldnames,
        [
            "cnn_rank",
            "rank",
            "selection_rank",
            "pose_rank",
        ],
    )

    pose_col = choose_column(
        fieldnames,
        [
            "pose",
            "pose_index",
            "pose_number",
            "mode",
            "mode_index",
        ],
    )

    seed_col = choose_column(
        fieldnames,
        [
            "seed",
            "random_seed",
        ],
    )

    if rmsd_col is None:
        return {
            "metrics_path": str(metrics_path),
            "status": "missing_rmsd_column",
            "available_columns": ",".join(fieldnames),
        }

    parsed = []
    for original_index, row in enumerate(rows, start=1):
        rmsd = parse_float(row.get(rmsd_col))
        if rmsd is None:
            continue

        cnn = parse_float(row.get(cnn_col)) if cnn_col else None
        explicit_rank = parse_float(row.get(rank_col)) if rank_col else None

        parsed.append(
            {
                "original_index": original_index,
                "row": row,
                "rmsd": rmsd,
                "cnn": cnn,
                "explicit_rank": explicit_rank,
                "pose": row.get(pose_col, "") if pose_col else "",
                "seed": row.get(seed_col, "") if seed_col else "",
            }
        )

    if not parsed:
        return {
            "metrics_path": str(metrics_path),
            "status": "no_numeric_rmsd_rows",
            "rmsd_column": rmsd_col,
            "available_columns": ",".join(fieldnames),
        }

    # CNN score is higher-is-better for GNINA CNN pose score.
    # If no CNN column exists, preserve file order as fallback.
    if cnn_col:
        cnn_sorted = sorted(
            parsed,
            key=lambda item: (
                item["cnn"] is None,
                -(item["cnn"] if item["cnn"] is not None else -1e9),
                item["original_index"],
            ),
        )
    elif rank_col:
        cnn_sorted = sorted(
            parsed,
            key=lambda item: (
                item["explicit_rank"] is None,
                item["explicit_rank"] if item["explicit_rank"] is not None else 1e9,
                item["original_index"],
            ),
        )
    else:
        cnn_sorted = sorted(parsed, key=lambda item: item["original_index"])

    for rank, item in enumerate(cnn_sorted, start=1):
        item["computed_cnn_rank"] = rank

    best_rmsd = min(parsed, key=lambda item: item["rmsd"])
    top_cnn = cnn_sorted[0]

    top_cnn_score = top_cnn["cnn"]
    best_cnn_score = best_rmsd["cnn"]

    result = {
        "metrics_path": str(metrics_path),
        "snapshot_run": metrics_path.parent.name,
        "status": "complete",
        "row_count": len(parsed),
        "rmsd_column": rmsd_col,
        "cnn_column": cnn_col or "",
        "rank_column": rank_col or "",
        "pose_column": pose_col or "",
        "seed_column": seed_col or "",
        "top_cnn_pose": top_cnn["pose"],
        "top_cnn_seed": top_cnn["seed"],
        "top_cnn_score": top_cnn_score if top_cnn_score is not None else "",
        "top_cnn_rmsd": top_cnn["rmsd"],
        "best_rmsd_pose": best_rmsd["pose"],
        "best_rmsd_seed": best_rmsd["seed"],
        "best_rmsd": best_rmsd["rmsd"],
        "best_rmsd_cnn_score": best_cnn_score if best_cnn_score is not None else "",
        "best_rmsd_cnn_rank": best_rmsd["computed_cnn_rank"],
        "best_rmsd_pass": best_rmsd["rmsd"] <= threshold,
        "ranking_pass": top_cnn["rmsd"] <= threshold,
        "cnn_score_gap_top_minus_best_rmsd": (
            top_cnn_score - best_cnn_score
            if top_cnn_score is not None and best_cnn_score is not None
            else ""
        ),
    }

    for n in top_ns:
        top_n = cnn_sorted[:n]
        best_in_top_n = min(top_n, key=lambda item: item["rmsd"])
        result[f"best_rmsd_top_{n}"] = best_in_top_n["rmsd"]
        result[f"pass_in_top_{n}"] = best_in_top_n["rmsd"] <= threshold
        result[f"best_pose_in_top_{n}"] = best_in_top_n["pose"]
        result[f"best_seed_in_top_{n}"] = best_in_top_n["seed"]

    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--runroot",
        required=True,
        type=Path,
        help="Completed pc_overnight_md_ensemble_* result directory.",
    )
    parser.add_argument("--threshold", type=float, default=2.0)
    parser.add_argument("--top-n", type=int, nargs="+", default=[3, 5, 10, 20])
    args = parser.parse_args()

    runroot = args.runroot
    metrics_files = sorted(runroot.glob("g39_md_snapshots/*/pose_set_recovery_metrics.csv"))

    if not metrics_files:
        raise SystemExit(f"No pose_set_recovery_metrics.csv files found under {runroot}")

    rows = [
        summarize_metrics(path, threshold=args.threshold, top_ns=args.top_n)
        for path in metrics_files
    ]

    out_csv = runroot / "pose_retention_opportunity_summary.csv"
    out_md = runroot / "pose_retention_opportunity_report.md"

    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    complete = [row for row in rows if row.get("status") == "complete"]

    best_overall = min(
        complete,
        key=lambda row: parse_float(row.get("best_rmsd")) or 1e9,
        default=None,
    )

    sampling_passes = [
        row for row in complete
        if truthy(row.get("best_rmsd_pass"))
    ]

    ranking_passes = [
        row for row in complete
        if truthy(row.get("ranking_pass"))
    ]

    retained_counts = {}
    for n in args.top_n:
        retained_counts[n] = sum(
            1 for row in complete
            if truthy(row.get(f"pass_in_top_{n}"))
        )

    md_lines = [
        "# Pose Retention Opportunity Report",
        "",
        f"- Run root: `{runroot}`",
        f"- RMSD threshold: {args.threshold:.3f} Å",
        f"- Snapshot metrics files: {len(metrics_files)}",
        "",
        "## Summary",
        "",
        f"- Complete snapshot analyses: {len(complete)}",
        f"- Sampling passes: {len(sampling_passes)} / {len(complete)}",
        f"- Ranking passes: {len(ranking_passes)} / {len(complete)}",
        "",
    ]

    if best_overall:
        md_lines.extend(
            [
                "## Best Overall Sampled Pose",
                "",
                f"- Snapshot run: `{best_overall.get('snapshot_run')}`",
                f"- Best RMSD: {best_overall.get('best_rmsd')} Å",
                f"- CNN rank of best-RMSD pose: {best_overall.get('best_rmsd_cnn_rank')}",
                f"- Top-CNN RMSD for same snapshot: {best_overall.get('top_cnn_rmsd')} Å",
                f"- CNN score gap, top minus best-RMSD pose: {best_overall.get('cnn_score_gap_top_minus_best_rmsd')}",
                "",
            ]
        )

    md_lines.extend(
        [
            "## Would Top-N Retention Help?",
            "",
            "| Top-N retained poses | Snapshots with sub-threshold pose retained |",
            "|---:|---:|",
        ]
    )

    for n in args.top_n:
        md_lines.append(f"| {n} | {retained_counts[n]} / {len(complete)} |")

    md_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "If sampling passes but ranking does not, the docking engine is finding a plausible cognate-like pose but not ranking it first.",
            "That supports retaining alternative hypotheses instead of treating the top CNN pose as the only output.",
            "",
            "Recommended policy: keep the top CNN pose as the primary hypothesis, but retain alternative poses when they are near-ranked, benchmark-supported, contact-supported, or recurrent across receptor snapshots.",
            "",
            f"CSV output: `{out_csv.name}`",
            "",
        ]
    )

    out_md.write_text("\n".join(md_lines), encoding="utf-8")

    print("Wrote:", out_csv)
    print("Wrote:", out_md)
    print()

    print("=== SUMMARY ===")
    print(f"Complete snapshots: {len(complete)}")
    print(f"Sampling passes: {len(sampling_passes)} / {len(complete)}")
    print(f"Ranking passes: {len(ranking_passes)} / {len(complete)}")
    if best_overall:
        print(
            "Best sampled:",
            best_overall.get("best_rmsd"),
            "Å in",
            best_overall.get("snapshot_run"),
            "| CNN rank:",
            best_overall.get("best_rmsd_cnn_rank"),
            "| top-CNN RMSD:",
            best_overall.get("top_cnn_rmsd"),
        )

    print()
    print("Top-N retention:")
    for n in args.top_n:
        print(f"  Top {n}: {retained_counts[n]} / {len(complete)} snapshots retain sub-{args.threshold} Å pose")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
