#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 RUNROOT"
  echo
  echo "Optional:"
  echo "  KNOWN_RESIDUES='ARG118:B,ASP151:B,...' $0 RUNROOT"
  exit 2
fi

RUNROOT="$1"
SNAPROOT="$RUNROOT/retrieved_on_top_g39_snapshots"
KNOWN_RESIDUES="${KNOWN_RESIDUES:-}"

if [ ! -d "$RUNROOT" ]; then
  echo "ERROR: run root does not exist: $RUNROOT" >&2
  exit 1
fi

if [ ! -d "$SNAPROOT" ]; then
  echo "ERROR: retrieved snapshot directory does not exist: $SNAPROOT" >&2
  exit 1
fi

echo "Run root: $RUNROOT"
echo "Retrieved snapshot root: $SNAPROOT"

if [ -n "$KNOWN_RESIDUES" ]; then
  echo "Known residues: $KNOWN_RESIDUES"
else
  echo "Known residues: none"
fi

echo
echo "=== Per-snapshot production pose retention ==="

count=0
failures=0

for OUTDIR in "$SNAPROOT"/*; do
  if [ ! -d "$OUTDIR" ]; then
    continue
  fi

  count=$((count + 1))
  echo
  echo "--- $OUTDIR ---"

  if [ -n "$KNOWN_RESIDUES" ]; then
    python -m compoundrank.pose_contact_retention \
      --output-dir "$OUTDIR" \
      --known-residues "$KNOWN_RESIDUES" || failures=$((failures + 1))
  else
    python -m compoundrank.pose_contact_retention \
      --output-dir "$OUTDIR" || failures=$((failures + 1))
  fi
done

echo
echo "Per-snapshot outputs attempted: $count"
echo "Per-snapshot failures: $failures"

if [ "$count" -eq 0 ]; then
  echo "ERROR: no retrieved snapshot output directories found." >&2
  exit 1
fi

if [ "$failures" -ne 0 ]; then
  echo "ERROR: one or more per-snapshot retention runs failed." >&2
  exit 1
fi

echo
echo "=== Ensemble production pose retention ==="
python -m compoundrank.pose_retention_ensemble \
  --runroot "$RUNROOT"

echo
echo "=== Regenerating top-level run report ==="
python -m compoundrank.run_report \
  --output-dir "$RUNROOT"

echo
echo "=== Retention artifact check ==="
find "$SNAPROOT" -name "production_pose_retention_candidates.csv" | wc -l
find "$SNAPROOT" -name "production_pose_retention_report.md" | wc -l

test -f "$RUNROOT/production_pose_retention_ensemble_summary.csv"
test -f "$RUNROOT/production_pose_retention_ensemble_report.md"
test -f "$RUNROOT/compoundrank_run_report.md"

echo
echo "=== Top-level report tier check ==="
grep -n "Compound-Level Confidence Tier Summary\|cleaner_support\|warning_dominated\|primary_supported\|strong_alternative\|physically_warned" \
  "$RUNROOT/compoundrank_run_report.md" || true

echo
echo "Production pose-retention suite complete."
