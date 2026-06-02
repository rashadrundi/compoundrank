#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="${1:-}"

if [[ -z "$RUN_DIR" ]]; then
  echo "Usage: 02_run_interproscan.sh <run_dir>" >&2
  exit 2
fi

IPS_CMD="${INTERPROSCAN_CMD:-interproscan.sh}"
CPU="${INTERPROSCAN_CPU:-2}"

INPUT="$RUN_DIR/input/protein.fasta"
OUT_DIR="$RUN_DIR/annotation/interpro"
OUT_FILE="$OUT_DIR/interpro.tsv"
LOG_DIR="$RUN_DIR/logs"

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

if [[ ! -f "$INPUT" ]]; then
  echo "[InterProScan] Missing input FASTA: $INPUT" >&2
  exit 1
fi

echo "[InterProScan] Input: $INPUT"
echo "[InterProScan] Output: $OUT_FILE"
echo "[InterProScan] Command: $IPS_CMD"

"$IPS_CMD" \
  -i "$INPUT" \
  -f TSV \
  -o "$OUT_FILE" \
  -iprlookup \
  -goterms \
  -pa \
  -cpu "$CPU" \
  > "$LOG_DIR/interproscan_stdout.log" \
  2> "$LOG_DIR/interproscan_stderr.log"

echo "[InterProScan] Done."
echo "[InterProScan] Wrote: $OUT_FILE"
