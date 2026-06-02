#!/usr/bin/env bash
set -euo pipefail
# ### EXTERNAL TOOL REQUIRED HERE: Pfam / HMMER ###
# Usage: bash external_scripts/run_tool_scripts/02_run_hmmscan_pfam.sh runs/YOUR_RUN /path/to/Pfam-A.hmm
RUN_DIR="${1:?Usage: $0 RUN_DIR PFAM_A_HMM}"
PFAM_HMM="${2:?Usage: $0 RUN_DIR PFAM_A_HMM}"
FASTA="$RUN_DIR/input/protein.fasta"
OUT="$RUN_DIR/annotation/pfam/pfam.tbl"
mkdir -p "$(dirname "$OUT")"
command -v hmmscan >/dev/null 2>&1 || { echo '[ERROR] hmmscan not found'; exit 1; }
[[ -f "$PFAM_HMM" ]] || { echo "[ERROR] Pfam HMM missing: $PFAM_HMM"; exit 1; }
[[ -f "$FASTA" ]] || { echo "[ERROR] FASTA missing: $FASTA"; exit 1; }
hmmscan --tblout "$OUT" "$PFAM_HMM" "$FASTA"
echo "[DONE] $OUT"
