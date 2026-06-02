#!/usr/bin/env bash
set -euo pipefail
# ### EXTERNAL TOOL REQUIRED HERE: MAFFT / MUSCLE ###
# Usage: bash external_scripts/run_tool_scripts/03_run_mafft_muscle.sh runs/YOUR_RUN
RUN_DIR="${1:?Usage: $0 RUN_DIR}"
HOMOLOGS="$RUN_DIR/homologs/homologs.fasta"
mkdir -p "$RUN_DIR/msa"
[[ -f "$HOMOLOGS" ]] || { echo "[ERROR] Homolog FASTA missing: $HOMOLOGS"; exit 1; }
if command -v mafft >/dev/null 2>&1; then mafft --auto "$HOMOLOGS" > "$RUN_DIR/msa/mafft_alignment.fasta"; echo '[DONE] MAFFT'; else echo '[MISSING] mafft'; fi
if command -v muscle >/dev/null 2>&1; then muscle -align "$HOMOLOGS" -output "$RUN_DIR/msa/muscle_alignment.fasta"; echo '[DONE] MUSCLE'; else echo '[MISSING] muscle'; fi
