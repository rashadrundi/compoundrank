#!/usr/bin/env bash
set -euo pipefail
# ### EXTERNAL TOOL REQUIRED HERE: fpocket ###
# Usage: bash external_scripts/run_tool_scripts/04_run_fpocket.sh runs/YOUR_RUN
RUN_DIR="${1:?Usage: $0 RUN_DIR}"
RECEPTOR="$RUN_DIR/structure/receptor.pdb"
OUT_PARENT="$RUN_DIR/pockets"
EXPECTED="$RUN_DIR/pockets/fpocket_out"
mkdir -p "$OUT_PARENT"
command -v fpocket >/dev/null 2>&1 || { echo '[ERROR] fpocket not found'; exit 1; }
[[ -f "$RECEPTOR" ]] || { echo "[ERROR] receptor missing: $RECEPTOR"; exit 1; }
pushd "$OUT_PARENT" >/dev/null
fpocket -f "../structure/receptor.pdb"
popd >/dev/null
FOUND="$(find "$OUT_PARENT" -maxdepth 1 -type d -name '*_out' | head -n 1 || true)"
if [[ -n "$FOUND" ]]; then rm -rf "$EXPECTED"; cp -R "$FOUND" "$EXPECTED"; echo "[DONE] $EXPECTED"; else echo '[WARN] No fpocket *_out directory found'; fi
