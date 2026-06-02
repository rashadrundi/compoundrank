#!/usr/bin/env bash
set -euo pipefail
# ### EXTERNAL TOOL REQUIRED HERE: GNINA ###
# Usage with env vars for box:
# CENTER_X=10 CENTER_Y=2 CENTER_Z=1 bash external_scripts/run_tool_scripts/06_run_gnina_template.sh runs/YOUR_RUN
RUN_DIR="${1:?Usage: $0 RUN_DIR}"
RECEPTOR="$RUN_DIR/receptor/receptor_prepared.pdbqt"
LIG_DIR="$RUN_DIR/ligands/prepared/pdbqt_ligands"
OUT_DIR="$RUN_DIR/docking/gnina_raw"
mkdir -p "$OUT_DIR"
command -v gnina >/dev/null 2>&1 || { echo '[ERROR] gnina not found'; exit 1; }
[[ -f "$RECEPTOR" ]] || { echo "[ERROR] receptor PDBQT missing: $RECEPTOR"; exit 1; }
[[ -d "$LIG_DIR" ]] || { echo "[ERROR] ligand dir missing: $LIG_DIR"; exit 1; }
CX="${CENTER_X:-0}"; CY="${CENTER_Y:-0}"; CZ="${CENTER_Z:-0}"
SX="${SIZE_X:-20}"; SY="${SIZE_Y:-20}"; SZ="${SIZE_Z:-20}"
for lig in "$LIG_DIR"/*.pdbqt; do
  [ -e "$lig" ] || continue
  name="$(basename "$lig" .pdbqt)"
  gnina --receptor "$RECEPTOR" --ligand "$lig" --center_x "$CX" --center_y "$CY" --center_z "$CZ" --size_x "$SX" --size_y "$SY" --size_z "$SZ" --out "$OUT_DIR/${name}_out.sdf" > "$OUT_DIR/${name}.log" 2>&1
  echo "[DONE] $name"
done
echo '[NEXT] Summarize into docking/gnina_scores.csv using insertion_scripts/insert_gnina_scores.py'
