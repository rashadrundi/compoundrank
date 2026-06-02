#!/usr/bin/env bash
set -euo pipefail
# ### EXTERNAL TOOL REQUIRED HERE: Open Babel / Meeko ligand prep ###
# Usage: bash external_scripts/run_tool_scripts/05_openbabel_meeko_ligand_prep.sh runs/YOUR_RUN
RUN_DIR="${1:?Usage: $0 RUN_DIR}"
CLEAN="$RUN_DIR/ligands/clean/clean_ligands.sdf"
WITH_H="$RUN_DIR/ligands/prepared/ligands_with_h.sdf"
PDBQT_DIR="$RUN_DIR/ligands/prepared/pdbqt_ligands"
mkdir -p "$RUN_DIR/ligands/prepared"
[[ -f "$CLEAN" ]] || { echo "[ERROR] clean_ligands.sdf missing. Run main pipeline after inserting ligand_candidates.csv first."; exit 1; }
command -v obabel >/dev/null 2>&1 || { echo '[ERROR] obabel not found'; exit 1; }
command -v mk_prepare_ligand.py >/dev/null 2>&1 || { echo '[ERROR] mk_prepare_ligand.py not found'; exit 1; }
obabel "$CLEAN" -O "$WITH_H" -h
mkdir -p "$PDBQT_DIR"
mk_prepare_ligand.py -i "$WITH_H" --multimol_outdir "$PDBQT_DIR"
echo "[DONE] $PDBQT_DIR"
echo '[REMINDER] Manually inspect top ligands for protonation, tautomer, stereochemistry, and charge.'
