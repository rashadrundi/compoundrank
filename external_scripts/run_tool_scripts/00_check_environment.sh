#!/usr/bin/env bash
set -euo pipefail
echo 'Checking common external programs...'
for tool in interproscan.sh hmmscan mafft muscle fpocket gnina obabel mk_prepare_ligand.py mk_prepare_receptor.py; do
  if command -v "$tool" >/dev/null 2>&1; then echo "[FOUND] $tool -> $(command -v $tool)"; else echo "[MISSING] $tool"; fi
done
