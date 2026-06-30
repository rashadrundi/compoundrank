#!/bin/bash
#SBATCH --job-name=exorcist-g39-mdpose
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-50%5
#SBATCH --output=logs/g39_mdpose_%A_%a.out
#SBATCH --error=logs/g39_mdpose_%A_%a.err

set -euo pipefail

# EDIT THESE FOR AMAREL
PROJECT="$HOME/compoundrank"
DATA="$HOME/compoundrank-data"
BENCH="$DATA/benchmarks/influenza_neuraminidase_2HU0"
SNAPLIST="$BENCH/md_snapshot_receptors.txt"
REF="$BENCH/2HU0_chainB_bound_ligand_reference.sdf"
MANIFEST_G39="$BENCH/docking_manifest_g39_reference_only.csv"
BOX="$BENCH/reference_box_chainB_bound_ligand_tight18.json"
RESULT_ROOT="$DATA/results/md_ensemble_g39_pose_recovery"

mkdir -p logs
mkdir -p "$RESULT_ROOT"

cd "$PROJECT"

# Activate your Amarel compoundrank environment.
# EDIT if your venv path differs.
source "$HOME/.venvs/compoundrank-docking/bin/activate"

SNAPSHOT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SNAPLIST")"

if [ -z "${SNAPSHOT:-}" ]; then
  echo "No snapshot for task ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 0
fi

SNAP_BASENAME="$(basename "$SNAPSHOT" .pdb)"
OUT="$RESULT_ROOT/${SNAP_BASENAME}_g39_seed10_modes20_exh24"

echo "Task: $SLURM_ARRAY_TASK_ID"
echo "Snapshot: $SNAPSHOT"
echo "Output: $OUT"

python -m compoundrank \
  --receptor "$SNAPSHOT" \
  --data-root "$DATA" \
  --output-dir "$OUT" \
  --ligand-manifest "$MANIFEST_G39" \
  --box-json "$BOX" \
  --reference-ligand "$REF" \
  --seeds 10 \
  --num-modes 20 \
  --exhaustiveness 24 \
  --skip-validity \
  --keep-workdir \
  --overwrite

grep -R "POSE_RECOVERY\|RMSD\|Sampling pass\|Ranking pass" "$OUT" -n || true
