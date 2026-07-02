#!/usr/bin/env bash
set -uo pipefail

IMAGE="${IMAGE:-exorcist:latest}"
DATA="${DATA:-/mnt/c/Users/kausr/OneDrive/Desktop/compoundrank-data}"
RUN_NAME="${RUN_NAME:-colabfold_smoke_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="$DATA/results/$RUN_NAME"
INPUT_DIR="$DATA/inputs/colabfold_smoke"
CACHE_DIR="$DATA/cache"
CONTAINER_NAME="${CONTAINER_NAME:-colabfold_smoke_${RUN_NAME}}"

mkdir -p "$RUN_DIR" "$INPUT_DIR" "$CACHE_DIR"

cat > "$INPUT_DIR/gb1.faa" <<'FASTA'
>gb1_tiny_fold_smoke
MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE
FASTA

docker rm -f "$CONTAINER_NAME" >/dev/null 2>&1 || true

echo "IMAGE=$IMAGE"
echo "DATA=$DATA"
echo "RUN_DIR=$RUN_DIR"
echo "CONTAINER_NAME=$CONTAINER_NAME"

set +e
docker run \
  --name "$CONTAINER_NAME" \
  --gpus all \
  -v "$DATA:/data" \
  -v "$CACHE_DIR:/cache" \
  -e RUN_NAME="$RUN_NAME" \
  -e XLA_PYTHON_CLIENT_PREALLOCATE=false \
  -e XLA_PYTHON_CLIENT_MEM_FRACTION=0.60 \
  -e TF_FORCE_UNIFIED_MEMORY=0 \
  --entrypoint /bin/bash \
  "$IMAGE" -lc '
set -uo pipefail

RUN_DIR="/data/results/$RUN_NAME"

{
  echo "=== DATE ==="
  date

  echo
  echo "=== GPU ==="
  nvidia-smi || true

  echo
  echo "=== COLABFOLD BIN ==="
  which colabfold_batch || true
  colabfold_batch --help 2>/dev/null | head -30 || true

  echo
  echo "=== PYTHON/JAX ENVIRONMENT ==="
  python - <<PY
import sys
import importlib.metadata as md

print("python:", sys.version)

for pkg in [
    "colabfold",
    "jax",
    "jaxlib",
    "numpy",
    "tensorflow",
    "tensorflow-cpu",
    "absl-py",
]:
    try:
        print(f"{pkg}:", md.version(pkg))
    except Exception as exc:
        print(f"{pkg}: not found ({exc})")

try:
    import jax
    import jax.numpy as jnp
    print("jax devices:", jax.devices())
    x = jnp.ones((256, 256))
    y = (x @ x).block_until_ready()
    print("jax matmul ok:", float(y[0, 0]))
except BaseException as exc:
    print("jax smoke failed:", repr(exc))
PY
  echo "python_probe_exit=$?"

  echo
  echo "=== CUDA LIBS ==="
  ldconfig -p 2>/dev/null | grep -E "libcuda|cudnn|cublas|cusolver" | head -80 || true
} 2>&1 | tee "$RUN_DIR/environment_probe.log"

echo
echo "=== RUN COLABFOLD ===" | tee "$RUN_DIR/colabfold_console.log"

colabfold_batch \
  --msa-mode single_sequence \
  --num-models 1 \
  --num-recycle 1 \
  --num-seeds 1 \
  --model-type alphafold2_ptm \
  --disable-unified-memory \
  --overwrite-existing-results \
  /data/inputs/colabfold_smoke/gb1.faa \
  "$RUN_DIR" \
  2>&1 | tee -a "$RUN_DIR/colabfold_console.log"

COLABFOLD_EXIT=$?
echo "$COLABFOLD_EXIT" > "$RUN_DIR/colabfold_exit_code.txt"

echo
echo "=== PDB OUTPUTS ===" | tee -a "$RUN_DIR/colabfold_console.log"
find "$RUN_DIR" -type f -name "*.pdb" -print | tee "$RUN_DIR/pdb_outputs.txt"

exit "$COLABFOLD_EXIT"
'
DOCKER_EXIT=$?
set -e

docker logs "$CONTAINER_NAME" > "$RUN_DIR/docker.log" 2>&1 || true
docker inspect "$CONTAINER_NAME" > "$RUN_DIR/docker_inspect.json" 2>&1 || true
docker stats --no-stream "$CONTAINER_NAME" > "$RUN_DIR/docker_stats.txt" 2>&1 || true

PDB_COUNT="$(find "$RUN_DIR" -type f -name '*.pdb' | wc -l | tr -d ' ')"

python - <<PY
import json
from pathlib import Path

run_dir = Path("$RUN_DIR")
payload = {
    "image": "$IMAGE",
    "run_name": "$RUN_NAME",
    "run_dir": "$RUN_DIR",
    "container_name": "$CONTAINER_NAME",
    "docker_exit_code": $DOCKER_EXIT,
    "pdb_count": int("$PDB_COUNT"),
    "status": "pass" if int("$PDB_COUNT") > 0 else "fail",
    "interpretation": (
        "ColabFold produced at least one PDB."
        if int("$PDB_COUNT") > 0
        else "ColabFold did not produce a PDB. Inspect environment_probe.log and colabfold_console.log."
    ),
}
(run_dir / "diagnostic_summary.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

BUNDLE="$DATA/results/${RUN_NAME}_diagnostic_bundle.tgz"
tar -czf "$BUNDLE" -C "$RUN_DIR" .

echo
echo "=== SUMMARY ==="
python -m json.tool "$RUN_DIR/diagnostic_summary.json"

echo
echo "Diagnostic bundle:"
ls -lh "$BUNDLE"

echo
echo "Container kept for inspection:"
echo "docker logs $CONTAINER_NAME | tail -120"
echo "docker rm $CONTAINER_NAME"
