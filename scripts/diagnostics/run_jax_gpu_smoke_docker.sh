#!/usr/bin/env bash
set -uo pipefail

IMAGE="${IMAGE:-exorcist:latest}"
DATA="${DATA:-/mnt/c/Users/kausr/OneDrive/Desktop/compoundrank-data}"
RUN_NAME="${RUN_NAME:-jax_gpu_smoke_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="$DATA/results/$RUN_NAME"
CACHE_DIR="$DATA/cache"
CONTAINER_NAME="${CONTAINER_NAME:-jax_gpu_smoke_${RUN_NAME}}"

mkdir -p "$RUN_DIR" "$CACHE_DIR"

cat > "$RUN_DIR/jax_probe.py" <<'PY'
import os
import sys
import importlib.metadata as md

print("python:", sys.version)
print("executable:", sys.executable)

for key in [
    "CUDA_VISIBLE_DEVICES",
    "JAX_PLATFORMS",
    "JAX_PLATFORM_NAME",
    "XLA_PYTHON_CLIENT_PREALLOCATE",
    "XLA_PYTHON_CLIENT_MEM_FRACTION",
    "LD_LIBRARY_PATH",
]:
    print(f"{key}={os.environ.get(key)}")

for pkg in [
    "colabfold",
    "jax",
    "jaxlib",
    "numpy",
    "tensorflow",
    "tensorflow-cpu",
]:
    try:
        print(f"{pkg}: {md.version(pkg)}")
    except Exception as exc:
        print(f"{pkg}: not found ({exc})")

import jax
import jax.numpy as jnp

print("jax backend:", jax.default_backend())
print("jax devices:", jax.devices())

x = jnp.ones((256, 256), dtype=jnp.float32)
y = (x @ x).block_until_ready()

print("jax matmul ok:", float(y[0, 0]))
PY

echo "IMAGE=$IMAGE"
echo "DATA=$DATA"
echo "RUN_DIR=$RUN_DIR"
echo "CONTAINER_NAME=$CONTAINER_NAME"

echo
echo "=== CPU-FORCED JAX PROBE ==="

set +e
docker run --rm \
  -v "$DATA:/data" \
  -v "$CACHE_DIR:/cache" \
  -e JAX_PLATFORMS=cpu \
  -e JAX_PLATFORM_NAME=cpu \
  -e CUDA_VISIBLE_DEVICES="" \
  --entrypoint /bin/bash \
  "$IMAGE" -lc "python /data/results/$RUN_NAME/jax_probe.py" \
  > "$RUN_DIR/jax_cpu_probe.log" 2>&1
CPU_EXIT=$?
set -e

echo "$CPU_EXIT" > "$RUN_DIR/jax_cpu_exit_code.txt"
echo "CPU probe exit: $CPU_EXIT"

echo
echo "=== GPU JAX PROBE ==="

docker rm -f "$CONTAINER_NAME" >/dev/null 2>&1 || true

set +e
docker run \
  --name "$CONTAINER_NAME" \
  --gpus all \
  -v "$DATA:/data" \
  -v "$CACHE_DIR:/cache" \
  -e XLA_PYTHON_CLIENT_PREALLOCATE=false \
  -e XLA_PYTHON_CLIENT_MEM_FRACTION=0.60 \
  -e TF_FORCE_UNIFIED_MEMORY=0 \
  --entrypoint /bin/bash \
  "$IMAGE" -lc "
set -uo pipefail
nvidia-smi > /data/results/$RUN_NAME/nvidia_smi_inside_container.txt 2>&1 || true
python /data/results/$RUN_NAME/jax_probe.py
" \
  > "$RUN_DIR/jax_gpu_probe.log" 2>&1
GPU_EXIT=$?
set -e

echo "$GPU_EXIT" > "$RUN_DIR/jax_gpu_exit_code.txt"
echo "GPU probe exit: $GPU_EXIT"

docker logs "$CONTAINER_NAME" > "$RUN_DIR/docker.log" 2>&1 || true
docker inspect "$CONTAINER_NAME" > "$RUN_DIR/docker_inspect.json" 2>&1 || true
docker stats --no-stream "$CONTAINER_NAME" > "$RUN_DIR/docker_stats.txt" 2>&1 || true

python - <<PY
import json
from pathlib import Path

cpu_exit = int("$CPU_EXIT")
gpu_exit = int("$GPU_EXIT")

if gpu_exit == 0:
    status = "pass"
    reason_code = "jax_gpu_smoke_passed"
elif gpu_exit == 139:
    status = "fail"
    reason_code = "jax_gpu_segmentation_fault"
elif cpu_exit != 0:
    status = "fail"
    reason_code = "jax_environment_failed"
else:
    status = "fail"
    reason_code = "jax_gpu_smoke_failed"

payload = {
    "image": "$IMAGE",
    "run_name": "$RUN_NAME",
    "run_dir": "$RUN_DIR",
    "container_name": "$CONTAINER_NAME",
    "cpu_exit_code": cpu_exit,
    "gpu_exit_code": gpu_exit,
    "status": status,
    "reason_code": reason_code,
    "interpretation": (
        "JAX GPU smoke passed."
        if gpu_exit == 0
        else "JAX GPU smoke failed before ColabFold should be trusted."
    ),
}

Path("$RUN_DIR/diagnostic_summary.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\\n",
    encoding="utf-8",
)
PY

BUNDLE="$DATA/results/${RUN_NAME}_diagnostic_bundle.tgz"
tar -czf "$BUNDLE" -C "$RUN_DIR" .

echo
echo "=== SUMMARY ==="
python -m json.tool "$RUN_DIR/diagnostic_summary.json"

echo
echo "=== CPU LOG TAIL ==="
tail -80 "$RUN_DIR/jax_cpu_probe.log" || true

echo
echo "=== GPU LOG TAIL ==="
tail -120 "$RUN_DIR/jax_gpu_probe.log" || true

echo
echo "Diagnostic bundle:"
ls -lh "$BUNDLE"

echo
echo "Container kept for inspection:"
echo "docker logs $CONTAINER_NAME | tail -120"
echo "docker rm $CONTAINER_NAME"
