#!/usr/bin/env bash
set -uo pipefail

###############################################################################
# EXORCIST / CompoundRank PC overnight MD-ensemble benchmark
#
# What this does:
# 1. Validates local PC paths.
# 2. Builds/repairs md_snapshot_receptors.txt using absolute PDB paths.
# 3. Creates G39 reference manifest.
# 4. Creates oseltamivir/peramivir filtered manifest if needed.
# 5. Runs static G39 cognate redocking baseline.
# 6. Runs G39 pose recovery across every MD snapshot.
# 7. Aggregates best RMSD per snapshot.
# 8. Docks oseltamivir/peramivir on the top G39-performing snapshots.
#
# IMPORTANT:
# This consumes existing MD snapshot receptor PDBs.
# Put receptor-only snapshot PDBs in:
#   $BENCH/md_snapshots/
###############################################################################

REPO="/mnt/c/Users/kausr/OneDrive/Desktop/compoundrank"
DATA="/mnt/c/Users/kausr/OneDrive/Desktop/compoundrank-data"
BENCH="$DATA/benchmarks/influenza_neuraminidase_2HU0"
PREV="$DATA/results/auto_retrieve_2hu0_pc_structure_fetch_20260629_193117"

VENV="$HOME/.venvs/compoundrank-docking/bin/activate"

CHAINB_RECEPTOR="$BENCH/2HU0_chainB_receptor_clean.pdb"
REF="$BENCH/2HU0_chainB_bound_ligand_reference.sdf"
BOX="$BENCH/reference_box_chainB_bound_ligand_tight18.json"

SNAPDIR="$BENCH/md_snapshots"
SNAPLIST="$BENCH/md_snapshot_receptors.txt"

MANIFEST_G39="$BENCH/docking_manifest_g39_reference_only.csv"
SOURCE_MANIFEST="$PREV/stage4a_compound_retrieval/docking_manifest.csv"
MANIFEST_RETRIEVED="$PREV/stage4a_compound_retrieval/docking_manifest_oseltamivir_peramivir.csv"

STAMP="$(date +%Y%m%d_%H%M%S)"
RUNROOT="$DATA/results/pc_overnight_md_ensemble_$STAMP"
LOGDIR="$RUNROOT/logs"

# Tune these if needed.
MAX_SNAPSHOTS="${MAX_SNAPSHOTS:-9999}"
TOP_SNAPSHOTS_FOR_RETRIEVED="${TOP_SNAPSHOTS_FOR_RETRIEVED:-5}"

STATIC_SEEDS="${STATIC_SEEDS:-10}"
STATIC_MODES="${STATIC_MODES:-20}"
STATIC_EXHAUSTIVENESS="${STATIC_EXHAUSTIVENESS:-24}"

G39_SEEDS="${G39_SEEDS:-10}"
G39_MODES="${G39_MODES:-20}"
G39_EXHAUSTIVENESS="${G39_EXHAUSTIVENESS:-24}"

RETRIEVED_SEEDS="${RETRIEVED_SEEDS:-5}"
RETRIEVED_MODES="${RETRIEVED_MODES:-10}"
RETRIEVED_EXHAUSTIVENESS="${RETRIEVED_EXHAUSTIVENESS:-16}"

mkdir -p "$RUNROOT" "$LOGDIR"

MASTER_LOG="$RUNROOT/overnight_master.log"
exec > >(tee -a "$MASTER_LOG") 2>&1

echo "========================================================================"
echo "EXORCIST PC OVERNIGHT MD-ENSEMBLE BENCHMARK"
echo "Started: $(date)"
echo "Run root: $RUNROOT"
echo "========================================================================"

require_file() {
  local label="$1"
  local path="$2"
  if [ ! -f "$path" ]; then
    echo "ERROR: Missing $label:"
    echo "  $path"
    exit 2
  fi
}

require_dir() {
  local label="$1"
  local path="$2"
  if [ ! -d "$path" ]; then
    echo "ERROR: Missing $label:"
    echo "  $path"
    exit 2
  fi
}

require_file "repo module" "$REPO/compoundrank/__main__.py"
require_file "chain-B receptor" "$CHAINB_RECEPTOR"
require_file "G39 reference ligand SDF" "$REF"
require_file "tight18 box JSON" "$BOX"
require_file "source Stage 4A docking manifest" "$SOURCE_MANIFEST"

if [ ! -f "$VENV" ]; then
  echo "ERROR: Missing venv activation file:"
  echo "  $VENV"
  exit 2
fi

cd "$REPO"
source "$VENV"

echo
echo "=== ENVIRONMENT ==="
echo "Repo: $REPO"
echo "Data: $DATA"
echo "Python: $(which python)"
python -c "import compoundrank; print('compoundrank import OK')"
command -v gnina || echo "WARNING: gnina not found on PATH"
command -v obabel || echo "WARNING: obabel not found on PATH"

echo
echo "=== CREATE G39 MANIFEST ==="
cat > "$MANIFEST_G39" <<MANIFEST
name,source_type,value
G39_REFERENCE,file,$REF
MANIFEST
cat "$MANIFEST_G39"

echo
echo "=== CREATE RETRIEVED OS/PER MANIFEST ==="
python - <<PY
from pathlib import Path
import csv

src = Path("$SOURCE_MANIFEST")
dst = Path("$MANIFEST_RETRIEVED")

wanted = {"OSELTAMIVIR", "PERAMIVIR ANHYDROUS", "peramivir_anhydrous", "oseltamivir"}

with src.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    rows = [
        row for row in reader
        if row.get("name") in wanted or row.get("compound_name") in wanted
    ]

if not rows:
    raise SystemExit("No oseltamivir/peramivir rows found in source manifest")

with dst.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(dst)
print("rows:", len(rows))
PY

echo
echo "=== BUILD ABSOLUTE MD SNAPSHOT LIST ==="
mkdir -p "$SNAPDIR"

# Only real files ending in .pdb. This avoids the previous '.' bug.
find "$SNAPDIR" -maxdepth 1 -type f -name "*.pdb" -exec readlink -f {} \; | sort > "$SNAPLIST"

SNAPCOUNT="$(wc -l < "$SNAPLIST" | tr -d ' ')"
echo "Snapshot directory: $SNAPDIR"
echo "Snapshot list: $SNAPLIST"
echo "Snapshot count: $SNAPCOUNT"

if [ "$SNAPCOUNT" -eq 0 ]; then
  echo
  echo "ERROR: No MD snapshot PDBs found."
  echo "Put receptor-only MD snapshot PDB files here:"
  echo "  $SNAPDIR"
  echo
  echo "Example filenames:"
  echo "  snapshot_000.pdb"
  echo "  snapshot_001.pdb"
  echo
  echo "Stopping instead of pretending this is an MD ensemble."
  exit 3
fi

echo
echo "First snapshots:"
nl -ba "$SNAPLIST" | head -20

echo
echo "=== RUN STATIC G39 BASELINE ==="
STATIC_OUT="$RUNROOT/00_static_chainB_g39_reference"

if [ ! -f "$STATIC_OUT/pose_set_recovery_summary.json" ]; then
  mkdir -p "$STATIC_OUT"
  (
    python -m compoundrank \
      --receptor "$CHAINB_RECEPTOR" \
      --data-root "$DATA" \
      --output-dir "$STATIC_OUT" \
      --ligand-manifest "$MANIFEST_G39" \
      --box-json "$BOX" \
      --reference-ligand "$REF" \
      --seeds "$STATIC_SEEDS" \
      --num-modes "$STATIC_MODES" \
      --exhaustiveness "$STATIC_EXHAUSTIVENESS" \
      --skip-validity \
      --keep-workdir \
      --overwrite
  ) > "$LOGDIR/00_static_g39.log" 2>&1
  STATUS=$?
  echo "Static baseline exit status: $STATUS"
else
  echo "Static baseline already complete; skipping."
fi

grep -R "POSE_RECOVERY\|RMSD\|Sampling pass\|Ranking pass" "$STATIC_OUT" -n || true

echo
echo "=== RUN G39 POSE RECOVERY ACROSS MD SNAPSHOTS ==="

MAP_CSV="$RUNROOT/snapshot_run_map.csv"
echo "index,snapshot_path,out_dir,status" > "$MAP_CSV"

i=0
while IFS= read -r SNAPSHOT; do
  i=$((i + 1))

  if [ "$i" -gt "$MAX_SNAPSHOTS" ]; then
    echo "Reached MAX_SNAPSHOTS=$MAX_SNAPSHOTS; stopping snapshot loop."
    break
  fi

  if [ ! -f "$SNAPSHOT" ]; then
    echo "$i,$SNAPSHOT,,missing_snapshot" >> "$MAP_CSV"
    echo "[$i] Missing snapshot: $SNAPSHOT"
    continue
  fi

  SNAP_BASENAME="$(basename "$SNAPSHOT" .pdb)"
  OUT="$RUNROOT/g39_md_snapshots/${i}__${SNAP_BASENAME}_g39_seed${G39_SEEDS}_modes${G39_MODES}_exh${G39_EXHAUSTIVENESS}"
  LOG="$LOGDIR/g39_snapshot_${i}__${SNAP_BASENAME}.log"
  mkdir -p "$OUT"

  echo
  echo "[$i/$SNAPCOUNT] G39 snapshot:"
  echo "  Snapshot: $SNAPSHOT"
  echo "  Output:   $OUT"
  echo "  Log:      $LOG"

  (
    python -m compoundrank \
      --receptor "$SNAPSHOT" \
      --data-root "$DATA" \
      --output-dir "$OUT" \
      --ligand-manifest "$MANIFEST_G39" \
      --box-json "$BOX" \
      --reference-ligand "$REF" \
      --seeds "$G39_SEEDS" \
      --num-modes "$G39_MODES" \
      --exhaustiveness "$G39_EXHAUSTIVENESS" \
      --skip-validity \
      --keep-workdir \
      --overwrite
  ) > "$LOG" 2>&1

  STATUS=$?
  if [ "$STATUS" -eq 0 ]; then
    echo "$i,$SNAPSHOT,$OUT,complete" >> "$MAP_CSV"
    echo "  Status: complete"
    grep -E "POSE_RECOVERY|RMSD|Sampling pass|Ranking pass" "$LOG" || true
  else
    echo "$i,$SNAPSHOT,$OUT,failed_exit_${STATUS}" >> "$MAP_CSV"
    echo "  Status: failed_exit_${STATUS}"
    tail -80 "$LOG" || true
  fi

done < "$SNAPLIST"

echo
echo "=== AGGREGATE G39 POSE RECOVERY ==="

python - <<PY
from pathlib import Path
import csv
import json
import re
import math

runroot = Path("$RUNROOT")
map_csv = runroot / "snapshot_run_map.csv"
summary_csv = runroot / "g39_md_pose_recovery_summary.csv"
top_txt = runroot / "top_g39_snapshots.txt"

def parse_float(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    m = re.search(r"[-+]?[0-9]*\\.?[0-9]+", str(value))
    return float(m.group(0)) if m else None

def deep_find(obj, needles):
    if isinstance(obj, dict):
        for k, v in obj.items():
            kl = str(k).lower()
            if all(n in kl for n in needles):
                fv = parse_float(v)
                if fv is not None:
                    return fv
        for v in obj.values():
            found = deep_find(v, needles)
            if found is not None:
                return found
    elif isinstance(obj, list):
        for v in obj:
            found = deep_find(v, needles)
            if found is not None:
                return found
    return None

def parse_report(path):
    result = {}
    if not path.exists():
        return result
    text = path.read_text(errors="replace")
    patterns = {
        "top_cnn_pose_rmsd": r"Top CNN pose RMSD\\s*\\|\\s*([0-9.]+)",
        "best_sampled_rmsd": r"Best sampled(?: pose)? RMSD\\s*\\|\\s*([0-9.]+)",
        "sampling_pass_text": r"Sampling pass\\s*\\|\\s*([^|\\n]+)",
        "ranking_pass_text": r"Ranking pass\\s*\\|\\s*([^|\\n]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text, re.I)
        if m:
            result[key] = m.group(1).strip()
    return result

rows = []

with map_csv.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    for mrow in reader:
        out_dir = Path(mrow.get("out_dir") or "")
        row = {
            "index": mrow.get("index"),
            "snapshot_path": mrow.get("snapshot_path"),
            "out_dir": str(out_dir) if str(out_dir) != "." else "",
            "status": mrow.get("status"),
            "sampling_pass": "",
            "ranking_pass": "",
            "best_sampled_rmsd": "",
            "top_cnn_pose_rmsd": "",
            "rmsd_threshold": "",
        }

        if out_dir and out_dir.exists():
            summary_json = out_dir / "pose_set_recovery_summary.json"
            if summary_json.exists():
                try:
                    data = json.loads(summary_json.read_text())
                    row["sampling_pass"] = data.get("sampling_pass", "")
                    row["ranking_pass"] = data.get("ranking_pass", "")
                    row["best_sampled_rmsd"] = (
                        data.get("best_sampled_rmsd")
                        or data.get("best_sampled_pose_rmsd")
                        or deep_find(data, ["best", "rmsd"])
                        or ""
                    )
                    row["top_cnn_pose_rmsd"] = (
                        data.get("top_cnn_pose_rmsd")
                        or deep_find(data, ["top", "rmsd"])
                        or ""
                    )
                    row["rmsd_threshold"] = (
                        data.get("cognate_rmsd_threshold")
                        or data.get("rmsd_threshold")
                        or deep_find(data, ["threshold"])
                        or ""
                    )
                except Exception as e:
                    row["status"] = f"{row['status']};json_error:{e}"

            report_data = parse_report(out_dir / "pose_set_recovery_report.md")
            for key, val in report_data.items():
                if key == "sampling_pass_text" and not row["sampling_pass"]:
                    row["sampling_pass"] = val
                elif key == "ranking_pass_text" and not row["ranking_pass"]:
                    row["ranking_pass"] = val
                elif key in row and not row[key]:
                    row[key] = val

        rows.append(row)

fieldnames = [
    "index",
    "snapshot_path",
    "status",
    "sampling_pass",
    "ranking_pass",
    "best_sampled_rmsd",
    "top_cnn_pose_rmsd",
    "rmsd_threshold",
    "out_dir",
]

with summary_csv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

valid = []
for row in rows:
    rmsd = parse_float(row.get("best_sampled_rmsd"))
    if rmsd is not None and math.isfinite(rmsd):
        valid.append((rmsd, row))

valid.sort(key=lambda item: item[0])

with top_txt.open("w", encoding="utf-8") as handle:
    for rmsd, row in valid[:int("$TOP_SNAPSHOTS_FOR_RETRIEVED")]:
        handle.write(row["snapshot_path"] + "\\n")

print("Summary CSV:", summary_csv)
print("Top snapshots:", top_txt)
print()
print("Top 10 by best sampled RMSD:")
for rmsd, row in valid[:10]:
    print(
        f"RMSD={rmsd:.3f}",
        "top=", row.get("top_cnn_pose_rmsd"),
        "sampling=", row.get("sampling_pass"),
        "ranking=", row.get("ranking_pass"),
        "snapshot=", Path(row.get("snapshot_path", "")).name,
    )
PY

echo
echo "=== DOCK RETRIEVED LIGANDS ON TOP G39 SNAPSHOTS ==="

TOP_SNAPSHOTS="$RUNROOT/top_g39_snapshots.txt"

if [ ! -s "$TOP_SNAPSHOTS" ]; then
  echo "No top snapshots found. Skipping retrieved-ligand docking."
else
  j=0
  while IFS= read -r SNAPSHOT; do
    j=$((j + 1))

    if [ ! -f "$SNAPSHOT" ]; then
      echo "Top snapshot missing: $SNAPSHOT"
      continue
    fi

    SNAP_BASENAME="$(basename "$SNAPSHOT" .pdb)"
    OUT="$RUNROOT/retrieved_on_top_g39_snapshots/${j}__${SNAP_BASENAME}_retrieved_seed${RETRIEVED_SEEDS}_modes${RETRIEVED_MODES}_exh${RETRIEVED_EXHAUSTIVENESS}"
    LOG="$LOGDIR/retrieved_top_${j}__${SNAP_BASENAME}.log"
    mkdir -p "$OUT"

    echo
    echo "Retrieved docking on top snapshot $j:"
    echo "  Snapshot: $SNAPSHOT"
    echo "  Output:   $OUT"

    (
      python -m compoundrank \
        --receptor "$SNAPSHOT" \
        --data-root "$DATA" \
        --output-dir "$OUT" \
        --ligand-manifest "$MANIFEST_RETRIEVED" \
        --box-json "$BOX" \
        --seeds "$RETRIEVED_SEEDS" \
        --num-modes "$RETRIEVED_MODES" \
        --exhaustiveness "$RETRIEVED_EXHAUSTIVENESS" \
        --skip-validity \
        --keep-workdir \
        --overwrite
    ) > "$LOG" 2>&1

    STATUS=$?
    echo "  Status: $STATUS"
    if [ "$STATUS" -ne 0 ]; then
      tail -80 "$LOG" || true
    else
      grep -E "COMPOUND PRIORITY|CNNscore|HYPOTHESES|Final PDB" "$LOG" || true
    fi

  done < "$TOP_SNAPSHOTS"
fi

echo
echo "========================================================================"
echo "FINISHED PC OVERNIGHT MD-ENSEMBLE BENCHMARK"
echo "Finished: $(date)"
echo "Run root: $RUNROOT"
echo "Master log: $MASTER_LOG"
echo "G39 summary: $RUNROOT/g39_md_pose_recovery_summary.csv"
echo "Top snapshots: $RUNROOT/top_g39_snapshots.txt"
echo "========================================================================"
