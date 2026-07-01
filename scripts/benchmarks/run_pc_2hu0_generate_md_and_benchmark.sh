#!/usr/bin/env bash
set -euo pipefail

REPO="/mnt/c/Users/kausr/OneDrive/Desktop/compoundrank"
DATA="/mnt/c/Users/kausr/OneDrive/Desktop/compoundrank-data"
BENCH="$DATA/benchmarks/influenza_neuraminidase_2HU0"

VENV="$HOME/.venvs/compoundrank-docking/bin/activate"

REFERENCE_RECEPTOR="$BENCH/2HU0_chainB_receptor_clean.pdb"
PREPARED_OPENMM_RECEPTOR="$BENCH/derived/openmm_md_ensemble_2HU0/prepared/receptor_pdb2pqr.pdb"
TARGET_SNAPDIR="$BENCH/md_snapshots"

STAMP="$(date +%Y%m%d_%H%M%S)"
RUNROOT="$DATA/results/pc_2hu0_openmm_md_generation_$STAMP"
RAW_SNAPDIR="$RUNROOT/raw_snapshots"
ALIGNED_SNAPDIR="$RUNROOT/aligned_snapshots"
LOG="$RUNROOT/openmm_md_generation.log"

N_SNAPSHOTS="${N_SNAPSHOTS:-20}"
EQUIL_PS="${EQUIL_PS:-25}"
PRODUCTION_PS="${PRODUCTION_PS:-250}"
TIMESTEP_FS="${TIMESTEP_FS:-2}"
TEMPERATURE_K="${TEMPERATURE_K:-300}"

mkdir -p "$RUNROOT" "$RAW_SNAPDIR" "$ALIGNED_SNAPDIR" "$TARGET_SNAPDIR"

exec > >(tee -a "$LOG") 2>&1

echo "========================================================================"
echo "PC 2HU0 OPENMM MD GENERATION + BENCHMARK"
echo "Started: $(date)"
echo "RUNROOT=$RUNROOT"
echo "N_SNAPSHOTS=$N_SNAPSHOTS"
echo "EQUIL_PS=$EQUIL_PS"
echo "PRODUCTION_PS=$PRODUCTION_PS"
echo "TIMESTEP_FS=$TIMESTEP_FS"
echo "TEMPERATURE_K=$TEMPERATURE_K"
echo "========================================================================"

cd "$REPO"
source "$VENV"

echo
echo "=== ENVIRONMENT ==="
python -c "import sys; print(sys.executable)"
python -c "import compoundrank; print('compoundrank import OK')"

echo
echo "=== OPENMM CHECK ==="
python - <<'PY'
try:
    import openmm
    import openmm.app
    print("OpenMM import OK")
except Exception as error:
    raise SystemExit(
        "OpenMM is not available in this environment. Try:\n"
        "  python -m pip install openmm\n\n"
        f"Import error: {error}"
    )
PY

echo
echo "=== SELECT OPENMM INPUT RECEPTOR ==="

if [ -f "$PREPARED_OPENMM_RECEPTOR" ]; then
  OPENMM_INPUT_RECEPTOR="$PREPARED_OPENMM_RECEPTOR"
  echo "Using prepared OpenMM receptor:"
  echo "$OPENMM_INPUT_RECEPTOR"
else
  OPENMM_INPUT_RECEPTOR="$REFERENCE_RECEPTOR"
  echo "Prepared OpenMM receptor not found; using cleaned receptor:"
  echo "$OPENMM_INPUT_RECEPTOR"
fi

if [ ! -f "$OPENMM_INPUT_RECEPTOR" ]; then
  echo "ERROR: Missing OpenMM input receptor:"
  echo "$OPENMM_INPUT_RECEPTOR"
  exit 2
fi

echo
echo "=== GENERATE RAW OPENMM SNAPSHOTS ==="

python - <<PY
from pathlib import Path

from openmm import unit
from openmm import LangevinMiddleIntegrator, Platform
from openmm.app import (
    PDBFile,
    Modeller,
    ForceField,
    Simulation,
    NoCutoff,
    HBonds,
)

input_receptor = Path("$OPENMM_INPUT_RECEPTOR")
raw_snapdir = Path("$RAW_SNAPDIR")
runroot = Path("$RUNROOT")

n_snapshots = int("$N_SNAPSHOTS")
equil_ps = float("$EQUIL_PS")
production_ps = float("$PRODUCTION_PS")
timestep_fs = float("$TIMESTEP_FS")
temperature_k = float("$TEMPERATURE_K")

raw_snapdir.mkdir(parents=True, exist_ok=True)

print("Input receptor:", input_receptor)
print("Raw snapshot dir:", raw_snapdir)

pdb = PDBFile(str(input_receptor))
forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")

modeller = Modeller(pdb.topology, pdb.positions)

topology = None
positions = None

print("Trying addHydrogens...")
try:
    modeller.addHydrogens(forcefield, pH=7.4)
    topology = modeller.topology
    positions = modeller.positions
    print("addHydrogens OK")
except Exception as error:
    print("addHydrogens failed; trying existing topology as-is.")
    print("addHydrogens error:", error)
    topology = pdb.topology
    positions = pdb.positions

print("Creating OpenMM system...")
try:
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )
except Exception as first_error:
    print("System creation failed with selected topology.")
    print("Error:", first_error)

    if topology is not pdb.topology:
        print("Retrying system creation with original PDB topology.")
        topology = pdb.topology
        positions = pdb.positions
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
        )
    else:
        raise

integrator = LangevinMiddleIntegrator(
    temperature_k * unit.kelvin,
    1.0 / unit.picosecond,
    timestep_fs * unit.femtosecond,
)

platform = None
platform_name = None
for candidate in ["CUDA", "OpenCL", "CPU"]:
    try:
        platform = Platform.getPlatformByName(candidate)
        platform_name = candidate
        break
    except Exception:
        pass

if platform is None:
    raise SystemExit("No usable OpenMM platform found.")

print("Using OpenMM platform:", platform_name)

if platform_name == "CUDA":
    simulation = Simulation(
        topology,
        system,
        integrator,
        platform,
        {"Precision": "mixed"},
    )
else:
    simulation = Simulation(topology, system, integrator, platform)

simulation.context.setPositions(positions)

print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=1000)

timestep_ps = timestep_fs / 1000.0
equil_steps = int(round(equil_ps / timestep_ps))
production_steps = int(round(production_ps / timestep_ps))
snapshot_interval = max(1, production_steps // n_snapshots)

print("Equilibration steps:", equil_steps)
print("Production steps:", production_steps)
print("Snapshot interval:", snapshot_interval)

if equil_steps > 0:
    print("Equilibrating...")
    simulation.step(equil_steps)

state = simulation.context.getState(getPositions=True)
equilibrated = runroot / "receptor_equilibrated_openmm.pdb"
with equilibrated.open("w") as handle:
    PDBFile.writeFile(simulation.topology, state.getPositions(), handle)

print("Equilibrated receptor:", equilibrated)
print("Production...")

for idx in range(1, n_snapshots + 1):
    simulation.step(snapshot_interval)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    outfile = raw_snapdir / f"snapshot_{idx:04d}.pdb"
    with outfile.open("w") as handle:
        PDBFile.writeFile(simulation.topology, state.getPositions(), handle)
    print(f"snapshot {idx:04d}: {outfile} | potential={state.getPotentialEnergy()}")

print("Snapshots written:", n_snapshots)
PY

echo
echo "=== ALIGN SNAPSHOTS BACK TO REFERENCE FRAME ==="

python - <<PY
from pathlib import Path
import numpy as np

reference = Path("$REFERENCE_RECEPTOR")
raw_dir = Path("$RAW_SNAPDIR")
aligned_dir = Path("$ALIGNED_SNAPDIR")
target_dir = Path("$TARGET_SNAPDIR")

aligned_dir.mkdir(parents=True, exist_ok=True)
target_dir.mkdir(parents=True, exist_ok=True)

BACKBONE = {"N", "CA", "C", "O"}

def parse_atoms(path):
    atoms = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            atom = {
                "line": line,
                "atom_name": line[12:16].strip(),
                "chain": line[21].strip(),
                "resseq": line[22:26].strip(),
                "icode": line[26].strip(),
                "coord": np.array(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ],
                    dtype=float,
                ),
            }
        except Exception:
            continue
        atom["key"] = (
            atom["chain"],
            atom["resseq"],
            atom["icode"],
            atom["atom_name"],
        )
        atoms.append(atom)
    return atoms

def kabsch_mobile_to_ref(mobile, ref):
    mobile_centroid = mobile.mean(axis=0)
    ref_centroid = ref.mean(axis=0)
    mobile_centered = mobile - mobile_centroid
    ref_centered = ref - ref_centroid

    h = mobile_centered.T @ ref_centered
    u, s, vt = np.linalg.svd(h)
    r = vt.T @ u.T

    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T

    t = ref_centroid - mobile_centroid @ r
    return r, t

def transform_coord(coord, r, t):
    return coord @ r + t

def rewrite_pdb(path, atoms, r, t, outpath):
    atom_idx = 0
    output = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")) and atom_idx < len(atoms):
            atom = atoms[atom_idx]
            atom_idx += 1
            new_coord = transform_coord(atom["coord"], r, t)
            line = (
                f"{line[:30]}"
                f"{new_coord[0]:8.3f}"
                f"{new_coord[1]:8.3f}"
                f"{new_coord[2]:8.3f}"
                f"{line[54:]}"
            )
        output.append(line)

    if not output or output[-1].strip() != "END":
        output.append("END")

    outpath.write_text("\\n".join(output) + "\\n", encoding="utf-8")

ref_atoms = parse_atoms(reference)
ref_by_key = {a["key"]: a for a in ref_atoms if a["atom_name"] in BACKBONE}
ref_ca = [a for a in ref_atoms if a["atom_name"] == "CA"]

if len(ref_ca) < 20:
    raise SystemExit("Reference has too few CA atoms for fallback alignment.")

for raw in sorted(raw_dir.glob("snapshot_*.pdb")):
    mob_atoms = parse_atoms(raw)
    mob_by_key = {a["key"]: a for a in mob_atoms if a["atom_name"] in BACKBONE}

    common = sorted(set(ref_by_key) & set(mob_by_key))

    if len(common) >= 20:
        ref_coords = np.array([ref_by_key[k]["coord"] for k in common])
        mob_coords = np.array([mob_by_key[k]["coord"] for k in common])
        align_label = f"keyed backbone atoms={len(common)}"
    else:
        mob_ca = [a for a in mob_atoms if a["atom_name"] == "CA"]
        n = min(len(ref_ca), len(mob_ca))
        if n < 20:
            raise SystemExit(f"Too few atoms to align {raw}: keyed={len(common)}, CA={n}")
        ref_coords = np.array([a["coord"] for a in ref_ca[:n]])
        mob_coords = np.array([a["coord"] for a in mob_ca[:n]])
        align_label = f"sequential CA atoms={n}"

    r, t = kabsch_mobile_to_ref(mob_coords, ref_coords)
    aligned = aligned_dir / raw.name
    rewrite_pdb(raw, mob_atoms, r, t, aligned)

    target = target_dir / raw.name
    target.write_text(aligned.read_text(encoding="utf-8"), encoding="utf-8")

    aligned_atoms = parse_atoms(aligned)
    aligned_ca = [a for a in aligned_atoms if a["atom_name"] == "CA"]
    n = min(len(ref_ca), len(aligned_ca))
    rmsd = np.sqrt(
        ((np.array([a["coord"] for a in aligned_ca[:n]]) -
          np.array([a["coord"] for a in ref_ca[:n]])) ** 2).sum(axis=1).mean()
    )

    print(f"{raw.name}: aligned RMSD={rmsd:.3f} Å using {align_label}")

print("Aligned snapshots copied to:", target_dir)
PY

echo
echo "=== ACTIVE SNAPSHOT LIST ==="
find "$TARGET_SNAPDIR" -maxdepth 1 -type f -name "snapshot_*.pdb" -exec readlink -f {} \; | sort > "$BENCH/md_snapshot_receptors.txt"
nl -ba "$BENCH/md_snapshot_receptors.txt" | head -50
echo "Snapshot count: $(wc -l < "$BENCH/md_snapshot_receptors.txt")"

echo
echo "=== RUN EXORCIST MD-ENSEMBLE BENCHMARK ==="
cd "$REPO"
./scripts/benchmarks/run_pc_2hu0_md_ensemble_benchmark.sh

echo
echo "========================================================================"
echo "FINISHED PC 2HU0 OPENMM MD GENERATION + BENCHMARK"
echo "Finished: $(date)"
echo "MD generation root: $RUNROOT"
echo "Benchmark results are under the newest pc_overnight_md_ensemble_* folder."
echo "========================================================================"
