# CompoundRank local docking hypotheses

This production reset starts from an externally stored receptor PDB and one or
more externally stored or retrieved ligands. It deliberately does **not** use a
custom pose-ranking formula.

The active workflow is:

```text
receptor PDB + ligand source
-> receptor/ligand preparation in external cache
-> explicit, autobox, or fpocket pocket definition
-> GNINA docking across several seeds
-> PoseBusters physical-validity rejection
-> geometric clustering into distinct binding hypotheses
-> generic distance-based interaction evidence
-> compound ordering by GNINA CNN score only
-> uncertainty and alternate hypotheses
-> final receptor-ligand PDB files
```

All runtime data must be outside the repository. The standard sibling layout is:

```text
compoundrank-data/
├── inputs/
├── preserved/pre-reset-output/
├── results/
├── work/
└── cache/
```

## Runtime dependencies

Python environment:

```bash
source ~/.venvs/compoundrank-docking/bin/activate
python -m pip install -e . --no-deps
```

Required commands:

```text
gnina
fpocket
obabel
pdb2pqr
mk_prepare_receptor.py
mk_prepare_ligand.py
~/.venvs/posebusters/bin/bust
```

The first five should already be available in the docking environment or host.
PoseBusters remains isolated in its separate environment and is called by its
absolute executable path.

## Basic run using the preserved HIV-1 fold

Locate an existing preserved receptor without moving or refolding it:

```bash
DATA=/mnt/c/Users/Kausr/OneDrive/Desktop/compoundrank-data
find "$DATA/preserved/pre-reset-output/hiv1_protease_dimer" \
  -type f -name '*rank_001*.pdb' -print
```

Copy or place the ligand source under the external input directory, then run:

```bash
REPO=/mnt/c/Users/Kausr/OneDrive/Desktop/compoundrank
DATA=/mnt/c/Users/Kausr/OneDrive/Desktop/compoundrank-data
RECEPTOR=/absolute/path/to/the/preserved/receptor.pdb
LIGAND=/absolute/path/to/ligand.sdf

cd "$REPO"
source ~/.venvs/compoundrank-docking/bin/activate

python -m compoundrank \
  --receptor "$RECEPTOR" \
  --ligand-file "$LIGAND" \
  --data-root "$DATA" \
  --run-name first_clean_run \
  --seeds 2026 3101 4202 \
  --exhaustiveness 32 \
  --num-modes 20 \
  --max-hypotheses 3
```

With no pocket coordinates supplied, the program runs fpocket inside a temporary
external work directory and uses the highest fpocket score.

## Explicit pocket

```bash
python -m compoundrank \
  --receptor "$RECEPTOR" \
  --ligand-file "$LIGAND" \
  --data-root "$DATA" \
  --run-name explicit_box_test \
  --center-x 9.675 --center-y 2.414 --center-z -0.713 \
  --size-x 20.265 --size-y 19.865 --size-z 21.682
```

## Reference-ligand autobox

```bash
python -m compoundrank \
  --receptor "$RECEPTOR" \
  --ligand-file "$LIGAND" \
  --autobox-ligand /absolute/path/to/reference_ligand.sdf \
  --data-root "$DATA" \
  --run-name autobox_test
```

## Compound retrieval and multiple compounds

Repeatable PubChem retrieval:

```bash
python -m compoundrank \
  --receptor "$RECEPTOR" \
  --ligand-cid 213039 \
  --ligand-cid 441243 \
  --data-root "$DATA" \
  --run-name cid_screen
```

Manifest format:

```csv
name,source_type,value
darunavir,file,/absolute/path/to/darunavir.sdf
saquinavir,cid,441243
ethanol,smiles,CCO
```

Run:

```bash
python -m compoundrank \
  --receptor "$RECEPTOR" \
  --ligand-manifest /absolute/path/to/ligands.csv \
  --data-root "$DATA" \
  --run-name manifest_screen
```

## Output policy

The results directory contains only PDB files such as:

```text
01__darunavir__hypothesis_01.pdb
01__darunavir__hypothesis_02.pdb
02__saquinavir__hypothesis_01.pdb
```

Each PDB embeds GNINA scores, validity status, uncertainty, cluster recurrence,
pocket source, and interaction evidence in `REMARK 900` records. Temporary SDF,
PDBQT, logs, CSV files, and fpocket output remain in the external work directory
only during execution and are deleted after success. Use `--keep-workdir` only
for debugging.

## Important interpretation

- GNINA provides pose generation and pose ordering.
- PoseBusters rejects physically invalid poses; it does not establish the true
  biological pose.
- Clustering removes redundant poses and exposes alternate binding hypotheses;
  it does not create a new score.
- Compound priority is GNINA CNN score only.
- Interaction evidence and uncertainty are descriptive, not proof of binding,
  inhibition, antiviral efficacy, selectivity, or safety.
