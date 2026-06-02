# Runbook Commands

## 1. Test integrated example

```bash
python run_pipeline.py --config example_project/pipeline_config.yaml
```

## 2. Make a new run

```bash
mkdir -p runs
cp -r example_project runs/first_real_test
perl -pi -e 's#example_project#runs/first_real_test#g' runs/first_real_test/pipeline_config.yaml
```

## 3. Put protein sequence in

```text
runs/first_real_test/input/protein.fasta
```

## 4. Create empty expected files

```bash
python external_scripts/insertion_scripts/make_empty_expected_files.py --run-dir runs/first_real_test
```

## 5. Check installed programs

```bash
bash external_scripts/run_tool_scripts/00_check_environment.sh
```

## 6. Tool runs / insertions

```bash
bash external_scripts/run_tool_scripts/01_run_interproscan.sh runs/first_real_test
bash external_scripts/run_tool_scripts/02_run_hmmscan_pfam.sh runs/first_real_test /path/to/Pfam-A.hmm
python external_scripts/insertion_scripts/insert_cdd_results.py --run-dir runs/first_real_test --input raw_outputs/cdd_results.csv
python external_scripts/insertion_scripts/insert_vogdb_hits.py --run-dir runs/first_real_test --input raw_outputs/vogdb_hits.csv
bash external_scripts/run_tool_scripts/03_run_mafft_muscle.sh runs/first_real_test
python external_scripts/insertion_scripts/copy_structure_pdb.py --run-dir runs/first_real_test --input raw_outputs/receptor.pdb
python external_scripts/insertion_scripts/insert_consurf_scores.py --run-dir runs/first_real_test --input raw_outputs/consurf_scores.csv
python external_scripts/insertion_scripts/insert_ring_edges.py --run-dir runs/first_real_test --input raw_outputs/ring_edges.csv
bash external_scripts/run_tool_scripts/04_run_fpocket.sh runs/first_real_test
python external_scripts/insertion_scripts/insert_p2rank_predictions.py --run-dir runs/first_real_test --input raw_outputs/p2rank_predictions.csv
python external_scripts/insertion_scripts/insert_dogsite_pockets.py --run-dir runs/first_real_test --input raw_outputs/dogsite_pockets.csv
python external_scripts/insertion_scripts/insert_ligand_candidates.py --run-dir runs/first_real_test --input raw_outputs/ligand_candidates.csv
```

## 7. Run pipeline once to create clean ligand SDF

```bash
python run_pipeline.py --config runs/first_real_test/pipeline_config.yaml
```

## 8. Prepare ligands

```bash
bash external_scripts/run_tool_scripts/05_openbabel_meeko_ligand_prep.sh runs/first_real_test
```

## 9. Record receptor prep and docking box

```bash
python external_scripts/insertion_scripts/insert_receptor_prep_status.py \
  --run-dir runs/first_real_test \
  --center-x 10 --center-y 2 --center-z 1 \
  --size-x 20 --size-y 20 --size-z 20
```

## 10. Run GNINA template

```bash
CENTER_X=10 CENTER_Y=2 CENTER_Z=1 SIZE_X=20 SIZE_Y=20 SIZE_Z=20 \
bash external_scripts/run_tool_scripts/06_run_gnina_template.sh runs/first_real_test
```

Then summarize GNINA results into:

```text
runs/first_real_test/docking/gnina_scores.csv
```

or insert an existing summary:

```bash
python external_scripts/insertion_scripts/insert_gnina_scores.py --run-dir runs/first_real_test --input raw_outputs/gnina_scores.csv
```

## 11. Validate and final run

```bash
python external_scripts/insertion_scripts/validate_required_files.py --run-dir runs/first_real_test
python run_pipeline.py --config runs/first_real_test/pipeline_config.yaml
```
