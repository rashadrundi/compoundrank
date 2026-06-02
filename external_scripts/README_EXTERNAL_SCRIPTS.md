# External Run + Insert Scripts Add-on

Drop this folder next to the existing `antiviral_stage1_pipeline` code. It provides one script per external-tool/output-insertion point.

## Basic pattern

```text
run external tool or download result
        ↓
use an insertion script to place/normalize output into the framework
        ↓
run the master pipeline
```

Example:

```bash
bash external_scripts/run_tool_scripts/01_run_interproscan.sh runs/sars_mpro_test
python external_scripts/insertion_scripts/insert_cdd_results.py --run-dir runs/sars_mpro_test --input raw_cdd.csv
bash external_scripts/run_tool_scripts/03_run_mafft_muscle.sh runs/sars_mpro_test
python run_pipeline.py --config runs/sars_mpro_test/pipeline_config.yaml
```

## What this add-on contains

- `run_tool_scripts/`: shell scripts for tools that can be run locally.
- `insertion_scripts/`: Python scripts that normalize/copy outputs into the exact folders expected by the framework.
- `templates/`: blank CSV templates for manual insertion.

## Expected framework locations

| Stage | Expected file/folder |
|---|---|
| InterProScan | `annotation/interpro/interpro.tsv` |
| CDD | `annotation/cdd/cdd_results.csv` |
| Pfam/HMMER | `annotation/pfam/pfam.tbl` |
| VOGDB | `homologs/vogdb_hits.csv` |
| Homolog FASTA | `homologs/homologs.fasta` |
| MAFFT | `msa/mafft_alignment.fasta` |
| MUSCLE | `msa/muscle_alignment.fasta` |
| Structure | `structure/receptor.pdb` |
| ConSurf | `conservation/consurf_scores.csv` |
| RING/PSN | `psn/ring_edges.csv` |
| fpocket | `pockets/fpocket_out/` |
| P2Rank | `pockets/p2rank_predictions.csv` |
| DoGSiteScorer | `pockets/dogsite_pockets.csv` |
| Receptor prep | `receptor/receptor_prep_status.yaml` |
| Ligands | `ligands/raw/ligand_candidates.csv` |
| GNINA | `docking/gnina_scores.csv` |
