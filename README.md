# Antiviral Stage 1 Pipeline

This is a defensive computational pipeline scaffold for creating a **protein–pocket–ligand evidence package** that can be handed to a biophysics simulation workflow.

## Goal

```text
viral protein FASTA
→ annotation evidence
→ homolog/MSA evidence
→ structure-quality evidence
→ conserved pocket prioritization
→ ligand discovery/filtering
→ receptor/ligand preparation tracking
→ docking evidence
→ biophysics handoff report
```

The pipeline does **not** prove a compound is an antiviral. It creates a structured hypothesis package for later simulation and validation.

## Safety scope

This codebase is only for therapeutic/defensive antiviral target and ligand prioritization. It must not be used to enhance viral fitness, infectivity, immune evasion, resistance, tropism, or pathogenicity.

## External-tool markers included

The package includes clearly marked points for:

```text
### EXTERNAL TOOL REQUIRED HERE: InterProScan ###
### EXTERNAL TOOL REQUIRED HERE: CDD / CD-Search ###
### EXTERNAL TOOL REQUIRED HERE: Pfam / HMMER ###
### EXTERNAL TOOL REQUIRED HERE: VOGDB ###
### EXTERNAL TOOL REQUIRED HERE: MAFFT ###
### EXTERNAL TOOL REQUIRED HERE: MUSCLE ###
### EXTERNAL TOOL REQUIRED HERE: AlphaFold / ColabFold ###
### EXTERNAL TOOL REQUIRED HERE: ConSurf ###
### EXTERNAL TOOL REQUIRED HERE: RING / Cytoscape / RINalyzer ###
### EXTERNAL TOOL REQUIRED HERE: fpocket ###
### EXTERNAL TOOL REQUIRED HERE: P2Rank / PrankWeb ###
### EXTERNAL TOOL REQUIRED HERE: DoGSiteScorer ###
### EXTERNAL TOOL REQUIRED HERE: PDBFixer ###
### EXTERNAL TOOL REQUIRED HERE: ChimeraX ###
### EXTERNAL TOOL REQUIRED HERE: Meeko receptor prep ###
### EXTERNAL TOOL REQUIRED HERE: ChEMBL ###
### EXTERNAL TOOL REQUIRED HERE: BindingDB ###
### EXTERNAL TOOL REQUIRED HERE: PubChem Structure Search ###
### EXTERNAL TOOL REQUIRED HERE: ZINC-22 ###
### EXTERNAL TOOL REQUIRED HERE: Open Babel ###
### EXTERNAL TOOL REQUIRED HERE: Meeko ligand prep ###
### EXTERNAL TOOL REQUIRED HERE: GNINA ###
```

## Install

Recommended:

```bash
conda env create -f environment.yml
conda activate antiviral-stage1
```

Or:

```bash
pip install -r requirements.txt
```

## Run example

```bash
python run_pipeline.py --config example_project/pipeline_config.yaml
```

Expected final outputs:

```text
example_project/reports/final_handoff_report.md
example_project/reports/biophysics_handoff_manifest.json
example_project/reports/candidate_table.csv
```

## Confidence flags

| Flag | Meaning |
|---|---|
| HIGH | Multiple evidence layers agree and no major blocking problems detected. |
| MODERATE | Useful hypothesis, but some evidence is missing or incomplete. |
| LOW | Exploratory only; evidence is sparse or inconsistent. |
| FAIL | Not ready for handoff without fixing missing/invalid evidence. |


## Python file key

See `PY_FILE_KEY.md` for a plain-English explanation of every `.py` file.
