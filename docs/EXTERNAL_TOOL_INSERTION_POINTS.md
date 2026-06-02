# External Tool Insertion Points

This file shows exactly where your partner inserts external program outputs.

## InterProScan
`annotation/interpro/interpro.tsv`

## CDD / CD-Search
`annotation/cdd/cdd_results.csv`

Columns:
`query,hit_id,short_name,description,start,end,evalue,bitscore,superfamily,sites`

## Pfam / HMMER
`annotation/pfam/pfam.tbl`

## VOGDB
`homologs/vogdb_hits.csv` and `homologs/homologs.fasta`

## MAFFT / MUSCLE
`msa/mafft_alignment.fasta` and `msa/muscle_alignment.fasta`

## AlphaFold / ColabFold
`structure/receptor.pdb` and optional `structure/pae.json`

## MolProbity / QMEAN
`structure/molprobity.csv` and `structure/qmean.csv`

## ConSurf
`conservation/consurf_scores.csv`

Columns:
`chain,residue_number,residue_name,consurf_grade,conservation_score,exposure`

## RING / Cytoscape / RINalyzer
`psn/ring_edges.csv`

## fpocket / P2Rank / DoGSiteScorer
`pockets/fpocket_out/`, `pockets/p2rank_predictions.csv`, `pockets/dogsite_pockets.csv`

## PDBFixer / ChimeraX / Meeko receptor prep
`receptor/receptor_prep_status.yaml`

## ChEMBL / BindingDB / PubChem Structure Search / ZINC-22
`ligands/raw/ligand_candidates.csv`

## RDKit / Open Babel / Meeko ligand prep
RDKit runs inside the code for initial filtering and clean SDF generation.
Open Babel and Meeko are marked external prep steps.

## GNINA
`docking/gnina_scores.csv`
