#!/usr/bin/env python3
import argparse
from pathlib import Path
FILES=['input/protein.fasta','annotation/interpro/interpro.tsv','annotation/cdd/cdd_results.csv','annotation/pfam/pfam.tbl','homologs/vogdb_hits.csv','homologs/homologs.fasta','msa/mafft_alignment.fasta','msa/muscle_alignment.fasta','structure/receptor.pdb','conservation/consurf_scores.csv','psn/ring_edges.csv','pockets/p2rank_predictions.csv','pockets/dogsite_pockets.csv','receptor/receptor_prep_status.yaml','ligands/raw/ligand_candidates.csv','docking/gnina_scores.csv']
ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); a=ap.parse_args(); run=Path(a.run_dir); missing=[]
for rel in FILES:
 p=run/rel
 ok=p.exists() and p.stat().st_size>0
 print(('[OK]      ' if ok else '[MISSING] ')+rel)
 if not ok: missing.append(rel)
fp=run/'pockets/fpocket_out'; ok=fp.exists() and any(fp.iterdir())
print(('[OK]      ' if ok else '[MISSING] ')+'pockets/fpocket_out/')
print('\nMissing/empty count:', len(missing)+(0 if ok else 1))
