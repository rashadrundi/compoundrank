#!/usr/bin/env python3
import argparse
from pathlib import Path
TEMPLATES={
'annotation/cdd/cdd_results.csv':'query,hit_id,short_name,description,start,end,evalue,bitscore,superfamily,sites\n',
'homologs/vogdb_hits.csv':'query,vog_id,vfam_id,vfold_id,annotation,evalue,score,virus_group\n',
'conservation/consurf_scores.csv':'chain,residue_number,residue_name,consurf_grade,conservation_score,exposure\n',
'conservation/functional_residues.csv':'chain,residue_number,label\n',
'psn/ring_edges.csv':'res1_chain,res1_number,res1_name,res2_chain,res2_number,res2_name,interaction_type,weight\n',
'pockets/p2rank_predictions.csv':'rank,name,score,probability,center_x,center_y,center_z,residue_ids\n',
'pockets/dogsite_pockets.csv':'pocket_id,score,druggability_score,volume,surface,depth,center_x,center_y,center_z\n',
'ligands/raw/ligand_candidates.csv':'ligand_id,name,smiles,source,target_class,evidence_type,evidence_notes,known_activity_value,known_activity_units\n',
'docking/gnina_scores.csv':'ligand_id,name,pocket_id,cnn_score,cnn_affinity,affinity,pose_file,contacts,notes\n'}
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); a=ap.parse_args(); run=Path(a.run_dir)
    for rel,header in TEMPLATES.items():
        p=run/rel; p.parent.mkdir(parents=True,exist_ok=True)
        if not p.exists(): p.write_text(header); print('[CREATE]',p)
        else: print('[EXISTS]',p)
if __name__=='__main__': main()
