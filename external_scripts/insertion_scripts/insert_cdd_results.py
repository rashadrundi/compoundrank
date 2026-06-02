#!/usr/bin/env python3
import argparse, pandas as pd
from pathlib import Path
from _common import read_table_auto, write_csv
EXPECTED=['query','hit_id','short_name','description','start','end','evalue','bitscore','superfamily','sites']
ALIASES={'query':['query','query_id','protein','protein_id'],'hit_id':['hit_id','accession','acc','pssm_id','cdd_id'],'short_name':['short_name','shortname','name','domain','title'],'description':['description','desc','definition','full_name'],'start':['start','from','query_start','q_start'],'end':['end','to','query_end','q_end'],'evalue':['evalue','e-value','expect','eval'],'bitscore':['bitscore','bit_score','score'],'superfamily':['superfamily','clade','superfamily_id'],'sites':['sites','features','functional_sites','site_annotations']}
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); ap.add_argument('--input',required=True); a=ap.parse_args(); df=read_table_auto(a.input); lower={c.lower():c for c in df.columns}; out={}
 for col,names in ALIASES.items(): out[col]=df[next((lower[n] for n in names if n in lower), None)] if any(n in lower for n in names) else ''
 write_csv(pd.DataFrame(out), Path(a.run_dir)/'annotation/cdd/cdd_results.csv', EXPECTED)
if __name__=='__main__': main()
