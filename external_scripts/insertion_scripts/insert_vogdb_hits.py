#!/usr/bin/env python3
import argparse
from pathlib import Path
from _common import read_table_auto, write_csv
EXPECTED=['query', 'vog_id', 'vfam_id', 'vfold_id', 'annotation', 'evalue', 'score', 'virus_group']
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); ap.add_argument('--input',required=True); a=ap.parse_args()
    df=read_table_auto(a.input)
    for c in EXPECTED:
        if c not in df.columns: df[c]=''
    write_csv(df[EXPECTED], Path(a.run_dir)/'homologs/vogdb_hits.csv', EXPECTED)
if __name__=='__main__': main()
