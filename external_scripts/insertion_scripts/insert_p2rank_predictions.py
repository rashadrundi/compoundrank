#!/usr/bin/env python3
import argparse
from pathlib import Path
from _common import read_table_auto, write_csv
EXPECTED=['rank', 'name', 'score', 'probability', 'center_x', 'center_y', 'center_z', 'residue_ids']
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); ap.add_argument('--input',required=True); a=ap.parse_args()
    df=read_table_auto(a.input)
    for c in EXPECTED:
        if c not in df.columns: df[c]=''
    write_csv(df[EXPECTED], Path(a.run_dir)/'pockets/p2rank_predictions.csv', EXPECTED)
if __name__=='__main__': main()
