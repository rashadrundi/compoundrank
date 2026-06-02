#!/usr/bin/env python3
import argparse
from pathlib import Path
from _common import copy_file
import subprocess, sys
ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); ap.add_argument('--hits',required=True); ap.add_argument('--homolog-fasta',required=True); a=ap.parse_args()
subprocess.check_call([sys.executable, str(Path(__file__).with_name('insert_vogdb_hits.py')), '--run-dir', a.run_dir, '--input', a.hits])
copy_file(a.homolog_fasta, Path(a.run_dir)/'homologs/homologs.fasta')
