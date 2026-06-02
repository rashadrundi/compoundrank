#!/usr/bin/env python3
import argparse
from pathlib import Path
from _common import copy_file
ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); ap.add_argument('--input',required=True); a=ap.parse_args()
copy_file(a.input, Path(a.run_dir)/'structure/receptor.pdb')
