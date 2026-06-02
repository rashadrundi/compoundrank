#!/usr/bin/env python3
import argparse
from pathlib import Path
T="""pdbfixer:
  status: {pdbfixer_status}
  output_file: "{pdbfixer_output}"
  notes: "Update after actual PDBFixer repair."
chimerax:
  status: {chimerax_status}
  output_file: "{chimerax_output}"
  notes: "Update after actual ChimeraX inspection."
meeko_receptor:
  status: {meeko_status}
  output_file: "{meeko_output}"
  notes: "Update after actual Meeko receptor prep."
docking_box:
  status: {box_status}
  center_x: {center_x}
  center_y: {center_y}
  center_z: {center_z}
  size_x: {size_x}
  size_y: {size_y}
  size_z: {size_z}
  notes: "Docking box centered on prioritized pocket."
"""
ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); ap.add_argument('--center-x',default='0'); ap.add_argument('--center-y',default='0'); ap.add_argument('--center-z',default='0'); ap.add_argument('--size-x',default='20'); ap.add_argument('--size-y',default='20'); ap.add_argument('--size-z',default='20'); ap.add_argument('--pdbfixer-status',default='done'); ap.add_argument('--chimerax-status',default='done'); ap.add_argument('--meeko-status',default='done'); ap.add_argument('--box-status',default='done'); a=ap.parse_args(); run=Path(a.run_dir); out=run/'receptor/receptor_prep_status.yaml'; out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(T.format(pdbfixer_status=a.pdbfixer_status,pdbfixer_output=str(run/'receptor/protein_fixed.pdb'),chimerax_status=a.chimerax_status,chimerax_output=str(run/'receptor/protein_clean_checked.pdb'),meeko_status=a.meeko_status,meeko_output=str(run/'receptor/receptor_prepared.pdbqt'),box_status=a.box_status,center_x=a.center_x,center_y=a.center_y,center_z=a.center_z,size_x=a.size_x,size_y=a.size_y,size_z=a.size_z)); print('[WRITE]',out)
