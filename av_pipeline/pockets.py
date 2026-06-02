from pathlib import Path
import pandas as pd
import re
from .utils import safe_read_csv, confidence_from_score

def parse_residue_string(value):
    if value is None or pd.isna(value):
        return set()
    residues = set()
    for token in re.split(r"[;, ]+", str(value)):
        m = re.match(r"([A-Za-z_])[:]?([0-9]+)", token.strip())
        if m:
            residues.add((m.group(1), int(m.group(2))))
    return residues

def parse_pdb_residues(path):
    residues = set()
    p = Path(path)
    if not p.exists():
        return residues
    for line in p.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                residues.add((line[21].strip() or "_", int(line[22:26])))
            except Exception:
                pass
    return residues

def parse_fpocket_dir(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: fpocket ###
    Parses permissive fpocket info files.
    """
    d = Path(path)
    if not d.exists():
        return pd.DataFrame()
    rows = []
    for f in list(d.rglob("*info*.txt")) + list(d.rglob("*.info")):
        text = f.read_text(errors="ignore")
        data = {"source": "fpocket", "pocket_id": f.stem, "file": str(f)}
        for key, pat in {
            "fpocket_score": r"Score\s*:\s*([-+0-9.eE]+)",
            "druggability_score": r"Druggability Score\s*:\s*([-+0-9.eE]+)",
            "volume": r"Volume\s*:\s*([-+0-9.eE]+)",
            "polarity_score": r"Polarity score\s*:\s*([-+0-9.eE]+)",
            "hydrophobicity_score": r"Hydrophobicity score\s*:\s*([-+0-9.eE]+)",
            "num_alpha_spheres": r"Number of Alpha Spheres\s*:\s*([-+0-9.eE]+)",
        }.items():
            m = re.search(pat, text, re.I)
            if m:
                try: data[key] = float(m.group(1))
                except Exception: data[key] = m.group(1)
        residues = set()
        for pdb in f.parent.glob(f"{f.stem}*.pdb"):
            residues |= parse_pdb_residues(pdb)
        if residues:
            data["residues"] = ";".join(f"{c}:{n}" for c,n in sorted(residues))
        rows.append(data)
    return pd.DataFrame(rows)

def parse_p2rank(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: P2Rank / PrankWeb ###
    Expected flexible columns: name,rank,score,probability,center_x,center_y,center_z,residue_ids
    """
    df = safe_read_csv(path)
    if df.empty: return df
    df = df.copy()
    df["source"] = "P2Rank"
    if "pocket_id" not in df.columns:
        df["pocket_id"] = df["name"].astype(str) if "name" in df.columns else [f"p2rank_{i+1}" for i in range(len(df))]
    return df

def parse_dogsite(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: DoGSiteScorer ###
    Expected flexible columns: pocket_id,score,druggability_score,volume,surface,depth,center_x,center_y,center_z
    """
    df = safe_read_csv(path)
    if df.empty: return df
    df = df.copy()
    df["source"] = "DoGSiteScorer"
    if "pocket_id" not in df.columns:
        df["pocket_id"] = [f"dogsite_{i+1}" for i in range(len(df))]
    return df

def prioritize_pockets(fpocket, p2rank, dogsite, conserved=set(), central=set(), functional=set()):
    primary = fpocket if fpocket is not None and not fpocket.empty else (p2rank if p2rank is not None and not p2rank.empty else dogsite)
    if primary is None or primary.empty:
        return pd.DataFrame()
    rows = []
    for i, r in primary.iterrows():
        pid = str(r.get("pocket_id", f"pocket_{i+1}"))
        residues = parse_residue_string(r.get("residues", r.get("residue_ids", "")))
        score, notes = 0, []
        try:
            drug = float(r.get("druggability_score"))
            score += min(20, drug*20 if drug <= 1 else drug)
            notes.append("Druggability score parsed.")
        except Exception: pass
        try:
            fp = float(r.get("fpocket_score", r.get("score")))
            score += min(15, max(0, fp))
            notes.append("Pocket score parsed.")
        except Exception: pass
        try:
            vol = float(r.get("volume"))
            if 150 <= vol <= 900:
                score += 15
                notes.append("Volume broadly ligand-compatible.")
            elif 80 <= vol < 150 or 900 < vol <= 1500:
                score += 8
                notes.append("Volume marginal.")
        except Exception: pass
        if p2rank is not None and not p2rank.empty:
            score += 15
            notes.append("P2Rank output available.")
        if dogsite is not None and not dogsite.empty:
            score += 10
            notes.append("DoGSiteScorer output available.")
        if residues:
            nc = len(residues & conserved)
            ncent = len(residues & central)
            nf = len(residues & functional)
            if nc:
                score += min(20, 5 + 3*nc); notes.append(f"{nc} conserved residues in pocket.")
            if ncent:
                score += min(10, 3 + 2*ncent); notes.append(f"{ncent} PSN-central residues in pocket.")
            if nf:
                score += min(15, 5 + 5*nf); notes.append(f"{nf} functional residues in pocket.")
        score = round(min(100, score), 2)
        rows.append({"pocket_id": pid, "pocket_priority_score": score, "pocket_confidence": confidence_from_score(score), "residues": r.get("residues", r.get("residue_ids", "")), "volume": r.get("volume", ""), "druggability_score": r.get("druggability_score", ""), "notes": " ".join(notes)})
    return pd.DataFrame(rows).sort_values("pocket_priority_score", ascending=False)
