from pathlib import Path
import pandas as pd
import statistics, json
from .utils import safe_read_csv, confidence_from_score

def parse_plddt_from_pdb(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: AlphaFold / ColabFold ###
    pLDDT is often stored in the PDB B-factor column.
    """
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    residues = {}
    for line in p.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            chain = line[21].strip() or "_"
            resi = int(line[22:26])
            resn = line[17:20].strip()
            b = float(line[60:66])
        except Exception:
            continue
        residues.setdefault((chain, resi, resn), []).append(b)
    rows = [{"chain": c, "residue_number": n, "residue_name": r, "plddt": statistics.mean(vals)} for (c,n,r), vals in residues.items()]
    return pd.DataFrame(rows)

def summarize_plddt(df):
    if df is None or df.empty:
        return {"structure_quality_score": 0, "structure_confidence": "FAIL", "notes": "No pLDDT values parsed."}
    vals = df["plddt"].astype(float)
    mean = vals.mean()
    pct70 = 100 * (vals >= 70).sum() / len(vals)
    pct90 = 100 * (vals >= 90).sum() / len(vals)
    score = 0
    score += 45 if mean >= 90 else 35 if mean >= 70 else 20 if mean >= 50 else 5
    score += 35 if pct70 >= 90 else 25 if pct70 >= 70 else 15 if pct70 >= 50 else 0
    score += 20 if pct90 >= 50 else 10 if pct90 >= 25 else 0
    score = min(100, round(score, 2))
    return {"structure_quality_score": score, "structure_confidence": confidence_from_score(score), "mean_plddt": round(mean, 2), "pct_residues_over_70": round(pct70, 2), "pct_residues_over_90": round(pct90, 2), "notes": "Structure-quality estimate from PDB B-factors/pLDDT."}

def parse_pae_summary(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {"pae_available": False, "notes": "No PAE JSON found."}
    try:
        obj = json.loads(p.read_text())
    except Exception as e:
        return {"pae_available": False, "notes": f"Unreadable PAE JSON: {e}"}
    pae = None
    if isinstance(obj, dict):
        pae = obj.get("predicted_aligned_error") or obj.get("pae")
    elif isinstance(obj, list) and obj:
        pae = obj[0].get("predicted_aligned_error") if isinstance(obj[0], dict) else obj
    vals = []
    if pae:
        for row in pae:
            for x in row:
                try: vals.append(float(x))
                except Exception: pass
    if not vals:
        return {"pae_available": False, "notes": "No numeric PAE matrix parsed."}
    return {"pae_available": True, "mean_pae": round(sum(vals)/len(vals), 2), "notes": "Lower PAE is better; inspect pocket/domain-specific PAE manually."}

def parse_molprobity(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: MolProbity ###
    Optional manual CSV columns: clashscore,ramachandran_outliers_pct,rotamer_outliers_pct,molprobity_score
    """
    df = safe_read_csv(path)
    return {"molprobity_available": False} if df.empty else {"molprobity_available": True, **df.iloc[0].to_dict()}

def parse_qmean(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: QMEAN / QMEANDisCo ###
    Optional manual CSV columns: qmean_score,qmean_zscore,local_quality_notes
    """
    df = safe_read_csv(path)
    return {"qmean_available": False} if df.empty else {"qmean_available": True, **df.iloc[0].to_dict()}
