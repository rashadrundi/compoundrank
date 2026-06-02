import pandas as pd
from .utils import confidence_from_score

def weighted_score(ann, hom, struct, cons, psn, rec, pocket_score, ligand_score, docking_score):
    vals = {
        "annotation": ann.get("annotation_score", 0) or 0,
        "homolog": hom.get("homolog_score", 0) or 0,
        "structure": struct.get("structure_quality_score", 0) or 0,
        "conservation": cons.get("conservation_score", 0) or 0,
        "psn": psn.get("psn_score", 0) or 0,
        "receptor": rec.get("receptor_prep_score", 0) or 0,
        "pocket": pocket_score or 0,
        "ligand": ligand_score or 0,
        "docking": docking_score or 0,
    }
    weights = {"annotation":0.12,"homolog":0.10,"structure":0.14,"conservation":0.12,"psn":0.08,"receptor":0.10,"pocket":0.14,"ligand":0.10,"docking":0.10}
    total = round(sum(float(vals[k]) * weights[k] for k in vals), 2)
    conf = confidence_from_score(total)
    if vals["docking"] == 0 and conf == "HIGH":
        conf = "MODERATE"
    out = {"overall_score": total, "overall_confidence": conf}
    out.update({f"{k}_component": round(float(v),2) for k, v in vals.items()})
    return out

def make_candidates(ann, hom, struct, cons, psn, rec, pockets, ligands, docking):
    rows = []
    if docking is not None and not docking.empty:
        for _, d in docking.head(25).iterrows():
            pid = str(d.get("pocket_id",""))
            lname = d.get("name", d.get("ligand_id",""))
            pscore = 0
            if pockets is not None and not pockets.empty and pid:
                m = pockets[pockets["pocket_id"].astype(str) == pid]
                if not m.empty: pscore = m.iloc[0].get("pocket_priority_score", 0)
            lscore = 0
            if ligands is not None and not ligands.empty:
                lid = str(d.get("ligand_id","")); nm = str(d.get("name",""))
                mask = False
                if "ligand_id" in ligands.columns: mask = (ligands["ligand_id"].astype(str) == lid)
                if "name" in ligands.columns: mask = mask | (ligands["name"].astype(str) == nm) if hasattr(mask, "__or__") else (ligands["name"].astype(str) == nm)
                try:
                    m = ligands[mask]
                    if not m.empty: lscore = m.iloc[0].get("ligand_score", 0)
                except Exception: pass
            base = weighted_score(ann, hom, struct, cons, psn, rec, pscore, lscore, d.get("docking_score",0))
            rows.append({"candidate_type":"docked_pair","pocket_id":pid,"ligand":lname,**base,"notes":"Candidate from GNINA docking table."})
    elif pockets is not None and not pockets.empty and ligands is not None and not ligands.empty:
        for _, p in pockets.head(5).iterrows():
            for _, l in ligands.head(10).iterrows():
                base = weighted_score(ann, hom, struct, cons, psn, rec, p.get("pocket_priority_score",0), l.get("ligand_score",0), 0)
                rows.append({"candidate_type":"pre_docking_pair","pocket_id":p.get("pocket_id",""),"ligand":l.get("name",l.get("ligand_id","")),**base,"notes":"Pre-docking candidate; GNINA score unavailable."})
    return pd.DataFrame(rows).sort_values("overall_score", ascending=False) if rows else pd.DataFrame()
