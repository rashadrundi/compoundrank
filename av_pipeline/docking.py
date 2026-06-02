import pandas as pd
from .utils import safe_read_csv, confidence_from_score

def parse_gnina(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: GNINA ###
    Expected columns:
    ligand_id,name,pocket_id,cnn_score,cnn_affinity,affinity,pose_file,contacts,notes
    """
    df = safe_read_csv(path)
    for c in ["cnn_score", "cnn_affinity", "affinity"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def score_docking(df):
    if df is None or df.empty:
        return pd.DataFrame()
    rows = []
    for _, r in df.iterrows():
        score, notes = 0, []
        if pd.notna(r.get("cnn_score", None)):
            cnn = float(r.get("cnn_score")); score += min(40, cnn*40 if cnn <= 1 else cnn); notes.append("CNN score parsed.")
        aff = None
        for col in ["cnn_affinity", "affinity"]:
            if pd.notna(r.get(col, None)):
                aff = float(r.get(col)); break
        if aff is not None:
            score += 30 if aff <= -9 else 22 if aff <= -7 else 12 if aff <= -5 else 5
            notes.append("Docking affinity parsed.")
        contacts = str(r.get("contacts","")).lower()
        if any(k in contacts for k in ["hbond","hydrogen","salt","pi","hydrophobic","metal"]):
            score += 15; notes.append("Contact annotations present.")
        if str(r.get("pose_file","")).strip(): score += 10
        if str(r.get("pocket_id","")).strip(): score += 5
        score = min(100, score)
        row = r.to_dict()
        row.update({"docking_score": score, "docking_confidence": confidence_from_score(score), "docking_notes": " ".join(notes)})
        rows.append(row)
    return pd.DataFrame(rows).sort_values("docking_score", ascending=False)
