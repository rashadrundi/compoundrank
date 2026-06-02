from pathlib import Path
import yaml
from .utils import confidence_from_score

def load_status(path):
    """
    ### EXTERNAL TOOL REQUIRED HERE: PDBFixer ###
    ### EXTERNAL TOOL REQUIRED HERE: ChimeraX ###
    ### EXTERNAL TOOL REQUIRED HERE: Meeko receptor prep ###
    """
    p = Path(path)
    if not p.exists():
        return {
            "pdbfixer": {"status": "missing"},
            "chimerax": {"status": "missing"},
            "meeko_receptor": {"status": "missing"},
            "docking_box": {"status": "missing"},
        }
    return yaml.safe_load(p.read_text()) or {}

def score_status(status):
    score, notes = 0, []
    for key, label in [("pdbfixer","PDBFixer repair"),("chimerax","ChimeraX inspection"),("meeko_receptor","Meeko receptor prep"),("docking_box","Docking box")]:
        st = str(status.get(key, {}).get("status", "")).lower()
        if st in {"done","complete","completed","pass","passed","ok"}:
            score += 25; notes.append(f"{label} complete.")
        elif st in {"partial","manual","warning"}:
            score += 12; notes.append(f"{label} partial/manual.")
        else:
            notes.append(f"{label} missing.")
    return {"receptor_prep_score": score, "receptor_prep_confidence": confidence_from_score(score), "notes": " ".join(notes)}
