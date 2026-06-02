from pathlib import Path
import shutil, subprocess, json
from datetime import datetime
import pandas as pd

def now_iso():
    return datetime.utcnow().replace(microsecond=0).isoformat() + "Z"

def ensure_dir(path):
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

def ensure_parent(path):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p

def safe_read_csv(path, **kwargs):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(p, **kwargs)

def command_exists(cmd):
    return shutil.which(cmd) is not None

def run_command(cmd, cwd=None):
    print("[RUN]", " ".join(cmd))
    return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)

def write_json(path, obj):
    p = ensure_parent(path)
    p.write_text(json.dumps(obj, indent=2), encoding="utf-8")

def confidence_from_score(score, high=75, moderate=50, low=25):
    try:
        score = float(score)
    except Exception:
        return "FAIL"
    if score >= high:
        return "HIGH"
    if score >= moderate:
        return "MODERATE"
    if score >= low:
        return "LOW"
    return "FAIL"

def df_records(df, max_rows=25):
    if df is None or df.empty:
        return []
    return df.head(max_rows).fillna("").to_dict(orient="records")
