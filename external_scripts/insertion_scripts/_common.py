from pathlib import Path
import pandas as pd, shutil

def ensure_parent(path):
    p=Path(path); p.parent.mkdir(parents=True, exist_ok=True); return p

def copy_file(src,dst):
    src=Path(src); dst=ensure_parent(dst)
    if not src.exists(): raise FileNotFoundError(src)
    shutil.copy2(src,dst); print(f'[COPY] {src} -> {dst}')

def read_table_auto(path):
    p=Path(path)
    if not p.exists(): raise FileNotFoundError(p)
    if p.suffix.lower() in ['.tsv','.tab']:
        return pd.read_csv(p, sep='\t')
    if p.suffix.lower()=='.csv':
        return pd.read_csv(p)
    try: return pd.read_csv(p, sep='\t')
    except Exception: return pd.read_csv(p)

def write_csv(df,path,columns=None):
    p=ensure_parent(path)
    if columns:
        for c in columns:
            if c not in df.columns: df[c]=''
        df=df[columns]
    df.to_csv(p,index=False); print(f'[WRITE] {p} ({len(df)} rows)')
