from pathlib import Path
from Bio import SeqIO
import re, hashlib

AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUOJ\-*]+$", re.I)

def read_protein_fasta(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"FASTA not found: {p}")
    records = list(SeqIO.parse(str(p), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {p}")
    out = []
    for rec in records:
        seq = str(rec.seq).replace(" ", "").replace("\n", "").upper()
        if not AA_PATTERN.match(seq):
            raise ValueError(f"Record {rec.id} contains non-protein characters.")
        out.append((rec.id, seq))
    return out

def summarize_fasta(path):
    records = read_protein_fasta(path)
    return {
        "num_records": len(records),
        "record_ids": [x[0] for x in records],
        "lengths": [len(x[1].replace("-", "").replace("*", "")) for x in records],
        "sequence_hashes": {rid: hashlib.sha256(seq.encode()).hexdigest()[:16] for rid, seq in records},
    }
