import csv
import os
import subprocess
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parent


def _nonempty(path):
    path = Path(path)
    return path.exists() and path.stat().st_size > 0


def _run(cmd, cwd=None):
    print("[VOGDB] Running:")
    print(" ".join(str(x) for x in cmd))

    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.stdout:
        print("[VOGDB] STDOUT:")
        print(result.stdout)

    if result.stderr:
        print("[VOGDB] STDERR:")
        print(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            "[VOGDB] Command failed\n"
            f"Exit code: {result.returncode}\n\n"
            f"Command:\n{' '.join(str(x) for x in cmd)}\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}\n"
        )

    return result


def _default_vogdb_data_dir():
    return Path(
        os.getenv(
            "VOGDB_DATA_ROOT",
            "/opt/compoundrank/data/vogdb_data",
        )
    ).resolve()


def parse_vogdb_tblout(tblout_path, out_csv, max_hits=10):
    tblout_path = Path(tblout_path)
    out_csv = Path(out_csv)

    rows = []
    per_query_counts = {}

    with open(tblout_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#"):
                continue

            parts = line.split(maxsplit=18)

            if len(parts) < 18:
                continue

            target_name = parts[0]
            target_accession = parts[1]
            query_name = parts[2]
            query_accession = parts[3]
            full_evalue = parts[4]
            full_score = parts[5]
            full_bias = parts[6]
            best_domain_evalue = parts[7]
            best_domain_score = parts[8]
            best_domain_bias = parts[9]
            description = parts[18] if len(parts) >= 19 else ""

            current_count = per_query_counts.get(query_name, 0)

            if current_count >= int(max_hits):
                continue

            per_query_counts[query_name] = current_count + 1

            notes = (
                f"query={query_name}; query_accession={query_accession}; "
                f"full_bias={full_bias}; "
                f"best_domain_evalue={best_domain_evalue}; "
                f"best_domain_score={best_domain_score}; "
                f"best_domain_bias={best_domain_bias}; "
                f"description={description}"
            )

            rows.append(
                {
                    "source": "VOGDB",
                    "name": target_name.replace(",", ";"),
                    "accession": target_accession.replace(",", ";"),
                    "start": "",
                    "end": "",
                    "evalue": full_evalue,
                    "score": full_score,
                    "notes": notes.replace(",", ";"),
                }
            )

    fieldnames = [
        "source",
        "name",
        "accession",
        "start",
        "end",
        "evalue",
        "score",
        "notes",
    ]

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)

    print(f"[VOGDB] Parsed {len(rows)} hit(s) to: {out_csv}")

    return out_csv


def run_vogdb_tool(
    fasta_path,
    output_dir,
    hmm_db=None,
    vogdb_data_dir=None,
    force=False,
    cpu=2,
    max_hits=10,
):
    """
    Runs hmmscan against VOGDB.

    Expected FastAPI job layout:
    output_dir = /opt/compoundrank/jobs/<job_id>/homologs

    Writes:
    output_dir/vogdb_hits.csv
    output_dir/vogdb_raw/
    """

    fasta_path = Path(fasta_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_dir = output_dir / "vogdb_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    out_csv = output_dir / "vogdb_hits.csv"

    if _nonempty(out_csv) and not force:
        print(f"[VOGDB] Existing output found. Skipping: {out_csv}")
        return out_csv

    if not fasta_path.exists():
        raise FileNotFoundError(f"[VOGDB] Missing FASTA file: {fasta_path}")

    if vogdb_data_dir is None:
        vogdb_data_dir = _default_vogdb_data_dir()
    else:
        vogdb_data_dir = Path(vogdb_data_dir).resolve()

    if hmm_db is None:
        hmm_db = vogdb_data_dir / "vog224.hmm"
    else:
        hmm_db = Path(hmm_db).resolve()

    if not hmm_db.exists():
        raise FileNotFoundError(
            "[VOGDB] Missing VOGDB HMM database:\n"
            + str(hmm_db)
            + "\n\nExpected default location:\n"
            + str(vogdb_data_dir / "vog224.hmm")
        )

    tblout = raw_dir / "vogdb.tblout"
    domtblout = raw_dir / "vogdb.domtblout"
    hmmscan_txt = raw_dir / "vogdb.hmmscan.txt"

    pressed_files = [
        Path(str(hmm_db) + ".h3f"),
        Path(str(hmm_db) + ".h3i"),
        Path(str(hmm_db) + ".h3m"),
        Path(str(hmm_db) + ".h3p"),
    ]

    if not all(p.exists() for p in pressed_files):
        print(f"[VOGDB] HMMER index files missing. Running hmmpress on: {hmm_db}")
        _run(["hmmpress", str(hmm_db)], cwd=vogdb_data_dir)

    print("[VOGDB] Running hmmscan.")
    print(f"[VOGDB] HMM DB: {hmm_db}")
    print(f"[VOGDB] FASTA: {fasta_path}")
    print(f"[VOGDB] Output CSV: {out_csv}")

    cmd = [
        "hmmscan",
        "--tblout",
        str(tblout),
        "--domtblout",
        str(domtblout),
        "--cpu",
        str(cpu),
        str(hmm_db),
        str(fasta_path),
    ]

    with open(hmmscan_txt, "w", encoding="utf-8") as stdout_file:
        result = subprocess.run(
            cmd,
            cwd=str(output_dir),
            stdout=stdout_file,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    if result.stderr:
        print("[VOGDB] hmmscan STDERR:")
        print(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            "[VOGDB] hmmscan failed\n"
            f"Exit code: {result.returncode}\n\n"
            f"Command:\n{' '.join(cmd)}\n\n"
            f"STDERR:\n{result.stderr}\n"
        )

    if not _nonempty(tblout):
        raise RuntimeError(f"[VOGDB] Expected tblout missing or empty: {tblout}")

    parse_vogdb_tblout(
        tblout_path=tblout,
        out_csv=out_csv,
        max_hits=max_hits,
    )

    if not _nonempty(out_csv):
        raise RuntimeError(f"[VOGDB] Expected output CSV missing or empty: {out_csv}")

    print(f"[VOGDB] VOGDB hits ready: {out_csv}")

    return out_csv