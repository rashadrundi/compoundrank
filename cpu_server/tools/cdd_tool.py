import csv
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


BASE_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"


def _nonempty(path):
    path = Path(path)
    return path.exists() and path.stat().st_size > 0


def read_text(path):
    with open(path, "r", encoding="utf-8") as f:
        return f.read().strip() + "\n\n"


def fetch_get(params):
    url = BASE_URL + "?" + urllib.parse.urlencode(params)

    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "compoundrank-cdd-tool/0.1 educational-use"
        },
    )

    with urllib.request.urlopen(req, timeout=90) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_post(params):
    data = urllib.parse.urlencode(params).encode("utf-8")

    req = urllib.request.Request(
        BASE_URL,
        data=data,
        headers={
            "User-Agent": "compoundrank-cdd-tool/0.1 educational-use",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )

    with urllib.request.urlopen(req, timeout=90) as response:
        return response.read().decode("utf-8", errors="replace")


def extract_search_id(text):
    patterns = [
        r"#cdsid\s+([A-Za-z0-9_.\-]+)",
        r"Search-ID:\s*([A-Za-z0-9_.\-]+)",
        r"Search ID:\s*([A-Za-z0-9_.\-]+)",
        r"cdsid=([A-Za-z0-9_.\-]+)",
        r'name=["\']cdsid["\']\s+value=["\']([^"\']+)["\']',
    ]

    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1)

    return None


def looks_like_hit_table(text):
    header = "Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort name\tSuperfamily"
    return header in text


def parse_cdd_status(text):
    match = re.search(r"#status\s+(\d+)", text)

    if match:
        return int(match.group(1))

    lower = text.lower()

    if "job is still running" in lower:
        return 3

    if "query\thit type\tpssm-id" in lower:
        return 0

    return None


def parse_hits(text):
    rows = []

    for line in text.splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        if line.startswith("Query\tHit type"):
            continue

        parts = line.split("\t")

        if len(parts) < 10:
            continue

        query = parts[0]
        hit_type = parts[1]
        pssm_id = parts[2]
        start = parts[3]
        end = parts[4]
        evalue = parts[5]
        bitscore = parts[6]
        accession = parts[7]
        short_name = parts[8]
        superfamily = parts[9]

        notes = (
            f"hit_type={hit_type}; pssm_id={pssm_id}; "
            f"query={query}; superfamily={superfamily}"
        )

        rows.append(
            {
                "source": "CDD",
                "name": short_name.replace(",", ";"),
                "accession": accession.replace(",", ";"),
                "start": start,
                "end": end,
                "evalue": evalue,
                "score": bitscore,
                "notes": notes.replace(",", ";"),
            }
        )

    return rows


def write_cdd_csv(rows, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

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

    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def submit_cdd_job(fasta_text):
    params = {
        "queries": fasta_text,
        "db": "cdd",
        "smode": "auto",
        "tdata": "hits",
        "maxhit": "500",
    }

    print("[CDD] Submitting FASTA to NCBI Batch CD-Search...")

    try:
        return fetch_post(params)
    except urllib.error.HTTPError:
        print("[CDD] POST failed. Trying GET instead...")
        return fetch_get(params)


def retrieve_hits(search_id, cdd_dir, max_attempts=30, sleep_seconds=70):
    print(f"[CDD] Search ID: {search_id}")
    print("[CDD] Waiting before first result check...")

    time.sleep(sleep_seconds)

    params = {
        "cdsid": search_id,
        "tdata": "hits",
    }

    for attempt in range(1, max_attempts + 1):
        print(f"[CDD] Checking results attempt {attempt}/{max_attempts}")

        text = fetch_get(params)
        status = parse_cdd_status(text)

        print(f"[CDD] Current status: {status}")

        if status == 0 or looks_like_hit_table(text):
            print("[CDD] Job finished.")
            return text

        if status == 3:
            print(f"[CDD] Job still running. Waiting {sleep_seconds} seconds...")
            time.sleep(sleep_seconds)
            continue

        debug_path = Path(cdd_dir) / "cdd_unknown_status_debug.txt"

        with open(debug_path, "w", encoding="utf-8") as f:
            f.write(text)

        print(f"[CDD] Unknown status. Saved debug output to: {debug_path}")
        time.sleep(sleep_seconds)

    raise RuntimeError("[CDD] Job did not finish before timeout.")


def run_cdd_tool(
    fasta_path,
    output_dir,
    force=False,
    max_attempts=30,
    sleep_seconds=70,
):
    """
    Runs NCBI Batch CD-Search.

    Expected FastAPI job layout:
    output_dir = /opt/compoundrank/jobs/<job_id>/annotation/cdd

    Writes:
    output_dir/cdd_results.csv
    output_dir/cdd_submit_response.txt
    output_dir/cdd_raw_hits.txt
    """

    fasta_path = Path(fasta_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    out_csv = output_dir / "cdd_results.csv"

    if _nonempty(out_csv) and not force:
        print(f"[CDD] Existing output found. Skipping: {out_csv}")
        return out_csv

    if not fasta_path.exists():
        raise FileNotFoundError(f"[CDD] Missing FASTA file: {fasta_path}")

    fasta_text = read_text(fasta_path)

    initial_text = submit_cdd_job(fasta_text)

    submit_response_path = output_dir / "cdd_submit_response.txt"

    with open(submit_response_path, "w", encoding="utf-8") as f:
        f.write(initial_text)

    if looks_like_hit_table(initial_text):
        result_text = initial_text
    else:
        search_id = extract_search_id(initial_text)

        if not search_id:
            debug_path = output_dir / "cdd_submit_debug.html"

            with open(debug_path, "w", encoding="utf-8") as f:
                f.write(initial_text)

            raise RuntimeError(
                "[CDD] Could not find Search-ID. Debug saved to: "
                + str(debug_path)
            )

        result_text = retrieve_hits(
            search_id=search_id,
            cdd_dir=output_dir,
            max_attempts=max_attempts,
            sleep_seconds=sleep_seconds,
        )

    raw_hits_path = output_dir / "cdd_raw_hits.txt"

    with open(raw_hits_path, "w", encoding="utf-8") as f:
        f.write(result_text)

    rows = parse_hits(result_text)

    if not rows:
        debug_path = output_dir / "cdd_results_debug.txt"

        with open(debug_path, "w", encoding="utf-8") as f:
            f.write(result_text)

        raise RuntimeError(
            "[CDD] Finished, but no hits were parsed. Debug saved to: "
            + str(debug_path)
        )

    write_cdd_csv(rows, out_csv)

    print(f"[CDD] Wrote {len(rows)} hit(s) to: {out_csv}")

    return out_csv