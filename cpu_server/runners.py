import contextlib
import csv
import io
import json
import os
import shutil
import subprocess
import traceback
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional
from uuid import uuid4

from .settings import (
    DATA_ROOT,
    DEPLOY_ROOT,
    INTERPRO_DATA_ROOT,
    JOBS_ROOT,
    REPO_ROOT,
    VOGDB_DATA_ROOT,
    apply_runtime_env,
    ensure_runtime_dirs,
    runtime_env,
)

from .tools.cdd_tool import run_cdd_tool
from .tools.interpro_tool import run_interpro_tool
from .tools.vogdb_tool import run_vogdb_tool


def create_job_dir() -> tuple[str, Path]:
    """
    Creates a new API job folder outside the repo.

    Expected droplet output:
    /opt/compoundrank/jobs/<job_id>/
    """

    ensure_runtime_dirs()

    job_id = str(uuid4())
    job_dir = JOBS_ROOT / job_id

    (job_dir / "input").mkdir(parents=True, exist_ok=True)
    (job_dir / "annotation" / "cdd").mkdir(parents=True, exist_ok=True)
    (job_dir / "annotation" / "interpro").mkdir(parents=True, exist_ok=True)
    (job_dir / "homologs").mkdir(parents=True, exist_ok=True)
    (job_dir / "logs").mkdir(parents=True, exist_ok=True)

    return job_id, job_dir


def save_uploaded_fasta(upload_file, job_dir: Path) -> Path:
    fasta_path = job_dir / "input" / "protein.fasta"

    with fasta_path.open("wb") as out:
        shutil.copyfileobj(upload_file.file, out)

    return fasta_path


def run_command(
    command: list[str],
    cwd: Optional[Path] = None,
    timeout: int = 7200,
) -> tuple[bool, str]:
    """
    Keeps subprocess support for any code that still relies on command-based tools.

    Important:
    New CDD/InterPro/VOGDB wrapper functions below do not use subprocess.
    They call run_*_tool() directly through run_python_tool().
    """

    try:
        env = {
            **os.environ,
            **runtime_env(),
        }

        completed = subprocess.run(
            command,
            cwd=str(cwd or REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
            env=env,
        )

        output = f"COMMAND: {' '.join(command)}\n"
        output += f"CWD: {str(cwd or REPO_ROOT)}\n\n"

        if completed.stdout:
            output += completed.stdout

        if completed.stderr:
            output += "\nSTDERR:\n" + completed.stderr

        return completed.returncode == 0, output

    except subprocess.TimeoutExpired as exc:
        return False, f"Command timed out after {timeout} seconds: {exc}"

    except Exception as exc:
        return False, str(exc)


def run_python_tool(
    command: list[str],
    tool_func: Callable[..., Any],
    log_path: Path,
    **kwargs,
) -> tuple[bool, str]:
    """
    Runs one of your Python tool wrappers directly.

    This replaces:
        ok, logs = run_command(command)

    with:
        ok, logs = run_python_tool(command, run_cdd_tool, ...)

    Why keep command?
    Because logs/manifests/debugging may still expect a command-like item.
    Here it is a descriptive label, not something executed by subprocess.
    """

    apply_runtime_env()

    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()

    try:
        with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
            result = tool_func(**kwargs)

        logs = ""
        logs += f"COMMAND_LABEL: {' '.join(command)}\n"
        logs += f"REPO_ROOT: {REPO_ROOT}\n"
        logs += f"DEPLOY_ROOT: {DEPLOY_ROOT}\n"
        logs += f"DATA_ROOT: {DATA_ROOT}\n"
        logs += f"JOBS_ROOT: {JOBS_ROOT}\n"
        logs += f"INTERPRO_DATA_ROOT: {INTERPRO_DATA_ROOT}\n"
        logs += f"VOGDB_DATA_ROOT: {VOGDB_DATA_ROOT}\n\n"

        stdout = stdout_buffer.getvalue()
        stderr = stderr_buffer.getvalue()

        if stdout:
            logs += "STDOUT:\n" + stdout + "\n"

        if stderr:
            logs += "STDERR:\n" + stderr + "\n"

        if result is not None:
            logs += f"RETURN_VALUE: {repr(result)}\n"

        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(logs, encoding="utf-8")

        return True, logs

    except Exception:
        logs = ""
        logs += f"COMMAND_LABEL: {' '.join(command)}\n"
        logs += f"REPO_ROOT: {REPO_ROOT}\n"
        logs += f"DEPLOY_ROOT: {DEPLOY_ROOT}\n"
        logs += f"DATA_ROOT: {DATA_ROOT}\n"
        logs += f"JOBS_ROOT: {JOBS_ROOT}\n"
        logs += f"INTERPRO_DATA_ROOT: {INTERPRO_DATA_ROOT}\n"
        logs += f"VOGDB_DATA_ROOT: {VOGDB_DATA_ROOT}\n\n"
        logs += "EXCEPTION:\n"
        logs += traceback.format_exc()

        stdout = stdout_buffer.getvalue()
        stderr = stderr_buffer.getvalue()

        if stdout:
            logs += "\nSTDOUT BEFORE ERROR:\n" + stdout

        if stderr:
            logs += "\nSTDERR BEFORE ERROR:\n" + stderr

        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(logs, encoding="utf-8")

        return False, logs


INTERPRO_TSV_COLUMNS = [
    "protein_accession",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "end",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_annotations",
    "pathways",
]


def read_table(path: Path, max_rows: int = 200) -> List[Dict[str, Any]]:
    if not path.exists():
        return []

    rows: List[Dict[str, Any]] = []

    # Special case: InterProScan TSV usually has no header row.
    if path.name == "interpro.tsv":
        with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            for idx, line in enumerate(f):
                if idx >= max_rows:
                    break

                line = line.rstrip("\n")

                if not line:
                    continue

                parts = line.split("\t")

                row = {}

                for col_idx, col_name in enumerate(INTERPRO_TSV_COLUMNS):
                    row[col_name] = parts[col_idx] if col_idx < len(parts) else ""

                # Preserve any extra fields just in case InterPro adds columns.
                if len(parts) > len(INTERPRO_TSV_COLUMNS):
                    row["extra_fields"] = parts[len(INTERPRO_TSV_COLUMNS):]

                rows.append(row)

        return rows

    delimiter = "\t" if path.suffix.lower() in [".tsv", ".gff3"] else ","

    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)

        for idx, row in enumerate(reader):
            if idx >= max_rows:
                break

            rows.append(dict(row))

    return rows


def find_existing_output(candidates: list[Path]) -> Optional[Path]:
    for path in candidates:
        if path.exists():
            return path

    return None


def tool_response(
    status: str,
    output_file: Optional[Path] = None,
    error: Optional[str] = None,
    command: Optional[list[str]] = None,
    log_file: Optional[Path] = None,
) -> Dict[str, Any]:
    rows = read_table(output_file) if output_file else []

    return {
        "status": status,
        "command": command or [],
        "output_file": str(output_file) if output_file else None,
        "log_file": str(log_file) if log_file else None,
        "row_count_returned": len(rows),
        "rows": rows,
        "error": error,
    }

def validate_fasta(fasta_path: Path) -> None:
    text = fasta_path.read_text(encoding="utf-8", errors="replace").strip()

    if not text:
        raise ValueError("Uploaded FASTA file is empty.")

    if not text.startswith(">"):
        raise ValueError("Uploaded file does not look like FASTA. It must start with '>'.")

    sequence_lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ]

    sequence = "".join(sequence_lines)

    if not sequence:
        raise ValueError("FASTA contains a header but no sequence.")

    allowed = set("ABCDEFGHIKLMNPQRSTVWXYZ*-.abcdefghiklmnpqrstvwxyz")
    invalid = sorted(set(sequence) - allowed)

    if invalid:
        raise ValueError(f"FASTA contains invalid characters: {invalid}")


def run_cdd(job_dir: Path, fasta_path: Path) -> Dict[str, Any]:
    """
    CPU step:
    FASTA → CDD results

    Expected output:
    /opt/compoundrank/jobs/<job_id>/annotation/cdd/cdd_results.csv
    """

    output_dir = job_dir / "annotation" / "cdd"
    expected_output = output_dir / "cdd_results.csv"
    log_file = output_dir / "cdd.log"

    command = [
        "python-call",
        "run_cdd_tool",
        "--fasta",
        str(fasta_path),
        "--outdir",
        str(output_dir),
    ]

    ok, logs = run_python_tool(
        command=command,
        tool_func=run_cdd_tool,
        log_path=log_file,
        fasta_path=fasta_path,
        output_dir=output_dir,
    )

    output_file = find_existing_output(
        [
            expected_output,
            output_dir / "cdd.csv",
            output_dir / "cdd_results.tsv",
            output_dir / "cdd.tsv",
        ]
    )

    if not ok:
        return tool_response(
            status="failed",
            output_file=output_file or expected_output,
            error=logs,
            command=command,
            log_file=log_file,
        )

    if output_file is None:
        return tool_response(
            status="failed",
            output_file=expected_output,
            error="CDD tool completed but no expected output file was found.",
            command=command,
            log_file=log_file,
        )

    return tool_response(
        status="complete",
        output_file=output_file,
        command=command,
        log_file=log_file,
    )


def run_interpro(job_dir: Path, fasta_path: Path) -> Dict[str, Any]:
    """
    CPU step:
    FASTA → InterPro results

    Expected output:
    /opt/compoundrank/jobs/<job_id>/annotation/interpro/interpro.tsv
    """

    output_dir = job_dir / "annotation" / "interpro"
    expected_output = output_dir / "interpro.tsv"
    log_file = output_dir / "interpro.log"

    command = [
        "python-call",
        "run_interpro_tool",
        "--fasta",
        str(fasta_path),
        "--outdir",
        str(output_dir),
        "--interpro-data",
        str(INTERPRO_DATA_ROOT),
    ]

    ok, logs = run_python_tool(
        command=command,
        tool_func=run_interpro_tool,
        log_path=log_file,
        fasta_path=fasta_path,
        output_dir=output_dir,
        interpro_data_dir=INTERPRO_DATA_ROOT,
    )

    output_file = find_existing_output(
        [
            expected_output,
            output_dir / "protein.fasta.tsv",
            output_dir / "protein.tsv",
            output_dir / "interproscan.tsv",
            output_dir / "nextflow_out" / "protein.fasta.tsv",
        ]
    )

    if not ok:
        return tool_response(
            status="failed",
            output_file=output_file or expected_output,
            error=logs,
            command=command,
            log_file=log_file,
        )

    if output_file is None:
        return tool_response(
            status="failed",
            output_file=expected_output,
            error="InterPro tool completed but no expected output file was found.",
            command=command,
            log_file=log_file,
        )

    return tool_response(
        status="complete",
        output_file=output_file,
        command=command,
        log_file=log_file,
    )


def run_vogdb(job_dir: Path, fasta_path: Path) -> Dict[str, Any]:
    """
    CPU step:
    FASTA → VOGDB results

    Expected output:
    /opt/compoundrank/jobs/<job_id>/homologs/vogdb_hits.csv
    """

    output_dir = job_dir / "homologs"
    expected_output = output_dir / "vogdb_hits.csv"
    log_file = output_dir / "vogdb.log"

    command = [
        "python-call",
        "run_vogdb_tool",
        "--fasta",
        str(fasta_path),
        "--outdir",
        str(output_dir),
        "--vogdb-data",
        str(VOGDB_DATA_ROOT),
    ]

    ok, logs = run_python_tool(
        command=command,
        tool_func=run_vogdb_tool,
        log_path=log_file,
        fasta_path=fasta_path,
        output_dir=output_dir,
        vogdb_data_dir=VOGDB_DATA_ROOT,
    )

    output_file = find_existing_output(
        [
            expected_output,
            output_dir / "vogdb.csv",
            output_dir / "vog_hits.csv",
            output_dir / "vogdb_hits.tsv",
            output_dir / "vog_hits.tsv",
        ]
    )

    if not ok:
        return tool_response(
            status="failed",
            output_file=output_file or expected_output,
            error=logs,
            command=command,
            log_file=log_file,
        )

    if output_file is None:
        return tool_response(
            status="failed",
            output_file=expected_output,
            error="VOGDB tool completed but no expected output file was found.",
            command=command,
            log_file=log_file,
        )

    return tool_response(
        status="complete",
        output_file=output_file,
        command=command,
        log_file=log_file,
    )


def run_cpu_analysis(upload_file) -> Dict[str, Any]:
    """
    Main FastAPI entrypoint.

    Upload FASTA
    → validate FASTA
    → run CDD
    → run InterPro
    → run VOGDB
    → write CPU manifest and local-workflow handoff
    → return JSON response
    """

    apply_runtime_env()
    ensure_runtime_dirs()

    job_id, job_dir = create_job_dir()
    fasta_path = save_uploaded_fasta(upload_file, job_dir)

    validate_fasta(fasta_path)

    cdd_result = run_cdd(job_dir, fasta_path)
    interpro_result = run_interpro(job_dir, fasta_path)
    vogdb_result = run_vogdb(job_dir, fasta_path)

    tool_results = [
        cdd_result,
        interpro_result,
        vogdb_result,
    ]

    manifest_path = job_dir / "cpu_manifest.json"
    handoff_path = job_dir / "cpu_to_local_handoff.json"

    file_map = {
        "input_fasta": str(fasta_path),
        "manifest": str(manifest_path),
        "local_handoff_json": str(handoff_path),

        "cdd_csv": cdd_result.get("output_file"),
        "cdd_log": cdd_result.get("log_file"),

        "interpro_tsv": interpro_result.get("output_file"),
        "interpro_log": interpro_result.get("log_file"),
        "interpro_gff3": str(
            job_dir
            / "annotation"
            / "interpro"
            / "nextflow_out"
            / "protein.fasta.gff3"
        ),
        "interpro_json": str(
            job_dir
            / "annotation"
            / "interpro"
            / "nextflow_out"
            / "protein.fasta.json"
        ),
        "interpro_jsonl": str(
            job_dir
            / "annotation"
            / "interpro"
            / "nextflow_out"
            / "protein.fasta.jsonl"
        ),
        "interpro_xml": str(
            job_dir
            / "annotation"
            / "interpro"
            / "nextflow_out"
            / "protein.fasta.xml"
        ),

        "vogdb_csv": vogdb_result.get("output_file"),
        "vogdb_log": vogdb_result.get("log_file"),
        "vogdb_tblout": str(
            job_dir
            / "homologs"
            / "vogdb_raw"
            / "vogdb.tblout"
        ),
        "vogdb_domtblout": str(
            job_dir
            / "homologs"
            / "vogdb_raw"
            / "vogdb.domtblout"
        ),
        "vogdb_hmmscan_txt": str(
            job_dir
            / "homologs"
            / "vogdb_raw"
            / "vogdb.hmmscan.txt"
        ),
    }

    result_counts = {
        "cdd": len(cdd_result.get("rows", [])),
        "interpro": len(interpro_result.get("rows", [])),
        "vogdb": len(vogdb_result.get("rows", [])),
    }

    response = {
        "job_id": job_id,
        "status": (
            "complete"
            if all(
                result["status"] == "complete"
                for result in tool_results
            )
            else "partial_or_failed"
        ),
        "paths": {
            "repo_root": str(REPO_ROOT),
            "deploy_root": str(DEPLOY_ROOT),
            "data_root": str(DATA_ROOT),
            "jobs_root": str(JOBS_ROOT),
            "interpro_data_root": str(INTERPRO_DATA_ROOT),
            "vogdb_data_root": str(VOGDB_DATA_ROOT),
        },
        "fasta_file": str(fasta_path),
        "job_dir": str(job_dir),
        "manifest_file": str(manifest_path),

        "cdd": cdd_result,
        "interpro": interpro_result,
        "vogdb": vogdb_result,

        "results": {
            "cdd": cdd_result.get("rows", []),
            "interpro": interpro_result.get("rows", []),
            "vogdb": vogdb_result.get("rows", []),
        },

        "result_counts": result_counts,
        "files": file_map,
    }

    handoff = {
        "job_id": job_id,
        "status": response["status"],
        "purpose": "local_gpu_workflow_input",
        "requested_local_steps": [
            "colabfold",
            "fpocket",
            "gnina",
        ],
        "input": {
            "fasta_file": str(fasta_path),
        },
        "annotations": response["results"],
        "result_counts": result_counts,
        "cpu_files": file_map,
        "notes": {
            "fpocket_location": "local_gpu_workstation",
            "fpocket_reason": (
                "fpocket requires the predicted or supplied PDB "
                "structure and should run after ColabFold and "
                "before GNINA."
            ),
            "cpu_server_steps_completed": [
                "cdd",
                "interpro",
                "vogdb",
            ],
            "gpu_local_steps_pending": [
                "colabfold",
                "fpocket",
                "gnina",
            ],
        },
    }

    response["local_handoff"] = handoff

    handoff_path.write_text(
        json.dumps(handoff, indent=2),
        encoding="utf-8",
    )

    manifest_path.write_text(
        json.dumps(response, indent=2),
        encoding="utf-8",
    )

    return response
