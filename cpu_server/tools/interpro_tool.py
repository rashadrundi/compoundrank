import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parent


def _nonempty(path):
    path = Path(path)
    return path.exists() and path.stat().st_size > 0


def _run(cmd, cwd=None):
    print("[InterPro] Running:")
    print(" ".join(str(x) for x in cmd))

    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.stdout:
        print("[InterPro] STDOUT:")
        print(result.stdout)

    if result.stderr:
        print("[InterPro] STDERR:")
        print(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            "[InterPro] Command failed\n"
            f"Exit code: {result.returncode}\n\n"
            f"Command:\n{' '.join(str(x) for x in cmd)}\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}\n"
        )

    return result


def _win_to_wsl(path):
    path = Path(path).resolve()
    s = str(path)

    if len(s) >= 3 and s[1] == ":":
        drive = s[0].lower()
        rest = s[2:].replace("\\", "/")
        return f"/mnt/{drive}{rest}"

    return s.replace("\\", "/")


def _default_interpro_data_dir():
    return Path(
        os.getenv(
            "INTERPRO_DATA_ROOT",
            "/opt/compoundrank/data/interpro_data",
        )
    ).resolve()


class InterProDatabasePreflightError(RuntimeError):
    """Raised before Nextflow when required InterPro data are unusable."""


def _parse_applications(applications):
    if applications is None:
        return []

    if isinstance(applications, str):
        raw_values = applications.split(",")
    else:
        raw_values = applications

    normalized = []

    for value in raw_values:
        application = str(value).strip().lower()

        if application and application not in normalized:
            normalized.append(application)

    return normalized


def _matches_application_hmm(application, path):
    name = Path(path).name.lower()
    normalized_name = name.replace("-", "_")

    if application == "ncbifam":
        return normalized_name.startswith("ncbifam") and name.endswith(".hmm")

    if application == "pfam":
        return (
            normalized_name.startswith("pfam_a")
            or normalized_name == "pfam.hmm"
        ) and name.endswith(".hmm")

    return False


def _candidate_version_key(path, datadir):
    path = Path(path)
    datadir = Path(datadir)

    try:
        relative = str(path.relative_to(datadir))
    except ValueError:
        relative = str(path)

    version_numbers = tuple(
        int(value)
        for value in re.findall(r"\d+", relative)
    )

    return (
        version_numbers,
        -len(relative),
        relative.lower(),
    )


def _find_application_hmms(datadir, application):
    datadir = Path(datadir)

    if not datadir.exists():
        return []

    candidates = [
        path.resolve()
        for path in datadir.rglob("*.hmm")
        if _matches_application_hmm(
            application,
            path,
        )
    ]

    candidates = sorted(
        set(candidates),
        key=lambda candidate: _candidate_version_key(
            candidate,
            datadir,
        ),
        reverse=True,
    )

    return candidates


def _limited_text(value, maximum_length=4000):
    text = str(value or "").strip()

    if len(text) <= maximum_length:
        return text

    return text[: maximum_length - 3] + "..."


def _write_preflight_report(report_path, report):
    report_path = Path(report_path)
    report_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def _run_hmmstat(
    hmm_path,
    *,
    mode,
    wsl_distro,
):
    hmm_path = Path(hmm_path).resolve()

    if mode == "local":
        executable = shutil.which("hmmstat")

        if executable is None:
            return {
                "status": "failed",
                "error_type": "executable_missing",
                "error": (
                    "hmmstat was not found on PATH. "
                    "Install HMMER before running InterProScan."
                ),
                "command": ["hmmstat", str(hmm_path)],
                "return_code": None,
                "stdout": "",
                "stderr": "",
            }

        command = [
            executable,
            str(hmm_path),
        ]

    elif mode == "wsl":
        executable = shutil.which("wsl")

        if executable is None:
            return {
                "status": "failed",
                "error_type": "executable_missing",
                "error": (
                    "wsl was not found on PATH; "
                    "the InterPro HMM database could "
                    "not be validated inside WSL."
                ),
                "command": [
                    "wsl",
                    "-d",
                    str(wsl_distro),
                    "--",
                    "hmmstat",
                    _win_to_wsl(hmm_path),
                ],
                "return_code": None,
                "stdout": "",
                "stderr": "",
            }

        command = [
            executable,
            "-d",
            str(wsl_distro),
            "--",
            "hmmstat",
            _win_to_wsl(hmm_path),
        ]

    else:
        return {
            "status": "failed",
            "error_type": "unsupported_mode",
            "error": (
                "Unsupported InterPro preflight "
                f"execution mode: {mode}"
            ),
            "command": [],
            "return_code": None,
            "stdout": "",
            "stderr": "",
        }

    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )

    result = {
        "status": (
            "passed"
            if completed.returncode == 0
            else "failed"
        ),
        "error_type": (
            None
            if completed.returncode == 0
            else "database_invalid"
        ),
        "error": (
            None
            if completed.returncode == 0
            else (
                "hmmstat rejected the HMM database file."
            )
        ),
        "command": [
            str(value)
            for value in command
        ],
        "return_code": completed.returncode,
        "stdout": _limited_text(
            completed.stdout
        ),
        "stderr": _limited_text(
            completed.stderr
        ),
    }

    return result


def preflight_interpro_databases(
    datadir,
    output_dir,
    applications,
    *,
    mode="local",
    wsl_distro="Ubuntu",
):
    """
    Validate required HMM databases before starting Nextflow.

    A machine-readable report is always written to:
        output_dir/interpro_database_preflight.json
    """

    datadir = Path(datadir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_path = (
        output_dir
        / "interpro_database_preflight.json"
    )

    requested_applications = _parse_applications(
        applications
    )

    report = {
        "schema_version": (
            "interpro_database_preflight.v1"
        ),
        "status": "passed",
        "mode": mode,
        "datadir": str(datadir),
        "applications_requested": (
            requested_applications
        ),
        "checks": [],
    }

    failures = []

    if not requested_applications:
        failure = {
            "application": None,
            "status": "failed",
            "error_type": "invalid_configuration",
            "error": (
                "No InterPro applications were requested."
            ),
            "selected_path": None,
            "candidates_found": [],
        }

        report["checks"].append(failure)
        failures.append(failure)

    for application in requested_applications:
        if application not in {
            "pfam",
            "ncbifam",
        }:
            report["checks"].append(
                {
                    "application": application,
                    "status": "not_checked",
                    "error_type": None,
                    "error": (
                        "No HMM-file preflight rule "
                        "is defined for this application."
                    ),
                    "selected_path": None,
                    "candidates_found": [],
                }
            )
            continue

        candidates = _find_application_hmms(
            datadir,
            application,
        )

        candidate_strings = [
            str(candidate)
            for candidate in candidates
        ]

        if not candidates:
            failure = {
                "application": application,
                "status": "failed",
                "error_type": "database_missing",
                "error": (
                    "No matching non-index HMM database "
                    f"file was found for {application} "
                    f"under {datadir}."
                ),
                "selected_path": None,
                "candidates_found": [],
            }

            report["checks"].append(failure)
            failures.append(failure)
            continue

        selected_path = candidates[0]

        check = {
            "application": application,
            "status": "passed",
            "error_type": None,
            "error": None,
            "selected_path": str(selected_path),
            "selected_size_bytes": (
                selected_path.stat().st_size
                if selected_path.exists()
                else 0
            ),
            "candidates_found": candidate_strings,
        }

        if (
            not selected_path.exists()
            or not selected_path.is_file()
            or selected_path.stat().st_size <= 0
        ):
            check.update(
                {
                    "status": "failed",
                    "error_type": "database_missing",
                    "error": (
                        "The selected HMM database file "
                        "does not exist or is empty."
                    ),
                }
            )

            report["checks"].append(check)
            failures.append(check)
            continue

        hmmstat_result = _run_hmmstat(
            selected_path,
            mode=mode,
            wsl_distro=wsl_distro,
        )

        check["hmmstat"] = hmmstat_result

        if hmmstat_result["status"] != "passed":
            check.update(
                {
                    "status": "failed",
                    "error_type": (
                        hmmstat_result.get(
                            "error_type"
                        )
                        or "database_invalid"
                    ),
                    "error": (
                        hmmstat_result.get("error")
                        or "HMM database validation failed."
                    ),
                }
            )

            failures.append(check)

        report["checks"].append(check)

    if failures:
        report["status"] = "failed"

    _write_preflight_report(
        report_path,
        report,
    )

    if failures:
        failure_messages = []

        for failure in failures:
            application = (
                failure.get("application")
                or "configuration"
            )
            error_type = (
                failure.get("error_type")
                or "preflight_failed"
            )
            error = (
                failure.get("error")
                or "Unknown preflight failure."
            )

            hmmstat = failure.get(
                "hmmstat",
                {},
            )

            stderr = (
                hmmstat.get("stderr")
                if isinstance(hmmstat, dict)
                else ""
            )

            detail = (
                f"{application}: "
                f"{error_type}: {error}"
            )

            if stderr:
                detail += (
                    f" hmmstat stderr: {stderr}"
                )

            failure_messages.append(detail)

        raise InterProDatabasePreflightError(
            "[InterPro] Database preflight failed "
            "before Nextflow execution.\n"
            + "\n".join(failure_messages)
            + "\nPreflight report: "
            + str(report_path)
        )

    print(
        "[InterPro] Database preflight passed: "
        f"{report_path}"
    )

    return report


def _copy_first_tsv(nextflow_out, final_tsv):
    nextflow_out = Path(nextflow_out)
    final_tsv = Path(final_tsv)

    candidates = []
    candidates.extend(nextflow_out.rglob("*.tsv"))
    candidates.extend(nextflow_out.rglob("*.tsv.gz"))

    if not candidates:
        print("[InterPro] Files found in output directory:")

        if nextflow_out.exists():
            for p in nextflow_out.rglob("*"):
                print(f"  {p}")

        raise RuntimeError(f"[InterPro] No TSV output found in: {nextflow_out}")

    source_tsv = candidates[0]
    print(f"[InterPro] Found TSV output: {source_tsv}")

    final_tsv.parent.mkdir(parents=True, exist_ok=True)

    if str(source_tsv).lower().endswith(".gz"):
        with gzip.open(source_tsv, "rb") as src, open(final_tsv, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(source_tsv, final_tsv)

    print(f"[InterPro] Copied InterPro TSV to: {final_tsv}")

    return final_tsv


def run_interpro_local(
    fasta_path,
    output_dir,
    force=False,
    nxf_ver="25.04.6",
    revision="6.0.0",
    datadir=None,
    interpro_data_dir=None,
    workdir=None,
    applications="pfam,ncbifam",
    goterms=False,
    pathways=False,
):
    """
    Runs InterProScan 6 through Nextflow/Docker.

    Expected FastAPI job layout:
    output_dir = /opt/compoundrank/jobs/<job_id>/annotation/interpro

    Writes:
    output_dir/interpro.tsv
    output_dir/nextflow_out/
    output_dir/interpro_work/
    """

    fasta_path = Path(fasta_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    nextflow_out = output_dir / "nextflow_out"
    final_tsv = output_dir / "interpro.tsv"

    if _nonempty(final_tsv) and not force:
        print(f"[InterPro] Existing output found. Skipping: {final_tsv}")
        return final_tsv

    if not fasta_path.exists():
        raise FileNotFoundError(f"[InterPro] Missing FASTA file: {fasta_path}")

    if interpro_data_dir is not None and datadir is None:
        datadir = interpro_data_dir

    if datadir is None:
        datadir = _default_interpro_data_dir()
    else:
        datadir = Path(datadir).resolve()

    if workdir is None:
        workdir = output_dir / "interpro_work"
    else:
        workdir = Path(workdir).resolve()

    preflight_interpro_databases(
        datadir=datadir,
        output_dir=output_dir,
        applications=applications,
        mode="local",
    )

    nextflow_out.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    extra_flags = []

    if goterms:
        extra_flags.append("--goterms")

    if pathways:
        extra_flags.append("--pathways")

    extra = " ".join(extra_flags)

    bash_command = f"""
set -e

mkdir -p {shlex.quote(str(nextflow_out))}
mkdir -p {shlex.quote(str(datadir))}
mkdir -p {shlex.quote(str(workdir))}

NXF_VER={shlex.quote(str(nxf_ver))} nextflow run ebi-pf-team/interproscan6 \\
  -r {shlex.quote(str(revision))} \\
  -profile docker \\
  --input {shlex.quote(str(fasta_path))} \\
  --datadir {shlex.quote(str(datadir))} \\
  --interpro latest \\
  --outdir {shlex.quote(str(nextflow_out))} \\
  --applications {shlex.quote(str(applications))} \\
  {extra} \\
  -w {shlex.quote(str(workdir))}
"""

    print("[InterPro] Running local Nextflow InterProScan.")
    print(f"[InterPro] FASTA: {fasta_path}")
    print(f"[InterPro] Data dir: {datadir}")
    print(f"[InterPro] Output dir: {output_dir}")
    print(f"[InterPro] Work dir: {workdir}")

    _run(["bash", "-lc", bash_command], cwd=MODULE_DIR)

    return _copy_first_tsv(nextflow_out, final_tsv)


def run_interpro_wsl(
    fasta_path,
    output_dir,
    force=False,
    nxf_ver="25.04.6",
    revision="6.0.0",
    datadir=None,
    interpro_data_dir=None,
    workdir=None,
    applications="pfam,ncbifam,superfamily",
    goterms=False,
    pathways=False,
    wsl_distro="Ubuntu",
):
    fasta_path = Path(fasta_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    nextflow_out = output_dir / "nextflow_out"
    final_tsv = output_dir / "interpro.tsv"

    if _nonempty(final_tsv) and not force:
        print(f"[InterPro] Existing output found. Skipping: {final_tsv}")
        return final_tsv

    if not fasta_path.exists():
        raise FileNotFoundError(f"[InterPro] Missing FASTA file: {fasta_path}")

    if interpro_data_dir is not None and datadir is None:
        datadir = interpro_data_dir

    if datadir is None:
        datadir = _default_interpro_data_dir()
    else:
        datadir = Path(datadir).resolve()

    if workdir is None:
        workdir = output_dir / "interpro_work"
    else:
        workdir = Path(workdir).resolve()

    preflight_interpro_databases(
        datadir=datadir,
        output_dir=output_dir,
        applications=applications,
        mode="wsl",
        wsl_distro=wsl_distro,
    )

    nextflow_out.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    fasta_wsl = _win_to_wsl(fasta_path)
    nextflow_out_wsl = _win_to_wsl(nextflow_out)
    datadir_wsl = _win_to_wsl(datadir)
    workdir_wsl = _win_to_wsl(workdir)
    module_dir_wsl = _win_to_wsl(MODULE_DIR)

    extra_flags = []

    if goterms:
        extra_flags.append("--goterms")

    if pathways:
        extra_flags.append("--pathways")

    extra = " ".join(extra_flags)

    bash_command = f"""
set -e

cd {shlex.quote(module_dir_wsl)}

mkdir -p {shlex.quote(nextflow_out_wsl)}
mkdir -p {shlex.quote(datadir_wsl)}
mkdir -p {shlex.quote(workdir_wsl)}

NXF_VER={shlex.quote(str(nxf_ver))} nextflow run ebi-pf-team/interproscan6 \\
  -r {shlex.quote(str(revision))} \\
  -profile docker \\
  --input {shlex.quote(fasta_wsl)} \\
  --datadir {shlex.quote(datadir_wsl)} \\
  --interpro latest \\
  --outdir {shlex.quote(nextflow_out_wsl)} \\
  --applications {shlex.quote(str(applications))} \\
  {extra} \\
  -w {shlex.quote(workdir_wsl)}
"""

    print("[InterPro] Running WSL Nextflow InterProScan.")
    print(f"[InterPro] FASTA: {fasta_wsl}")
    print(f"[InterPro] Data dir: {datadir_wsl}")
    print(f"[InterPro] Output dir: {nextflow_out_wsl}")
    print(f"[InterPro] Work dir: {workdir_wsl}")

    _run(
        ["wsl", "-d", wsl_distro, "--", "bash", "-lc", bash_command],
        cwd=MODULE_DIR,
    )

    return _copy_first_tsv(nextflow_out, final_tsv)


def run_interpro_tool(
    fasta_path,
    output_dir,
    mode="local",
    force=False,
    nxf_ver="25.04.6",
    revision="6.0.0",
    datadir=None,
    interpro_data_dir=None,
    workdir=None,
    applications="pfam,ncbifam",
    goterms=False,
    pathways=False,
    wsl_distro="Ubuntu",
):
    if interpro_data_dir is not None and datadir is None:
        datadir = interpro_data_dir

    if datadir is None:
        datadir = _default_interpro_data_dir()

    if mode == "local":
        return run_interpro_local(
            fasta_path=fasta_path,
            output_dir=output_dir,
            force=force,
            nxf_ver=nxf_ver,
            revision=revision,
            datadir=datadir,
            workdir=workdir,
            applications=applications,
            goterms=goterms,
            pathways=pathways,
        )

    if mode == "wsl":
        return run_interpro_wsl(
            fasta_path=fasta_path,
            output_dir=output_dir,
            force=force,
            nxf_ver=nxf_ver,
            revision=revision,
            datadir=datadir,
            workdir=workdir,
            applications=applications,
            goterms=goterms,
            pathways=pathways,
            wsl_distro=wsl_distro,
        )

    raise ValueError(f"[InterPro] Unsupported mode: {mode}")