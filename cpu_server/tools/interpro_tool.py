import gzip
import os
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

    nextflow_out.mkdir(parents=True, exist_ok=True)
    datadir.mkdir(parents=True, exist_ok=True)
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

    nextflow_out.mkdir(parents=True, exist_ok=True)
    datadir.mkdir(parents=True, exist_ok=True)
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