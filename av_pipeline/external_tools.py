import os
import subprocess
import sys
import shutil
import gzip
from pathlib import Path


def _nonempty(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def _run(cmd, cwd=None):
    print("[AUTO] Running:")
    print(" ".join(str(x) for x in cmd))
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}")


def _win_to_wsl(path):
    p = Path(path).resolve()
    s = str(p)

    # Windows path like C:\Users\...
    if len(s) >= 3 and s[1] == ":":
        drive = s[0].lower()
        rest = s[2:].replace("\\", "/")
        return f"/mnt/{drive}{rest}"

    # Already Linux-like
    return s.replace("\\", "/")


def _run_cdd(project_root, run_dir, force=False):
    out_csv = Path(run_dir) / "annotation" / "cdd" / "cdd_results.csv"

    if _nonempty(out_csv) and not force:
        print(f"[AUTO:CDD] Existing CDD output found. Skipping: {out_csv}")
        return

    script = Path(project_root) / "external_scripts" / "insertion_scripts" / "run_cdd_auto.py"

    if not script.exists():
        raise FileNotFoundError(f"Missing CDD script: {script}")

    _run([
        sys.executable,
        str(script),
        "--run-dir",
        str(run_dir),
    ], cwd=project_root)


def _run_interpro_wsl(project_root, run_dir, settings):
    out_tsv = Path(run_dir) / "annotation" / "interpro" / "interpro.tsv"
    out_dir_win = Path(run_dir) / "annotation" / "interpro" / "nextflow_out"
    force = bool(settings.get("force", False))

    if _nonempty(out_tsv) and not force:
        print(f"[AUTO:InterPro] Existing InterPro output found. Skipping: {out_tsv}")
        return

    out_dir_win.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    project_wsl = _win_to_wsl(project_root)
    run_wsl = _win_to_wsl(run_dir)
    outdir_wsl = f"{run_wsl}/annotation/interpro/nextflow_out"

    nxf_ver = str(settings.get("nxf_ver", "25.04.6"))
    revision = str(settings.get("revision", "6.0.0"))
    datadir = str(settings.get("datadir", "interpro_data"))
    workdir = str(settings.get("workdir", "interpro_work"))
    applications = str(settings.get("applications", "pfam,ncbifam,superfamily"))
    wsl_distro = str(settings.get("wsl_distro", "Ubuntu"))

    extra_flags = []
    if bool(settings.get("goterms", False)):
        extra_flags.append("--goterms")
    if bool(settings.get("pathways", False)):
        extra_flags.append("--pathways")

    extra = " ".join(extra_flags)

    bash_command = f'''
set -e
cd "{project_wsl}"

mkdir -p "{outdir_wsl}"
mkdir -p "{datadir}"
mkdir -p "{workdir}"
mkdir -p "{run_wsl}/annotation/interpro"

NXF_VER="{nxf_ver}" nextflow run ebi-pf-team/interproscan6 \\
  -r "{revision}" \\
  -profile docker \\
  --input "{run_wsl}/input/protein.fasta" \\
  --datadir "{datadir}" \\
  --interpro latest \\
  --outdir "{outdir_wsl}" \\
  --applications "{applications}" \\
  {extra} \\
  -w "{workdir}"
'''

    _run(["wsl", "-d", wsl_distro, "--", "bash", "-lc", bash_command], cwd=project_root)

    print(f"[AUTO:InterPro] Looking for TSV output in: {out_dir_win}")

    candidates = []
    candidates.extend(out_dir_win.rglob("*.tsv"))
    candidates.extend(out_dir_win.rglob("*.tsv.gz"))

    if not candidates:
        print("[AUTO:InterPro] Files found in output directory:")
        if out_dir_win.exists():
            for p in out_dir_win.rglob("*"):
                print(f"  {p}")
        raise RuntimeError(f"No InterPro TSV output found in: {out_dir_win}")

    source_tsv = candidates[0]
    print(f"[AUTO:InterPro] Found TSV output: {source_tsv}")

    if str(source_tsv).lower().endswith(".gz"):
        with gzip.open(source_tsv, "rb") as src, open(out_tsv, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(source_tsv, out_tsv)

    print(f"[AUTO:InterPro] Copied InterPro TSV to: {out_tsv}")


def _run_interpro_local(project_root, run_dir, settings):
    out_tsv = Path(run_dir) / "annotation" / "interpro" / "nextflow_out" / "protein.fasta.tsv"
    out_dir = Path(run_dir) / "annotation" / "interpro" / "nextflow_out"
    force = bool(settings.get("force", False))

    if _nonempty(out_tsv) and not force:
        print(f"[AUTO:InterPro] Existing InterPro output found. Skipping: {out_tsv}")
        return

    out_dir.mkdir(parents=True, exist_ok=True)

    nxf_ver = str(settings.get("nxf_ver", "25.04.6"))
    revision = str(settings.get("revision", "6.0.0"))
    datadir = str(settings.get("datadir", "interpro_data"))
    workdir = str(settings.get("workdir", "interpro_work"))
    applications = str(settings.get("applications", "pfam,ncbifam,superfamily"))

    extra_flags = []
    if bool(settings.get("goterms", False)):
        extra_flags.append("--goterms")
    if bool(settings.get("pathways", False)):
        extra_flags.append("--pathways")

    extra = " ".join(extra_flags)

    bash_command = f'''
    set -e
    cd "{project_root}"

    mkdir -p "{out_dir}"
    mkdir -p "{datadir}"
    mkdir -p "{workdir}"

    NXF_VER="{nxf_ver}" nextflow run ebi-pf-team/interproscan6 \\
    -r "{revision}" \\
    -profile docker \\
    --input "{Path(run_dir) / "input" / "protein.fasta"}" \\
    --datadir "{datadir}" \\
    --interpro latest \\
    --outdir "{out_dir}" \\
    --applications "{applications}" \\
    {extra} \\
    -w "{workdir}"
    '''

    _run(["bash", "-lc", bash_command], cwd=project_root)

    if not _nonempty(out_tsv):
        print("[AUTO:InterPro] Files found in InterPro output directory:")
        for p in out_dir.rglob("*"):
            print(f"  {p}")
        raise RuntimeError(f"Expected InterPro TSV not found: {out_tsv}")

    print(f"[AUTO:InterPro] InterPro TSV ready: {out_tsv}")

def maybe_run_external_tools(config):
    external = config.get("external_tools", {})

    if not external.get("auto_run", False):
        print("[INFO] External auto-run disabled. Use external_markers/ and insert output files manually.")
        return

    project_root = Path.cwd()
    run_dir = Path(config.get("run_dir", ""))

    if not str(run_dir):
        raise ValueError("Config is missing run_dir.")

    if not run_dir.is_absolute():
        run_dir = project_root / run_dir

    print(f"[AUTO] External auto-run enabled for: {run_dir}")

    cdd_settings = external.get("cdd", {})
    if cdd_settings.get("enabled", False):
        _run_cdd(
            project_root=project_root,
            run_dir=run_dir,
            force=bool(cdd_settings.get("force", False)),
        )

    interpro_settings = external.get("interpro", {})
    if interpro_settings.get("enabled", False):
        mode = interpro_settings.get("mode", "local")

        if mode == "wsl":
            _run_interpro_wsl(
                project_root=project_root,
                run_dir=run_dir,
                settings=interpro_settings,
            )
        elif mode == "local":
            _run_interpro_local(
                project_root=project_root,
                run_dir=run_dir,
                settings=interpro_settings,
            )
        else:
            raise ValueError(f"Unsupported InterPro mode: {mode}")