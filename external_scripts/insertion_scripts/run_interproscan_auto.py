import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


def run_command(cmd, cwd=None):
    print("[InterProScan] Running command:")
    print(" ".join(str(x) for x in cmd))
    result = subprocess.run(cmd, cwd=cwd)

    if result.returncode != 0:
        raise RuntimeError(f"InterProScan command failed with exit code {result.returncode}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument(
        "--mode",
        choices=["auto", "bash", "direct"],
        default="auto",
        help="auto tries direct command first, then bash shell script.",
    )
    parser.add_argument(
        "--interproscan-cmd",
        default=os.environ.get("INTERPROSCAN_CMD", "interproscan.sh"),
        help="Path/name of InterProScan command. Example: /opt/interproscan/interproscan.sh",
    )
    parser.add_argument(
        "--cpu",
        default=os.environ.get("INTERPROSCAN_CPU", "2"),
        help="CPU count for InterProScan.",
    )
    args = parser.parse_args()

    project_root = Path.cwd()
    run_dir = Path(args.run_dir)
    fasta = run_dir / "input" / "protein.fasta"
    out_dir = run_dir / "annotation" / "interpro"
    out_file = out_dir / "interpro.tsv"
    log_dir = run_dir / "logs"

    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    if not fasta.exists():
        raise FileNotFoundError(f"Missing FASTA file: {fasta}")

    # If output already exists, remove it so we know this run generated a fresh file.
    if out_file.exists():
        out_file.unlink()

    if args.mode in ["auto", "direct"]:
        cmd_path = shutil.which(args.interproscan_cmd) or args.interproscan_cmd

        direct_cmd = [
            cmd_path,
            "-i", str(fasta),
            "-f", "TSV",
            "-o", str(out_file),
            "-iprlookup",
            "-goterms",
            "-pa",
            "-cpu", str(args.cpu),
        ]

        if args.mode == "direct":
            run_command(direct_cmd)
        else:
            try:
                run_command(direct_cmd)
            except Exception as e:
                print(f"[InterProScan] Direct mode failed: {e}")
                print("[InterProScan] Falling back to bash mode...")
                args.mode = "bash"

    if args.mode == "bash":
        shell_script = project_root / "external_scripts" / "run_tool_scripts" / "02_run_interproscan.sh"

        if not shell_script.exists():
            raise FileNotFoundError(f"Missing shell script: {shell_script}")

        bash_cmd = shutil.which("bash")

        if not bash_cmd:
            raise RuntimeError(
                "Could not find bash. Use WSL/Git Bash, or run with --mode direct and provide --interproscan-cmd."
            )

        env = os.environ.copy()
        env["INTERPROSCAN_CMD"] = args.interproscan_cmd
        env["INTERPROSCAN_CPU"] = str(args.cpu)

        print("[InterProScan] Running through bash shell script...")
        result = subprocess.run(
            [bash_cmd, str(shell_script), str(run_dir)],
            env=env,
        )

        if result.returncode != 0:
            raise RuntimeError(f"InterProScan bash script failed with exit code {result.returncode}")

    if not out_file.exists() or out_file.stat().st_size == 0:
        raise RuntimeError(f"InterProScan did not produce a non-empty TSV: {out_file}")

    print("[InterProScan] Success.")
    print(f"[InterProScan] Output written to: {out_file}")


if __name__ == "__main__":
    main()