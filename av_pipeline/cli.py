import argparse
from .config import load_config
from .pipeline import run_pipeline
from .setup_run import create_run_dirs, write_external_markers

def main():
    parser = argparse.ArgumentParser(description="Run Antiviral Stage 1 pipeline.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--init-only", action="store_true")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    if args.strict:
        config["strict"] = True

    if args.init_only:
        create_run_dirs(config["run_dir"])
        write_external_markers(config["run_dir"])
        print(f"[INIT] Created folders and external markers under {config['run_dir']}")
        return

    run_pipeline(config)
