from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from .ligand import combine_requests, read_manifest
from .paths import (
    require_absolute_external_dir,
    require_absolute_external_file,
    sanitize_name,
)
from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "CompoundRank: GNINA-scored docking hypotheses with physical-validity "
            "filtering, clustering, interaction evidence, and uncertainty."
        )
    )
    parser.add_argument("--receptor", required=True, help="Absolute receptor PDB path")
    parser.add_argument("--data-root", required=True, help="Absolute external data root")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--run-name", default=None)

    parser.add_argument("--ligand-file", action="append", default=[])
    parser.add_argument("--ligand-cid", action="append", default=[])
    parser.add_argument(
        "--ligand-smiles",
        action="append",
        default=[],
        help='Repeatable NAME=SMILES value',
    )
    parser.add_argument("--ligand-manifest", default=None)

    parser.add_argument("--center-x", type=float)
    parser.add_argument("--center-y", type=float)
    parser.add_argument("--center-z", type=float)
    parser.add_argument("--size-x", type=float)
    parser.add_argument("--size-y", type=float)
    parser.add_argument("--size-z", type=float)
    parser.add_argument("--autobox-ligand", default=None)
    parser.add_argument("--fpocket-padding", type=float, default=4.0)
    parser.add_argument("--fpocket-pocket", type=int, default=None)
    parser.add_argument(
        "--fpocket-top-n",
        type=int,
        default=1,
        help=(
            "When no explicit box or autobox ligand is provided, dock the top N "
            "fpocket pockets instead of only the highest-scoring pocket."
        ),
    )

    parser.add_argument("--seeds", nargs="+", type=int, default=[2026, 3101, 4202])
    parser.add_argument("--exhaustiveness", type=int, default=32)
    parser.add_argument("--num-modes", type=int, default=20)
    parser.add_argument(
        "--cnn-scoring",
        choices=["none", "rescore", "refinement", "all"],
        default="refinement",
    )
    parser.add_argument("--cluster-threshold", type=float, default=2.0)
    parser.add_argument("--max-hypotheses", type=int, default=3)
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--cpu", type=int, default=None)
    parser.add_argument("--device", type=int, default=None)

    parser.add_argument("--gnina-bin", default="gnina")
    parser.add_argument("--fpocket-bin", default="fpocket")
    parser.add_argument("--obabel-bin", default="obabel")
    parser.add_argument("--pdb2pqr-bin", default="pdb2pqr")
    parser.add_argument("--meeko-receptor-bin", default="mk_prepare_receptor.py")
    parser.add_argument("--meeko-ligand-bin", default="mk_prepare_ligand.py")
    parser.add_argument(
        "--posebusters-bin",
        default=os.path.expanduser("~/.venvs/posebusters/bin/bust"),
    )
    parser.add_argument("--skip-validity", action="store_true")
    parser.add_argument("--keep-workdir", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    receptor = require_absolute_external_file(args.receptor, "Receptor PDB")
    if receptor.suffix.lower() != ".pdb":
        raise ValueError("--receptor must be a PDB file")
    data_root = require_absolute_external_dir(args.data_root, "Data root", create=True)

    if args.output_dir:
        output_dir = require_absolute_external_dir(
            args.output_dir,
            "Output directory",
            create=True,
        )
    else:
        run_name = sanitize_name(args.run_name or receptor.stem)
        output_dir = require_absolute_external_dir(
            data_root / "results" / run_name,
            "Output directory",
            create=True,
        )

    manifest_requests = []
    if args.ligand_manifest:
        manifest_path = require_absolute_external_file(
            args.ligand_manifest,
            "Ligand manifest",
        )
        manifest_requests = read_manifest(manifest_path)

    # Validate file ligand paths as external before passing them onward.
    ligand_files = [
        str(require_absolute_external_file(value, "Ligand file"))
        for value in args.ligand_file
    ]
    requests = combine_requests(
        ligand_files,
        args.ligand_cid,
        args.ligand_smiles,
        manifest_requests,
    )
    normalized_requests = []
    for request in requests:
        if request.source_type == "file":
            external_path = require_absolute_external_file(
                request.value,
                f"Ligand file ({request.name})",
            )
            normalized_requests.append(
                type(request)(request.name, request.source_type, str(external_path))
            )
        else:
            normalized_requests.append(request)
    requests = normalized_requests

    autobox_ligand = None
    if args.autobox_ligand:
        autobox_ligand = require_absolute_external_file(
            args.autobox_ligand,
            "Autobox ligand",
        )

    run_pipeline(
        receptor_pdb=receptor,
        ligand_requests=requests,
        data_root=data_root,
        output_dir=output_dir,
        seeds=args.seeds,
        center_x=args.center_x,
        center_y=args.center_y,
        center_z=args.center_z,
        size_x=args.size_x,
        size_y=args.size_y,
        size_z=args.size_z,
        autobox_ligand=autobox_ligand,
        fpocket_padding=args.fpocket_padding,
        fpocket_pocket=args.fpocket_pocket,
        fpocket_top_n=args.fpocket_top_n,
        max_hypotheses=args.max_hypotheses,
        cluster_threshold=args.cluster_threshold,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        cnn_scoring=args.cnn_scoring,
        ph=args.ph,
        gnina_bin=args.gnina_bin,
        fpocket_bin=args.fpocket_bin,
        obabel_bin=args.obabel_bin,
        pdb2pqr_bin=args.pdb2pqr_bin,
        meeko_receptor_bin=args.meeko_receptor_bin,
        meeko_ligand_bin=args.meeko_ligand_bin,
        posebusters_bin=args.posebusters_bin,
        skip_validity=args.skip_validity,
        keep_workdir=args.keep_workdir,
        overwrite=args.overwrite,
        cpu=args.cpu,
        device=args.device,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
