from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path
from typing import Any, Sequence

from compoundrank.functional_site_transfer import (
    read_single_fasta,
    write_homolog_pocket_evidence,
)
from compoundrank.uniprot_acquisition import (
    build_uniprot_acquisition,
    fetch_uniprot_entry,
    write_acquisition_outputs,
)


from .ramachandran import run_ramachandran_validation

WORKFLOW_SCHEMA_VERSION = (
    "reference_evidence_workflow.v0.1"
)


def choose_pdb_candidate(
    candidates: Sequence[
        dict[str, Any]
    ],
    *,
    requested_pdb_id: str | None = None,
) -> dict[str, Any]:
    usable = [
        candidate
        for candidate in candidates
        if isinstance(candidate, dict)
        and str(
            candidate.get("pdb_id")
            or ""
        ).strip()
    ]

    if not usable:
        raise ValueError(
            "UniProt entry has no usable "
            "PDB structure candidates."
        )

    if requested_pdb_id is None:
        return dict(usable[0])

    normalized = str(
        requested_pdb_id
    ).strip().upper()

    for candidate in usable:
        if (
            str(
                candidate.get("pdb_id")
                or ""
            ).strip().upper()
            == normalized
        ):
            return dict(candidate)

    raise ValueError(
        "Requested PDB structure "
        f"{normalized!r} was not present "
        "in the UniProt candidate list."
    )


def choose_reference_chain(
    candidate: dict[str, Any],
    *,
    functional_positions: Sequence[int],
    requested_chain: str | None = None,
) -> str:
    chain_ranges = candidate.get(
        "chain_ranges"
    )

    if not isinstance(
        chain_ranges,
        list,
    ):
        chain_ranges = []

    normalized_requested = (
        str(requested_chain)
        .strip()
        .upper()
        if requested_chain
        else None
    )

    coverage_by_chain: dict[
        str,
        set[int],
    ] = {}

    mapped_length_by_chain: dict[
        str,
        int,
    ] = {}

    for mapping in chain_ranges:
        if not isinstance(mapping, dict):
            continue

        chain = str(
            mapping.get("chain")
            or ""
        ).strip().upper()

        if not chain:
            continue

        try:
            start = int(
                mapping[
                    "uniprot_start"
                ]
            )

            end = int(
                mapping[
                    "uniprot_end"
                ]
            )

            mapped_length = int(
                mapping.get(
                    "mapped_length"
                )
                or end - start + 1
            )
        except (
            KeyError,
            TypeError,
            ValueError,
        ):
            continue

        covered = {
            int(position)
            for position
            in functional_positions
            if start
            <= int(position)
            <= end
        }

        coverage_by_chain.setdefault(
            chain,
            set(),
        ).update(covered)

        mapped_length_by_chain[
            chain
        ] = (
            mapped_length_by_chain.get(
                chain,
                0,
            )
            + mapped_length
        )

    if normalized_requested is not None:
        if (
            normalized_requested
            not in coverage_by_chain
        ):
            raise ValueError(
                "Requested chain "
                f"{normalized_requested!r} "
                "was not present in the "
                "selected PDB candidate."
            )

        return normalized_requested

    if not coverage_by_chain:
        raise ValueError(
            "Selected PDB candidate has no "
            "parseable chain mappings."
        )

    ranked_chains = sorted(
        coverage_by_chain,
        key=lambda chain: (
            -len(
                coverage_by_chain[
                    chain
                ]
            ),
            -mapped_length_by_chain.get(
                chain,
                0,
            ),
            chain,
        ),
    )

    return ranked_chains[0]


def _validate_pdb_bytes(
    data: bytes,
    *,
    pdb_id: str,
) -> None:
    if len(data) < 100:
        raise ValueError(
            f"Downloaded PDB {pdb_id} "
            "was unexpectedly small."
        )

    if (
        b"\nATOM  " not in data
        and not data.startswith(
            b"ATOM  "
        )
        and b"\nHETATM" not in data
        and not data.startswith(
            b"HETATM"
        )
    ):
        raise ValueError(
            f"Downloaded PDB {pdb_id} "
            "contains no coordinate records."
        )


def download_pdb_structure(
    pdb_id: str,
    output_path: Path,
    *,
    timeout_seconds: float = 60.0,
    force: bool = False,
) -> dict[str, Any]:
    normalized = str(
        pdb_id or ""
    ).strip().upper()

    if not normalized:
        raise ValueError(
            "PDB identifier is empty."
        )

    destination = Path(output_path)

    if (
        destination.is_file()
        and not force
    ):
        data = destination.read_bytes()

        _validate_pdb_bytes(
            data,
            pdb_id=normalized,
        )

        return {
            "pdb_id": normalized,
            "path": str(destination),
            "cached": True,
            "size_bytes": len(data),
            "url": (
                "https://files.rcsb.org/"
                f"download/{normalized}.pdb"
            ),
        }

    url = (
        "https://files.rcsb.org/"
        f"download/{normalized}.pdb"
    )

    request = urllib.request.Request(
        url,
        headers={
            "Accept": "chemical/x-pdb,"
            "text/plain,*/*",
            "User-Agent": (
                "CompoundRank-EXORCIST/0.1"
            ),
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=timeout_seconds,
    ) as response:
        data = response.read()

        metadata = {
            "pdb_id": normalized,
            "url": response.geturl(),
            "status": getattr(
                response,
                "status",
                None,
            ),
            "content_type": (
                response.headers.get(
                    "Content-Type"
                )
            ),
            "etag": response.headers.get(
                "ETag"
            ),
            "last_modified": (
                response.headers.get(
                    "Last-Modified"
                )
            ),
            "cached": False,
            "size_bytes": len(data),
        }

    _validate_pdb_bytes(
        data,
        pdb_id=normalized,
    )

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    partial_path = destination.with_name(
        destination.name + ".part"
    )

    partial_path.write_bytes(data)
    partial_path.replace(destination)

    metadata["path"] = str(
        destination
    )

    return metadata


def run_reference_evidence_workflow(
    *,
    payload: dict[str, Any],
    submitted_fasta: Path,
    receptor_pdb: Path,
    output_dir: Path,
    response_metadata: (
        dict[str, Any] | None
    ) = None,
    requested_pdb_id: str | None = None,
    reference_chain_id: str | None = None,
    receptor_chain_id: str | None = None,
    reference_pdb: Path | None = None,
    timeout_seconds: float = 60.0,
    force_download: bool = False,
    requested_selection_mode: str = (
        "prioritize_supported"
    ),
    minimum_mapping_fraction: float = 0.50,
    minimum_mapped_residues: int = 2,
    minimum_structure_identity: float = 0.90,
    minimum_structure_coverage: float = 0.80,
) -> dict[str, Any]:
    output = Path(output_dir)

    output.mkdir(
        parents=True,
        exist_ok=True,
    )

    acquisition = (
        build_uniprot_acquisition(
            payload,
            requested_selection_mode=(
                requested_selection_mode
            ),
        )
    )

    acquisition_output = (
        output / "acquisition"
    )

    acquisition_paths = (
        write_acquisition_outputs(
            acquisition_output,
            payload=payload,
            response_metadata=(
                response_metadata
            ),
            requested_selection_mode=(
                requested_selection_mode
            ),
        )
    )

    reference_record = acquisition[
        "reference_record"
    ]

    functional_positions = [
        int(
            residue[
                "sequence_position"
            ]
        )
        for residue in reference_record[
            "residues"
        ]
    ]

    candidate = choose_pdb_candidate(
        acquisition[
            "pdb_candidates"
        ],
        requested_pdb_id=(
            requested_pdb_id
        ),
    )

    selected_pdb_id = str(
        candidate["pdb_id"]
    ).upper()

    selected_chain = (
        choose_reference_chain(
            candidate,
            functional_positions=(
                functional_positions
            ),
            requested_chain=(
                reference_chain_id
            ),
        )
    )

    if reference_pdb is None:
        reference_pdb_path = (
            output
            / "reference_structure"
            / f"{selected_pdb_id}.pdb"
        )

        download_metadata = (
            download_pdb_structure(
                selected_pdb_id,
                reference_pdb_path,
                timeout_seconds=(
                    timeout_seconds
                ),
                force=force_download,
            )
        )
    else:
        reference_pdb_path = Path(
            reference_pdb
        )

        if not reference_pdb_path.is_file():
            raise FileNotFoundError(
                reference_pdb_path
            )

        reference_data = (
            reference_pdb_path
            .read_bytes()
        )

        _validate_pdb_bytes(
            reference_data,
            pdb_id=selected_pdb_id,
        )

        download_metadata = {
            "pdb_id": selected_pdb_id,
            "path": str(
                reference_pdb_path
            ),
            "cached": True,
            "provided_locally": True,
            "size_bytes": len(
                reference_data
            ),
        }

    reference_validation = (
        run_ramachandran_validation(
            reference_pdb_path,
            (
                output
                / "structure_validation"
                / "reference"
            ),
            chain_id=selected_chain,
            continue_on_error=True,
        )
    )

    reference_validation_report = (
        reference_validation[
            "report"
        ]
    )

    reference_validation_outputs = (
        reference_validation[
            "outputs"
        ]
    )

    reference_validation_summary = (
        reference_validation_report.get(
            "summary",
            {},
        )
    )

    if not isinstance(
        reference_validation_summary,
        dict,
    ):
        reference_validation_summary = {}

    print(
        "[RAMACHANDRAN] "
        "reference structure: "
        f"status="
        f"{reference_validation_report.get('status')}; "
        f"evaluable="
        f"{reference_validation_report.get('evaluable_residues')}; "
        f"flag="
        f"{reference_validation_summary.get('screening_flag')}"
    )

    receptor_path = Path(
        receptor_pdb
    )

    if not receptor_path.is_file():
        raise FileNotFoundError(
            receptor_path
        )

    submitted_record = (
        read_single_fasta(
            Path(submitted_fasta)
        )
    )

    pocket_evidence_path = (
        output / "pocket_evidence.json"
    )

    pocket_evidence = (
        write_homolog_pocket_evidence(
            pocket_evidence_path,
            reference_record=(
                reference_record
            ),
            reference_sequence=(
                acquisition["sequence"]
            ),
            submitted_sequence=(
                submitted_record[
                    "sequence"
                ]
            ),
            reference_pdb=(
                reference_pdb_path
            ),
            receptor_pdb=(
                receptor_path
            ),
            reference_chain_id=(
                selected_chain
            ),
            receptor_chain_id=(
                receptor_chain_id
            ),
            minimum_mapping_fraction=(
                minimum_mapping_fraction
            ),
            minimum_mapped_residues=(
                minimum_mapped_residues
            ),
            minimum_structure_identity=(
                minimum_structure_identity
            ),
            minimum_structure_coverage=(
                minimum_structure_coverage
            ),
        )
    )

    workflow_summary = {
        "schema_version": (
            WORKFLOW_SCHEMA_VERSION
        ),
        "status": "complete",
        "primary_accession": (
            acquisition[
                "summary"
            ]["primary_accession"]
        ),
        "submitted_fasta": str(
            submitted_fasta
        ),
        "receptor_pdb": str(
            receptor_path
        ),
        "receptor_chain": (
            receptor_chain_id
        ),
        "selected_reference": {
            "pdb_id": (
                selected_pdb_id
            ),
            "chain": selected_chain,
            "candidate": candidate,
            "structure_path": str(
                reference_pdb_path
            ),
            "download": (
                download_metadata
            ),
        },
        "reference_structure_validation": {
            "status": (
                reference_validation_report.get(
                    "status"
                )
            ),
            "selection_mode": (
                reference_validation_report.get(
                    "selection_mode"
                )
            ),
            "chain": selected_chain,
            "evaluable_residues": (
                reference_validation_report.get(
                    "evaluable_residues"
                )
            ),
            "summary": (
                reference_validation_report.get(
                    "summary"
                )
            ),
            "outputs": (
                reference_validation_outputs
            ),
        },
        "acquisition_outputs": (
            acquisition_paths
        ),
        "functional_site_positions": (
            functional_positions
        ),
        "pocket_evidence_path": str(
            pocket_evidence_path
        ),
        "selection_mode": (
            pocket_evidence[
                "selection_mode"
            ]
        ),
        "confidence": (
            pocket_evidence[
                "confidence"
            ]
        ),
        "mapped_residues": (
            pocket_evidence[
                "residues"
            ]
        ),
        "transfer_summary": (
            pocket_evidence[
                "transfer_summary"
            ]
        ),
    }

    summary_path = (
        output / "workflow_summary.json"
    )

    summary_path.write_text(
        json.dumps(
            workflow_summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    workflow_summary[
        "workflow_summary_path"
    ] = str(summary_path)

    return workflow_summary


def build_cli_parser() -> (
    argparse.ArgumentParser
):
    parser = argparse.ArgumentParser(
        description=(
            "Acquire UniProt functional-site "
            "annotations, select and download "
            "a PDB reference structure, and "
            "transfer the residues onto a "
            "submitted receptor."
        )
    )

    source_group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    source_group.add_argument(
        "--accession",
        default=None,
    )

    source_group.add_argument(
        "--input-json",
        type=Path,
        default=None,
    )

    parser.add_argument(
        "--submitted-fasta",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--receptor-pdb",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--pdb-id",
        default=None,
    )

    parser.add_argument(
        "--reference-chain",
        default=None,
    )

    parser.add_argument(
        "--receptor-chain",
        default=None,
    )

    parser.add_argument(
        "--reference-pdb",
        type=Path,
        default=None,
        help=(
            "Optional local PDB override. "
            "The selected UniProt-linked PDB "
            "identifier is still recorded."
        ),
    )

    parser.add_argument(
        "--selection-mode",
        choices=(
            "report_only",
            "prioritize_supported",
        ),
        default="prioritize_supported",
    )

    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=60.0,
    )

    parser.add_argument(
        "--force-download",
        action="store_true",
    )

    parser.add_argument(
        "--minimum-mapping-fraction",
        type=float,
        default=0.50,
    )

    parser.add_argument(
        "--minimum-mapped-residues",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--minimum-structure-identity",
        type=float,
        default=0.90,
    )

    parser.add_argument(
        "--minimum-structure-coverage",
        type=float,
        default=0.80,
    )

    return parser


def main(
    argv: Sequence[str] | None = None,
) -> int:
    parser = build_cli_parser()

    arguments = parser.parse_args(
        argv
    )

    if not (
        0.0
        <= arguments.minimum_mapping_fraction
        <= 1.0
    ):
        parser.error(
            "--minimum-mapping-fraction "
            "must be between 0 and 1."
        )

    if (
        arguments.minimum_mapped_residues
        < 1
    ):
        parser.error(
            "--minimum-mapped-residues "
            "must be at least 1."
        )

    if not (
        0.0
        <= arguments.minimum_structure_identity
        <= 1.0
    ):
        parser.error(
            "--minimum-structure-identity "
            "must be between 0 and 1."
        )

    if not (
        0.0
        <= arguments.minimum_structure_coverage
        <= 1.0
    ):
        parser.error(
            "--minimum-structure-coverage "
            "must be between 0 and 1."
        )

    if arguments.input_json is not None:
        payload = json.loads(
            arguments.input_json.read_text(
                encoding="utf-8"
            )
        )

        response_metadata = {
            "source": "input_json",
            "input_path": str(
                arguments.input_json
            ),
        }
    else:
        payload, response_metadata = (
            fetch_uniprot_entry(
                arguments.accession,
                timeout_seconds=(
                    arguments.timeout_seconds
                ),
            )
        )

    result = (
        run_reference_evidence_workflow(
            payload=payload,
            submitted_fasta=(
                arguments.submitted_fasta
            ),
            receptor_pdb=(
                arguments.receptor_pdb
            ),
            output_dir=(
                arguments.output_dir
            ),
            response_metadata=(
                response_metadata
            ),
            requested_pdb_id=(
                arguments.pdb_id
            ),
            reference_chain_id=(
                arguments.reference_chain
            ),
            receptor_chain_id=(
                arguments.receptor_chain
            ),
            reference_pdb=(
                arguments.reference_pdb
            ),
            timeout_seconds=(
                arguments.timeout_seconds
            ),
            force_download=(
                arguments.force_download
            ),
            requested_selection_mode=(
                arguments.selection_mode
            ),
            minimum_mapping_fraction=(
                arguments.minimum_mapping_fraction
            ),
            minimum_mapped_residues=(
                arguments.minimum_mapped_residues
            ),
            minimum_structure_identity=(
                arguments.minimum_structure_identity
            ),
            minimum_structure_coverage=(
                arguments.minimum_structure_coverage
            ),
        )
    )

    print(
        json.dumps(
            {
                "status": result["status"],
                "primary_accession": (
                    result[
                        "primary_accession"
                    ]
                ),
                "reference_pdb_id": (
                    result[
                        "selected_reference"
                    ]["pdb_id"]
                ),
                "reference_chain": (
                    result[
                        "selected_reference"
                    ]["chain"]
                ),
                "selection_mode": (
                    result[
                        "selection_mode"
                    ]
                ),
                "confidence": (
                    result["confidence"]
                ),
                "mapped_residue_count": (
                    len(
                        result[
                            "mapped_residues"
                        ]
                    )
                ),
                "mapped_residues": (
                    result[
                        "mapped_residues"
                    ]
                ),
                "mapping_fraction": (
                    result[
                        "transfer_summary"
                    ]["mapping_fraction"]
                ),
                "selection_checks": (
                    result[
                        "transfer_summary"
                    ]["selection_checks"]
                ),
                "pocket_evidence": (
                    result[
                        "pocket_evidence_path"
                    ]
                ),
                "workflow_summary": (
                    result[
                        "workflow_summary_path"
                    ]
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
