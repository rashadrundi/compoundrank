"""Generic target-to-ligand search utilities.

This module must not contain target-specific compound names.

Its first responsibility is to convert normalized target evidence into
general target, domain, enzyme, and target-class search queries. Later stages
will submit these queries to target-aware bioactivity databases such as ChEMBL.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class GenericLigandQuery:
    """One general ligand-discovery query derived from target evidence."""

    query: str
    retrieval_route: str
    specificity: int
    basis: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _clean_text(value: Any) -> str:
    if value is None:
        return ""

    text = " ".join(str(value).strip().split())
    if text.lower() in {"", "unknown", "none", "null", "n/a"}:
        return ""

    return text


def generate_generic_queries(
    target: dict[str, Any],
    *,
    max_queries: int = 16,
) -> list[dict[str, Any]]:
    """Generate target-independent ligand search queries.

    The function receives the normalized dictionary produced by
    compound_retrieval.target_context(). It intentionally does not consult
    LIGAND_RETRIEVAL_RULES and contains no target-specific drug names.
    """

    queries: list[GenericLigandQuery] = []
    seen: set[str] = set()

    def add(
        query: str,
        *,
        retrieval_route: str,
        specificity: int,
        basis: str,
    ) -> None:
        cleaned_query = _clean_text(query)
        cleaned_basis = _clean_text(basis)

        if not cleaned_query:
            return

        key = cleaned_query.casefold()
        if key in seen:
            return

        seen.add(key)
        queries.append(
            GenericLigandQuery(
                query=cleaned_query,
                retrieval_route=retrieval_route,
                specificity=specificity,
                basis=cleaned_basis or cleaned_query,
            )
        )

    # Terms already recommended by the target-evidence engine receive the
    # highest priority because they are produced from combined annotations.
    for term in target.get("query_terms") or []:
        add(
            str(term),
            retrieval_route="generic_evidence_query",
            specificity=100,
            basis="target_evidence.query_terms",
        )

    target_name = _clean_text(target.get("target_name"))
    domain_label = _clean_text(target.get("special_domain_label"))
    domain_accession = _clean_text(target.get("special_domain_accession"))
    enzyme_class = _clean_text(target.get("enzyme_class"))
    target_class = _clean_text(target.get("target_class"))

    if target_name:
        add(
            f"{target_name} inhibitor",
            retrieval_route="generic_exact_target_search",
            specificity=95,
            basis=target_name,
        )
        add(
            f"{target_name} ligand",
            retrieval_route="generic_exact_target_search",
            specificity=90,
            basis=target_name,
        )

    if domain_label:
        add(
            f"{domain_label} inhibitor",
            retrieval_route="generic_domain_search",
            specificity=85,
            basis=domain_label,
        )
        add(
            f"{domain_label} ligand",
            retrieval_route="generic_domain_search",
            specificity=80,
            basis=domain_label,
        )

    if domain_accession:
        add(
            domain_accession,
            retrieval_route="generic_domain_accession_search",
            specificity=82,
            basis=domain_accession,
        )

    if enzyme_class:
        add(
            f"{enzyme_class} inhibitor",
            retrieval_route="generic_enzyme_search",
            specificity=70,
            basis=enzyme_class,
        )
        add(
            f"{enzyme_class} ligand",
            retrieval_route="generic_enzyme_search",
            specificity=65,
            basis=enzyme_class,
        )

    if target_class:
        add(
            f"{target_class} inhibitor",
            retrieval_route="generic_target_class_search",
            specificity=55,
            basis=target_class,
        )
        add(
            f"{target_class} ligand",
            retrieval_route="generic_target_class_search",
            specificity=50,
            basis=target_class,
        )

    if max_queries > 0:
        queries = queries[:max_queries]

    return [query.to_dict() for query in queries]
