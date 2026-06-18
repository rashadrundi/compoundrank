"""Rule registry for evidence-grounded ligand retrieval.

This is intentionally small and explainable for Stage 4A.

The goal is not to be a full ligand database yet. The goal is to convert
target evidence into reasonable known-ligand candidates while leaving room for
future mutation-adapted and novel-ligand hypothesis generation.
"""

from __future__ import annotations

from typing import Any


def compound(
    name: str,
    *,
    design_status: str = "known_inhibitor",
    evidence_level: str = "moderate",
    notes: str = "",
) -> dict[str, Any]:
    """Create a seed compound definition for a retrieval rule."""
    return {
        "name": name,
        "design_status": design_status,
        "evidence_level": evidence_level,
        "notes": notes,
    }


LIGAND_RETRIEVAL_RULES: list[dict[str, Any]] = [
    {
        "rule_id": "retroviral_aspartyl_protease",
        "rule_label": "Retroviral aspartyl protease inhibitors",
        "match_terms": [
            "RVP",
            "pfam00077",
            "Retropepsins",
            "retropepsin",
            "retroviral aspartyl protease",
            "HIV-like retroviral aspartyl protease",
            "viral aspartyl protease",
            "HIV protease",
        ],
        "target_class": "viral protease",
        "enzyme_class": "aspartyl protease",
        "target_family_basis": "retroviral aspartyl protease",
        "retrieval_terms": [
            "HIV protease inhibitor",
            "retroviral protease inhibitor",
            "viral aspartyl protease inhibitor",
            "protease inhibitor antiviral",
        ],
        "seed_compounds": [
            compound("darunavir", evidence_level="strong"),
            compound("saquinavir", evidence_level="strong"),
            compound("ritonavir", evidence_level="strong"),
            compound("lopinavir", evidence_level="strong"),
            compound("atazanavir", evidence_level="strong"),
            compound("indinavir", evidence_level="strong"),
            compound("nelfinavir", evidence_level="strong"),
            compound("amprenavir", evidence_level="strong"),
            compound("tipranavir", evidence_level="strong"),
        ],
        "rule_evidence_level": "strong",
        "limitations": [
            "Retrieved compounds are known inhibitors of related retroviral protease targets.",
            "They are candidates for docking, not confirmed inhibitors of the submitted target.",
            "Mutation-aware adaptation is not evaluated in Stage 4A.",
        ],
    },
    {
        "rule_id": "viral_neuraminidase",
        "rule_label": "Viral neuraminidase inhibitors",
        "match_terms": [
            "neuraminidase",
            "sialidase",
            "influenza neuraminidase",
        ],
        "target_class": "viral neuraminidase",
        "enzyme_class": "glycosidase",
        "target_family_basis": "viral neuraminidase",
        "retrieval_terms": [
            "neuraminidase inhibitor",
            "influenza neuraminidase inhibitor",
        ],
        "seed_compounds": [
            compound("oseltamivir", evidence_level="strong"),
            compound("zanamivir", evidence_level="strong"),
            compound("peramivir", evidence_level="strong"),
            compound("laninamivir", evidence_level="moderate"),
        ],
        "rule_evidence_level": "strong",
        "limitations": [
            "Retrieved compounds are neuraminidase inhibitor candidates.",
            "Docking and target-specific validation are required.",
        ],
    },
    {
        "rule_id": "viral_rna_dependent_rna_polymerase",
        "rule_label": "Viral RNA-dependent RNA polymerase inhibitor candidates",
        "match_terms": [
            "RNA-dependent RNA polymerase",
            "RdRp",
            "viral polymerase",
            "RNA polymerase",
            "polymerase",
        ],
        "target_class": "viral polymerase",
        "enzyme_class": "polymerase",
        "target_family_basis": "viral RNA-dependent RNA polymerase or related viral polymerase",
        "retrieval_terms": [
            "viral polymerase inhibitor",
            "RNA-dependent RNA polymerase inhibitor",
            "nucleoside analog antiviral",
            "non-nucleoside polymerase inhibitor",
        ],
        "seed_compounds": [
            compound("remdesivir", evidence_level="moderate"),
            compound("molnupiravir", evidence_level="moderate"),
            compound("favipiravir", evidence_level="moderate"),
            compound("ribavirin", evidence_level="moderate"),
            compound("sofosbuvir", evidence_level="moderate"),
        ],
        "rule_evidence_level": "moderate",
        "limitations": [
            "Polymerase inhibitors can be target- and virus-specific.",
            "Some compounds are prodrugs or require biological activation not represented by docking alone.",
        ],
    },
    {
        "rule_id": "viral_integrase",
        "rule_label": "Viral integrase inhibitor candidates",
        "match_terms": [
            "integrase",
            "retroviral integrase",
            "viral integrase",
        ],
        "target_class": "viral integrase",
        "enzyme_class": "integrase",
        "target_family_basis": "viral integrase",
        "retrieval_terms": [
            "viral integrase inhibitor",
            "HIV integrase inhibitor",
            "integrase strand transfer inhibitor",
        ],
        "seed_compounds": [
            compound("raltegravir", evidence_level="strong"),
            compound("dolutegravir", evidence_level="strong"),
            compound("bictegravir", evidence_level="strong"),
            compound("elvitegravir", evidence_level="strong"),
        ],
        "rule_evidence_level": "moderate",
        "limitations": [
            "Known HIV integrase inhibitors may not transfer to all viral integrases.",
            "Metal coordination and active-site geometry need later validation.",
        ],
    },
    {
        "rule_id": "coronavirus_like_cysteine_protease",
        "rule_label": "Coronavirus-like cysteine protease inhibitor candidates",
        "match_terms": [
            "3CLpro",
            "Mpro",
            "main protease",
            "cysteine protease",
            "chymotrypsin-like protease",
            "peptidase C30",
        ],
        "target_class": "viral protease",
        "enzyme_class": "cysteine protease",
        "target_family_basis": "coronavirus-like cysteine protease",
        "retrieval_terms": [
            "viral cysteine protease inhibitor",
            "main protease inhibitor",
            "3CLpro inhibitor",
        ],
        "seed_compounds": [
            compound("nirmatrelvir", evidence_level="strong"),
            compound("boceprevir", evidence_level="moderate"),
            compound("GC376", evidence_level="moderate"),
        ],
        "rule_evidence_level": "moderate",
        "limitations": [
            "Some seed compounds are target-family candidates rather than general antiviral solutions.",
            "Covalent mechanisms are not fully represented by standard noncovalent docking.",
        ],
    },
]
