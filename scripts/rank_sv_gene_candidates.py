#!/usr/bin/env python3
"""Prioritize genes intersected by the genome-wide SV callset.

This is a research-prioritization layer, not a pathogenicity classifier.

The ranking keeps three ideas separate:
1. Phenotype relevance: does the gene match the optic-neuropathy/HPO context?
2. Gene-disease evidence: how strong is the known human gene-disease evidence?
3. SV evidence: how many master SVs affect the gene?

GenCC classifications are converted into a small transparent numerical score
only to help order candidates. The score is not a probability of pathogenicity
and does not replace clinical variant interpretation.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


# Internal prioritization weights for GenCC terms.
# These values are deliberately simple and monotonic. They are NOT an
# official GenCC pathogenicity scale and must not be interpreted as a
# probability that a gene or variant is pathogenic.
GENCC_SCORES = {
    "DEFINITIVE": 4.0,
    "STRONG": 3.5,
    "MODERATE": 2.5,
    "SUPPORTIVE": 1.5,
    "LIMITED": 0.5,
    "ANIMAL MODEL ONLY": 0.25,
    "DISPUTED": 0.0,
    "REFUTED": 0.0,
    "NO KNOWN DISEASE RELATIONSHIP": 0.0,
}

POSITIVE_GENCC = {
    "DEFINITIVE",
    "STRONG",
    "MODERATE",
    "SUPPORTIVE",
    "LIMITED",
}

CONTRADICTORY_GENCC = {
    "DISPUTED",
    "REFUTED",
    "NO KNOWN DISEASE RELATIONSHIP",
}


def read_list(path: str) -> set[str]:
    out: set[str] = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            value = line.strip()
            if value and not value.startswith("#"):
                out.add(value.upper())
    return out


def first(row: dict, names: list[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value not in ("", ".", None):
            return str(value)
    lower = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        value = lower.get(name.lower(), "")
        if value not in ("", ".", None):
            return str(value)
    return ""


def split_values(value: str) -> list[str]:
    if value in ("", ".", None):
        return []
    return [
        x.strip()
        for x in re.split(r"[;|]", str(value))
        if x.strip()
    ]


def is_positive_flag(value: str) -> bool:
    return str(value or "").strip().upper() in {
        "YES",
        "TRUE",
        "1",
        "Y",
    }


def normalize_gencc_term(value: str) -> str | None:
    """Map raw GenCC text to one harmonized term when possible."""
    text = str(value or "").strip().upper()
    if not text:
        return None

    if "NO KNOWN" in text:
        return "NO KNOWN DISEASE RELATIONSHIP"
    if "ANIMAL" in text:
        return "ANIMAL MODEL ONLY"
    if "DEFINITIVE" in text:
        return "DEFINITIVE"
    if "STRONG" in text:
        return "STRONG"
    if "MODERATE" in text:
        return "MODERATE"
    if "SUPPORTIVE" in text:
        return "SUPPORTIVE"
    if "LIMITED" in text:
        return "LIMITED"
    if "DISPUTED" in text:
        return "DISPUTED"
    if "REFUTED" in text:
        return "REFUTED"
    return None


def summarize_gencc(values: list[str], has_omim: bool) -> dict[str, str | float]:
    terms: set[str] = set()
    for raw in values:
        pieces = split_values(raw) or [raw]
        for piece in pieces:
            term = normalize_gencc_term(piece)
            if term:
                terms.add(term)

    positive = sorted(
        (term for term in terms if term in POSITIVE_GENCC),
        key=lambda term: GENCC_SCORES[term],
        reverse=True,
    )
    contradictory = terms & CONTRADICTORY_GENCC
    conflict = bool(positive and contradictory)

    if positive:
        best = positive[0]
        score = GENCC_SCORES[best]
        source = "GenCC"
    elif "ANIMAL MODEL ONLY" in terms:
        best = "ANIMAL MODEL ONLY"
        score = GENCC_SCORES[best]
        source = "GenCC"
    elif terms:
        # Disputed/refuted/no-known should not contribute positive evidence.
        best = sorted(terms)[0]
        score = 0.0
        source = "GenCC"
    elif has_omim:
        # OMIM disease evidence is useful for prioritization, but it is not
        # equivalent to a GenCC Definitive/Strong/Moderate curation.
        best = "OMIM_DISEASE_EVIDENCE"
        score = 1.0
        source = "OMIM"
    else:
        best = "NO_CURATED_EVIDENCE"
        score = 0.0
        source = "NONE"

    return {
        "best": best,
        "score": score,
        "source": source,
        "conflict": "YES" if conflict else "NO",
        "all_terms": ";".join(sorted(terms)) if terms else ".",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotsv", required=True)
    parser.add_argument("--genes", required=True, help="All genes intersected by master SVs")
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    genes = read_list(args.genes)
    panel = read_list(args.panel)

    pheno = defaultdict(
        lambda: {
            "hpo_count": 0,
            "anchor_count": 0,
            "hpos": set(),
            "sources": set(),
        }
    )

    with open(args.phenotypes, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row.get("gene_symbol", "").upper()
            if gene not in genes:
                continue
            hpo_id = row.get("hpo_id", "")
            if hpo_id.startswith("HP:"):
                pheno[gene]["hpo_count"] += 1
                pheno[gene]["hpos"].add(hpo_id)
            if row.get("optic_neuropathy_anchor", "0") == "1":
                pheno[gene]["anchor_count"] += 1
            source = row.get("source", "")
            if source:
                pheno[gene]["sources"].add(source)

    annotsv = defaultdict(
        lambda: {
            "sv_ids": set(),
            "rank": [],
            "ranking_score": [],
            "ranking_criteria": [],
            "acmg_class": [],
            "omim": [],
            "omim_morbid": [],
            "gencc_classification": [],
            "gencc_disease": [],
            "gencc_moi": [],
            "clinvar": [],
            "constraint": [],
            "hi": [],
            "ts": [],
        }
    )

    with open(args.annotsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_field = first(
                row,
                ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL"],
            )
            if not gene_field:
                continue

            row_genes = {
                x.strip().upper()
                for x in re.split(r"[;,|]", gene_field)
                if x.strip()
            }
            sv_id = first(row, ["SV_ID", "AnnotSV_ID", "ID"])

            for gene in row_genes:
                if gene not in genes:
                    continue

                data = annotsv[gene]
                if sv_id:
                    data["sv_ids"].add(sv_id)

                for key, aliases in {
                    "rank": ["AnnotSV ranking", "AnnotSV_rank", "ACMG_class"],
                    "ranking_score": ["AnnotSV_ranking_score", "AnnotSV ranking score"],
                    "ranking_criteria": ["AnnotSV_ranking_criteria", "AnnotSV ranking criteria"],
                    "acmg_class": ["ACMG_class", "ACMG class"],
                    "omim": [
                        "OMIM_phenotype",
                        "OMIM",
                        "AnnotSV_OMIM",
                        "AnnotSV_OMIM_evidence",
                    ],
                    "omim_morbid": ["OMIM_morbid", "OMIM_morbid_candidate"],
                    "gencc_classification": [
                        "GenCC_classification",
                        "GENCC_classification",
                        "GENCC",
                        "GenCC",
                    ],
                    "gencc_disease": ["GenCC_disease", "GENCC_disease"],
                    "gencc_moi": ["GenCC_moi", "GENCC_moi"],
                    "clinvar": ["ClinVar", "ClinVar_SV"],
                    "constraint": ["HI", "TS", "pLI", "LOEUF", "LOEUF_bin", "GnomAD_pLI"],
                    "hi": ["HI"],
                    "ts": ["TS"],
                }.items():
                    value = first(row, aliases)
                    if value:
                        data[key].append(value)

    rows = []
    for gene in sorted(genes):
        p = pheno[gene]
        a = annotsv[gene]
        sv_count = len(a["sv_ids"])

        # Phenotype remains the strongest discovery component.
        # Optic-neuropathy anchor terms receive more weight than broad HPO
        # associations because they are closer to the phenotype of interest.
        phenotype_component = min(10.0, p["anchor_count"] * 5.0) + min(
            3.0, p["hpo_count"] / 20.0
        )

        # Repeated SV calls in the same gene are weak supporting evidence, not
        # a proxy for pathogenicity. Their weight is therefore deliberately
        # small and capped at 2 points.
        sv_component = min(2.0, float(sv_count) * 0.5)

        has_omim = bool(a["omim"]) or any(
            is_positive_flag(value)
            for value in a["omim_morbid"]
        )
        gencc_summary = summarize_gencc(a["gencc_classification"], has_omim)
        gene_disease_component = float(gencc_summary["score"])

        # Research-priority score only:
        # phenotype (0-13) + curated gene-disease evidence (0-4)
        # + limited SV-count contribution (0-2).
        integrated_discovery_score = (
            phenotype_component
            + gene_disease_component
            + sv_component
        )

        has_disease_evidence = (
            gene_disease_component > 0
            or has_omim
            or bool(a["clinvar"])
        )

        if gene in panel:
            candidate_group = "PANEL_GENE"
        elif p["anchor_count"] > 0 and has_disease_evidence:
            candidate_group = "NONPANEL_HPO_AND_DISEASE_EVIDENCE"
        elif p["anchor_count"] > 0:
            candidate_group = "NONPANEL_HPO_OVERLAP"
        elif has_disease_evidence:
            candidate_group = "NONPANEL_HUMAN_DISEASE_GENE"
        else:
            candidate_group = "OTHER_NONPANEL_CANDIDATE"

        rows.append(
            {
                "gene": gene,
                "panel_gene": "YES" if gene in panel else "NO",
                "SV_count": sv_count,
                "human_HPO_count": p["hpo_count"],
                "optic_neuropathy_anchor_HPO_count": p["anchor_count"],
                "phenotype_score": round(phenotype_component, 3),
                "gene_disease_evidence_score": round(gene_disease_component, 3),
                "gene_disease_evidence_level": gencc_summary["best"],
                "gene_disease_evidence_source": gencc_summary["source"],
                "gene_disease_evidence_conflict": gencc_summary["conflict"],
                "GenCC_classifications": gencc_summary["all_terms"],
                "GenCC_disease": ";".join(sorted(set(a["gencc_disease"]))),
                "GenCC_moi": ";".join(sorted(set(a["gencc_moi"]))),
                "ClinGen_HI": ";".join(sorted(set(a["hi"]))),
                "ClinGen_TS": ";".join(sorted(set(a["ts"]))),
                "SV_evidence_score": round(sv_component, 3),
                "integrated_discovery_score": round(integrated_discovery_score, 3),
                "AnnotSV_ranking_scores": ";".join(sorted(set(a["ranking_score"]))),
                "AnnotSV_ranking_criteria": ";".join(sorted(set(a["ranking_criteria"]))),
                "AnnotSV_ACMG_classes": ";".join(sorted(set(a["acmg_class"]))),
                "AnnotSV_OMIM_evidence": ";".join(sorted(set(a["omim"]))),
                "AnnotSV_GENCC_evidence": ";".join(sorted(set(a["gencc_classification"]))),
                "AnnotSV_ClinVar_evidence": ";".join(sorted(set(a["clinvar"]))),
                "AnnotSV_constraint_evidence": ";".join(sorted(set(a["constraint"]))),

                # Explicit interpretation columns retained in the gene-ranking
                # output so the most useful biological context is visible
                # without reopening the full AnnotSV table.
                "PANEL_STATUS": "PANEL_GENE" if gene in panel else "NONPANEL_GENE",
                "CANDIDATE_CLASS": candidate_group,
                "OMIM": ";".join(sorted(set(a["omim"]))) or ".",
                "GENCC": gencc_summary["all_terms"],
                "ANNOTSV_GENERAL_CLASSIFICATION": (
                    ";".join(sorted(set(a["acmg_class"]))) or "."
                ),
                "candidate_group": candidate_group,
                "classification": candidate_group,
                "ranking_model": "phenotype13_geneDisease4_sv2_v2",
                "interpretation": (
                    "Research-priority score only. The gene-disease component "
                    "summarizes curated evidence and is not a gene pathogenicity "
                    "probability. Review the specific SV, inheritance, dosage "
                    "mechanism, population frequency, phenotype fit and clinical "
                    "evidence before pathogenicity assessment."
                ),
            }
        )

    rows.sort(
        key=lambda row: (
            -float(row["integrated_discovery_score"]),
            -float(row["gene_disease_evidence_score"]),
            -int(row["optic_neuropathy_anchor_HPO_count"]),
            -int(row["SV_count"]),
            row["gene"],
        )
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else [
        "gene",
        "panel_gene",
        "gene_disease_evidence_score",
        "gene_disease_evidence_level",
        "candidate_group",
        "interpretation",
    ]
    with output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] ranked_genes={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
