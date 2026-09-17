#!/usr/bin/env python3
"""Prioritize genes overlapped by genome-wide SVs for exploratory review.

This module is deliberately *not* a clinical pathogenicity classifier and it is
not a formal HPO semantic-similarity implementation.  The phenotype component
counts direct human Monarch gene→HPO associations and exact matches to the
optic-neuropathy anchor HPO terms supplied to the offline Monarch step.

The output therefore uses the terms ``phenotype_anchor_score`` and
``discovery_relevance_class``.  Backwards-compatible aliases
``phenotype_score`` and ``classification`` are retained because downstream
visualization/integration code may already expect them.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


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
        if name in row and row[name] not in ("", ".", None):
            return str(row[name])
    lower = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        value = lower.get(name.lower())
        if value not in ("", ".", None):
            return str(value)
    return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotsv", required=True)
    ap.add_argument("--genes", required=True, help="All genes overlapped by the genome-wide SV callset")
    ap.add_argument("--phenotypes", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genes = read_list(args.genes)
    panel = read_list(args.panel)

    phenotype = defaultdict(
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
            hpo = row.get("hpo_id", "")
            if hpo.startswith("HP:"):
                phenotype[gene]["hpo_count"] += 1
                phenotype[gene]["hpos"].add(hpo)
            if row.get("optic_neuropathy_anchor", "0") == "1":
                phenotype[gene]["anchor_count"] += 1
            source = row.get("source", "")
            if source:
                phenotype[gene]["sources"].add(source)

    annotsv = defaultdict(
        lambda: {
            "sv_count": 0,
            "rank": [],
            "omim": [],
            "gencc": [],
            "clinvar": [],
            "constraint": [],
        }
    )

    with open(args.annotsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            value = first(row, ["Gene_name", "Gene", "GENE", "Genes", "SYMBOL"])
            if not value:
                continue
            row_genes = {
                x.strip().upper()
                for x in re.split(r"[;,|]", value)
                if x.strip()
            }
            for gene in row_genes:
                if gene not in genes:
                    continue
                a = annotsv[gene]
                a["sv_count"] += 1

                rank = first(row, ["AnnotSV ranking", "Ranking", "ACMG_class"])
                if rank:
                    a["rank"].append(rank)

                for key, names in {
                    "omim": ["OMIM", "OMIM_inheritance"],
                    "gencc": ["GENCC", "GenCC", "GENCC_ID"],
                    "clinvar": ["ClinVar", "ClinVar_SV"],
                    "constraint": ["HI", "TS", "pLI", "LOEUF"],
                }.items():
                    value = first(row, names)
                    if value:
                        a[key].append(value)

    rows = []
    for gene in sorted(genes):
        p = phenotype[gene]
        a = annotsv[gene]

        # Exploratory relevance score only. Exact anchor-HPO associations are
        # deliberately weighted more heavily than general HPO breadth.
        phenotype_anchor_score = min(10.0, p["anchor_count"] * 5.0) + min(
            3.0, p["hpo_count"] / 20.0
        )
        sv_evidence_score = min(5, a["sv_count"])
        integrated_discovery_score = phenotype_anchor_score + sv_evidence_score

        has_disease_evidence = bool(a["omim"] or a["gencc"] or a["clinvar"])

        if gene in panel:
            discovery_class = "Established optic-neuropathy panel gene"
        elif p["anchor_count"] > 0 and has_disease_evidence:
            discovery_class = "Human disease gene with direct optic-neuropathy anchor-HPO evidence"
        elif p["anchor_count"] > 0:
            discovery_class = "Non-panel gene with direct optic-neuropathy anchor-HPO evidence"
        elif has_disease_evidence:
            discovery_class = "Known human disease gene without direct anchor-HPO match in this analysis"
        else:
            discovery_class = "Other SV-overlapped gene"

        rows.append(
            {
                "gene": gene,
                "panel_gene": "YES" if gene in panel else "NO",
                "SV_annotation_count": a["sv_count"],
                "human_HPO_count": p["hpo_count"],
                "optic_neuropathy_anchor_HPO_count": p["anchor_count"],
                "phenotype_anchor_score": round(phenotype_anchor_score, 3),
                # Backwards-compatible alias; interpretation is explicitly an
                # anchor score, not ontology semantic similarity.
                "phenotype_score": round(phenotype_anchor_score, 3),
                "SV_evidence_score": sv_evidence_score,
                "integrated_discovery_score": round(integrated_discovery_score, 3),
                "AnnotSV_OMIM_evidence": ";".join(sorted(set(a["omim"]))),
                "AnnotSV_GENCC_evidence": ";".join(sorted(set(a["gencc"]))),
                "AnnotSV_ClinVar_evidence": ";".join(sorted(set(a["clinvar"]))),
                "AnnotSV_constraint_evidence": ";".join(sorted(set(a["constraint"]))),
                "discovery_relevance_class": discovery_class,
                # Backwards-compatible alias.
                "classification": discovery_class,
                "interpretation": (
                    "Exploratory gene prioritization only; not an ACMG/ClinGen "
                    "pathogenicity classification and not a formal HPO semantic-similarity score."
                ),
            }
        )

    rows.sort(
        key=lambda r: (
            -float(r["integrated_discovery_score"]),
            -int(r["optic_neuropathy_anchor_HPO_count"]),
            -int(r["SV_annotation_count"]),
            r["gene"],
        )
    )

    fields = list(rows[0].keys()) if rows else [
        "gene",
        "discovery_relevance_class",
        "interpretation",
    ]
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] genes={len(rows)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
