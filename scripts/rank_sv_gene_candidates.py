#!/usr/bin/env python3
"""Prioritize genes intersected by the genome-wide SV callset.

This is a discovery-prioritization layer, not a pathogenicity classifier.
Panel membership, human phenotype associations and AnnotSV disease evidence are
kept as separate interpretable dimensions.  The numerical score is used only
to order candidates for review.
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
            "omim": [],
            "gencc": [],
            "clinvar": [],
            "constraint": [],
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

                value = first(row, ["AnnotSV ranking", "AnnotSV_rank", "ACMG_class"])
                if value:
                    data["rank"].append(value)

                for key, aliases in {
                    "omim": ["OMIM", "OMIM_inheritance"],
                    "gencc": ["GENCC", "GenCC", "GENCC_ID"],
                    "clinvar": ["ClinVar", "ClinVar_SV"],
                    "constraint": ["HI", "TS", "pLI", "LOEUF"],
                }.items():
                    value = first(row, aliases)
                    if value:
                        data[key].append(value)

    rows = []
    for gene in sorted(genes):
        p = pheno[gene]
        a = annotsv[gene]
        sv_count = len(a["sv_ids"])

        # Discovery-priority score only.  It does not encode pathogenicity.
        phenotype_component = min(10.0, p["anchor_count"] * 5.0) + min(
            3.0, p["hpo_count"] / 20.0
        )
        sv_component = min(5.0, float(sv_count))
        discovery_score = phenotype_component + sv_component

        has_disease_evidence = bool(a["omim"] or a["gencc"] or a["clinvar"])
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
                "SV_evidence_score": round(sv_component, 3),
                "integrated_discovery_score": round(discovery_score, 3),
                "AnnotSV_OMIM_evidence": ";".join(sorted(set(a["omim"]))),
                "AnnotSV_GENCC_evidence": ";".join(sorted(set(a["gencc"]))),
                "AnnotSV_ClinVar_evidence": ";".join(sorted(set(a["clinvar"]))),
                "AnnotSV_constraint_evidence": ";".join(sorted(set(a["constraint"]))),
                "candidate_group": candidate_group,
                "classification": candidate_group,
                "interpretation": (
                    "Discovery-priority category only; review the underlying SV, "
                    "inheritance, population frequency, phenotype fit and clinical "
                    "evidence before any pathogenicity assessment."
                ),
            }
        )

    rows.sort(
        key=lambda row: (
            -float(row["integrated_discovery_score"]),
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
