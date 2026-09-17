#!/usr/bin/env python3
"""Prioritize genes affected by genome-wide SVs using transparent discovery evidence.

This is a discovery/prioritization table, not a pathogenicity classifier.
Panel membership, human phenotype associations and AnnotSV disease evidence are
kept as separate interpretable dimensions. The numeric score is only a ranking
heuristic for review and must not be interpreted as ACMG/AMP evidence.
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

GENE_COLUMNS = ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL", "AnnotSV_Gene"]


def read_list(path: str) -> set[str]:
    out = set()
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
    lowered = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        value = lowered.get(name.lower())
        if value not in ("", ".", None):
            return str(value)
    return ""


def split_genes(value: str) -> list[str]:
    return sorted({
        x.strip().upper()
        for x in re.split(r"[;,|/]", value or "")
        if x.strip() and x.strip() not in {".", "NA", "N/A"}
    })


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotsv", required=True)
    ap.add_argument("--genes", required=True, help="All genes intersected by the genome-wide SV callset")
    ap.add_argument("--phenotypes", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genes = read_list(args.genes)
    panel = read_list(args.panel)

    pheno = defaultdict(lambda: {
        "hpo_count": 0,
        "anchor_count": 0,
        "hpos": set(),
        "sources": set(),
    })

    with open(args.phenotypes, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row.get("gene_symbol", "").upper()
            if not gene or gene not in genes:
                continue
            hpo_id = row.get("hpo_id", "")
            if hpo_id.startswith("HP:"):
                pheno[gene]["hpo_count"] += 1
                pheno[gene]["hpos"].add(hpo_id)
            if row.get("optic_neuropathy_anchor") == "1":
                pheno[gene]["anchor_count"] += 1
            if row.get("source"):
                pheno[gene]["sources"].add(row["source"])

    annotsv = defaultdict(lambda: {
        "sv_count": 0,
        "rank": [],
        "omim": [],
        "gencc": [],
        "clinvar": [],
        "constraint": [],
    })

    with open(args.annotsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        gene_col = next((c for c in GENE_COLUMNS if c in fieldnames), None)
        if gene_col is None:
            lowered = {c.lower(): c for c in fieldnames}
            gene_col = next((lowered[c.lower()] for c in GENE_COLUMNS if c.lower() in lowered), None)
        if gene_col is None:
            raise ValueError("No supported AnnotSV gene column found: " + ", ".join(fieldnames))

        for row in reader:
            for gene in split_genes(row.get(gene_col, "")):
                if gene not in genes:
                    continue
                evidence = annotsv[gene]
                evidence["sv_count"] += 1

                rank = first(row, ["AnnotSV_ranking", "AnnotSV ranking", "Ranking", "ACMG_class"])
                if rank:
                    evidence["rank"].append(rank)

                for key, names in {
                    "omim": ["OMIM", "OMIM_inheritance"],
                    "gencc": ["GENCC", "GenCC", "GENCC_ID"],
                    "clinvar": ["ClinVar", "ClinVar_SV"],
                    "constraint": ["HI", "TS", "pLI", "LOEUF"],
                }.items():
                    value = first(row, names)
                    if value:
                        evidence[key].append(value)

    rows = []
    for gene in sorted(genes):
        p = pheno[gene]
        a = annotsv[gene]

        # Transparent discovery score only. Exact matches to user-defined optic
        # neuropathy anchor HPO terms contribute more strongly than phenotype
        # breadth; SV annotation count gives only a capped secondary signal.
        phenotype_score = min(10.0, p["anchor_count"] * 5.0) + min(3.0, p["hpo_count"] / 20.0)
        sv_score = min(5, a["sv_count"])
        discovery_score = phenotype_score + sv_score

        has_disease_evidence = bool(a["omim"] or a["gencc"] or a["clinvar"])
        if gene in panel:
            category = "PANEL_GENE"
        elif p["anchor_count"] > 0 and has_disease_evidence:
            category = "NONPANEL_HPO_AND_DISEASE_EVIDENCE"
        elif p["anchor_count"] > 0:
            category = "NONPANEL_HPO_EVIDENCE"
        elif has_disease_evidence:
            category = "NONPANEL_DISEASE_EVIDENCE"
        else:
            category = "NONPANEL_OTHER"

        rows.append({
            "gene": gene,
            "panel_gene": "YES" if gene in panel else "NO",
            "SV_annotation_count": a["sv_count"],
            "human_HPO_count": p["hpo_count"],
            "optic_neuropathy_anchor_HPO_count": p["anchor_count"],
            "phenotype_score": round(phenotype_score, 3),
            "SV_evidence_score": sv_score,
            "integrated_discovery_score": round(discovery_score, 3),
            "AnnotSV_OMIM_evidence": ";".join(sorted(set(a["omim"]))),
            "AnnotSV_GENCC_evidence": ";".join(sorted(set(a["gencc"]))),
            "AnnotSV_ClinVar_evidence": ";".join(sorted(set(a["clinvar"]))),
            "AnnotSV_constraint_evidence": ";".join(sorted(set(a["constraint"]))),
            "candidate_category": category,
            "interpretation": "Discovery candidate only; not an ACMG/AMP or pathogenicity classification.",
        })

    rows.sort(
        key=lambda row: (
            -float(row["integrated_discovery_score"]),
            -int(row["optic_neuropathy_anchor_HPO_count"]),
            -int(row["SV_annotation_count"]),
            row["gene"],
        )
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else ["gene", "candidate_category"]
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] ranked_genes={len(rows)} output={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
