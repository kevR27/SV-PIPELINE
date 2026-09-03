#!/usr/bin/env python3
"""
Rank non-panel genes affected by genome-wide SVs.

The ranking is deliberately evidence-based and non-diagnostic:
- panel membership is retained as a flag but input genes are non-panel;
- exact matches to optic-neuropathy anchor HPO terms receive the strongest
  phenotype evidence;
- human HPO breadth provides additional support;
- AnnotSV's existing annotation/ranking columns are retained in the output
  where available;
- no score is converted into a pathogenic/disease-causing classification.

Classification is:
  Established panel gene
  Established/known human disease gene (if AnnotSV provides disease evidence)
  Phenotypically relevant non-panel candidate
  Other non-panel candidate
"""

import argparse
import csv
import re
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--annotsv", required=True)
parser.add_argument("--genes", required=True)
parser.add_argument("--phenotypes", required=True)
parser.add_argument("--panel", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

def read_list(path):
    out = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            x = line.strip()
            if x and not x.startswith("#"):
                out.add(x.upper())
    return out

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
        gene = row["gene_symbol"].upper()
        if gene not in genes:
            continue
        if row["hpo_id"].startswith("HP:"):
            pheno[gene]["hpo_count"] += 1
            pheno[gene]["hpos"].add(row["hpo_id"])
        if row["optic_neuropathy_anchor"] == "1":
            pheno[gene]["anchor_count"] += 1
        if row["source"]:
            pheno[gene]["sources"].add(row["source"])

# Collect useful AnnotSV evidence per gene.
annotsv = defaultdict(lambda: {
    "sv_count": 0,
    "gene_count": 0,
    "rank": [],
    "omim": [],
    "gencc": [],
    "clinvar": [],
    "lof": [],
})

with open(args.annotsv, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []

    def col(row, names):
        for n in names:
            if n in row and row[n] not in ("", "."):
                return row[n]
        return ""

    for row in reader:
        value = row.get("Gene_name", "")
        if not value or value == ".":
            continue
        row_genes = [x.strip().upper() for x in re.split(r"[;,|]", value) if x.strip()]
        for gene in row_genes:
            if gene not in genes:
                continue
            a = annotsv[gene]
            a["sv_count"] += 1

            rank = col(row, ["AnnotSV ranking", "Ranking", "ACMG_class"])
            if rank:
                a["rank"].append(rank)

            for key, names in {
                "omim": ["OMIM", "OMIM_inheritance"],
                "gencc": ["GENCC", "GENCC_ID"],
                "clinvar": ["ClinVar", "ClinVar_SV"],
                "lof": ["HI", "TS", "pLI", "LOEUF"],
            }.items():
                v = col(row, names)
                if v:
                    a[key].append(v)

rows = []

for gene in sorted(genes):
    p = pheno[gene]
    a = annotsv[gene]

    # Transparent discovery score, not a pathogenicity score.
    phenotype_score = min(10, p["anchor_count"] * 5) + min(3, p["hpo_count"] / 20.0)
    sv_score = min(5, a["sv_count"])
    evidence_score = phenotype_score + sv_score

    has_disease_evidence = bool(a["omim"] or a["gencc"] or a["clinvar"])

    if gene in panel:
        classification = "Established panel gene"
    elif p["anchor_count"] > 0 and has_disease_evidence:
        classification = "Known human disease gene with optic-neuropathy phenotype evidence"
    elif p["anchor_count"] > 0:
        classification = "Phenotypically relevant non-panel candidate"
    elif has_disease_evidence:
        classification = "Known human disease gene; optic-neuropathy link not established here"
    else:
        classification = "Other non-panel candidate"

    rows.append({
        "gene": gene,
        "panel_gene": "YES" if gene in panel else "NO",
        "SV_annotation_count": a["sv_count"],
        "human_HPO_count": p["hpo_count"],
        "optic_neuropathy_anchor_HPO_count": p["anchor_count"],
        "phenotype_score": round(phenotype_score, 3),
        "SV_evidence_score": sv_score,
        "integrated_discovery_score": round(evidence_score, 3),
        "AnnotSV_OMIM_evidence": ";".join(sorted(set(a["omim"]))),
        "AnnotSV_GENCC_evidence": ";".join(sorted(set(a["gencc"]))),
        "AnnotSV_ClinVar_evidence": ";".join(sorted(set(a["clinvar"]))),
        "AnnotSV_constraint_evidence": ";".join(sorted(set(a["lof"]))),
        "classification": classification,
        "interpretation": "Candidate only; not automatically classified as disease-causing.",
    })

rows.sort(
    key=lambda r: (
        -float(r["integrated_discovery_score"]),
        -int(r["optic_neuropathy_anchor_HPO_count"]),
        -int(r["SV_annotation_count"]),
        r["gene"],
    )
)

with open(args.output, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [
        "gene", "classification"
    ], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
