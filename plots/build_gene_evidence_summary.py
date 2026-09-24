#!/usr/bin/env python3
"""Build a gene-centric summary from the extended SV-gene evidence table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv


MISSING = {"", ".", "NA", "N/A", "nan", "NaN", "None"}


def parse_args():
    p = argparse.ArgumentParser(description="Build a gene-level multimodal evidence summary.")
    p.add_argument("--integrated", required=True)
    p.add_argument("--ranking", default=None)
    p.add_argument("--phenotypes", default=None)
    p.add_argument("--rare-af", type=float, default=0.01)
    p.add_argument("--output", required=True)
    return p.parse_args()


def split_genes(value):
    text = str(value).strip()
    if text in MISSING:
        return []
    return [x for x in re.split(r"[;,|\s]+", text) if x and x not in MISSING]


def yes(series):
    return series.fillna("").astype(str).str.upper().isin(["YES", "TRUE", "1"])


def join_values(series):
    values = []
    for value in series.fillna("").astype(str):
        value = value.strip()
        if value in MISSING:
            continue
        for token in re.split(r"[;|]+", value):
            token = token.strip()
            if token and token not in MISSING:
                values.append(token)
    return ";".join(sorted(set(values))) if values else "."


def main():
    args = parse_args()
    df = read_tsv(args.integrated)

    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "gene"])
    id_col = first_existing(df, ["SV_ID", "ID"])
    if gene_col is None or id_col is None:
        raise ValueError("Integrated table requires gene and SV_ID columns.")

    work = df.copy()
    work["_gene_list"] = work[gene_col].apply(split_genes)
    work = work.explode("_gene_list")
    work = work[work["_gene_list"].notna() & work["_gene_list"].ne("")].copy()
    work["gene"] = work["_gene_list"].astype(str)

    caller_count_col = first_existing(work, ["CALLER_COUNT", "SUPP"])
    af_col = first_existing(work, ["NEEDLR_AF"])
    pheno_col = first_existing(work, ["PHENOTYPE_SCORE"])
    panel_col = first_existing(work, ["PANEL_STATUS"])
    candidate_col = first_existing(work, ["CANDIDATE_CLASS"])

    work["_caller_count"] = numeric(work[caller_count_col]).fillna(0) if caller_count_col else 0
    work["_af"] = numeric(work[af_col]) if af_col else np.nan
    work["_phenotype"] = numeric(work[pheno_col]).fillna(0) if pheno_col else 0
    work["_panel"] = (work[panel_col].fillna("").astype(str).str.upper().str.contains("PANEL_GENE|^YES$", regex=True)) if panel_col else False

    for source in ["STRAGLR_MATCH", "TLDR_MATCH", "LONGPHASE_MATCH", "LONGPHASE_PHASED"]:
        if source in work.columns:
            work["_" + source.lower()] = yes(work[source])
        else:
            work["_" + source.lower()] = False

    rows = []
    for gene, group in work.groupby("gene", sort=False):
        unique = group.drop_duplicates(id_col)
        af = unique["_af"] if "_af" in unique else pd.Series(dtype=float)
        rows.append(
            {
                "gene": gene,
                "master_SV_count": int(unique[id_col].nunique()),
                "multicaller_SV_count": int((unique["_caller_count"] >= 2).sum()),
                "needLR_rare_SV_count": int(((af <= args.rare_af) & af.notna()).sum()) if len(af) else 0,
                "needLR_unobserved_SV_count": int((af.eq(0) & af.notna()).sum()) if len(af) else 0,
                "straglr_matched_SV_count": int(unique["_straglr_match"].sum()),
                "tldr_matched_SV_count": int(unique["_tldr_match"].sum()),
                "longphase_matched_SV_count": int(unique["_longphase_match"].sum()),
                "longphase_phased_SV_count": int(unique["_longphase_phased"].sum()),
                "max_phenotype_score": float(unique["_phenotype"].max()) if len(unique) else 0.0,
                "panel_gene": "YES" if bool(unique["_panel"].any()) else "NO",
                "candidate_class": join_values(group[candidate_col]) if candidate_col else ".",
            }
        )

    summary = pd.DataFrame(rows)

    if args.ranking:
        rank = read_tsv(args.ranking)
        rgene = first_existing(rank, ["gene", "Gene", "GENE"])
        if rgene:
            rank = rank.rename(columns={rgene: "gene"})
            keep = [
                c for c in [
                    "gene", "SV_count", "human_HPO_count",
                    "optic_neuropathy_anchor_HPO_count", "phenotype_score",
                    "SV_evidence_score", "integrated_discovery_score",
                    "AnnotSV_OMIM_evidence", "AnnotSV_GENCC_evidence",
                    "AnnotSV_ClinVar_evidence", "AnnotSV_constraint_evidence",
                    "candidate_group", "classification", "interpretation",
                ]
                if c in rank.columns
            ]
            summary = summary.merge(rank[keep].drop_duplicates("gene"), on="gene", how="left")

    if args.phenotypes:
        hpo = read_tsv(args.phenotypes)
        hgene = first_existing(hpo, ["gene_symbol", "gene", "Gene"])
        hpo_id = first_existing(hpo, ["hpo_id", "HPO_ID", "HPO"])
        anchor = first_existing(hpo, ["optic_neuropathy_anchor"])
        if hgene and hpo_id:
            hpo = hpo[hpo[hpo_id].fillna("").astype(str).str.startswith("HP:")].copy()
            hpo["_anchor"] = yes(hpo[anchor]) if anchor else False
            hsum = (
                hpo.groupby(hgene)[hpo_id]
                .nunique()
                .rename("retrieved_HPO_count")
                .reset_index()
                .rename(columns={hgene: "gene"})
            )
            anchor_sum = (
                hpo[hpo["_anchor"]]
                .groupby(hgene)[hpo_id]
                .nunique()
                .rename("retrieved_anchor_HPO_count")
                .reset_index()
                .rename(columns={hgene: "gene"})
            )
            hsum = hsum.merge(anchor_sum, on="gene", how="left")
            hsum["retrieved_anchor_HPO_count"] = hsum["retrieved_anchor_HPO_count"].fillna(0)
            summary = summary.merge(hsum, on="gene", how="left")

    sort_cols = [c for c in ["integrated_discovery_score", "max_phenotype_score", "master_SV_count"] if c in summary.columns]
    if sort_cols:
        for col in sort_cols:
            summary[col] = pd.to_numeric(summary[col], errors="coerce")
        summary = summary.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, sep="\t", index=False)
    print(f"[OK] genes={len(summary)} output={output}")


if __name__ == "__main__":
    main()
