#!/usr/bin/env python3
"""Detailed quantitative view of prioritized SV-gene associations.

One point represents one unique (master SV, overlapping gene) pair.  The plot
keeps technical support, biological prioritization and gene-disease context
separate; none of the displayed scores is interpreted as a pathogenicity
probability.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import (
    SVTYPE_COLORS,
    abs_svlen,
    first_existing,
    normalize_svtype,
    numeric,
    read_tsv,
    save_figure,
    set_thesis_style,
    style_axis,
)


def parse_args():
    p = argparse.ArgumentParser(description="Plot detailed SV-gene evidence relationships.")
    p.add_argument("--input", required=True, help="Integrated/multimodal SV-gene TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=30)
    p.add_argument("--title", default="Detailed SV-gene evidence relationships")
    return p.parse_args()


def present(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.upper()
    return ~text.isin(["", ".", "NA", "N/A", "NONE", "NAN", "UNKNOWN", "NOT_AVAILABLE"])


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Integrated SV-gene table is empty.")

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    type_col = first_existing(df, ["SVTYPE", "SV_type"])
    caller_col = first_existing(df, ["CALLER_COUNT", "SUPP"])
    score_col = first_existing(
        df,
        ["INTEGRATED_DISCOVERY_SCORE", "integrated_discovery_score", "PHENOTYPE_SCORE", "ALLELE_RESEARCH_SCORE"],
    )
    pheno_col = first_existing(df, ["PHENOTYPE_SCORE", "phenotype_score"])
    disease_col = first_existing(df, ["GENE_DISEASE_EVIDENCE_SCORE", "gene_disease_evidence_score"])
    panel_col = first_existing(df, ["PANEL_STATUS", "panel_gene"])
    af_col = first_existing(df, ["NEEDLR_AF"])

    if id_col is None or gene_col is None:
        raise ValueError("Input needs SV_ID and overlapping-gene columns.")

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    work["_caller_count"] = numeric(work[caller_col]).fillna(0) if caller_col else 0
    work["_priority"] = numeric(work[score_col]).fillna(0) if score_col else 0
    work["_phenotype"] = numeric(work[pheno_col]).fillna(0) if pheno_col else 0
    work["_disease"] = numeric(work[disease_col]).fillna(0) if disease_col else 0
    work["_svlen"] = abs_svlen(work).fillna(50).clip(lower=50)
    work["_svtype"] = normalize_svtype(work[type_col]) if type_col else "OTHER"
    work["_af"] = numeric(work[af_col]) if af_col else np.nan

    if panel_col:
        ptxt = work[panel_col].fillna("").astype(str).str.upper().str.strip()
        work["_panel"] = ptxt.isin(["PANEL_GENE", "YES", "TRUE", "1"])
    else:
        work["_panel"] = False

    pathogenic_col = first_existing(work, ["SV_PATHOGENIC_DB_SOURCE"])
    orthogonal_cols = [
        first_existing(work, ["STRAGLR_MATCH"]),
        first_existing(work, ["TLDR_MATCH"]),
        first_existing(work, ["LONGPHASE_MATCH"]),
    ]
    orthogonal_cols = [c for c in orthogonal_cols if c]
    work["_path_db"] = present(work[pathogenic_col]) if pathogenic_col else False
    if orthogonal_cols:
        work["_orthogonal"] = False
        for col in orthogonal_cols:
            work["_orthogonal"] = work["_orthogonal"] | work[col].fillna("").astype(str).str.upper().eq("YES")
    else:
        work["_orthogonal"] = False

    work = (
        work.sort_values(
            ["_priority", "_phenotype", "_disease", "_caller_count", id_col, "_gene"],
            ascending=[False, False, False, False, True, True],
        )
        .drop_duplicates([id_col, gene_col], keep="first")
        .head(args.top_n)
        .copy()
    )

    work["PAIR_LABEL"] = work[id_col].astype(str) + " | " + work["_gene"]
    work["PANEL_GROUP"] = np.where(work["_panel"], "PANEL", "NONPANEL")
    work["PATHOGENIC_DB_REPORTED"] = np.where(work["_path_db"], "YES", "NO")
    work["ORTHOGONAL_MATCH"] = np.where(work["_orthogonal"], "YES", "NO")

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    source_cols = [
        id_col, gene_col, "PAIR_LABEL", "PANEL_GROUP", "_svtype", "_svlen",
        "_caller_count", "_priority", "_phenotype", "_disease", "_af",
        "PATHOGENIC_DB_REPORTED", "ORTHOGONAL_MATCH",
    ]
    work[source_cols].to_csv(
        prefix.with_name(prefix.name + "_sv_gene_pairs.tsv"),
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(1, 2, figsize=(14.8, 7.4))
    ax1, ax2 = axes

    sizes = 35 + 24 * np.log10(work["_svlen"].clip(lower=50))
    edge = np.where(work["_orthogonal"], "black", "white")
    colors = [SVTYPE_COLORS.get(x, "#999999") for x in work["_svtype"]]

    ax1.scatter(
        work["_caller_count"],
        work["_priority"],
        s=sizes,
        c=colors,
        edgecolors=edge,
        linewidths=np.where(work["_orthogonal"], 1.3, 0.5),
        alpha=0.88,
    )
    label_rows = (
        work.sort_values(["_priority", "_phenotype", "_disease", "_caller_count"], ascending=False)
        .drop_duplicates("_gene")
        .head(min(8, len(work)))
    )
    offsets = [(7, 8), (7, -15), (12, 18), (12, -25), (18, 6), (18, -18), (24, 16), (24, -28)]
    for j, (_, row) in enumerate(label_rows.iterrows()):
        ax1.annotate(
            str(row["_gene"]),
            (row["_caller_count"], row["_priority"]),
            xytext=offsets[j % len(offsets)],
            textcoords="offset points",
            fontsize=8,
            arrowprops={"arrowstyle": "-", "linewidth": 0.5, "alpha": 0.55},
        )
    ax1.set_xlabel("Caller count")
    ax1.set_ylabel((score_col or "priority score").replace("_", " "))
    ax1.set_title("Technical support versus candidate priority")
    style_axis(ax1, "both")

    ax2.scatter(
        work["_phenotype"],
        work["_disease"],
        s=sizes,
        c=np.where(work["_panel"], "#0072B2", "#E69F00"),
        edgecolors=np.where(work["_path_db"], "black", "white"),
        linewidths=np.where(work["_path_db"], 1.3, 0.5),
        alpha=0.88,
    )
    for j, (_, row) in enumerate(label_rows.iterrows()):
        ax2.annotate(
            str(row["_gene"]),
            (row["_phenotype"], row["_disease"]),
            xytext=offsets[j % len(offsets)],
            textcoords="offset points",
            fontsize=8,
            arrowprops={"arrowstyle": "-", "linewidth": 0.5, "alpha": 0.55},
        )
    ax2.set_xlabel("Phenotype relevance score")
    ax2.set_ylabel("Gene-disease evidence score")
    ax2.set_title("Phenotype relevance versus known disease evidence")
    style_axis(ax2, "both")

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.012,
        "Each point is one master SV-gene pair. Point size scales with SV length. Left: color denotes SV type and black outline denotes an orthogonal match. Right: blue=panel gene, orange=non-panel gene, black outline=reported pathogenic-SV database evidence. Scores are research-prioritization variables, not pathogenicity probabilities.",
        ha="center",
        fontsize=8.5,
        wrap=True,
    )
    fig.tight_layout(rect=[0, 0.065, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
