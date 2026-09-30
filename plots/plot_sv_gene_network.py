#!/usr/bin/env python3
"""Bipartite map of prioritized merged SVs and their overlapping genes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import SVTYPE_COLORS, first_existing, normalize_svtype, numeric, read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description="Plot explicit SV-to-gene relationships.")
    p.add_argument("--input", required=True, help="Integrated/multimodal SV-gene TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n-pairs", type=int, default=35)
    p.add_argument("--title", default="Prioritized structural variants and overlapping genes")
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    type_col = first_existing(df, ["SVTYPE", "SV_type"])
    caller_col = first_existing(df, ["CALLER_COUNT", "SUPP"])
    score_col = first_existing(df, ["EVENT_GENE_RELEVANCE_SCORE", "INTEGRATED_DISCOVERY_SCORE", "integrated_discovery_score", "PHENOTYPE_SCORE", "ALLELE_RESEARCH_SCORE"])
    panel_col = first_existing(df, ["PANEL_STATUS", "panel_gene"])
    if id_col is None or gene_col is None:
        raise ValueError("Input needs SV_ID and overlapping-gene columns.")

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    work = work[~work["_gene"].isin(["", ".", "NA", "N/A", "nan", "None"])].copy()
    work["_priority"] = numeric(work[score_col]).fillna(0) if score_col else 0
    work["_caller_count"] = numeric(work[caller_col]).fillna(0) if caller_col else 0
    work["_svtype"] = normalize_svtype(work[type_col]) if type_col else "OTHER"
    relationship_col = first_existing(work, ["SV_GENE_RELATIONSHIP"])
    work["_interval_context_only"] = (
        work[relationship_col]
        .fillna("")
        .astype(str)
        .isin(["INTERVAL_CONTEXT_ONLY", "INVERSION_SPANS_INTACT_GENE"])
        if relationship_col
        else False
    )
    if panel_col:
        ptxt = work[panel_col].fillna("").astype(str).str.upper().str.strip()
        work["_panel"] = ptxt.isin(["PANEL_GENE", "YES", "TRUE", "1"])
    else:
        work["_panel"] = False

    work = (
        work.sort_values(
            ["_interval_context_only", "_priority", "_caller_count", id_col, "_gene"],
            ascending=[True, False, False, True, True],
        )
        .drop_duplicates([id_col, gene_col], keep="first")
        .head(args.top_n_pairs)
        .copy()
    )
    if work.empty:
        raise ValueError("No SV-gene pairs with resolved genes were available.")

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work[
        [id_col, gene_col, "_priority", "_caller_count", "_svtype", "_panel"]
    ].to_csv(prefix.with_name(prefix.name + "_edges.tsv"), sep="\t", index=False)

    sv_order = (
        work.groupby(id_col)["_priority"].max().sort_values(ascending=False).index.astype(str).tolist()
    )
    gene_order = (
        work.groupby("_gene")["_priority"].max().sort_values(ascending=False).index.astype(str).tolist()
    )
    sv_positions = np.linspace(len(sv_order) - 1, 0, len(sv_order))
    gene_positions = np.linspace(len(gene_order) - 1, 0, len(gene_order))
    sv_y = dict(zip(sv_order, sv_positions))
    gene_y = dict(zip(gene_order, gene_positions))

    height = max(7.5, 0.36 * max(len(sv_order), len(gene_order)) + 2.5)
    fig, ax = plt.subplots(figsize=(15.5, height))

    for _, row in work.iterrows():
        sv = str(row[id_col])
        gene = str(row["_gene"])
        width = 0.7 + 0.45 * min(float(row["_caller_count"]), 3.0)
        interval_only = bool(row["_interval_context_only"])
        ax.plot(
            [0.2, 0.8],
            [sv_y[sv], gene_y[gene]],
            linewidth=width,
            alpha=0.22 if interval_only else 0.42,
            color="#8A949E",
            linestyle="--" if interval_only else "-",
            zorder=1,
        )

    sv_meta = work.drop_duplicates(id_col).set_index(id_col)
    for sv in sv_order:
        row = sv_meta.loc[sv]
        ax.scatter(
            [0.2], [sv_y[sv]], s=95,
            color=SVTYPE_COLORS.get(str(row["_svtype"]), "#999999"),
            edgecolor="white", linewidth=0.6, zorder=3,
        )
        ax.text(0.185, sv_y[sv], sv, ha="right", va="center", fontsize=8.5)

    gene_meta = work.sort_values("_priority", ascending=False).drop_duplicates("_gene").set_index("_gene")
    for gene in gene_order:
        row = gene_meta.loc[gene]
        ax.scatter(
            [0.8], [gene_y[gene]], s=95,
            color="#0072B2" if bool(row["_panel"]) else "#E69F00",
            edgecolor="white", linewidth=0.6, zorder=3,
        )
        ax.text(0.815, gene_y[gene], gene, ha="left", va="center", fontsize=9)

    ax.text(0.2, max(list(sv_y.values())) + 1.0, "Merged SV", ha="center", fontweight="bold")
    ax.text(0.8, max(list(gene_y.values())) + 1.0, "Overlapping gene", ha="center", fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, max(max(sv_y.values()), max(gene_y.values())) + 1.5)
    ax.axis("off")

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.012,
        "Each line is one explicit (SV_ID, gene) association. Solid edges indicate dosage or direct/proximal breakpoint context; dashed edges indicate inversion-spanned or rearrangement interval context without a direct transcript breakpoint. Edge width scales with caller count. Gene nodes: blue=panel, orange=non-panel.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.965])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
