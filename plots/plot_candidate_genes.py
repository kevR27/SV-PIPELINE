#!/usr/bin/env python3
"""Plot genome-wide candidate gene prioritization from ranked_candidates.tsv."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import add_panel_label, first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description="Plot prioritized SV-associated genes.")
    p.add_argument("--input", required=True, help="*_ranked_candidates.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--title", default="Genome-wide SV-associated gene prioritization")
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Ranked candidate table is empty.")

    gene_col = first_existing(df, ["gene", "Gene", "GENE", "SYMBOL"])
    score_col = first_existing(df, ["integrated_discovery_score", "discovery_score", "phenotype_score", "PHENOTYPE_SCORE"])
    sv_count_col = first_existing(df, ["SV_count", "sv_count"])
    anchor_col = first_existing(df, ["optic_neuropathy_anchor_HPO_count", "anchor_HPO_count"])
    panel_col = first_existing(df, ["panel_gene", "PANEL_STATUS"])
    class_col = first_existing(df, ["candidate_group", "classification", "CANDIDATE_CLASS"])
    if gene_col is None or score_col is None:
        raise ValueError("Could not infer gene and prioritization score columns.")

    work = df.copy()
    work["_score"] = numeric(work[score_col]).fillna(0)
    work["_sv_count"] = numeric(work[sv_count_col]).fillna(0) if sv_count_col else 0
    work["_anchor"] = numeric(work[anchor_col]).fillna(0) if anchor_col else 0
    if panel_col:
        ptxt = work[panel_col].fillna("").astype(str).str.upper()
        work["_panel"] = ptxt.str.contains("YES|PANEL_GENE", regex=True)
    else:
        work["_panel"] = False

    work = work.sort_values(["_score", "_anchor", "_sv_count"], ascending=False).head(args.top_n).copy()
    work = work.iloc[::-1].reset_index(drop=True)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work.to_csv(prefix.with_name(prefix.name + "_top_candidates.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, max(6.5, 0.31 * len(work) + 2.2)), gridspec_kw={"width_ratios": [1.65, 1.0]})
    ax1, ax2 = axes
    y = np.arange(len(work))

    colors = np.where(work["_panel"], "#0072B2", "#E69F00")
    ax1.hlines(y, 0, work["_score"], color="#D0D0D0", linewidth=1.1)
    sizes = 40 + np.sqrt(work["_sv_count"].clip(lower=0)) * 35
    ax1.scatter(work["_score"], y, s=sizes, c=colors, edgecolor="white", linewidth=0.5, zorder=3)
    ax1.set_yticks(y)
    ax1.set_yticklabels(work[gene_col])
    ax1.set_xlabel(score_col.replace("_", " "))
    ax1.set_ylabel("Gene")
    style_axis(ax1, "x")
    add_panel_label(ax1, "A")

    x_sv = work["_sv_count"].values
    x_anchor = work["_anchor"].values
    ax2.scatter(x_sv, x_anchor, s=sizes, c=colors, edgecolor="white", linewidth=0.5)
    for i, gene in enumerate(work[gene_col]):
        ax2.annotate(str(gene), (x_sv[i], x_anchor[i]), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax2.set_xlabel("SV count")
    ax2.set_ylabel("Optic-neuropathy anchor HPO count")
    style_axis(ax2, "both")
    add_panel_label(ax2, "B")

    if class_col:
        groups = work[class_col].value_counts().to_dict()
        group_text = "Candidate groups: " + "; ".join(f"{k}={v}" for k, v in groups.items())
    else:
        group_text = ""

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.01, "Blue = panel gene; orange = non-panel gene. Bubble area scales with the number of intersecting SVs. Discovery score is prioritization, not pathogenicity.", ha="center", fontsize=8.2)
    if group_text:
        fig.text(0.5, 0.027, group_text, ha="center", fontsize=7.5)
    fig.tight_layout(rect=[0, 0.045, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
