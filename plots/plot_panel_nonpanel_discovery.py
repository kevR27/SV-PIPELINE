#!/usr/bin/env python3
"""Summarize panel versus non-panel SV-gene discovery evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


CATEGORY_ORDER = [
    "PANEL_GENE",
    "NONPANEL_HPO_AND_DISEASE",
    "NONPANEL_HPO_ONLY",
    "NONPANEL_DISEASE_ONLY",
    "NONPANEL_OTHER",
    "UNRESOLVED_GENE",
]


def parse_args():
    p = argparse.ArgumentParser(description="Plot panel and non-panel SV-gene discovery categories.")
    p.add_argument("--input", required=True, help="Integrated/multimodal SV-gene TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="Panel and non-panel SV-gene discovery")
    return p.parse_args()


def classify(row):
    gene = str(row["_gene"]).strip()
    if gene in {"", ".", "NA", "N/A", "nan", "None"}:
        return "UNRESOLVED_GENE"
    if bool(row["_panel"]):
        return "PANEL_GENE"
    # Generic HON semantic similarity is intentionally thresholded rather
    # than treating any tiny ontology overlap as meaningful phenotype support.
    # 0.25 corresponds to the current SUPPORTING boundary in ranking_common.py.
    phenotype_supported = row["_phenotype"] >= 0.25
    if phenotype_supported and row["_disease"] > 0:
        return "NONPANEL_HPO_AND_DISEASE"
    if phenotype_supported:
        return "NONPANEL_HPO_ONLY"
    if row["_disease"] > 0:
        return "NONPANEL_DISEASE_ONLY"
    return "NONPANEL_OTHER"


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Integrated SV-gene table is empty.")

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    panel_col = first_existing(df, ["PANEL_STATUS", "panel_gene"])
    pheno_col = first_existing(
        df,
        [
            "HPO_QUERY_RESNIK_NORMALIZED",
            "HON_SEMANTIC_SIMILARITY_NORMALIZED",
            "PHENOTYPE_SCORE",
            "phenotype_score",
        ],
    )
    disease_col = first_existing(df, ["GENE_DISEASE_EVIDENCE_SCORE", "gene_disease_evidence_score"])
    score_col = first_existing(
        df,
        [
            "FINAL_GENE_RELEVANCE_DISPLAY_SCORE",
            "GENE_RELEVANCE_DISPLAY_SCORE",
            "EVENT_GENE_RELEVANCE_SCORE",
            "INTEGRATED_DISCOVERY_SCORE",
            "integrated_discovery_score",
            "PHENOTYPE_SCORE",
            "ALLELE_RESEARCH_SCORE",
        ],
    )
    if id_col is None or gene_col is None:
        raise ValueError("Input needs SV_ID and overlapping-gene columns.")

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    work["_phenotype"] = numeric(work[pheno_col]).fillna(0) if pheno_col else 0
    work["_disease"] = numeric(work[disease_col]).fillna(0) if disease_col else 0
    work["_priority"] = numeric(work[score_col]).fillna(0) if score_col else 0
    if panel_col:
        ptxt = work[panel_col].fillna("").astype(str).str.upper().str.strip()
        work["_panel"] = ptxt.isin(["PANEL_GENE", "YES", "TRUE", "1"])
    else:
        work["_panel"] = False

    work = work.drop_duplicates([id_col, gene_col], keep="first").copy()
    work["DISCOVERY_CATEGORY"] = work.apply(classify, axis=1)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work[
        [id_col, gene_col, "DISCOVERY_CATEGORY", "_phenotype", "_disease", "_priority"]
    ].to_csv(
        prefix.with_name(prefix.name + "_classified_pairs.tsv"),
        sep="\t",
        index=False,
    )

    rows = []
    for category in CATEGORY_ORDER:
        subset = work[work["DISCOVERY_CATEGORY"].eq(category)]
        rows.append(
            {
                "category": category,
                "SV_gene_pairs": len(subset),
                "unique_SVs": subset[id_col].nunique(),
                "unique_genes": subset["_gene"].replace(".", np.nan).dropna().nunique(),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)

    plot = summary[(summary[["SV_gene_pairs", "unique_SVs", "unique_genes"]].sum(axis=1) > 0)].copy()
    labels = [
        {
            "PANEL_GENE": "Panel gene",
            "NONPANEL_HPO_AND_DISEASE": "Non-panel: HPO + disease",
            "NONPANEL_HPO_ONLY": "Non-panel: HPO only",
            "NONPANEL_DISEASE_ONLY": "Non-panel: disease only",
            "NONPANEL_OTHER": "Non-panel: other",
            "UNRESOLVED_GENE": "Unresolved gene",
        }[x]
        for x in plot["category"]
    ]

    label_map = dict(zip(plot["category"], labels))
    prioritized_categories = [
        "PANEL_GENE",
        "NONPANEL_HPO_AND_DISEASE",
        "NONPANEL_HPO_ONLY",
        "NONPANEL_DISEASE_ONLY",
    ]
    focused = plot[plot["category"].isin(prioritized_categories)].copy()
    background = plot[~plot["category"].isin(prioritized_categories)].copy()

    fig, axes = plt.subplots(1, 2, figsize=(14.6, 7.2))
    ax1, ax2 = axes

    y1 = np.arange(len(focused))
    width = 0.36
    ax1.barh(y1 - width / 2, focused["unique_SVs"], height=width, label="Unique SVs")
    ax1.barh(y1 + width / 2, focused["unique_genes"], height=width, label="Unique genes")
    ax1.set_yticks(y1)
    ax1.set_yticklabels([label_map[x] for x in focused["category"]])
    ax1.set_xlabel("Count")
    ax1.set_title("Prioritized SV-gene discovery categories")
    ax1.legend(frameon=False)
    style_axis(ax1, "x")
    for container in ax1.containers:
        ax1.bar_label(container, fmt="%d", padding=3, fontsize=8)

    total_pairs = max(int(summary["SV_gene_pairs"].sum()), 1)
    background = background.copy()
    background["pair_percent"] = 100.0 * background["SV_gene_pairs"] / total_pairs
    y2 = np.arange(len(background))
    bars = ax2.barh(y2, background["pair_percent"])
    ax2.set_yticks(y2)
    ax2.set_yticklabels([label_map[x] for x in background["category"]])
    ax2.set_xlabel("Share of all SV-gene rows (%)")
    ax2.set_title("Background / unresolved annotation context")
    style_axis(ax2, "x")
    for bar, (_, row) in zip(bars, background.iterrows()):
        ax2.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2,
            f"  {row['pair_percent']:.1f}%  (n={int(row['SV_gene_pairs']):,})",
            va="center",
            fontsize=8.5,
        )

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.012,
        "Left: categories with phenotype and/or curated disease relevance are shown on a linear count scale. Right: the large background/unresolved classes are shown as proportions so they do not visually suppress the prioritized categories. Categories do not establish causality.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.055, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
