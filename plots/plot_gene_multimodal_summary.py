#!/usr/bin/env python3
"""Plot a gene-level multimodal evidence overview."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description="Plot gene-centric multimodal evidence.")
    p.add_argument("--input", required=True, help="*_gene_multimodal_evidence_summary.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--title", default="Integrated multimodal evidence by candidate gene")
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Gene evidence summary is empty.")

    gene_col = first_existing(df, ["gene", "Gene"])
    score_col = first_existing(df, ["integrated_discovery_score", "max_phenotype_score"])
    if gene_col is None:
        raise ValueError("Gene summary needs a gene column.")

    work = df.copy()
    if score_col:
        work["_score"] = numeric(work[score_col]).fillna(0)
    else:
        work["_score"] = numeric(work.get("master_SV_count", pd.Series(0, index=work.index))).fillna(0)

    anchor_col = first_existing(work, ["retrieved_anchor_HPO_count", "optic_neuropathy_anchor_HPO_count"])
    if anchor_col:
        work["_anchor"] = numeric(work[anchor_col]).fillna(0)
    else:
        work["_anchor"] = 0

    work = work.sort_values(["_score", "_anchor"], ascending=False).head(args.top_n).copy()
    work = work.iloc[::-1].reset_index(drop=True)

    evidence_specs = [
        ("multicaller_SV_count", "Multicaller SV"),
        ("needLR_rare_SV_count", "Rare needLR SV"),
        ("needLR_unobserved_SV_count", "needLR AF=0"),
        ("straglr_matched_SV_count", "Straglr overlap"),
        ("tldr_matched_SV_count", "TLDR overlap"),
        ("longphase_phased_SV_count", "LongPhase phased"),
    ]

    matrix_cols = []
    labels = []
    for column, label in evidence_specs:
        if column in work.columns:
            matrix_cols.append((numeric(work[column]).fillna(0) > 0).astype(int))
            labels.append(label)

    if anchor_col:
        matrix_cols.append((work["_anchor"] > 0).astype(int))
        labels.append("Anchor HPO")

    if "panel_gene" in work.columns:
        matrix_cols.append(
            work["panel_gene"].fillna("").astype(str).str.upper().isin(["YES", "TRUE", "1"]).astype(int)
        )
        labels.append("Panel gene")

    if not matrix_cols:
        raise ValueError("No multimodal evidence columns were detected.")

    matrix = np.column_stack([series.to_numpy() for series in matrix_cols])
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    source = pd.DataFrame(matrix, columns=labels)
    source.insert(0, "gene", work[gene_col].values)
    source.to_csv(prefix.with_name(prefix.name + "_matrix.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(
        1, 2,
        figsize=(14.2, max(7.0, 0.38 * len(work) + 2.4)),
        gridspec_kw={"width_ratios": [0.8, 1.8]},
    )
    ax1, ax2 = axes
    y = np.arange(len(work))

    ax1.barh(y, work["_score"], color="#0072B2")
    ax1.set_yticks(y)
    ax1.set_yticklabels(work[gene_col], fontsize=9.5)
    ax1.set_xlabel(score_col.replace("_", " ") if score_col else "Priority")
    ax1.set_ylabel("Gene")
    style_axis(ax1, "x")

    cmap = ListedColormap(["#F3F4F4", "#0B6E69"])
    ax2.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    ax2.set_xticks(np.arange(len(labels)))
    ax2.set_xticklabels(labels, rotation=38, ha="right", fontsize=9)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("Evidence layer")
    ax2.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, len(work), 1), minor=True)
    ax2.grid(which="minor", color="white", linewidth=1.0)
    ax2.tick_params(which="minor", bottom=False, left=False)

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.008,
        "Dark teal indicates that at least one SV/evidence item for that gene satisfies the indicated layer; this is prioritization evidence, not pathogenicity.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
