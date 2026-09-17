#!/usr/bin/env python3
"""Plot the genome-wide LRS structural-variant landscape from the integrated table."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import (
    SUPPORT_COLORS,
    SVTYPE_COLORS,
    abs_svlen,
    add_panel_label,
    chromosome_sort_key,
    first_existing,
    normalize_svtype,
    read_tsv,
    save_figure,
    set_thesis_style,
    style_axis,
    unique_master_svs,
)


def parse_args():
    p = argparse.ArgumentParser(description="Plot the master LRS SV landscape.")
    p.add_argument("--input", required=True, help="*_integrated_SV_gene_analysis.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="Genome-wide structural-variant landscape")
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()

    raw = read_tsv(args.input)
    df = unique_master_svs(raw)
    if df.empty:
        raise ValueError("Integrated table contains no SVs.")

    svtype_col = first_existing(df, ["SVTYPE", "SV_type", "Type"])
    chrom_col = first_existing(df, ["CHROM", "Chr", "chrom"])
    if svtype_col is None or chrom_col is None:
        raise ValueError("Input needs SVTYPE and CHROM columns.")

    df["SVTYPE_NORM"] = normalize_svtype(df[svtype_col])
    df["ABS_SVLEN"] = abs_svlen(df)
    df["LOG10_SVLEN"] = np.log10(df["ABS_SVLEN"].where(df["ABS_SVLEN"] > 0))

    support_col = first_existing(df, ["CALLER_SUPPORT_CLASS"])
    if support_col is None:
        count_col = first_existing(df, ["CALLER_COUNT", "SUPP"])
        if count_col:
            count = pd.to_numeric(df[count_col], errors="coerce")
            df["CALLER_SUPPORT_CLASS"] = np.select(
                [count >= 3, count == 2, count == 1],
                ["MULTICALLER", "TWO_CALLER", "SINGLE_CALLER"],
                default="UNKNOWN",
            )
        else:
            df["CALLER_SUPPORT_CLASS"] = "UNKNOWN"
        support_col = "CALLER_SUPPORT_CLASS"

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    type_counts = df["SVTYPE_NORM"].value_counts().rename_axis("SVTYPE").reset_index(name="count")
    type_counts.to_csv(prefix.with_name(prefix.name + "_svtype_counts.tsv"), sep="\t", index=False)

    chr_counts = df.groupby(chrom_col).size().rename("count").reset_index()
    chr_counts = chr_counts.sort_values(chrom_col, key=lambda s: s.map(chromosome_sort_key))
    chr_counts.to_csv(prefix.with_name(prefix.name + "_chromosome_counts.tsv"), sep="\t", index=False)

    support_counts = (
        df.groupby(["SVTYPE_NORM", support_col]).size().rename("count").reset_index()
    )
    support_counts.to_csv(prefix.with_name(prefix.name + "_support_by_svtype.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.7))
    ax1, ax2, ax3, ax4 = axes.flatten()

    order = [x for x in ["DEL", "INS", "DUP", "INV", "BND", "CNV", "OTHER"] if x in set(type_counts["SVTYPE"])]
    vals = type_counts.set_index("SVTYPE").reindex(order)["count"].fillna(0)
    ax1.bar(order, vals, color=[SVTYPE_COLORS.get(x, "#999999") for x in order])
    ax1.set_ylabel("Number of unique master SVs")
    ax1.set_xlabel("SV type")
    style_axis(ax1, "y")
    add_panel_label(ax1, "A")
    for i, value in enumerate(vals):
        ax1.text(i, value, f"{int(value):,}", ha="center", va="bottom", fontsize=8)

    for svtype in order:
        x = df.loc[df["SVTYPE_NORM"].eq(svtype), "LOG10_SVLEN"].dropna()
        if x.empty:
            continue
        ax2.hist(
            x,
            bins=np.linspace(max(1.5, df["LOG10_SVLEN"].min()), min(8.5, df["LOG10_SVLEN"].max()), 40),
            histtype="step",
            linewidth=1.7,
            label=svtype,
            color=SVTYPE_COLORS.get(svtype, "#999999"),
        )
    ax2.set_xlabel(r"SV size, $\log_{10}$(bp)")
    ax2.set_ylabel("Number of SVs")
    ax2.legend(frameon=False, ncol=2)
    style_axis(ax2, "y")
    add_panel_label(ax2, "B")

    ax3.bar(chr_counts[chrom_col], chr_counts["count"], color="#6E6E6E")
    ax3.set_xlabel("Chromosome")
    ax3.set_ylabel("Number of unique master SVs")
    ax3.tick_params(axis="x", rotation=55)
    style_axis(ax3, "y")
    add_panel_label(ax3, "C")

    support_order = ["SINGLE_CALLER", "TWO_CALLER", "MULTICALLER", "UNKNOWN"]
    pivot = support_counts.pivot(index="SVTYPE_NORM", columns=support_col, values="count").fillna(0)
    pivot = pivot.reindex(order).fillna(0)
    bottom = np.zeros(len(pivot))
    for category in support_order:
        if category not in pivot.columns:
            continue
        values = pivot[category].values
        ax4.bar(
            pivot.index,
            values,
            bottom=bottom,
            label=category.replace("_", " ").title(),
            color=SUPPORT_COLORS.get(category, "#999999"),
        )
        bottom += values
    ax4.set_xlabel("SV type")
    ax4.set_ylabel("Number of unique master SVs")
    ax4.legend(frameon=False, fontsize=8)
    style_axis(ax4, "y")
    add_panel_label(ax4, "D")

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)

    print(f"[OK] unique_master_SVs={len(df)}")
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
