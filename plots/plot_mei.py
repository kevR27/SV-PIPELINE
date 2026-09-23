#!/usr/bin/env python3
"""Plot TLDR mobile-element insertion calls for one sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description="Plot TLDR mobile-element insertion evidence.")
    p.add_argument("--input", required=True, help="TLDR *.table.txt")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--title", default="Mobile-element insertion evidence")
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("TLDR table is empty.")

    # Exact TLDR output schema currently produced by the pipeline is supported first.
    family_col = first_existing(df, ["Family", "family", "FAMILY", "repeat_family", "REPEAT_FAMILY", "element", "ELEMENT"])
    subfamily_col = first_existing(df, ["Subfamily", "subfamily", "SUBFAMILY", "repeat_name", "REPEAT_NAME"])
    chrom_col = first_existing(df, ["Chrom", "chrom", "CHROM", "chr", "Chr"])
    pos_col = first_existing(df, ["Start", "start", "START", "pos", "POS", "position", "POSITION"])
    support_col = first_existing(df, ["UsedReads", "used_reads", "num_support", "support", "SUPPORT", "supporting_reads", "N_SUPPORT", "reads"])
    spanning_col = first_existing(df, ["SpanReads", "span_reads", "num_spanning", "spanning_reads", "SPANNING", "N_SPANNING"])
    len_col = first_existing(df, ["LengthIns", "length_ins", "insert_len", "insertion_length", "INS_LEN", "length", "LENGTH"])
    filter_col = first_existing(df, ["Filter", "filter", "FILTER", "status", "STATUS"])

    work = df.copy()
    if family_col is None:
        if subfamily_col:
            work["_family"] = (
                work[subfamily_col]
                .fillna("Unknown")
                .astype(str)
                .str.extract(r"^([A-Za-z0-9]+)", expand=False)
                .fillna("Unknown")
            )
        else:
            work["_family"] = "Unknown"
    else:
        work["_family"] = work[family_col].fillna("Unknown").astype(str)

    work["_support"] = numeric(work[support_col]).fillna(0) if support_col else 0
    work["_spanning"] = numeric(work[spanning_col]).fillna(0) if spanning_col else 0
    work["_length"] = numeric(work[len_col]).abs() if len_col else np.nan

    work["_label"] = np.arange(len(work)).astype(str)
    if chrom_col and pos_col:
        work["_label"] = work[chrom_col].astype(str) + ":" + work[pos_col].astype(str)
    if subfamily_col:
        work["_label"] = work["_label"] + " | " + work[subfamily_col].fillna("").astype(str)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    family_summary = (
        work.groupby("_family")
        .agg(
            n_insertions=("_family", "size"),
            mean_support=("_support", "mean"),
            median_support=("_support", "median"),
        )
        .reset_index()
        .sort_values("n_insertions", ascending=False)
    )
    family_summary.to_csv(prefix.with_name(prefix.name + "_family_summary.tsv"), sep="\t", index=False)

    top = work.sort_values(["_support", "_spanning"], ascending=False).head(args.top_n).copy().iloc[::-1]
    top.to_csv(prefix.with_name(prefix.name + "_top_insertions.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.3, max(6.0, 0.30 * len(top) + 2.0)))
    ax1, ax2 = axes

    fam = family_summary.iloc[::-1]
    ax1.barh(fam["_family"], fam["n_insertions"], color="#009E73")
    ax1.set_xlabel("Number of TLDR insertions")
    ax1.set_ylabel("Mobile-element family")
    style_axis(ax1, "x")

    y = np.arange(len(top))
    ax2.barh(y, top["_support"], color="#0072B2", label="Used reads")
    if top["_spanning"].gt(0).any():
        ax2.scatter(top["_spanning"], y, marker="D", s=28, color="#E69F00", label="Spanning reads", zorder=3)
    ax2.set_yticks(y)
    ax2.set_yticklabels(top["_label"], fontsize=7.2)
    ax2.set_xlabel("Read evidence")
    ax2.set_ylabel("Insertion")
    ax2.legend(frameon=False, fontsize=8)
    style_axis(ax2, "x")

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    note = (
        "TLDR calls are an orthogonal MEI layer. Candidate overlap with Jasmine INS calls "
        "is established by breakpoint-aware matching, not gene name alone."
    )
    if filter_col:
        note += f" TLDR filter/status column: {filter_col}."
    fig.text(0.5, 0.01, note, ha="center", fontsize=8.0)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
