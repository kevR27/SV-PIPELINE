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
    p.add_argument("--top-n", type=int, default=15)
    p.add_argument("--include-nonpass", action="store_true", help="Include non-PASS TLDR records in the main figure")
    p.add_argument("--title", default="Mobile-element insertion evidence")
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("TLDR table is empty.")

    family_col = first_existing(df, ["Family", "family", "FAMILY", "repeat_family", "element"])
    subfamily_col = first_existing(df, ["Subfamily", "subfamily", "SUBFAMILY", "repeat_name"])
    chrom_col = first_existing(df, ["Chrom", "chrom", "CHROM", "chr", "Chr"])
    pos_col = first_existing(df, ["Start", "start", "START", "pos", "POS"])
    support_col = first_existing(df, ["UsedReads", "used_reads", "support", "SUPPORT", "supporting_reads"])
    spanning_col = first_existing(df, ["SpanReads", "span_reads", "spanning_reads"])
    len_col = first_existing(df, ["LengthIns", "length_ins", "insertion_length", "length"])
    filter_col = first_existing(df, ["Filter", "filter", "FILTER", "status"])

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    if filter_col:
        filter_summary = (
            df[filter_col].fillna("NA").astype(str)
            .value_counts()
            .rename_axis("filter")
            .reset_index(name="count")
        )
        filter_summary.to_csv(prefix.with_name(prefix.name + "_filter_summary.tsv"), sep="\t", index=False)

    work = df.copy()
    if filter_col and not args.include_nonpass:
        work = work[work[filter_col].fillna("").astype(str).str.upper().eq("PASS")].copy()

    if work.empty:
        raise ValueError("No TLDR records remain after filter selection.")

    if family_col is None:
        work["_family"] = "Unknown"
    else:
        work["_family"] = work[family_col].fillna("Unknown").replace({"NA": "Unknown"}).astype(str)

    work["_support"] = numeric(work[support_col]).fillna(0) if support_col else 0
    work["_spanning"] = numeric(work[spanning_col]).fillna(0) if spanning_col else 0
    work["_length"] = numeric(work[len_col]).abs() if len_col else np.nan

    work["_label"] = np.arange(len(work)).astype(str)
    if chrom_col and pos_col:
        work["_label"] = work[chrom_col].astype(str) + ":" + work[pos_col].astype(str)
    if subfamily_col:
        sub = work[subfamily_col].fillna("").replace({"NA": ""}).astype(str)
        mask = sub.ne("")
        work.loc[mask, "_label"] = work.loc[mask, "_label"] + " | " + sub.loc[mask]

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

    fig, axes = plt.subplots(1, 2, figsize=(14.0, max(7.0, 0.42 * len(top) + 2.5)), gridspec_kw={"width_ratios": [0.9, 1.4]})
    ax1, ax2 = axes

    fam = family_summary.iloc[::-1]
    ax1.barh(fam["_family"], fam["n_insertions"], color="#0B6E69")
    ax1.set_xlabel("Number of TLDR insertions")
    ax1.set_ylabel("Mobile-element family")
    style_axis(ax1, "x")
    for y, value in enumerate(fam["n_insertions"]):
        ax1.text(value, y, f" {int(value):,}", va="center", fontsize=9)

    y = np.arange(len(top))
    ax2.barh(y, top["_support"], color="#0072B2", label="Used reads")
    if top["_spanning"].gt(0).any():
        ax2.scatter(top["_spanning"], y, marker="D", s=35, color="#E69F00", label="Spanning reads", zorder=3)
    ax2.set_yticks(y)
    ax2.set_yticklabels(top["_label"], fontsize=9)
    ax2.set_xlabel("Read evidence")
    ax2.set_ylabel("Insertion")
    ax2.legend(frameon=False)
    style_axis(ax2, "x")

    status = "all TLDR records" if args.include_nonpass else "PASS TLDR records only"
    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.008,
        f"Main figure uses {status}. TLDR is an orthogonal MEI layer; Jasmine overlap is established by breakpoint-aware matching.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print(f"[OK] plotted_records={len(work)}")
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
