#!/usr/bin/env python3
"""Show the SV size/type spectrum contributing to prioritized genes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


SIZE_GROUPS = [
    ("LT100KB", "<100 kb"),
    ("100KB_1MB", "100 kb–1 Mb"),
    ("1MB_10MB", "1–10 Mb"),
    ("GE10MB", "≥10 Mb"),
    ("BREAKEND", "Breakend"),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=20)
    return p.parse_args()


def group_size(value):
    text = str(value)
    if text.startswith("BREAKEND"):
        return "BREAKEND"
    if text in {"LT50_BP", "50BP_1KB", "1KB_10KB", "10KB_100KB"}:
        return "LT100KB"
    if text == "100KB_1MB":
        return "100KB_1MB"
    if text == "1MB_10MB":
        return "1MB_10MB"
    if text == "GE_10MB":
        return "GE10MB"
    return "OTHER"


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    score_col = first_existing(df, ["EVENT_GENE_RELEVANCE_SCORE", "INTEGRATED_DISCOVERY_SCORE"])
    if id_col is None or gene_col is None:
        raise ValueError("Event table needs SV_ID and gene columns.")
    if "SV_EVENT_SIZE_CLASS" not in df:
        raise ValueError("Input needs SV_EVENT_SIZE_CLASS from rank_sv_gene_events.py.")

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    work = work[~work["_gene"].isin(["", ".", "NA", "N/A", "nan", "None"])].copy()
    work["_score"] = numeric(work[score_col]).fillna(0) if score_col else 0
    work["_size_group"] = work["SV_EVENT_SIZE_CLASS"].map(group_size)
    work["_breakpoint"] = (
        work.get("BREAKPOINT_DEFINED_EVENT", "NO")
        .fillna("NO")
        .astype(str)
        .str.upper()
        .eq("YES")
    )

    unique = (
        work.sort_values("_score", ascending=False)
        .drop_duplicates([id_col, gene_col], keep="first")
        .copy()
    )

    gene_priority = (
        unique.groupby("_gene")["_score"]
        .max()
        .sort_values(ascending=False)
        .head(args.top_n)
        .index
        .tolist()
    )
    unique = unique[unique["_gene"].isin(gene_priority)].copy()

    rows = []
    for gene in gene_priority:
        sub = unique[unique["_gene"].eq(gene)]
        row = {
            "gene": gene,
            "gene_relevance_score": float(sub["_score"].max()),
            "unique_SVs": int(sub[id_col].nunique()),
            "breakpoint_defined_INV_BND": int(sub.loc[sub["_breakpoint"], id_col].nunique()),
        }
        for key, _label in SIZE_GROUPS:
            row[key] = int(sub.loc[sub["_size_group"].eq(key), id_col].nunique())
        rows.append(row)

    summary = pd.DataFrame(rows)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(
        prefix.with_name(prefix.name + "_summary.tsv"),
        sep="\t",
        index=False,
    )

    plot = summary.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15.5, max(7.5, 0.39 * len(plot) + 2.8)),
        gridspec_kw={"width_ratios": [1.45, 0.75]},
    )
    ax1, ax2 = axes

    left = np.zeros(len(plot))
    for key, label in SIZE_GROUPS:
        values = plot[key].to_numpy()
        ax1.barh(y, values, left=left, label=label)
        left += values

    ax1.set_yticks(y)
    ax1.set_yticklabels(plot["gene"])
    ax1.set_xlabel("Unique SVs overlapping gene")
    ax1.set_ylabel("Prioritized gene")
    ax1.set_title("SV size spectrum per prioritized gene")
    ax1.legend(frameon=False, fontsize=8, ncol=2)
    style_axis(ax1, "x")

    bars = ax2.barh(y, plot["breakpoint_defined_INV_BND"])
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("Unique INV/BND events")
    ax2.set_title("Breakpoint-defined events")
    style_axis(ax2, "x")
    for bar, count in zip(bars, plot["breakpoint_defined_INV_BND"]):
        if count:
            ax2.text(
                bar.get_width(),
                bar.get_y() + bar.get_height() / 2,
                f"  {int(count)}",
                va="center",
                fontsize=8,
            )

    fig.suptitle(
        "Structural-variant spectrum of prioritized genes",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.012,
        "Counts are unique master SVs, not annotation rows. Size/type composition is descriptive and does not itself increase pathogenicity priority. INV/BND counts identify breakpoint-defined events requiring breakpoint-aware interpretation.",
        ha="center",
        fontsize=8.8,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.965])
    outputs = save_figure(fig, prefix)
    plt.close(fig)

    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
