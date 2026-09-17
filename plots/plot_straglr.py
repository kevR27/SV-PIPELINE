#!/usr/bin/env python3
"""Plot annotated Straglr tandem-repeat loci."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, parse_semicolon_numeric, read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description="Plot Straglr locus-level repeat evidence.")
    p.add_argument("--input", required=True, help="*_straglr.annotated.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--title", default="Tandem-repeat expansion evidence")
    return p.parse_args()


def max_numeric(series: pd.Series) -> pd.Series:
    return series.apply(lambda x: max(parse_semicolon_numeric(x), default=np.nan))


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Straglr annotation table is empty.")

    locus = first_existing(df, ["locus"])
    gene = first_existing(df, ["overlapping_genes"])
    cn = first_existing(df, ["copy_number"])
    size = first_existing(df, ["size"])
    support = first_existing(df, ["supporting_reads"])
    coverage = first_existing(df, ["coverage"])
    if locus is None:
        raise ValueError("Could not identify Straglr locus column.")

    work = df.copy()
    work["_cn"] = max_numeric(work[cn]) if cn else np.nan
    work["_size"] = max_numeric(work[size]) if size else np.nan
    work["_support"] = pd.to_numeric(work[support], errors="coerce").fillna(0) if support else 0
    work["_coverage"] = pd.to_numeric(work[coverage], errors="coerce") if coverage else np.nan
    work["_label"] = work[locus].astype(str)
    if gene:
        g = work[gene].fillna("").astype(str).str.strip()
        work.loc[g.ne(""), "_label"] = work.loc[g.ne(""), locus].astype(str) + " | " + g[g.ne("")]

    work["_priority"] = work["_support"] + work["_cn"].fillna(0) / max(work["_cn"].max(skipna=True) or 1, 1)
    work = work.sort_values(["_priority", "_cn", "_size"], ascending=False).head(args.top_n).iloc[::-1].copy()

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work.to_csv(prefix.with_name(prefix.name + "_top_loci.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, max(6.2, 0.30 * len(work) + 2.0)))
    ax1, ax2 = axes
    y = np.arange(len(work))

    ax1.barh(y, work["_support"], color="#0072B2")
    ax1.set_yticks(y)
    ax1.set_yticklabels(work["_label"], fontsize=7.5)
    ax1.set_xlabel("Supporting reads")
    ax1.set_ylabel("Repeat locus")
    style_axis(ax1, "x")

    x = work["_cn"] if work["_cn"].notna().any() else work["_size"]
    xlabel = "Maximum copy number" if work["_cn"].notna().any() else "Maximum repeat size"
    ax2.scatter(x, work["_support"], s=45 + work["_coverage"].fillna(0).clip(lower=0) * 2, color="#E69F00", edgecolor="white", linewidth=0.5)
    for _, row in work.iterrows():
        ax2.annotate(str(row["_label"]).split(" | ")[0], (row["_cn"] if pd.notna(row["_cn"]) else row["_size"], row["_support"]), xytext=(3, 3), textcoords="offset points", fontsize=6.5)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Supporting reads")
    style_axis(ax2, "both")

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.01, "Straglr is displayed as a separate tandem-repeat evidence layer; loci are not treated as Jasmine caller confirmations.", ha="center", fontsize=8.2)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
