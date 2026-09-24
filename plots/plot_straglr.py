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
    p.add_argument("--top-n", type=int, default=20)
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
        mask = g.ne("")
        work.loc[mask, "_label"] = work.loc[mask, locus].astype(str) + " | " + g.loc[mask]

    cn_max = work["_cn"].max(skipna=True)
    if pd.isna(cn_max) or cn_max <= 0:
        cn_component = pd.Series(0.0, index=work.index)
    else:
        cn_component = work["_cn"].fillna(0) / cn_max
    work["_priority"] = work["_support"] + cn_component
    work = work.sort_values(["_priority", "_cn", "_size"], ascending=False, na_position="last").head(args.top_n).iloc[::-1].copy()

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work.to_csv(prefix.with_name(prefix.name + "_top_loci.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14.0, max(7.0, 0.36 * len(work) + 2.2)), gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax1, ax2 = axes
    y = np.arange(len(work))

    ax1.barh(y, work["_support"], color="#0072B2")
    ax1.set_yticks(y)
    ax1.set_yticklabels(work["_label"], fontsize=9)
    ax1.set_xlabel("Supporting reads")
    ax1.set_ylabel("Repeat locus")
    style_axis(ax1, "x")

    use_cn = work["_cn"].notna().any()
    x = work["_cn"] if use_cn else work["_size"]
    xlabel = "Maximum copy number" if use_cn else "Maximum repeat size"
    if x.notna().any():
        ax2.scatter(x, work["_support"], s=45 + work["_coverage"].fillna(0).clip(lower=0) * 2, color="#E69F00", edgecolor="white", linewidth=0.5)
        label_rows = work.sort_values(["_support", "_cn"], ascending=False).head(8)
        for _, row in label_rows.iterrows():
            xpos = row["_cn"] if use_cn else row["_size"]
            if pd.notna(xpos):
                ax2.annotate(
                    str(row["_label"]).split(" | ")[0],
                    (xpos, row["_support"]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                )
    else:
        ax2.text(0.5, 0.5, "No numeric copy-number/size values detected", transform=ax2.transAxes, ha="center", va="center")
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Supporting reads")
    style_axis(ax2, "both")

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(0.5, 0.01, "Straglr is displayed as a separate tandem-repeat evidence layer; loci are not treated as Jasmine caller confirmations.", ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
