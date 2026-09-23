#!/usr/bin/env python3
"""Plot LRS coverage QC from mosdepth outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from plot_utils import save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description="Plot mosdepth coverage QC.")
    p.add_argument("--summary", required=True, help="*.mosdepth.summary.txt")
    p.add_argument("--global-dist", default=None, help="Optional *.mosdepth.global.dist.txt")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="Long-read sequencing coverage QC")
    return p.parse_args()


def load_summary(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    lower = {c.lower(): c for c in df.columns}
    required = {"chrom", "length", "bases", "mean", "min", "max"}
    missing = required - set(lower)
    if missing:
        raise ValueError(f"Unexpected mosdepth summary columns; missing {sorted(missing)} from {list(df.columns)}")
    return df


def load_global(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "depth", "fraction"], comment="#")
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    df["fraction"] = pd.to_numeric(df["fraction"], errors="coerce")
    return df.dropna(subset=["depth", "fraction"])


def main():
    args = parse_args()
    set_thesis_style()
    summary = load_summary(args.summary)

    lower = {c.lower(): c for c in summary.columns}
    chrom_col = lower["chrom"]
    mean_col = lower["mean"]
    min_col = lower["min"]
    max_col = lower["max"]

    for col in [mean_col, min_col, max_col]:
        summary[col] = pd.to_numeric(summary[col], errors="coerce")

    genome_row = summary[summary[chrom_col].astype(str).str.lower().isin(["total", "genome", "all"])]
    autosomal = summary[summary[chrom_col].astype(str).str.match(r"^(chr)?([1-9]|1[0-9]|2[0-2])$")].copy()

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    autosomal[[chrom_col, mean_col, min_col, max_col]].to_csv(
        prefix.with_name(prefix.name + "_chromosome_depth.tsv"), sep="\t", index=False
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.9))
    ax1, ax2 = axes

    if args.global_dist and Path(args.global_dist).exists():
        dist = load_global(args.global_dist)
        total = dist[dist["chrom"].astype(str).str.lower().isin(["total", "genome", "all"])]
        plot_dist = total if not total.empty else dist.groupby("depth", as_index=False)["fraction"].mean()
        ax1.plot(plot_dist["depth"], plot_dist["fraction"], color="#0072B2", linewidth=1.8)
        ax1.set_xlabel("Depth")
        ax1.set_ylabel("Fraction of genome at ≥ depth")
        ax1.set_xlim(left=0)
        ax1.set_ylim(0, 1.02)
    elif not autosomal.empty:
        autosomal["chrom_sort"] = (
            autosomal[chrom_col].astype(str).str.replace("chr", "", regex=False).astype(int)
        )
        autosomal = autosomal.sort_values("chrom_sort")
        x = range(len(autosomal))
        ax1.errorbar(
            x,
            autosomal[mean_col],
            yerr=[
                (autosomal[mean_col] - autosomal[min_col]).clip(lower=0),
                (autosomal[max_col] - autosomal[mean_col]).clip(lower=0),
            ],
            fmt="o",
            markersize=4,
            capsize=2,
            linewidth=0.8,
        )
        ax1.set_xticks(list(x))
        ax1.set_xticklabels(autosomal[chrom_col], rotation=55)
        ax1.set_xlabel("Chromosome")
        ax1.set_ylabel("Depth (mean with min/max range)")
    else:
        ax1.text(0.5, 0.5, "Coverage distribution unavailable", transform=ax1.transAxes, ha="center", va="center")
    style_axis(ax1, "both")

    if not autosomal.empty:
        if "chrom_sort" not in autosomal:
            autosomal["chrom_sort"] = (
                autosomal[chrom_col].astype(str).str.replace("chr", "", regex=False).astype(int)
            )
            autosomal = autosomal.sort_values("chrom_sort")
        ax2.bar(autosomal[chrom_col], autosomal[mean_col], color="#6E6E6E")
        ax2.tick_params(axis="x", rotation=55)
    else:
        ax2.text(0.5, 0.5, "Autosomal rows not detected", transform=ax2.transAxes, ha="center", va="center")
    ax2.set_xlabel("Chromosome")
    ax2.set_ylabel("Mean depth")
    style_axis(ax2, "y")

    if not genome_row.empty:
        genome_mean = pd.to_numeric(genome_row.iloc[0][mean_col], errors="coerce")
        if pd.notna(genome_mean):
            ax2.axhline(
                genome_mean,
                linestyle="--",
                linewidth=1.0,
                color="#D55E00",
                label=f"Genome mean = {genome_mean:.1f}×",
            )
            ax2.legend(frameon=False, fontsize=8)

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
