#!/usr/bin/env python3
"""Descriptive comparison of custom SV-gene and Genomiser gene rankings.

This is NOT a sensitivity, specificity, diagnostic accuracy or truth benchmark:
the candidate universes differ and only independently confirmed cases can
provide diagnostic validation.
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import save_figure, set_thesis_style, style_axis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--details", required=True, help="*.ranking_comparison.tsv")
    parser.add_argument("--summary", required=True, help="*.ranking_comparison.json")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()
    if args.top_n < 1:
        parser.error("--top-n must be positive")

    data = pd.read_csv(args.details, sep="\t", dtype={"GENE": str})
    metrics = json.loads(Path(args.summary).read_text(encoding="utf-8"))
    required = {"GENE", "CUSTOM_GENE_RANK_GLOBAL", "GENOMISER_GENE_RANK"}
    if not required.issubset(data.columns):
        raise ValueError(f"Missing rank fields: {sorted(required - set(data.columns))}")

    data["custom_rank"] = pd.to_numeric(data["CUSTOM_GENE_RANK_GLOBAL"], errors="coerce")
    data["genomiser_rank"] = pd.to_numeric(data["GENOMISER_GENE_RANK"], errors="coerce")
    both = data.dropna(subset=["custom_rank", "genomiser_rank"]).copy()
    chosen = both[
        (both["custom_rank"] <= args.top_n)
        | (both["genomiser_rank"] <= args.top_n)
    ].copy()

    ks = []
    agreements = []
    overlaps = []
    for key, entry in metrics.get("top_k", {}).items():
        if not str(key).isdigit():
            continue
        value = entry.get("jaccard")
        if value is None:
            continue
        ks.append(int(key))
        agreements.append(float(value))
        overlaps.append(int(entry["shared_n"]))
    order = np.argsort(ks)
    ks = [ks[i] for i in order]
    agreements = [agreements[i] for i in order]
    overlaps = [overlaps[i] for i in order]

    set_thesis_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    ax, bx = axes
    if ks:
        bars = ax.bar([str(v) for v in ks], agreements, color="#0072B2")
        ax.set_ylim(0, 1)
        ax.set_xlabel("Top-N genes")
        ax.set_ylabel("Jaccard overlap")
        for bar, overlap in zip(bars, overlaps):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.02,
                    f"{overlap} shared", ha="center", fontsize=9)
    else:
        ax.text(0.5, 0.5, "Top-N overlap unavailable", ha="center", transform=ax.transAxes)
    style_axis(ax, "y")
    ax.set_title("Candidate-set agreement")

    if not chosen.empty:
        panel = chosen.get("CUSTOM_PANEL_STATUS", pd.Series("", index=chosen.index))
        colors = np.where(panel.astype(str).str.upper().isin(["PANEL_GENE", "YES"]),
                          "#0072B2", "#009E73")
        bx.scatter(chosen["custom_rank"], chosen["genomiser_rank"],
                   c=colors, alpha=0.8, s=45)
        limit = max(chosen["custom_rank"].max(), chosen["genomiser_rank"].max(), args.top_n)
        bx.plot([1, limit], [1, limit], "--", color="gray", linewidth=1)
        bx.set_xscale("log")
        bx.set_yscale("log")
        bx.set_xlabel("Custom global SV-gene rank (log)")
        bx.set_ylabel("Genomiser gene rank (log)")
        bx.set_xlim(0.8, limit * 1.2)
        bx.set_ylim(0.8, limit * 1.2)
        for _, row in chosen.nsmallest(min(7, len(chosen)), "custom_rank").iterrows():
            bx.annotate(str(row["GENE"]), (row["custom_rank"], row["genomiser_rank"]),
                        xytext=(3, 3), textcoords="offset points", fontsize=8)
    else:
        bx.text(0.5, 0.5, "No shared genes within the displayed top ranks",
                ha="center", va="center", transform=bx.transAxes)
        bx.set_xlabel("Custom SV-gene rank")
        bx.set_ylabel("Genomiser rank")
    bx.set_title("Shared genes near either top-N")
    fig.suptitle("Custom prioritization versus Genomiser", fontsize=15, fontweight="bold")
    fig.text(
        0.5, 0.015,
        "Descriptive agreement, not diagnostic accuracy. Custom ranking covers SV-associated genes; "
        "Genomiser covers a broader genome-wide variant set.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.07, 1, 0.94])
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *paths, sep="\n")


if __name__ == "__main__":
    main()
