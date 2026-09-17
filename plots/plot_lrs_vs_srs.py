#!/usr/bin/env python3
"""Plot LRS-vs-SRS SV concordance from a Truvari summary and optional match table.

This is a technology-concordance visualization, not an accuracy claim. When neither
technology is a truth set, precision/recall from Truvari are direction-dependent and
should be described as base/comparison recovery metrics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description="Plot Truvari-based LRS/SRS concordance.")
    p.add_argument("--summary", required=True, help="Truvari summary.json")
    p.add_argument("--match-tsv", default=None, help="Optional TSV of matched SVs with LRS/SRS size and breakpoint fields")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="Long-read versus short-read SV concordance")
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()
    with open(args.summary, encoding="utf-8") as fh:
        summary = json.load(fh)

    tp_base = int(summary.get("TP-base", summary.get("TP_base", 0)))
    tp_call = int(summary.get("TP-call", summary.get("TP_call", tp_base)))
    fp = int(summary.get("FP", 0))
    fn = int(summary.get("FN", 0))
    precision = float(summary.get("precision", np.nan))
    recall = float(summary.get("recall", np.nan))
    f1 = float(summary.get("f1", summary.get("f1-score", np.nan)))

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)

    match = None
    if args.match_tsv:
        match = read_tsv(args.match_tsv)

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.3))
    ax1, ax2, ax3, ax4 = axes.flatten()

    labels = ["Matched base", "Matched comparison", "Comparison-only", "Base-only"]
    values = [tp_base, tp_call, fp, fn]
    ax1.bar(labels, values, color=["#009E73", "#56B4E9", "#E69F00", "#D55E00"])
    ax1.set_ylabel("SV count")
    ax1.tick_params(axis="x", rotation=25)
    style_axis(ax1, "y")

    metrics = pd.Series({"Base recovery": recall, "Comparison precision": precision, "F1": f1}).dropna()
    if not metrics.empty:
        ax2.bar(metrics.index, metrics.values * 100.0, color=["#0072B2", "#E69F00", "#009E73"][: len(metrics)])
        ax2.set_ylim(0, 105)
        ax2.set_ylabel("Truvari metric (%)")
        ax2.tick_params(axis="x", rotation=20)
    else:
        ax2.text(0.5, 0.5, "No Truvari metrics found", transform=ax2.transAxes, ha="center", va="center")
    style_axis(ax2, "y")

    if match is not None and not match.empty:
        lrs_len = first_existing(match, ["LRS_SVLEN", "lrs_svlen", "base_svlen"])
        srs_len = first_existing(match, ["SRS_SVLEN", "srs_svlen", "comp_svlen"])
        bp_diff = first_existing(match, ["BREAKPOINT_DIFF", "breakpoint_diff", "start_diff"])
        type_col = first_existing(match, ["SVTYPE", "svtype", "LRS_SVTYPE"])

        if lrs_len and srs_len:
            x = numeric(match[lrs_len]).abs()
            y = numeric(match[srs_len]).abs()
            ok = x.gt(0) & y.gt(0)
            ax3.scatter(np.log10(x[ok]), np.log10(y[ok]), s=14, alpha=0.55, color="#0072B2")
            if ok.any():
                lo = min(np.log10(x[ok]).min(), np.log10(y[ok]).min())
                hi = max(np.log10(x[ok]).max(), np.log10(y[ok]).max())
                ax3.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0, color="#777777")
            ax3.set_xlabel(r"LRS SV size, $\log_{10}$(bp)")
            ax3.set_ylabel(r"SRS SV size, $\log_{10}$(bp)")
        else:
            ax3.text(0.5, 0.5, "Match TSV has no paired SV-length columns", transform=ax3.transAxes, ha="center", va="center")
        style_axis(ax3, "both")

        if bp_diff:
            vals = numeric(match[bp_diff]).abs().dropna()
            if not vals.empty:
                ax4.hist(np.log10(vals + 1), bins=35, color="#CC79A7")
                ax4.set_xlabel(r"Breakpoint difference, $\log_{10}$(bp + 1)")
                ax4.set_ylabel("Matched SV pairs")
        elif type_col:
            counts = match[type_col].value_counts()
            ax4.bar(counts.index, counts.values, color="#6E6E6E")
            ax4.set_xlabel("SV type")
            ax4.set_ylabel("Matched SV pairs")
        else:
            ax4.text(0.5, 0.5, "No breakpoint-difference/SVTYPE column", transform=ax4.transAxes, ha="center", va="center")
        style_axis(ax4, "y")
    else:
        for ax, msg in [(ax3, "Optional paired-match TSV not supplied"), (ax4, "Optional paired-match TSV not supplied")]:
            ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center")
            ax.axis("off")

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.008, "If neither LRS nor SRS is a validated truth set, interpret Truvari outputs as concordance/recovery rather than sensitivity, specificity or diagnostic accuracy.", ha="center", fontsize=8.1)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
