#!/usr/bin/env python3
"""Plot mechanism-aware large and breakpoint-defined SV-gene candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


DIRECT_BUCKETS = [
    "BREAKPOINT_GENE_CANDIDATE",
    "LARGE_CNV_GENE_CANDIDATE",
    "VERY_LARGE_CNV_GENE_CONTEXT",
]

CONTEXT_BUCKETS = [
    "LARGE_COMPLEX_INTERVAL_CONTEXT",
    "COMPLEX_INTERVAL_CONTEXT",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-per-bucket", type=int, default=6)
    return p.parse_args()


def size_label(row):
    span = pd.to_numeric(pd.Series([row.get("SV_EVENT_SPAN_BP")]), errors="coerce").iloc[0]
    svtype = str(row.get("SVTYPE", "."))
    if pd.isna(span):
        return svtype
    if span >= 1_000_000:
        return f"{svtype} {span / 1e6:.1f} Mb"
    if span >= 1_000:
        return f"{svtype} {span / 1e3:.1f} kb"
    return f"{svtype} {span:.0f} bp"


def compact_population(value):
    mapping = {
        "RARE_AF_LE_0.01": "rare",
        "COMMON_AF_GT_0.01": "common",
        "NOT_EVALUABLE_GE10MB": "AF n/e ≥10Mb",
        "NOT_EVALUABLE_BND": "AF n/e BND",
        "NO_MATCH": "AF no match",
        "UNKNOWN_OR_MISSING": "AF unknown",
    }
    return mapping.get(str(value), str(value))


def choose(df, buckets, n):
    parts = []
    for bucket in buckets:
        sub = df[df["EVENT_REVIEW_BUCKET"].eq(bucket)].copy()
        if sub.empty:
            continue
        sub["_score"] = numeric(sub["EVENT_GENE_RELEVANCE_SCORE"]).fillna(0)
        sub["_rank"] = numeric(sub["EVENT_RANK_WITHIN_BUCKET"]).fillna(np.inf)
        sub = sub.sort_values(["_rank", "_score"], ascending=[True, False]).head(n)
        parts.append(sub)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def plot_panel(ax, data, title):
    if data.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No candidates in this class", transform=ax.transAxes, ha="center")
        ax.set_title(title)
        return

    data = data.copy()
    data["_score"] = numeric(data["EVENT_GENE_RELEVANCE_SCORE"]).fillna(0)
    data = data.sort_values(
        ["EVENT_REVIEW_BUCKET", "_score"],
        ascending=[True, True],
    ).reset_index(drop=True)

    y = np.arange(len(data))
    bars = ax.barh(y, data["_score"])

    labels = []
    for _, row in data.iterrows():
        gene = str(row.get("GENES", "."))
        relation = str(row.get("SV_GENE_RELATIONSHIP", "."))
        labels.append(f"{gene} | {size_label(row)}\n{relation}")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.2)
    ax.set_xlabel("Gene relevance score (phenotype + curated disease evidence)")
    ax.set_title(title)
    style_axis(ax, "x")

    xmax = max(float(data["_score"].max()), 1.0)
    for bar, (_, row) in zip(bars, data.iterrows()):
        callers = str(row.get("CALLER_COUNT", "."))
        pop = compact_population(row.get("EVENT_POPULATION_TIER", "."))
        genes = str(row.get("SV_GENE_COUNT", "."))
        panel = "panel" if str(row.get("PANEL_STATUS", "")).upper() == "PANEL_GENE" else "non-panel"
        ax.text(
            bar.get_width() + 0.012 * xmax,
            bar.get_y() + bar.get_height() / 2,
            f"{callers} caller(s); {pop}; genes={genes}; {panel}",
            va="center",
            fontsize=7.8,
        )

    ax.set_xlim(0, xmax * 1.65)


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)

    required = {"EVENT_REVIEW_BUCKET", "EVENT_GENE_RELEVANCE_SCORE"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            "Input is not the mechanism-aware SV-gene event ranking; missing: "
            + ", ".join(sorted(missing))
        )

    direct = choose(df, DIRECT_BUCKETS, args.top_per_bucket)
    context = choose(df, CONTEXT_BUCKETS, args.top_per_bucket)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    selected = pd.concat(
        [
            direct.assign(PLOT_SECTION="gene-directed_or_dosage"),
            context.assign(PLOT_SECTION="interval_context_only"),
        ],
        ignore_index=True,
    )
    selected.to_csv(
        prefix.with_name(prefix.name + "_selected.tsv"),
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(17.0, max(8.0, 0.45 * max(len(direct), len(context), 10) + 3.2)),
        gridspec_kw={"width_ratios": [1.05, 1.0]},
    )

    plot_panel(
        axes[0],
        direct,
        "Large CNV and breakpoint-gene candidates",
    )
    plot_panel(
        axes[1],
        context,
        "Inversion/BND interval context only",
    )

    fig.suptitle(
        "Mechanism-aware large and breakpoint-defined SV review",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.012,
        "Left: events with dosage or breakpoint-level gene relevance. Right: genes lying inside inversion/BND intervals without breakpoint overlap; these are retained for context but are not treated as direct disruption. Size itself does not add pathogenicity points.",
        ha="center",
        fontsize=8.8,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.965])
    outputs = save_figure(fig, prefix)
    plt.close(fig)

    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
