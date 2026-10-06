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
    "INVERSION_SPANNED_GENE_CONTEXT",
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
    final_rank_col = first_existing(
        df,
        [
            "FINAL_EVENT_RANK_GLOBAL",
            "EVENT_RANK_GLOBAL",
            "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS",
            "EVENT_RANK_WITHIN_PANEL_STATUS",
        ],
    )

    for bucket in buckets:
        sub = df[df["EVENT_REVIEW_BUCKET"].eq(bucket)].copy()
        if sub.empty:
            continue

        sub["_score"] = numeric(
            sub["EVENT_GENE_RELEVANCE_SCORE"]
        ).fillna(0)
        sub["_rank"] = (
            numeric(sub[final_rank_col]).fillna(np.inf)
            if final_rank_col
            else numeric(sub["EVENT_RANK_WITHIN_BUCKET"]).fillna(np.inf)
        )

        if "PANEL_STATUS" in sub.columns:
            panel_mask = (
                sub["PANEL_STATUS"]
                .fillna("")
                .astype(str)
                .eq("PANEL_GENE")
            )
        else:
            panel_mask = pd.Series(False, index=sub.index)

        panel = sub[panel_mask].sort_values(
            ["_rank", "_score"],
            ascending=[True, False],
        )

        nonpanel = sub[~panel_mask].sort_values(
            ["_rank", "_score"],
            ascending=[True, False],
        )

        selected = []
        if not panel.empty:
            selected.append(panel.head(max(1, n // 2)))
        if not nonpanel.empty:
            selected.append(
                nonpanel.head(max(1, n - sum(len(x) for x in selected)))
            )

        chosen = (
            pd.concat(selected)
            if selected
            else sub.head(n)
        )

        if len(chosen) < n:
            remainder = sub.loc[
                ~sub.index.isin(chosen.index)
            ].sort_values(
                ["_rank", "_score"],
                ascending=[True, False],
            )
            chosen = pd.concat(
                [chosen, remainder.head(n - len(chosen))]
            )

        parts.append(chosen.head(n))

    return (
        pd.concat(parts, ignore_index=True)
        if parts
        else pd.DataFrame()
    )


def compact_gene(value, max_chars=18):
    text = str(value)
    return text if len(text) <= max_chars else text[: max_chars - 3] + "..."


def compact_relationship(value):
    mapping = {
        "WHOLE_GENE_DOSAGE_CONTEXT": "whole-gene dosage",
        "PARTIAL_GENE_OVERLAP": "partial overlap",
        "BREAKPOINT_WITHIN_TRANSCRIPT": "breakpoint in gene",
        "BREAKPOINT_PROXIMAL_TO_GENE": "breakpoint near gene",
        "INSERTION_WITHIN_TRANSCRIPT": "insertion in gene",
        "INSERTION_PROXIMAL_TO_GENE": "insertion near gene",
        "INVERSION_SPANS_INTACT_GENE": "INV spans intact gene",
        "INTERVAL_CONTEXT_ONLY": "interval context",
        "GENE_PROXIMAL_INTERVAL": "proximal interval",
    }
    text = str(value)
    return mapping.get(text, text.replace("_", " ").lower()[:30])


def plot_panel(ax, data, title):
    if data.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No candidates in this class", transform=ax.transAxes, ha="center")
        ax.set_title(title)
        return

    data = data.copy()
    data["_score"] = numeric(data["EVENT_GENE_RELEVANCE_SCORE"]).fillna(0)
    final_rank_col = first_existing(
        data,
        [
            "FINAL_EVENT_RANK_GLOBAL",
            "EVENT_RANK_GLOBAL",
            "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS",
            "EVENT_RANK_WITHIN_PANEL_STATUS",
        ],
    )
    data["_final_rank"] = (
        numeric(data[final_rank_col]).fillna(np.inf)
        if final_rank_col
        else np.inf
    )
    data = data.sort_values(
        ["EVENT_REVIEW_BUCKET", "_final_rank", "_score"],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    y = np.arange(len(data))
    bars = ax.barh(y, data["_score"])

    labels = []
    for _, row in data.iterrows():
        gene = compact_gene(row.get("GENES", "."))
        relation = compact_relationship(row.get("SV_GENE_RELATIONSHIP", "."))
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
        final_rank = (
            int(row["_final_rank"])
            if np.isfinite(row["_final_rank"])
            else "."
        )
        ax.text(
            bar.get_width() + 0.012 * xmax,
            bar.get_y() + bar.get_height() / 2,
            f"r{final_rank}; {callers} caller(s); {pop}; n={genes}; {panel}",
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
        "Inversion-spanned and rearrangement context",
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
        "Left: events with dosage or direct/proximal breakpoint-level gene relevance. Right: genes fully spanned by inversions or otherwise lying in rearranged intervals without a direct transcript breakpoint. These may warrant regulatory, position-effect or 3D-genome review but are not automatically considered non-functional. Size itself does not add pathogenicity points.",
        ha="center",
        fontsize=8.8,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.965])
    outputs = save_figure(fig, prefix)
    plt.close(fig)

    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
