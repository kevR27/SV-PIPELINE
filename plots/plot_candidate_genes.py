#!/usr/bin/env python3
"""Plot panel and non-panel gene rankings as separate thesis figures.

The input may be the early genome-wide gene ranking or the final patient-aware
ranking. Final patient-specific fields are preferred when present.

Panel membership is a display/grouping variable only. It does not contribute
ranking points.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import (
    first_existing,
    numeric,
    read_tsv,
    save_figure,
    set_thesis_style,
    style_axis,
)


TIER_ORDER = {
    "HIGH": 3,
    "MODERATE": 2,
    "SUPPORTING": 1,
    "LIMITED": 0,
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument(
        "--title",
        default="SV-associated gene prioritization",
    )
    return p.parse_args()


def panel_status(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.upper().str.strip()
    return np.where(
        text.isin(["PANEL_GENE", "YES", "TRUE", "1"]),
        "PANEL_GENE",
        "NONPANEL_GENE",
    )


def prepare(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    gene_col = first_existing(df, ["GENE", "gene", "Gene", "SYMBOL"])
    if gene_col is None:
        raise ValueError("Could not infer gene column.")

    panel_col = first_existing(df, ["PANEL_STATUS", "panel_gene"])
    if panel_col is None:
        raise ValueError("Candidate table requires PANEL_STATUS/panel_gene.")

    tier_col = first_existing(
        df,
        ["FINAL_GENE_RELEVANCE_TIER", "GENE_RELEVANCE_TIER"],
    )
    score_col = first_existing(
        df,
        [
            "FINAL_GENE_RELEVANCE_DISPLAY_SCORE",
            "GENE_RELEVANCE_DISPLAY_SCORE",
            "integrated_discovery_score",
            "GENE_RELEVANCE_SCORE",
            "phenotype_score",
        ],
    )
    rank_col = first_existing(
        df,
        [
            "FINAL_GENE_RANK_WITHIN_PANEL_STATUS",
            "GENE_RANK_WITHIN_PANEL_STATUS",
        ],
    )
    inheritance_col = first_existing(
        df,
        ["GENE_INHERITANCE_CLASS", "BEST_EVENT_INHERITANCE_MECHANISM_CLASS"],
    )
    mechanism_col = first_existing(
        df,
        [
            "BEST_EVENT_INHERITANCE_MECHANISM_CLASS",
            "INHERITANCE_MECHANISM_CLASS",
        ],
    )
    sv_count_col = first_existing(df, ["SV_COUNT", "SV_count", "master_SV_count"])

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    work["_panel_status"] = panel_status(work[panel_col])
    work["_tier"] = (
        work[tier_col].fillna("LIMITED").astype(str).str.upper()
        if tier_col
        else "LIMITED"
    )
    work["_tier_rank"] = (
        pd.Series(work["_tier"], index=work.index)
        .map(TIER_ORDER)
        .fillna(0)
    )
    work["_score"] = (
        numeric(work[score_col]).fillna(0)
        if score_col
        else pd.Series(0.0, index=work.index)
    )
    work["_sv_count"] = (
        numeric(work[sv_count_col]).fillna(0)
        if sv_count_col
        else pd.Series(0.0, index=work.index)
    )
    work["_rank"] = (
        numeric(work[rank_col])
        if rank_col
        else pd.Series(np.nan, index=work.index)
    )
    work["_inheritance"] = (
        work[inheritance_col].fillna("UNKNOWN").astype(str)
        if inheritance_col
        else "UNKNOWN"
    )
    work["_mechanism"] = (
        work[mechanism_col].fillna(".").astype(str)
        if mechanism_col
        else "."
    )

    work = work[
        ~work["_gene"].isin(["", ".", "NA", "N/A", "nan", "None"])
    ].copy()

    # Prefer the explicit within-panel rank. Older outputs fall back to the
    # same tier/score ordering used by the current gene ranker.
    work["_rank_missing"] = work["_rank"].isna().astype(int)
    work = work.sort_values(
        [
            "_panel_status",
            "_rank_missing",
            "_rank",
            "_tier_rank",
            "_score",
            "_gene",
        ],
        ascending=[True, True, True, False, False, True],
        na_position="last",
    )

    if rank_col is None:
        work["_rank"] = (
            work.groupby("_panel_status").cumcount() + 1
        )

    return work, score_col or "ranking score"


def plot_group(
    work: pd.DataFrame,
    group: str,
    label: str,
    score_label: str,
    prefix: Path,
    top_n: int,
):
    view = (
        work[work["_panel_status"].eq(group)]
        .head(top_n)
        .copy()
    )
    view = view.iloc[::-1].reset_index(drop=True)

    height = max(5.6, 0.42 * max(len(view), 1) + 2.4)
    fig, ax = plt.subplots(figsize=(13.5, height))

    if view.empty:
        ax.axis("off")
        ax.text(
            0.5,
            0.56,
            f"No {label.lower()} were available",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )
    else:
        y = np.arange(len(view))
        ax.hlines(
            y,
            0,
            view["_score"],
            linewidth=1.0,
            alpha=0.45,
        )
        sizes = 55 + np.sqrt(view["_sv_count"].clip(lower=0)) * 35
        ax.scatter(
            view["_score"],
            y,
            s=sizes,
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(
            [
                f"{gene}  [#{int(rank)} | {tier}]"
                for gene, rank, tier in zip(
                    view["_gene"],
                    view["_rank"],
                    view["_tier"],
                )
            ],
            fontsize=9,
        )
        ax.set_xlabel(score_label.replace("_", " "))
        ax.set_ylabel("Gene")
        style_axis(ax, "x")

        xmax = max(float(view["_score"].max()), 1.0)
        for i, (_, row) in enumerate(view.iterrows()):
            inheritance = str(row["_inheritance"]).replace("_", " ")
            mechanism = str(row["_mechanism"]).replace("_", " ")
            note = inheritance
            if mechanism not in {"", ".", "nan", "None"}:
                note += f" | {mechanism}"
            ax.text(
                row["_score"] + 0.02 * xmax,
                i,
                note,
                va="center",
                fontsize=7.6,
            )
        ax.set_xlim(0, xmax * 1.75)

    ax.set_title(label)
    fig.text(
        0.5,
        0.012,
        (
            "Rank is research prioritization, not pathogenicity. Gene relevance "
            "tier is primary; inheritance/mechanism is shown explicitly. "
            "Panel membership only determines which figure/rank group is used."
        ),
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    return outputs, view


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Ranked candidate table is empty.")

    work, score_label = prepare(df)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    all_outputs = []
    exported = []

    for group, suffix, label in [
        (
            "PANEL_GENE",
            "panel",
            "Optic-neuropathy panel genes affected by SVs",
        ),
        (
            "NONPANEL_GENE",
            "nonpanel",
            "Non-panel genes affected by SVs",
        ),
    ]:
        out_prefix = prefix.with_name(f"{prefix.name}_{suffix}")
        outputs, view = plot_group(
            work,
            group,
            label,
            score_label,
            out_prefix,
            args.top_n,
        )
        all_outputs.extend(outputs)
        if not view.empty:
            exported.append(view.assign(RANK_GROUP=group))

    if exported:
        pd.concat(exported, ignore_index=True).to_csv(
            prefix.with_name(prefix.name + "_panel_nonpanel_rankings.tsv"),
            sep="\t",
            index=False,
        )

    print("[OK]", *all_outputs, sep="\n")


if __name__ == "__main__":
    main()
