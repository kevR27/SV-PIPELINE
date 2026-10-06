#!/usr/bin/env python3
"""Plot mitochondrial gene rankings separated by encoding genome and panel status."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25)
    return p.parse_args()


def rank_column(df):
    if "MITO_RANK_WITHIN_ENCODING_PANEL_STATUS" in df.columns:
        return "MITO_RANK_WITHIN_ENCODING_PANEL_STATUS"
    return "MITO_RANK_WITHIN_ENCODING"


def plot_group(
    df,
    encoding,
    panel_status,
    title,
    out_prefix,
    top_n,
):
    rank_col = rank_column(df)
    group = df[
        df["ENCODING_GENOME"]
        .fillna("")
        .astype(str)
        .str.upper()
        .eq(encoding)
        & df["PANEL_STATUS"]
        .fillna("NONPANEL_GENE")
        .astype(str)
        .str.upper()
        .replace({"NON_PANEL": "NONPANEL_GENE"})
        .eq(panel_status)
    ].copy()

    group[rank_col] = pd.to_numeric(group[rank_col], errors="coerce")
    group = group.sort_values(rank_col).head(top_n)

    fig_height = max(5.5, 0.42 * max(len(group), 1) + 2.6)
    fig, ax = plt.subplots(figsize=(14.0, fig_height))

    if group.empty:
        ax.axis("off")
        ax.text(
            0.5,
            0.56,
            f"No genes in this ranking group",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )
    else:
        view = group.iloc[::-1].reset_index(drop=True)
        y = np.arange(len(view))
        bars = ax.barh(y, view["UNIQUE_SVS"].astype(float))
        ax.set_yticks(y)

        if encoding == "MTDNA" and "MTDNA_FUNCTION" in view.columns:
            labels = [
                (
                    f"{gene}  [{function.replace('_', ' ').lower()}] "
                    f"[rank {int(rank)}]"
                )
                for gene, function, rank in zip(
                    view["GENE"],
                    view["MTDNA_FUNCTION"].fillna(".").astype(str),
                    view[rank_col],
                )
            ]
        else:
            labels = [
                f"{gene}  [rank {int(rank)}]"
                for gene, rank in zip(view["GENE"], view[rank_col])
            ]

        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Unique structural variants linked to gene")
        style_axis(ax, "x")

        xmax = max(float(view["UNIQUE_SVS"].astype(float).max()), 1.0)
        for bar, (_, row) in zip(bars, view.iterrows()):
            mechanism = str(
                row.get(
                    "BEST_INHERITANCE_MECHANISM_CLASS",
                    row.get("MITO_PRIORITY_TIER", "."),
                )
            ).replace("_", " ")
            ax.text(
                bar.get_width() + 0.02 * xmax,
                bar.get_y() + bar.get_height() / 2,
                (
                    f"direct={row['DIRECT_SVS']}; "
                    f"proximal={row['PROXIMAL_SVS']}; "
                    f"context={row['CONTEXT_SVS']} | {mechanism}"
                ),
                va="center",
                fontsize=7.7,
            )
        ax.set_xlim(0, xmax * 1.8)

    ax.set_title(title)
    fig.text(
        0.5,
        0.01,
        (
            "Nuclear and mtDNA genes are ranked separately. Panel and non-panel "
            "genes also receive separate within-group ranks. Inheritance and "
            "SV mechanism are prioritization context, not pathogenicity calls."
        ),
        ha="center",
        fontsize=8.5,
    )
    fig.subplots_adjust(left=0.24, right=0.98, top=0.91, bottom=0.10)
    outputs = save_figure(fig, out_prefix)
    plt.close(fig)
    return outputs


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)

    required = {
        "GENE",
        "ENCODING_GENOME",
        "PANEL_STATUS",
        "MITO_RANK_WITHIN_ENCODING",
        "MITO_PRIORITY_TIER",
        "UNIQUE_SVS",
        "DIRECT_SVS",
        "PROXIMAL_SVS",
        "CONTEXT_SVS",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            "Mitochondrial ranking table is missing: " + ", ".join(missing)
        )

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    groups = [
        (
            "NUCLEAR",
            "PANEL_GENE",
            "nuclear_panel",
            "Panel: nuclear-encoded mitochondrial genes affected by SVs",
        ),
        (
            "NUCLEAR",
            "NONPANEL_GENE",
            "nuclear_nonpanel",
            "Non-panel: nuclear-encoded mitochondrial genes affected by SVs",
        ),
        (
            "MTDNA",
            "PANEL_GENE",
            "mtdna_panel",
            "Panel: mtDNA-encoded genes affected by SVs",
        ),
        (
            "MTDNA",
            "NONPANEL_GENE",
            "mtdna_nonpanel",
            "Non-panel: mtDNA-encoded genes affected by SVs",
        ),
    ]

    outputs = []
    for encoding, panel, suffix, title in groups:
        out_prefix = prefix.with_name(prefix.name + "_" + suffix)
        outputs.extend(
            plot_group(
                df,
                encoding,
                panel,
                title,
                out_prefix,
                args.top_n,
            )
        )

    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
