#!/usr/bin/env python3
"""Plot separate rankings for nuclear-encoded mitochondrial and mtDNA genes."""

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


def plot_group(df, encoding, title, out_prefix, top_n):
    group = df[
        df["ENCODING_GENOME"].fillna("").astype(str).str.upper().eq(encoding)
    ].copy()
    group["MITO_RANK_WITHIN_ENCODING"] = pd.to_numeric(
        group["MITO_RANK_WITHIN_ENCODING"], errors="coerce"
    )
    group = group.sort_values("MITO_RANK_WITHIN_ENCODING").head(top_n)

    fig_height = max(5.5, 0.38 * max(len(group), 1) + 2.6)
    fig, ax = plt.subplots(figsize=(12.5, fig_height))

    if group.empty:
        ax.axis("off")
        ax.text(
            0.5, 0.56,
            f"No {title.lower()} were linked to structural variants",
            transform=ax.transAxes,
            ha="center", va="center", fontsize=13, fontweight="bold",
        )
        ax.text(
            0.5, 0.42,
            "This is a negative result for the current SV table, not proof that the gene class is irrelevant.",
            transform=ax.transAxes,
            ha="center", va="center", fontsize=9.5,
        )
    else:
        view = group.iloc[::-1].reset_index(drop=True)
        y = np.arange(len(view))
        bars = ax.barh(y, view["UNIQUE_SVS"].astype(float))
        ax.set_yticks(y)
        ax.set_yticklabels(
            [
                f"{gene}  [rank {int(rank)}]"
                for gene, rank in zip(
                    view["GENE"], view["MITO_RANK_WITHIN_ENCODING"]
                )
            ],
            fontsize=9,
        )
        ax.set_xlabel("Unique structural variants linked to gene")
        ax.set_title(title)
        style_axis(ax, "x")

        xmax = max(float(view["UNIQUE_SVS"].astype(float).max()), 1.0)
        for bar, (_, row) in zip(bars, view.iterrows()):
            ax.text(
                bar.get_width() + 0.02 * xmax,
                bar.get_y() + bar.get_height() / 2,
                (
                    f"direct={row['DIRECT_SVS']}; proximal={row['PROXIMAL_SVS']}; "
                    f"context={row['CONTEXT_SVS']} | {row['MITO_PRIORITY_TIER']}"
                ),
                va="center",
                fontsize=8,
            )
        ax.set_xlim(0, xmax * 1.55)

    fig.subplots_adjust(left=0.22, right=0.98, top=0.91, bottom=0.12)
    outputs = save_figure(fig, out_prefix)
    plt.close(fig)
    return outputs


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)

    required = {
        "GENE", "ENCODING_GENOME", "MITO_RANK_WITHIN_ENCODING",
        "MITO_PRIORITY_TIER", "UNIQUE_SVS", "DIRECT_SVS",
        "PROXIMAL_SVS", "CONTEXT_SVS",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            "Mitochondrial ranking table is missing: " + ", ".join(missing)
        )

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    nuclear_prefix = prefix.with_name(prefix.name + "_nuclear")
    mtdna_prefix = prefix.with_name(prefix.name + "_mtdna")

    nuclear_outputs = plot_group(
        df,
        "NUCLEAR",
        "Ranked nuclear-encoded mitochondrial genes affected by SVs",
        nuclear_prefix,
        args.top_n,
    )
    mtdna_outputs = plot_group(
        df,
        "MTDNA",
        "Ranked mtDNA-encoded genes affected by SVs",
        mtdna_prefix,
        args.top_n,
    )

    print("[OK]", *nuclear_outputs, *mtdna_outputs, sep="\n")


if __name__ == "__main__":
    main()
