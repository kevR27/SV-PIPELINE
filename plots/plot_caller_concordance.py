#!/usr/bin/env python3
"""Plot caller concordance for a merged LRS or SRS sample."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import (
    CALLER_COLORS,
    add_panel_label,
    first_existing,
    read_tsv,
    save_figure,
    set_thesis_style,
    style_axis,
)

DEFAULT_CALLERS = ["Sniffles2", "cuteSV", "Delly"]


def parse_args():
    p = argparse.ArgumentParser(description="Plot Jasmine caller concordance.")
    p.add_argument("--input", required=True, help="*_caller_support_summary.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="Structural-variant caller concordance")
    p.add_argument("--caller-order", default=",".join(DEFAULT_CALLERS))
    return p.parse_args()


def infer_vector(row: pd.Series, callers: list[str]) -> str | None:
    if "SUPP_VEC" in row.index:
        vec = str(row.get("SUPP_VEC", "")).strip()
        if len(vec) == len(callers) and set(vec) <= {"0", "1"}:
            return vec

    caller_col = None
    for candidate in ["CALLERS", "CALLER"]:
        if candidate in row.index:
            caller_col = candidate
            break
    if caller_col is None:
        return None

    text = str(row.get(caller_col, ""))
    present = set()
    for token in text.replace("|", ";").replace(",", ";").split(";"):
        token = token.strip().lower()
        if not token:
            continue
        for caller in callers:
            if token == caller.lower():
                present.add(caller)
    return "".join("1" if c in present else "0" for c in callers)


def build_summary(df: pd.DataFrame, callers: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    vectors = []
    for _, row in df.iterrows():
        vec = infer_vector(row, callers)
        if vec and "1" in vec:
            vectors.append(vec)

    if not vectors:
        raise ValueError("No usable SUPP_VEC/CALLERS information found in the input table.")

    counts = Counter(vectors)
    rows = []
    for vec, n in counts.items():
        rows.append(
            {
                "SUPP_VEC": vec,
                "intersection_size": n,
                "caller_count": vec.count("1"),
                "callers": ";".join(c for c, bit in zip(callers, vec) if bit == "1"),
            }
        )
    intersections = pd.DataFrame(rows).sort_values(
        ["intersection_size", "caller_count", "SUPP_VEC"], ascending=[False, False, True]
    ).reset_index(drop=True)

    set_sizes = pd.DataFrame(
        {
            "caller": callers,
            "size": [sum(n for v, n in counts.items() if v[i] == "1") for i in range(len(callers))],
        }
    )
    return intersections, set_sizes


def plot_upset(intersections: pd.DataFrame, set_sizes: pd.DataFrame, callers: list[str], title: str):
    set_thesis_style()
    n = len(intersections)
    x = np.arange(n)

    fig = plt.figure(figsize=(13.5, 8.5))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.65, 6.8],
        height_ratios=[3.3, 2.25],
        hspace=0.06,
        wspace=0.10,
    )
    ax_empty = fig.add_subplot(gs[0, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax_sets = fig.add_subplot(gs[1, 0])
    ax_matrix = fig.add_subplot(gs[1, 1])
    ax_empty.axis("off")

    bar_colors = []
    for count in intersections["caller_count"]:
        if count == 3:
            bar_colors.append("#0072B2")
        elif count == 2:
            bar_colors.append("#56B4E9")
        else:
            bar_colors.append("#BDBDBD")

    bars = ax_top.bar(x, intersections["intersection_size"], color=bar_colors, width=0.72)
    ymax = max(float(intersections["intersection_size"].max()) * 1.17, 1.0)
    ax_top.set_ylim(0, ymax)
    ax_top.set_ylabel("Number of master SVs")
    ax_top.set_xticks([])
    style_axis(ax_top, "y")
    add_panel_label(ax_top, "A")

    for bar, value in zip(bars, intersections["intersection_size"]):
        ax_top.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ymax * 0.018,
            f"{int(value):,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    y = np.arange(len(callers))
    sizes = set_sizes.set_index("caller").loc[callers, "size"].values
    colors = [CALLER_COLORS.get(c, "#777777") for c in callers]
    ax_sets.barh(y, sizes, color=colors, height=0.58)
    ax_sets.invert_xaxis()
    ax_sets.invert_yaxis()
    ax_sets.set_yticks(y)
    ax_sets.set_yticklabels(callers)
    ax_sets.set_xlabel("Set size")
    style_axis(ax_sets, "x")
    add_panel_label(ax_sets, "B")

    xmax = max(float(max(sizes)), 1.0)
    ax_sets.set_xlim(xmax * 1.2, 0)
    for yi, value in zip(y, sizes):
        ax_sets.text(value * 0.98, yi, f"{int(value):,}", ha="left", va="center", fontsize=9, color="white", fontweight="bold")

    ax_matrix.set_xlim(-0.5, n - 0.5)
    ax_matrix.set_ylim(-0.5, len(callers) - 0.5)
    ax_matrix.invert_yaxis()
    ax_matrix.set_xticks(x)
    combo_labels = [
        " + ".join(c for c, bit in zip(callers, vec) if bit == "1")
        for vec in intersections["SUPP_VEC"]
    ]
    ax_matrix.set_xticklabels(combo_labels, rotation=35, ha="right", fontsize=8.5)
    ax_matrix.set_yticks(y)
    ax_matrix.set_yticklabels(callers)
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)

    for yi in y:
        if yi % 2:
            ax_matrix.axhspan(yi - 0.5, yi + 0.5, color="#F5F5F5", zorder=0)

    combination_labels = []
    for i, vec in enumerate(intersections["SUPP_VEC"]):
        active = [j for j, bit in enumerate(vec) if bit == "1"]
        ax_matrix.scatter([i] * len(callers), y, s=62, color="#D6D6D6", zorder=1)
        if len(active) > 1:
            ax_matrix.plot([i, i], [min(active), max(active)], color="#222222", lw=1.5, zorder=2)
        if active:
            ax_matrix.scatter([i] * len(active), active, s=78, color="#222222", zorder=3)
        names = [callers[j] for j in active]
        combination_labels.append(" + ".join(names))

    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels(combination_labels, rotation=38, ha="right", fontsize=8.5)
    ax_matrix.tick_params(axis="x", pad=5)
    add_panel_label(ax_matrix, "C")

    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.99)
    fig.text(
        0.5,
        0.01,
        "Bars show exact caller intersections; the complete merged callset remains the master SV universe.",
        ha="center",
        fontsize=9,
    )
    return fig


def main():
    args = parse_args()
    callers = [x.strip() for x in args.caller_order.split(",") if x.strip()]
    if len(callers) < 2:
        raise ValueError("--caller-order must contain at least two callers.")

    df = read_tsv(args.input)
    intersections, set_sizes = build_summary(df, callers)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    intersections.to_csv(prefix.with_name(prefix.name + "_intersections.tsv"), sep="\t", index=False)
    set_sizes.to_csv(prefix.with_name(prefix.name + "_set_sizes.tsv"), sep="\t", index=False)

    fig = plot_upset(intersections, set_sizes, callers, args.title)
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
