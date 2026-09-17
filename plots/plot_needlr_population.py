#!/usr/bin/env python3
"""Plot needLR population-frequency evidence for one sample."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import add_panel_label, first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis

ANCESTRIES = ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]


def parse_args():
    p = argparse.ArgumentParser(description="Plot needLR control population frequencies.")
    p.add_argument("--input", required=True, help="needLR *_RESULTS.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="needLR population-frequency evidence")
    return p.parse_args()


def infer_frequency_column(df: pd.DataFrame, ancestry: str) -> str | None:
    candidates = [
        f"Pop_Freq_{ancestry}", f"Population_Freq_{ancestry}", f"Allele_Freq_{ancestry}",
        f"Control_Allele_Freq_{ancestry}", f"AF_{ancestry}", f"AF_1KGP_{ancestry}",
        f"1KGP_AF_{ancestry}", f"Freq_{ancestry}",
    ]
    return first_existing(df, candidates)


def infer_overall_af(df: pd.DataFrame) -> tuple[str, pd.Series]:
    col = first_existing(df, [
        "Allele_Freq_ALL", "Pop_Freq_ALL", "AlleleFreqAll", "AF_ALL", "AF",
        "MAX_AF", "AF_MAX", "SV_AF", "AF_1KGP", "1KGP_AF",
    ])
    if col is None:
        raise ValueError("Could not identify an overall needLR population-frequency column.")
    return col, numeric(df[col])


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("needLR table is empty.")

    af_col, af = infer_overall_af(df)
    work = df.copy()
    work["OVERALL_AF"] = af
    work["AF_CLASS"] = pd.cut(
        work["OVERALL_AF"],
        bins=[-np.inf, 0, 0.001, 0.01, np.inf],
        labels=["Not observed", "AF ≤ 0.001", "0.001 < AF ≤ 0.01", "AF > 0.01"],
        include_lowest=True,
    ).astype(object)
    work.loc[work["OVERALL_AF"].isna(), "AF_CLASS"] = "Not evaluable"

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    class_order = ["Not observed", "AF ≤ 0.001", "0.001 < AF ≤ 0.01", "AF > 0.01", "Not evaluable"]
    burden = work["AF_CLASS"].value_counts().reindex(class_order).fillna(0).rename_axis("AF_class").reset_index(name="count")
    burden["percent"] = burden["count"] / max(len(work), 1) * 100.0
    burden.to_csv(prefix.with_name(prefix.name + "_frequency_classes.tsv"), sep="\t", index=False)

    ancestry_rows = []
    ancestry_values: dict[str, pd.Series] = {}
    for anc in ANCESTRIES:
        col = infer_frequency_column(df, anc)
        if col is None:
            continue
        values = numeric(df[col])
        ancestry_values[anc] = values
        positive = values[values > 0]
        ancestry_rows.append({
            "ancestry": anc, "column": col,
            "n_evaluable": int(values.notna().sum()),
            "n_present": int((values > 0).sum()),
            "percent_present": float((values > 0).sum() / max(values.notna().sum(), 1) * 100.0),
            "mean_AF_when_present": float(positive.mean()) if not positive.empty else np.nan,
            "median_AF_when_present": float(positive.median()) if not positive.empty else np.nan,
        })
    ancestry_summary = pd.DataFrame(ancestry_rows)
    ancestry_summary.to_csv(prefix.with_name(prefix.name + "_ancestry_summary.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.5))
    ax1, ax2, ax3, ax4 = axes.flatten()
    colors = ["#BDBDBD", "#56B4E9", "#0072B2", "#D55E00", "#777777"]
    ax1.bar(burden["AF_class"], burden["percent"], color=colors)
    ax1.set_ylabel("Variants (%)")
    ax1.tick_params(axis="x", rotation=25)
    style_axis(ax1, "y")
    add_panel_label(ax1, "A")

    positive = work.loc[work["OVERALL_AF"] > 0, "OVERALL_AF"]
    if positive.empty:
        ax2.text(0.5, 0.5, "No variants with AF > 0", transform=ax2.transAxes, ha="center", va="center")
    else:
        ax2.hist(np.log10(positive), bins=35, color="#0072B2", alpha=0.9)
    ax2.set_xlabel(r"Control AF, $\log_{10}$(AF)")
    ax2.set_ylabel("Variants")
    style_axis(ax2, "y")
    add_panel_label(ax2, "B")

    if not ancestry_summary.empty:
        tmp = ancestry_summary[ancestry_summary["ancestry"].ne("ALL")]
        ax3.bar(tmp["ancestry"], tmp["percent_present"], color="#009E73")
        ax3.set_ylabel("SVs present in control population (%)")
        ax3.set_xlabel("needLR control ancestry")
    else:
        ax3.text(0.5, 0.5, "Ancestry-specific columns not detected", transform=ax3.transAxes, ha="center", va="center")
    style_axis(ax3, "y")
    add_panel_label(ax3, "C")

    labels = [a for a in ANCESTRIES if a in ancestry_values]
    nonempty = []
    for anc in labels:
        vals = ancestry_values[anc]
        vals = vals[vals > 0]
        if not vals.empty:
            nonempty.append((anc, np.log10(vals)))
    if nonempty:
        ax4.boxplot([v for _, v in nonempty], tick_labels=[l for l, _ in nonempty], showfliers=False)
        ax4.set_ylabel(r"AF among observed SVs, $\log_{10}$(AF)")
        ax4.set_xlabel("Control population")
    else:
        ax4.text(0.5, 0.5, "No positive ancestry-specific AF values", transform=ax4.transAxes, ha="center", va="center")
    style_axis(ax4, "y")
    add_panel_label(ax4, "D")

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.008, f"Overall frequency source: {af_col}. AF=0 means not observed in the needLR control dataset; it is not a pathogenicity label.", ha="center", fontsize=8.3)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
