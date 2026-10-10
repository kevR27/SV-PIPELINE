#!/usr/bin/env python3
"""Plot the audited downstream filters without confusing genes with SVs."""

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from plot_utils import read_tsv, save_figure, set_thesis_style


CONTEXT_COLORS = {"EXONIC_OR_SPLICE": "#0072B2", "INTRONIC": "#E69F00", "EXONIC_AND_INTRONIC": "#009E73"}
CONTEXT_LABELS = {"EXONIC_OR_SPLICE": "Exonic / splice", "INTRONIC": "Intronic", "EXONIC_AND_INTRONIC": "Exonic and intronic"}


def export(fig, folder, name):
    folder.mkdir(parents=True, exist_ok=True)
    save_figure(fig, folder / name, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--filter-dir", required=True, help="Output of filter_sv_candidates_for_review.py")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    source, out = Path(args.filter_dir), Path(args.out_dir)
    settings = json.loads((source / "filter_run.json").read_text())
    data = read_tsv(source / "all_records.annotated.tsv", required=["DOWNSTREAM_SV_ID", "DOWNSTREAM_SIZE_GROUP", "DOWNSTREAM_DECISION", "DOWNSTREAM_CONTEXT", "DOWNSTREAM_CALLER_COUNT", "DOWNSTREAM_AF_STATUS"])
    set_thesis_style()
    sample = str(settings["sample"])
    safe_sample = re.sub(r"[^A-Za-z0-9_-]", "_", sample)
    cutoff = settings["size_cutoff_bp"]
    size_groups = ["GT_CUTOFF", "LT_CUTOFF", "AT_CUTOFF", "LENGTH_UNRESOLVED"]
    size_labels = [f">{cutoff:,} bp", f"<{cutoff:,} bp", f"Exactly {cutoff:,} bp", "Length unresolved / BND"]
    decisions = ["RETAINED", "REVIEW_REQUIRED", "EXCLUDED"]
    colors = ["#009E73", "#E69F00", "#999999"]

    # Every association belongs to one outcome and one size group.
    folder = out / "filter_outcomes"
    folder.mkdir(parents=True, exist_ok=True)
    counts = data.groupby(["DOWNSTREAM_SIZE_GROUP", "DOWNSTREAM_DECISION"], as_index=False).agg(SV_GENE_ASSOCIATIONS=("DOWNSTREAM_SV_ID", "size"), UNIQUE_SVS=("DOWNSTREAM_SV_ID", "nunique"))
    counts.to_csv(folder / f"{safe_sample}_filter_outcomes.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(11.5, 6))
    left = np.zeros(len(size_groups))
    for decision, color in zip(decisions, colors):
        sub = counts[counts["DOWNSTREAM_DECISION"].eq(decision)].set_index("DOWNSTREAM_SIZE_GROUP")
        values = sub["SV_GENE_ASSOCIATIONS"].reindex(size_groups, fill_value=0).to_numpy()
        ax.barh(size_labels, values, left=left, color=color, label=decision.replace("_", " ").title())
        left += values
    ax.invert_yaxis()
    ax.set_xlabel("SV–gene associations in the input table")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(f"{sample}: downstream filter outcomes", loc="left", fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.15), ncol=3, frameon=False)
    fig.text(0.5, 0.01, "Retained means passing non-size filters within that size group. Only the >cutoff group is the main large-SV view.\nA single SV can have several gene associations; unique-SV counts per outcome are in the TSV and may overlap.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.24, right=0.97, top=0.89, bottom=0.26)
    export(fig, folder, f"{safe_sample}_filter_outcomes")

    # Show small genic events alongside the main view, in separate figures.
    retained = data[data["DOWNSTREAM_DECISION"].eq("RETAINED")]
    folder = out / "retained_types_and_context"
    folder.mkdir(parents=True, exist_ok=True)
    for size, suffix, title in [("GT_CUTOFF", "large", f">{cutoff:,} bp"), ("LT_CUTOFF", "small", f"<{cutoff:,} bp")]:
        rows = retained[retained["DOWNSTREAM_SIZE_GROUP"].eq(size)]
        counts = rows.groupby(["DOWNSTREAM_SVTYPE", "DOWNSTREAM_CONTEXT"], as_index=False).agg(SV_GENE_ASSOCIATIONS=("DOWNSTREAM_SV_ID", "size"), UNIQUE_SVS=("DOWNSTREAM_SV_ID", "nunique"))
        counts.to_csv(folder / f"{safe_sample}_{suffix}_types.tsv", sep="\t", index=False)
        kinds = sorted(rows["DOWNSTREAM_SVTYPE"].unique())
        fig, ax = plt.subplots(figsize=(10.5, 5.5))
        bottom = np.zeros(len(kinds))
        for context, color in CONTEXT_COLORS.items():
            sub = counts[counts["DOWNSTREAM_CONTEXT"].eq(context)].set_index("DOWNSTREAM_SVTYPE")
            values = sub["SV_GENE_ASSOCIATIONS"].reindex(kinds, fill_value=0).to_numpy()
            ax.bar(kinds, values, bottom=bottom, color=color, label=CONTEXT_LABELS[context])
            bottom += values
        if not kinds:
            ax.text(0.5, 0.5, "No associations pass these filters", transform=ax.transAxes, ha="center")
        ax.set_ylabel("SV–gene associations")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(f"{sample}: retained {title} annotation context", loc="left", fontweight="bold")
        ax.legend(loc="upper left", bbox_to_anchor=(0, -0.12), ncol=3, frameon=False)
        fig.text(0.5, 0.01, "Annotated exon/splice/intron context is not proof of disruption or exact breakpoint location.\nMixed means both labels occur across the available transcript annotations. Small events can still matter biologically.", ha="center", fontsize=9)
        fig.subplots_adjust(left=0.12, right=0.97, top=0.89, bottom=0.27)
        export(fig, folder, f"{safe_sample}_{suffix}_types")

    # Event-wide evidence: count each retained large SV only once.
    events = retained[retained["DOWNSTREAM_SIZE_GROUP"].eq("GT_CUTOFF")].drop_duplicates("DOWNSTREAM_SV_ID").copy()
    events["CALLERS_DISPLAY"] = events["DOWNSTREAM_CALLER_COUNT"].fillna("Unknown").replace(".", "Unknown")
    folder = out / "callers_and_population"
    folder.mkdir(parents=True, exist_ok=True)
    counts = events.groupby(["CALLERS_DISPLAY", "DOWNSTREAM_AF_STATUS"], as_index=False).agg(UNIQUE_SVS=("DOWNSTREAM_SV_ID", "size"))
    counts.to_csv(folder / f"{safe_sample}_callers_and_population.tsv", sep="\t", index=False)
    caller_groups = sorted(events["CALLERS_DISPLAY"].unique(), key=lambda x: -float(x) if x != "Unknown" else float("inf"))
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    x = np.arange(len(caller_groups))
    bottom = np.zeros(len(caller_groups))
    for i, (state, label, color) in enumerate([("MEASURED_AT_OR_BELOW_CUTOFF", f"Measured AF ≤{settings['max_population_af']*100:g}%", "#009E73"), ("AF_UNKNOWN", "AF unknown", "#8064A2")]):
        sub = counts[counts["DOWNSTREAM_AF_STATUS"].eq(state)].set_index("CALLERS_DISPLAY")
        values = sub["UNIQUE_SVS"].reindex(caller_groups, fill_value=0).to_numpy()
        ax.bar(x, values, bottom=bottom, width=0.65, label=label, color=color)
        bottom += values
    ax.set_xticks(x, [f"{c} callers" for c in caller_groups])
    ax.set_ylabel("Unique retained SVs")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(f"{sample}: combined low-AF or missing-AF SVs by caller support", loc="left", fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.15), ncol=2, frameon=False)
    if events.empty:
        ax.text(0.5, 0.5, "No large SVs pass these filters", transform=ax.transAxes, ha="center")
    fig.text(0.5, 0.01, "Each bar combines measured low AF and missing AF; color preserves the evidence source.\nMissing AF permits review for possible novelty; it does not prove absence from a database. Each SV is counted once.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.89, bottom=0.27)
    export(fig, folder, f"{safe_sample}_callers_and_population")
    # Inversions need their own simple support view, with unique event counts.
    inversions = events[events["DOWNSTREAM_SVTYPE"].eq("INV")]
    folder = out / "inversion_support"
    counts = inversions.groupby(["CALLERS_DISPLAY", "DOWNSTREAM_INV_READ_GROUP"], as_index=False).agg(UNIQUE_SVS=("DOWNSTREAM_SV_ID", "size"))
    folder.mkdir(parents=True, exist_ok=True)
    counts.to_csv(folder / f"{safe_sample}_inversion_support.tsv", sep="\t", index=False)
    groups = sorted(inversions["CALLERS_DISPLAY"].unique(), key=lambda c: -float(c) if c != "Unknown" else float("inf"))
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    bottom = np.zeros(len(groups))
    for state, label, color in [
        ("MINIMUM_TO_STRONG", f"{settings['inv_min_reads']}–{settings['inv_strong_reads']-1} reported reads", "#E69F00"),
        ("AT_OR_ABOVE_STRONG", f"≥{settings['inv_strong_reads']} reported reads", "#0072B2"),
    ]:
        sub = counts[counts["DOWNSTREAM_INV_READ_GROUP"].eq(state)].set_index("CALLERS_DISPLAY")
        values = sub["UNIQUE_SVS"].reindex(groups, fill_value=0).to_numpy()
        ax.bar(groups, values, bottom=bottom, color=color, label=label)
        bottom += values
    if not groups:
        ax.text(0.5, 0.5, "No large inversions pass these filters", transform=ax.transAxes, ha="center")
    ax.set_xlabel("Number of callers")
    ax.set_ylabel("Unique retained inversions")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(f"{sample}: inversion caller and read support", loc="left", fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.18), ncol=2, frameon=False)
    fig.text(0.5, 0.01, "Support is the largest reported count from one caller. These are review thresholds, not validation.\nAn inversion spanning an intact gene does not establish disruption of that gene.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.89, bottom=0.28)
    export(fig, folder, f"{safe_sample}_inversion_support")
    print(f"[OK] five separate downstream filter figures: {out}")


if __name__ == "__main__":
    main()
