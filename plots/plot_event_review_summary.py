#!/usr/bin/env python3
"""Evidence availability and predicted mechanisms for the full candidate table."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from candidate_review import DOMAINS, GROUP_LABELS, PANEL_GROUPS, STATE_COLORS, STATE_LABELS, evidence_table, normalise_events, text
from plot_utils import read_tsv, save_figure, set_thesis_style

EVENT_DOMAINS = {"Call support", "Population", "Local depth", "SV phase", "Straglr", "TLDR", "Methylation"}


def availability_counts(observations):
    rows = []
    for domain in DOMAINS:
        sub = observations[observations["DOMAIN"].eq(domain)].copy()
        unit = "UNIQUE_SV" if domain in EVENT_DOMAINS else "SV_GENE_ASSOCIATION"
        if unit == "UNIQUE_SV":
            # Event-wide evidence must agree across gene rows. Conflicting
            # propagated observations stay unresolved rather than picking one.
            sub = sub.groupby("SV_ID", as_index=False).agg(STATE=("STATE", lambda x: x.iloc[0] if x.nunique() == 1 else "UNRESOLVED"))
        total = len(sub)
        for state in STATE_COLORS:
            count = int(sub["STATE"].eq(state).sum())
            rows.append({"DOMAIN": domain, "COUNTING_UNIT": unit, "STATE": state, "COUNT": count, "TOTAL": total, "PERCENT": 100 * count / total if total else 0})
    return pd.DataFrame(rows)


def mechanism_category(row):
    effect = text(row, ["SV_GENE_EFFECT", "SV_GENE_RELATIONSHIP"]).upper()
    svtype = text(row, ["SVTYPE"]).upper()
    if effect in {"WHOLE_GENE_DELETION", "WHOLE_GENE_DOSAGE_CONTEXT"} and svtype == "DEL":
        return "Whole-gene loss"
    if effect in {"WHOLE_GENE_DUPLICATION", "WHOLE_GENE_DOSAGE_CONTEXT"} and svtype == "DUP":
        return "Whole-gene copy gain"
    if "PARTIAL_GENE" in effect:
        return "Partial gene loss" if svtype == "DEL" else "Partial gene duplication" if svtype == "DUP" else "Partial gene overlap"
    if "SPANNED_BY_INVERSION" in effect or "SPANS_INTACT_GENE" in effect:
        return "Gene spanned by inversion"
    if "BREAKPOINT_WITHIN" in effect or "BREAKPOINT_IN_GENE" in effect:
        return "Breakpoint within gene"
    if "INSERTION_WITHIN" in effect or effect == "INSERTION_IN_GENE":
        return "Insertion within gene"
    if "INTRON" in effect or text(row, ["VEP_TRANSCRIPT_REGION_CLASS"]).upper() == "INTRONIC_ONLY":
        return "Intronic context"
    if "PROXIMAL" in effect or "NEAR_GENE" in effect:
        return "Proximal context"
    if "CONTEXT" in effect or effect == "GENE_FULLY_SPANNED":
        return "Other interval context"
    return "Unresolved / other"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample", required=True)
    args = p.parse_args()
    set_thesis_style()
    work = normalise_events(read_tsv(args.input))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    availability_out = out / "evidence_availability"
    mechanism_out = out / "mechanisms"
    availability_out.mkdir(exist_ok=True)
    mechanism_out.mkdir(exist_ok=True)
    counts = availability_counts(evidence_table(work))
    counts.to_csv(availability_out / f"{args.sample}_evidence_availability.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(13, 8.5))
    left = np.zeros(len(DOMAINS))
    for state in STATE_COLORS:
        values = counts[counts["STATE"].eq(state)].set_index("DOMAIN")["PERCENT"].reindex(DOMAINS).to_numpy()
        ax.barh(DOMAINS, values, left=left, label=STATE_LABELS[state], color=STATE_COLORS[state], edgecolor="white", height=0.7)
        left += values
    labels = [d + (" [SV]" if d in EVENT_DOMAINS else " [SV–gene]") for d in DOMAINS]
    ax.set_yticks(range(len(DOMAINS)), labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Percentage of the stated counting unit")
    ax.set_title(f"{args.sample}: evidence availability in the candidate table", loc="left", fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.22, -0.15), ncol=3, frameon=False, fontsize=9)
    fig.text(0.5, 0.012, "SV-wide domains count each SV once. Gene-specific domains count SV–gene associations.\nUnavailable means absent or unevaluable in this table; it is not evidence against the event. Conflicting event-wide rows are unresolved.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.23, right=0.97, top=0.92, bottom=0.26)
    save_figure(fig, availability_out / f"{args.sample}_evidence_availability", dpi=300)
    plt.close(fig)

    work["PREDICTED_MECHANISM"] = [mechanism_category(row) for _, row in work.iterrows()]
    summary = work.groupby(["_panel_group", "PREDICTED_MECHANISM"], as_index=False).agg(SV_GENE_ASSOCIATIONS=("SV_ID", "size"), UNIQUE_SVS=("SV_ID", "nunique"), UNIQUE_GENES=("GENE", "nunique"))
    summary = summary.rename(columns={"_panel_group": "PANEL_STATUS"})
    summary.to_csv(mechanism_out / f"{args.sample}_mechanism_summary.tsv", sep="\t", index=False)
    categories = list(dict.fromkeys(summary.sort_values("SV_GENE_ASSOCIATIONS", ascending=False)["PREDICTED_MECHANISM"]))
    fig, ax = plt.subplots(figsize=(13, max(5.3, len(categories)*0.48+2)))
    left = np.zeros(len(categories))
    for group, color in zip(PANEL_GROUPS, ["#0072B2", "#E69F00", "#999999"]):
        sub = summary[summary["PANEL_STATUS"].eq(group)].set_index("PREDICTED_MECHANISM")
        values = sub["SV_GENE_ASSOCIATIONS"].reindex(categories, fill_value=0).to_numpy()
        ax.barh(categories, values, left=left, label=GROUP_LABELS[group], color=color, height=0.7)
        left += values
    for i, value in enumerate(left):
        ax.text(value, i, f"  {int(value):,}", va="center", fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, max(left, default=0)*1.15+1)
    ax.set_xlabel("Unique SV–gene associations")
    ax.set_title(f"{args.sample}: predicted SV–gene mechanisms", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper right")
    fig.text(0.5, 0.012, "Categories describe annotated geometry, not demonstrated loss of function. One SV may affect several genes and contribute to several categories.\nPartial gene loss does not establish coding-exon loss. Unique SV counts per category are provided separately in the TSV and are not additive.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.29, right=0.97, top=0.9, bottom=0.19)
    save_figure(fig, mechanism_out / f"{args.sample}_mechanism_summary", dpi=300)
    plt.close(fig)
    print(f"[OK] {len(work)} unique SV-gene associations; review summaries: {out}")


if __name__ == "__main__":
    main()
