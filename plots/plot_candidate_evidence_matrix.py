#!/usr/bin/env python3
"""Plot separate, annotated evidence matrices for panel/non-panel SV-gene rows."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np

from candidate_review import (
    DOMAINS, domains_for, PANEL_GROUPS, GROUP_LABELS, GROUP_SUFFIXES, STATE_COLORS,
    STATE_LABELS, enrich_events, event_evidence, evidence_table, select_events, text,
)
from plot_utils import read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--evidence-table", help="Optional detailed integrated table, joined by SV_ID and gene")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25, help="Rows per panel group")
    p.add_argument("--rare-af", type=float, default=0.01, help="AF review threshold; not a pathogenicity threshold")
    p.add_argument("--platform", choices=["lrs", "srs"], default="lrs")
    p.add_argument("--title", default="Candidate event evidence")
    return p.parse_args()


def cell_text(label):
    replacements = {
        "direct effect recessive second allele required": "Second allele required",
        "direct effect moi unknown": "Inheritance unknown",
        "mixed ad ar disease model review": "Mixed AD/AR model",
        "gene fully spanned by inversion": "Inversion spans gene",
        "inversion spans intact gene": "Inversion spans gene",
        "breakpoint within transcript": "Transcript breakpoint",
        "whole gene dosage context": "Whole-gene dosage",
    }
    label = replacements.get(label, label)
    lines = []
    for line in label.splitlines():
        lines.extend(textwrap.wrap(line, width=14) or [""])
    if len(lines) > 4:
        lines = lines[:4]
        lines[-1] = lines[-1][:11] + "..."
    return "\n".join(lines)


def render_matrix(work, prefix, title, rare_af=0.01, platform="lrs"):
    if work.empty:
        return []
    states = list(STATE_COLORS)
    domains = domains_for(platform)
    observations = [event_evidence(row, rare_af, platform) for _, row in work.iterrows()]
    matrix = np.array([[states.index(obs[d][0]) for d in domains] for obs in observations])
    height = max(3.3, len(work) * 0.67 + 2.25)
    fig, ax = plt.subplots(figsize=(17 if platform == "srs" else 19, height))
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=ListedColormap(list(STATE_COLORS.values())), vmin=-0.5, vmax=len(states)-0.5)
    for i, obs in enumerate(observations):
        for j, domain in enumerate(domains):
            ax.text(j, i, cell_text(obs[domain][1]), ha="center", va="center", fontsize=9, color="#25364A")
    labels = []
    for _, row in work.iterrows():
        locus = text(row, ["CHROM"]) + ":" + text(row, ["START", "POS"])
        rank = row["_rank"]
        rank_label = str(int(rank)) if np.isfinite(rank) else "?"
        labels.append(f"{rank_label}. {row['GENE']} | {text(row, ['SVTYPE'])}\n{locus} | {row['SV_ID'][:32]}")
    ax.set_yticks(range(len(work)), labels, fontsize=10)
    ax.set_xticks(range(len(domains)), [d.replace(" ", "\n", 1) for d in domains], fontsize=11)
    ax.tick_params(axis="both", length=0)
    ax.set_xticks(np.arange(-0.5, len(domains)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(work)), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)
    ax.set_ylabel("SV–gene association", fontsize=11)
    ax.set_title(title, loc="left", pad=72, fontsize=15, fontweight="bold")
    ax.legend(handles=[Patch(facecolor=STATE_COLORS[s], edgecolor="#FFFFFF", label=STATE_LABELS[s]) for s in states], loc="lower left", bbox_to_anchor=(0, 1.01), ncol=4, frameon=False, fontsize=9)
    context_note = "Repeat-locus overlap is context, not proof of an expansion." if platform == "srs" else "Gene relevance and methylation are context."
    fig.text(0.5, 0.015, f"Colours describe evidence states, not pathogenicity. {context_note} Depth and caller agreement use the same sequencing data.\nPopulation cells show the largest explicit event AF across available sources. Full observations and identifiers are in the accompanying TSV.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.26, right=0.985, bottom=0.8/height, top=1-1.45/height)
    result = save_figure(fig, prefix, dpi=300)
    plt.close(fig)
    return result


def main():
    args = parse_args()
    set_thesis_style()
    work = enrich_events(select_events(read_tsv(args.input), args.top_n), args.evidence_table)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    observations = evidence_table(work, args.rare_af, args.platform)
    observations.to_csv(str(prefix) + "_observations.tsv", sep="\t", index=False)
    # Keep earlier exported caller fields, but do not show inapplicable callers.
    work["LongPhase match"] = [1 if text(r, ["LONGPHASE", "LONGPHASE_MATCH"]).upper() == "YES" else 0 if text(r, ["LONGPHASE", "LONGPHASE_MATCH"]).upper() == "NO" else np.nan for _, r in work.iterrows()]
    callers = work.get("CALLERS")
    if callers is not None:
        for caller in ["Sniffles2", "cuteSV", "Delly", "Manta"]:
            if callers.fillna("").astype(str).str.contains(caller, case=False, regex=False).any():
                work[caller] = callers.fillna("").astype(str).str.contains(caller, case=False, regex=False).astype(int)
    work.drop(columns=["_source_order"], errors="ignore").to_csv(str(prefix) + "_matrix.tsv", sep="\t", index=False)
    outputs = []
    for group in PANEL_GROUPS:
        subset = work[work["_panel_group"].eq(group)].copy()
        path = prefix.with_name(prefix.name + "_" + GROUP_SUFFIXES[group])
        outputs.extend(render_matrix(subset, path, args.title + ": " + GROUP_LABELS[group], args.rare_af, args.platform))
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()

