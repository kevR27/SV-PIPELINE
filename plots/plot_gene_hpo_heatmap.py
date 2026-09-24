#!/usr/bin/env python3
"""Plot human gene-HPO associations retrieved by the offline Monarch branch."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd

from plot_utils import first_existing, read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description="Plot gene-HPO association heatmap.")
    p.add_argument("--input", required=True, help="*_human_gene_phenotypes.tsv")
    p.add_argument("--ranking", default=None, help="Optional *_ranked_candidates.tsv used to select/order genes")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-genes", type=int, default=20)
    p.add_argument("--top-hpo", type=int, default=24)
    p.add_argument("--title", default="Human gene–phenotype associations for SV-overlapping genes")
    return p.parse_args()


def make_heatmap(matrix, xlabels, anchor_hpos, out_prefix, title, hpos):
    cmap = ListedColormap(["#F3F4F4", "#0B6E69"])
    fig_w = max(11.5, 0.55 * len(hpos) + 4.5)
    fig_h = max(6.5, 0.38 * len(matrix.index) + 2.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(matrix.values, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)

    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=10)
    ax.set_xticks(np.arange(len(hpos)))
    ax.set_xticklabels(xlabels, rotation=50, ha="right", fontsize=8.5)
    ax.set_xlabel("Human Phenotype Ontology term")
    ax.set_ylabel("SV-overlapping gene")
    ax.set_title(title, fontweight="bold", fontsize=14)

    for j, hpo in enumerate(hpos):
        if hpo in anchor_hpos:
            ax.get_xticklabels()[j].set_fontweight("bold")

    ax.set_xticks(np.arange(-0.5, len(hpos), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.text(
        0.5,
        0.008,
        "Dark teal = documented human gene–HPO association. Bold labels = optic-neuropathy anchor terms.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    outputs = save_figure(fig, out_prefix)
    plt.close(fig)
    return outputs


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    gene_col = first_existing(df, ["gene_symbol", "gene", "Gene"])
    hpo_col = first_existing(df, ["hpo_id", "HPO", "HPO_ID"])
    label_col = first_existing(df, ["hpo_label", "HPO_label", "phenotype"])
    anchor_col = first_existing(df, ["optic_neuropathy_anchor"])
    if gene_col is None or hpo_col is None:
        raise ValueError("Phenotype table requires gene_symbol and hpo_id columns.")

    work = df.copy()
    work = work[work[hpo_col].fillna("").astype(str).str.startswith("HP:")]
    if work.empty:
        raise ValueError("No human HPO associations found in the phenotype table.")

    if args.ranking:
        rank = read_tsv(args.ranking)
        rgene = first_existing(rank, ["gene", "Gene", "GENE"])
        score = first_existing(rank, ["integrated_discovery_score", "phenotype_score"])
        if rgene and score:
            rank["_score"] = pd.to_numeric(rank[score], errors="coerce").fillna(0)
            genes = rank.sort_values("_score", ascending=False)[rgene].astype(str).head(args.top_genes).tolist()
        else:
            genes = work[gene_col].value_counts().head(args.top_genes).index.tolist()
    else:
        genes = work[gene_col].value_counts().head(args.top_genes).index.tolist()

    work = work[work[gene_col].isin(genes)].copy()

    anchor_hpos = set()
    if anchor_col:
        anchor_mask = work[anchor_col].fillna("0").astype(str).str.upper().isin(["1", "YES", "TRUE"])
        anchor_hpos = set(work.loc[anchor_mask, hpo_col])

    counts = work[hpo_col].value_counts()
    ordered_anchor = [h for h in counts.index if h in anchor_hpos]
    ordered_other = [h for h in counts.index if h not in anchor_hpos]
    hpos = (ordered_anchor + ordered_other)[: args.top_hpo]
    work_main = work[work[hpo_col].isin(hpos)].copy()

    matrix = pd.crosstab(work_main[gene_col], work_main[hpo_col]).clip(upper=1)
    matrix = matrix.reindex(index=[g for g in genes if g in matrix.index], columns=hpos, fill_value=0)

    label_map = {}
    if label_col:
        label_map = (
            work.drop_duplicates(hpo_col)
            .set_index(hpo_col)[label_col]
            .fillna("")
            .astype(str)
            .to_dict()
        )

    def label_for(h):
        label = label_map.get(h, "")
        wrapped = "\n".join(textwrap.wrap(label, width=20)) if label else ""
        return f"{h}\n{wrapped}" if wrapped else h

    xlabels = [label_for(h) for h in hpos]

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(prefix.with_name(prefix.name + "_matrix.tsv"), sep="\t")

    outputs = make_heatmap(matrix, xlabels, anchor_hpos, prefix, args.title, hpos)

    anchor_present = [h for h in ordered_anchor if h in matrix.columns]
    if anchor_present:
        anchor_matrix = matrix.reindex(columns=anchor_present, fill_value=0)
        anchor_prefix = prefix.with_name(prefix.name + "_anchors")
        anchor_matrix.to_csv(anchor_prefix.with_name(anchor_prefix.name + "_matrix.tsv"), sep="\t")
        make_heatmap(
            anchor_matrix,
            [label_for(h) for h in anchor_present],
            set(anchor_present),
            anchor_prefix,
            "Optic-neuropathy anchor HPO associations",
            anchor_present,
        )

    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
