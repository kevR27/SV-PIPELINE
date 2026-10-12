#!/usr/bin/env python3
"""Plot gene-HPO associations and the separate custom/Exomiser comparison."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from plot_utils import first_existing, read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description="Plot gene-HPO association heatmap.")
    p.add_argument("--input", help="*_human_gene_phenotypes.tsv")
    p.add_argument("--ranking", default=None, help="Optional *_ranked_candidates.tsv used to select/order genes")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-genes", type=int, default=25)
    p.add_argument("--top-hpo", type=int, default=18)
    p.add_argument("--title", default="Human gene–phenotype associations for SV-overlapping genes")
    p.add_argument("--exomiser-comparison", help="Separate *_exomiser_comparison.tsv")
    p.add_argument("--comparison-only", action="store_true", help="Plot the comparison without an HPO heatmap")
    return p.parse_args()


def make_exomiser_comparison(table, out_prefix, top_genes=25):
    """Show rank changes with two colours; never average the two scores."""
    required = {"GENE", "COMPARISON_STATUS", "CUSTOM_QUERY_RANK", "EXOMISER_QUERY_RANK", "QUERY_SCOPE"}
    if not required.issubset(table.columns):
        raise ValueError("Exomiser comparison table is missing required columns.")
    prefix = Path(out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work = table[table.COMPARISON_STATUS.eq("PAIRED")].copy()
    for column in ["CUSTOM_QUERY_RANK", "EXOMISER_QUERY_RANK"]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna(subset=["CUSTOM_QUERY_RANK", "EXOMISER_QUERY_RANK"])
    work["_best_rank"] = work[["CUSTOM_QUERY_RANK", "EXOMISER_QUERY_RANK"]].min(axis=1)
    selected = work.sort_values(["_best_rank", "CUSTOM_QUERY_RANK", "GENE"]).head(top_genes)
    shown = set(selected.GENE)
    source = table.copy()
    source["SHOWN_IN_PLOT"] = source.GENE.isin(shown).map({True: "YES", False: "NO"})
    source.to_csv(prefix.with_name(prefix.name + "_source.tsv"), sep="\t", index=False)
    scope_labels = {"SHARED_HON_PROFILE": "Shared HON profile",
                    "PATIENT_HPO_PROFILE": "Patient HPO profile",
                    "EXOMISER_DISABLED": "Exomiser disabled"}
    scope = ", ".join(scope_labels.get(value, value) for value in sorted(set(table.QUERY_SCOPE))) if not table.empty else "No candidate genes"
    fig, ax = plt.subplots(figsize=(11.5, max(5, 0.36 * len(selected) + 2.5)))
    if selected.empty:
        explanations = []
        for label, column in [("Exomiser", "EXOMISER_STATUS"), ("Custom method", "CUSTOM_QUERY_STATUS")]:
            statuses = ", ".join(sorted(set(table.get(column, pd.Series(dtype=str)))))
            if statuses:
                explanations.append(textwrap.fill(label + ": " + statuses.replace("_", " ").lower(), 70))
        ax.axis("off")
        ax.text(0.5, 0.55, "No comparable gene–phenotype scores", ha="center", va="center", fontsize=15)
        ax.text(0.5, 0.40, "\n".join(explanations) or "No candidate genes", ha="center", va="center", fontsize=10)
    else:
        y = np.arange(len(selected))
        # A small vertical offset keeps both methods visible when ranks agree.
        for i, (_, row) in enumerate(selected.iterrows()):
            ax.plot([row.CUSTOM_QUERY_RANK, row.EXOMISER_QUERY_RANK], [i - 0.08, i + 0.08], color="#C7CCCF", linewidth=1.3, zorder=1)
        ax.scatter(selected.CUSTOM_QUERY_RANK, y - 0.08, color="#0B6E69", s=50, label="Custom phenotype method", zorder=2)
        backend = set(selected.get("EXOMISER_BACKEND", pd.Series(dtype=str)).dropna())
        label = "Exomiser offline CLI" if backend == {"OFFLINE_CLI_PHENOTYPE_ONLY"} else "Exomiser phenotype"
        ax.scatter(selected.EXOMISER_QUERY_RANK, y + 0.08, color="#775DA6", marker="D", s=38, label=label, zorder=3)
        labels = [str(row.GENE) + (" [panel]" if row.get("PANEL_STATUS") == "PANEL_GENE" else "")
                  for _, row in selected.iterrows()]
        ax.set_yticks(y, labels)
        ax.set_ylim(len(selected) - 0.5, -0.5)
        ax.set_xlim(0.5, max(selected.CUSTOM_QUERY_RANK.max(), selected.EXOMISER_QUERY_RANK.max()) + 0.5)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlabel("Phenotype rank among the same comparable genes (1 = highest)")
        ax.set_ylabel("SV-overlapping gene")
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=9)
        ax.grid(axis="x", alpha=0.18)
    ax.set_title("Custom HON method compared with Exomiser", fontsize=14, fontweight="bold", pad=44 if not selected.empty else 12)
    fig.text(0.5, 0.01, f"{scope} · {len(work)} comparable genes; showing up to {top_genes}.\n"
             "Same HPO query; phenotype-only ranks. The original integrated ranking is unchanged.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.065, 1, 1])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    return outputs


def make_heatmap(matrix, xlabels, anchor_hpos, out_prefix, title, hpos):
    cmap = ListedColormap(["#F3F4F4", "#0B6E69"])
    fig_w = max(12.0, 0.72 * len(hpos) + 4.5)
    fig_h = max(6.5, 0.38 * len(matrix.index) + 2.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(matrix.values, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)

    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=10)
    ax.set_xticks(np.arange(len(hpos)))
    ax.set_xticklabels(xlabels, rotation=55, ha="right", fontsize=8.5)
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
    if args.top_genes < 1:
        raise ValueError("--top-genes must be at least 1.")
    if args.exomiser_comparison:
        prefix = args.out_prefix if args.comparison_only else args.out_prefix + "_exomiser_comparison"
        outputs = make_exomiser_comparison(read_tsv(args.exomiser_comparison), prefix, args.top_genes)
        print("[OK]", *outputs, sep="\n")
    if args.comparison_only:
        if not args.exomiser_comparison:
            raise ValueError("--comparison-only needs --exomiser-comparison.")
        return
    if not args.input:
        raise ValueError("--input is required for the gene-HPO heatmap.")
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

    panel_genes = []
    nonpanel_genes = []
    if args.ranking:
        rank = read_tsv(args.ranking)
        rgene = first_existing(rank, ["gene", "Gene", "GENE"])
        score = first_existing(
            rank,
            [
                "FINAL_GENE_RELEVANCE_DISPLAY_SCORE",
                "GENE_RELEVANCE_DISPLAY_SCORE",
                "integrated_discovery_score",
                "phenotype_score",
            ],
        )
        panel_col = first_existing(rank, ["PANEL_STATUS", "panel_gene"])
        rank_col = first_existing(
            rank,
            [
                "FINAL_GENE_RANK_WITHIN_PANEL_STATUS",
                "GENE_RANK_WITHIN_PANEL_STATUS",
            ],
        )
        if rgene and score and panel_col:
            rank["_score"] = pd.to_numeric(
                rank[score], errors="coerce"
            ).fillna(0)
            rank["_panel_group"] = np.where(
                rank[panel_col]
                .fillna("")
                .astype(str)
                .str.upper()
                .isin(["PANEL_GENE", "YES", "TRUE", "1"]),
                "PANEL_GENE",
                "NONPANEL_GENE",
            )
            rank["_rank"] = (
                pd.to_numeric(rank[rank_col], errors="coerce")
                if rank_col
                else np.nan
            )
            rank = rank.sort_values(
                ["_panel_group", "_rank", "_score", rgene],
                ascending=[True, True, False, True],
                na_position="last",
            )
            panel_genes = (
                rank[rank["_panel_group"].eq("PANEL_GENE")][rgene]
                .astype(str)
                .head(args.top_genes)
                .tolist()
            )
            nonpanel_genes = (
                rank[rank["_panel_group"].eq("NONPANEL_GENE")][rgene]
                .astype(str)
                .head(args.top_genes)
                .tolist()
            )
            genes = panel_genes + nonpanel_genes
        else:
            genes = (
                work[gene_col]
                .value_counts()
                .head(args.top_genes)
                .index.tolist()
            )
    else:
        genes = (
            work[gene_col]
            .value_counts()
            .head(args.top_genes)
            .index.tolist()
        )

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

    outputs = make_heatmap(
        matrix,
        xlabels,
        anchor_hpos,
        prefix,
        args.title,
        hpos,
    )

    # Separate panel/non-panel HPO figures when a ranking table supplies both
    # groups. The combined matrix is kept as a compact overview.
    for label, selected in [
        ("panel", panel_genes),
        ("nonpanel", nonpanel_genes),
    ]:
        selected = [gene for gene in selected if gene in matrix.index]
        if not selected:
            continue
        group_matrix = matrix.reindex(index=selected)
        group_prefix = prefix.with_name(prefix.name + "_" + label)
        group_matrix.to_csv(
            group_prefix.with_name(group_prefix.name + "_matrix.tsv"),
            sep="\t",
        )
        make_heatmap(
            group_matrix,
            xlabels,
            anchor_hpos,
            group_prefix,
            (
                "Panel gene–HPO associations"
                if label == "panel"
                else "Non-panel gene–HPO associations"
            ),
            hpos,
        )

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
