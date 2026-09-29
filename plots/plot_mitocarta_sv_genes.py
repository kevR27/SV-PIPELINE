#!/usr/bin/env python3
"""Plot nuclear-encoded mitochondrial genes intersected by master SVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="MitoCarta-annotated integrated TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=20)
    return p.parse_args()


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    score_col = first_existing(df, ["ALLELE_RESEARCH_SCORE", "INTEGRATED_DISCOVERY_SCORE", "PHENOTYPE_SCORE"])
    panel_col = first_existing(df, ["PANEL_STATUS", "panel_gene"])
    if id_col is None or gene_col is None or "MITOCARTA_ENCODING" not in df:
        raise ValueError("Input is not a MitoCarta-annotated integrated table.")

    work = df[df["MITOCARTA_ENCODING"].eq("NUCLEAR_MITOCHONDRIAL_GENE")].copy()

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    if work.empty:
        pd.DataFrame(columns=["SV_ID", "gene"]).to_csv(
            prefix.with_name(prefix.name + "_sv_gene_pairs.tsv"), sep="\t", index=False
        )
        pd.DataFrame(columns=["gene", "unique_SVs"]).to_csv(
            prefix.with_name(prefix.name + "_genes.tsv"), sep="\t", index=False
        )
        pd.DataFrame(columns=["pathway", "unique_SVs", "unique_genes"]).to_csv(
            prefix.with_name(prefix.name + "_pathways.tsv"), sep="\t", index=False
        )
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.axis("off")
        ax.text(
            0.5, 0.55,
            "No nuclear-encoded MitoCarta3.0 genes overlap the master SV-gene set",
            ha="center", va="center", fontsize=13,
        )
        ax.text(
            0.5, 0.42,
            "This is a valid negative result, not a missing-data state.",
            ha="center", va="center", fontsize=10,
        )
        fig.suptitle("MitoCarta3.0 context for SV-overlapping genes", fontsize=16, fontweight="bold")
        outputs = save_figure(fig, prefix)
        plt.close(fig)
        print("[OK] no nuclear MitoCarta SV-gene overlaps")
        print("[OK]", *outputs, sep="\n")
        return

    work["_priority"] = numeric(work[score_col]).fillna(0) if score_col else 0
    if panel_col:
        ptxt = work[panel_col].fillna("").astype(str).str.upper().str.strip()
        work["_panel"] = ptxt.isin(["PANEL_GENE", "YES", "TRUE", "1"])
    else:
        work["_panel"] = False

    pairs = work.drop_duplicates([id_col, gene_col], keep="first")
    gene_summary = (
        pairs.groupby(gene_col)
        .agg(
            unique_SVs=(id_col, "nunique"),
            priority=("_priority", "max"),
            ON_context=("MITO_ON_CONTEXT", lambda s: int((s == "YES").any())),
            panel=("_panel", "max"),
            pathways=("MITOCARTA_MITOPATHWAYS", lambda s: ";".join(sorted({x for x in s if str(x) not in {".", "nan"}}))),
            subcompartment=("MITOCARTA_SUBCOMPARTMENT", lambda s: ";".join(sorted({x for x in s if str(x) not in {".", "nan"}}))),
        )
        .reset_index()
        .sort_values(["ON_context", "priority", "unique_SVs"], ascending=False)
        .head(args.top_n)
    )

    pathway_rows = []
    for _, row in pairs.iterrows():
        raw = str(row.get("MITOCARTA_TOP_LEVEL_PATHWAYS", "."))
        for pathway in [x.strip() for x in raw.split(";") if x.strip() and x.strip() != "."]:
            pathway_rows.append((pathway, row[id_col], row[gene_col]))
    pathway_df = pd.DataFrame(pathway_rows, columns=["pathway", "SV_ID", "gene"])
    if not pathway_df.empty:
        pathway_summary = (
            pathway_df.drop_duplicates()
            .groupby("pathway")
            .agg(unique_SVs=("SV_ID", "nunique"), unique_genes=("gene", "nunique"))
            .reset_index()
            .sort_values(["unique_genes", "unique_SVs"], ascending=False)
        )
    else:
        pathway_summary = pd.DataFrame(columns=["pathway", "unique_SVs", "unique_genes"])

    pairs.to_csv(prefix.with_name(prefix.name + "_sv_gene_pairs.tsv"), sep="\t", index=False)
    gene_summary.to_csv(prefix.with_name(prefix.name + "_genes.tsv"), sep="\t", index=False)
    pathway_summary.to_csv(prefix.with_name(prefix.name + "_pathways.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(15.0, max(7.2, 0.38 * len(gene_summary) + 2.6)))
    ax1, ax2 = axes

    gs = gene_summary.iloc[::-1].reset_index(drop=True)
    colors = np.where(gs["panel"], "#0072B2", "#E69F00")
    ax1.barh(np.arange(len(gs)), gs["priority"], color=colors)
    ax1.set_yticks(np.arange(len(gs)))
    labels = [
        f"{g}  (SV={n}; ON={'yes' if on else 'no'})"
        for g, n, on in zip(gs[gene_col], gs["unique_SVs"], gs["ON_context"])
    ]
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel((score_col or "priority").replace("_", " "))
    ax1.set_title("Nuclear mitochondrial SV-overlapping genes")
    style_axis(ax1, "x")

    ps = pathway_summary.head(12).iloc[::-1]
    ax2.barh(np.arange(len(ps)), ps["unique_genes"])
    ax2.set_yticks(np.arange(len(ps)))
    ax2.set_yticklabels(ps["pathway"], fontsize=9)
    ax2.set_xlabel("Unique SV-overlapping nuclear mitochondrial genes")
    ax2.set_title("MitoCarta top-level pathways")
    style_axis(ax2, "x")

    fig.suptitle("MitoCarta3.0 context for SV-overlapping genes", fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5, 0.01,
        "Blue=optic-neuropathy panel gene; orange=non-panel gene. ON context means the existing pipeline reports positive optic-neuropathy phenotype evidence for that gene. MitoCarta membership/localization is independent evidence and does not establish SV causality.",
        ha="center", fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
