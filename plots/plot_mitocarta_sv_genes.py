#!/usr/bin/env python3
"""Plot MitoCarta context for genes intersected by structural variants."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, numeric, read_tsv, save_figure, set_thesis_style, style_axis


MISSING = {"", ".", "NA", "N/A", "nan", "None"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="MitoCarta-annotated integrated TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25)
    return p.parse_args()


def split_values(series):
    values = set()
    for raw in series.fillna(".").astype(str):
        for item in raw.split(";"):
            item = item.strip()
            if item not in MISSING:
                values.add(item)
    return sorted(values)


def relevance_score(df):
    """Return a transparent gene-relevance value for display only."""
    direct = first_existing(
        df,
        [
            "GENE_RELEVANCE_SCORE",
            "EVENT_GENE_RELEVANCE_SCORE",
            "INTEGRATED_DISCOVERY_SCORE",
            "integrated_discovery_score",
        ],
    )
    if direct:
        return numeric(df[direct]).fillna(0), direct

    phenotype = first_existing(df, ["PHENOTYPE_SCORE"])
    disease = first_existing(df, ["GENE_DISEASE_EVIDENCE_SCORE"])

    phenotype_values = (
        numeric(df[phenotype]).fillna(0)
        if phenotype
        else pd.Series(0.0, index=df.index)
    )
    disease_values = (
        numeric(df[disease]).fillna(0)
        if disease
        else pd.Series(0.0, index=df.index)
    )
    return phenotype_values + disease_values, "phenotype + gene-disease evidence"


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    panel_col = first_existing(df, ["PANEL_STATUS", "panel_gene"])

    required = {"MITOCARTA_STATUS", "MITOCARTA_ENCODING"}
    if id_col is None or gene_col is None or not required.issubset(df.columns):
        raise ValueError("Input is not a MitoCarta-annotated SV-gene table.")

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    mito_mask = df["MITOCARTA_STATUS"].fillna("NO").astype(str).str.upper().eq("YES")
    work = df[mito_mask].copy()

    # The thesis analysis currently focuses on nuclear-encoded mitochondrial
    # genes. mtDNA-encoded genes remain identifiable in the diagnostic summary.
    nuclear = work[
        work["MITOCARTA_ENCODING"]
        .fillna("")
        .astype(str)
        .str.upper()
        .eq("NUCLEAR_MITOCHONDRIAL_GENE")
    ].copy()

    diagnostic = {
        "input_rows": len(df),
        "mitocarta_rows": len(work),
        "nuclear_mito_rows": len(nuclear),
        "unique_mitocarta_genes": int(work[gene_col].nunique()) if not work.empty else 0,
        "unique_nuclear_mito_genes": int(nuclear[gene_col].nunique()) if not nuclear.empty else 0,
        "rows_with_pathway": (
            int(
                ~work["MITOCARTA_MITOPATHWAYS"]
                .fillna(".")
                .astype(str)
                .isin([".", "", "nan", "None"])
            .sum()
            )
            if "MITOCARTA_MITOPATHWAYS" in work
            else 0
        ),
        "rows_with_subcompartment": (
            int(
                ~work["MITOCARTA_SUBCOMPARTMENT"]
                .fillna(".")
                .astype(str)
                .isin([".", "", "nan", "None"])
            .sum()
            )
            if "MITOCARTA_SUBCOMPARTMENT" in work
            else 0
        ),
    }
    pd.DataFrame([diagnostic]).to_csv(
        prefix.with_name(prefix.name + "_diagnostics.tsv"),
        sep="\t",
        index=False,
    )

    if nuclear.empty:
        pd.DataFrame(columns=["SV_ID", "gene"]).to_csv(
            prefix.with_name(prefix.name + "_sv_gene_pairs.tsv"),
            sep="\t",
            index=False,
        )
        fig, ax = plt.subplots(figsize=(11, 5.2))
        ax.axis("off")
        ax.text(
            0.5,
            0.62,
            "No nuclear MitoCarta genes were recovered from this SV-gene table",
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.47,
            (
                f"Input rows: {diagnostic['input_rows']:,} | "
                f"MitoCarta rows: {diagnostic['mitocarta_rows']:,} | "
                f"Unique MitoCarta genes: {diagnostic['unique_mitocarta_genes']}"
            ),
            ha="center",
            va="center",
            fontsize=10,
        )
        ax.text(
            0.5,
            0.34,
            "Check the companion *_diagnostics.tsv before interpreting this as a biological negative.",
            ha="center",
            va="center",
            fontsize=9.5,
        )
        fig.suptitle(
            "MitoCarta3.0 context for SV-overlapping genes",
            fontsize=16,
            fontweight="bold",
        )
        outputs = save_figure(fig, prefix)
        plt.close(fig)
        print(
            "[MitoCarta plot] no nuclear mitochondrial rows; "
            f"diagnostics={prefix.with_name(prefix.name + '_diagnostics.tsv')}"
        )
        print("[OK]", *outputs, sep="\n")
        return

    nuclear["_gene_relevance"], score_name = relevance_score(nuclear)

    if panel_col:
        panel_text = nuclear[panel_col].fillna("").astype(str).str.upper().str.strip()
        nuclear["_panel"] = panel_text.isin(["PANEL_GENE", "YES", "TRUE", "1"])
    else:
        nuclear["_panel"] = False

    pairs = nuclear.drop_duplicates([id_col, gene_col], keep="first").copy()

    gene_rows = []
    for gene, group in pairs.groupby(gene_col, sort=False):
        pathways = split_values(group["MITOCARTA_MITOPATHWAYS"]) if "MITOCARTA_MITOPATHWAYS" in group else []
        compartments = split_values(group["MITOCARTA_SUBCOMPARTMENT"]) if "MITOCARTA_SUBCOMPARTMENT" in group else []
        gene_rows.append(
            {
                "gene": gene,
                "unique_SVs": int(group[id_col].nunique()),
                "gene_relevance": float(group["_gene_relevance"].max()),
                "ON_context": int(
                    group.get("MITO_ON_CONTEXT", pd.Series("NO", index=group.index))
                    .fillna("NO")
                    .astype(str)
                    .str.upper()
                    .eq("YES")
                    .any()
                ),
                "panel": bool(group["_panel"].any()),
                "pathways": ";".join(pathways) if pathways else ".",
                "subcompartment": ";".join(compartments) if compartments else ".",
            }
        )

    gene_summary = pd.DataFrame(gene_rows).sort_values(
        ["ON_context", "gene_relevance", "unique_SVs", "gene"],
        ascending=[False, False, False, True],
    ).head(args.top_n)

    pathway_rows = []
    for _, row in pairs.iterrows():
        raw_top = str(row.get("MITOCARTA_TOP_LEVEL_PATHWAYS", "."))
        pathways = [
            item.strip()
            for item in raw_top.split(";")
            if item.strip() not in MISSING
        ]

        if not pathways:
            pathways = [
                item.strip()
                for item in str(row.get("MITOCARTA_MITOPATHWAYS", ".")).split(";")
                if item.strip() not in MISSING
            ]

        for pathway in pathways:
            pathway_rows.append((pathway, row[id_col], row[gene_col]))

    pathway_df = pd.DataFrame(pathway_rows, columns=["pathway", "SV_ID", "gene"])
    if not pathway_df.empty:
        pathway_summary = (
            pathway_df.drop_duplicates()
            .groupby("pathway")
            .agg(unique_SVs=("SV_ID", "nunique"), unique_genes=("gene", "nunique"))
            .reset_index()
            .sort_values(["unique_genes", "unique_SVs", "pathway"], ascending=[False, False, True])
        )
    else:
        pathway_summary = pd.DataFrame(columns=["pathway", "unique_SVs", "unique_genes"])

    compartment_rows = []
    for _, row in pairs.iterrows():
        for compartment in [
            item.strip()
            for item in str(row.get("MITOCARTA_SUBCOMPARTMENT", ".")).split(";")
            if item.strip() not in MISSING
        ]:
            compartment_rows.append((compartment, row[id_col], row[gene_col]))

    compartment_df = pd.DataFrame(
        compartment_rows,
        columns=["subcompartment", "SV_ID", "gene"],
    )
    if not compartment_df.empty:
        compartment_summary = (
            compartment_df.drop_duplicates()
            .groupby("subcompartment")
            .agg(unique_SVs=("SV_ID", "nunique"), unique_genes=("gene", "nunique"))
            .reset_index()
            .sort_values(["unique_genes", "unique_SVs", "subcompartment"], ascending=[False, False, True])
        )
    else:
        compartment_summary = pd.DataFrame(
            columns=["subcompartment", "unique_SVs", "unique_genes"]
        )

    pairs.to_csv(
        prefix.with_name(prefix.name + "_sv_gene_pairs.tsv"),
        sep="\t",
        index=False,
    )
    gene_summary.to_csv(
        prefix.with_name(prefix.name + "_genes.tsv"),
        sep="\t",
        index=False,
    )
    pathway_summary.to_csv(
        prefix.with_name(prefix.name + "_pathways.tsv"),
        sep="\t",
        index=False,
    )
    compartment_summary.to_csv(
        prefix.with_name(prefix.name + "_subcompartments.tsv"),
        sep="\t",
        index=False,
    )

    fig_height = max(8.0, 0.34 * len(gene_summary) + 3.0)
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(17.0, fig_height),
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )
    ax1, ax2 = axes

    gs = gene_summary.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(gs))

    # Always plot an integer count so the panel cannot disappear when a score
    # column is missing or zero. The relevance score is shown as text.
    bars = ax1.barh(y, gs["unique_SVs"])
    ax1.set_yticks(y)
    labels = [
        f"{gene}  ({'panel' if panel else 'non-panel'}; ON={'yes' if on else 'no'})"
        for gene, panel, on in zip(gs["gene"], gs["panel"], gs["ON_context"])
    ]
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel("Unique structural variants overlapping gene")
    ax1.set_title("Nuclear mitochondrial genes affected by SVs")
    style_axis(ax1, "x")

    xmax = max(float(gs["unique_SVs"].max()), 1.0)
    for bar, (_, row) in zip(bars, gs.iterrows()):
        ax1.text(
            bar.get_width() + 0.02 * xmax,
            bar.get_y() + bar.get_height() / 2,
            f"score={row['gene_relevance']:.1f}",
            va="center",
            fontsize=8,
        )
    ax1.set_xlim(0, xmax * 1.35)

    ps = pathway_summary.head(15).iloc[::-1]
    if not ps.empty:
        ax2.barh(np.arange(len(ps)), ps["unique_genes"])
        ax2.set_yticks(np.arange(len(ps)))
        ax2.set_yticklabels(ps["pathway"], fontsize=8.5)
        ax2.set_xlabel("Unique mitochondrial genes")
        ax2.set_title("MitoCarta pathway context")
        style_axis(ax2, "x")
    else:
        cs = compartment_summary.head(15).iloc[::-1]
        if not cs.empty:
            ax2.barh(np.arange(len(cs)), cs["unique_genes"])
            ax2.set_yticks(np.arange(len(cs)))
            ax2.set_yticklabels(cs["subcompartment"], fontsize=8.5)
            ax2.set_xlabel("Unique mitochondrial genes")
            ax2.set_title("MitoCarta sub-mitochondrial localization")
            style_axis(ax2, "x")
        else:
            ax2.axis("off")
            ax2.text(
                0.5,
                0.58,
                "MitoCarta genes were found,",
                transform=ax2.transAxes,
                ha="center",
                fontsize=12,
                fontweight="bold",
            )
            ax2.text(
                0.5,
                0.47,
                "but no pathway or subcompartment annotations were parsed.",
                transform=ax2.transAxes,
                ha="center",
                fontsize=10,
            )
            ax2.text(
                0.5,
                0.36,
                "Review *_diagnostics.tsv and the configured MitoCarta reference files.",
                transform=ax2.transAxes,
                ha="center",
                fontsize=9,
            )

    fig.suptitle(
        "MitoCarta3.0 context for SV-overlapping genes",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.012,
        (
            "Left: number of unique SVs per nuclear mitochondrial gene; the text value is the "
            f"gene-relevance score ({score_name}). Right: MitoCarta pathway context, with "
            "sub-mitochondrial localization used only when pathway assignments are unavailable. "
            "MitoCarta membership does not establish SV causality."
        ),
        ha="center",
        fontsize=8.5,
        wrap=True,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)

    print(
        "[MitoCarta plot] "
        f"nuclear_genes={diagnostic['unique_nuclear_mito_genes']} "
        f"pathway_rows={diagnostic['rows_with_pathway']} "
        f"subcompartment_rows={diagnostic['rows_with_subcompartment']}"
    )
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
