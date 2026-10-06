#!/usr/bin/env python3
"""Create an integrated technical/biological evidence matrix for top SV-gene candidates."""

from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from plot_utils import add_panel_label, first_existing, numeric, read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description="Plot integrated evidence for prioritized SV-gene pairs.")
    p.add_argument("--input", required=True, help="Integrated or extended SV-gene analysis TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--rare-af", type=float, default=0.01)
    p.add_argument("--title", default="Integrated evidence for prioritized SV-gene candidates")
    return p.parse_args()


def text_present(series: pd.Series) -> pd.Series:
    s = series.fillna("").astype(str).str.strip()
    return ~s.isin(["", ".", "NA", "N/A", "None", "nan", "NaN"])


def caller_present(series: pd.Series, caller: str) -> pd.Series:
    return series.fillna("").astype(str).str.contains(re.escape(caller), case=False, regex=True)


def yes_flag(series: pd.Series) -> pd.Series:
    values = series.fillna("").astype(str).str.upper().str.strip()
    return values.map({"YES": 1.0, "TRUE": 1.0, "1": 1.0, "NO": 0.0, "FALSE": 0.0, "0": 0.0})


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Integrated SV-gene table is empty.")

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "Gene", "GENE", "ANNotsv_Gene"])
    chrom_col = first_existing(df, ["CHROM", "chrom"])
    start_col = first_existing(df, ["START", "POS"])
    type_col = first_existing(df, ["SVTYPE", "SV_type"])
    if id_col is None or gene_col is None:
        raise ValueError("Input needs an SV ID and gene column.")

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    if chrom_col and start_col and type_col:
        start_mb = numeric(work[start_col]) / 1e6
        work["_label"] = (
            work["_gene"]
            + " | "
            + work[type_col].astype(str)
            + " "
            + work[chrom_col].astype(str)
            + ":"
            + start_mb.map(lambda x: f"{x:.2f} Mb" if pd.notna(x) else "?")
        )
    else:
        work["_label"] = work["_gene"] + " | " + work[id_col].fillna(".").astype(str)

    callers_col = first_existing(work, ["CALLERS"])
    caller_count_col = first_existing(work, ["CALLER_COUNT", "SUPP"])
    af_col = first_existing(work, ["NEEDLR_AF"])
    panel_col = first_existing(work, ["PANEL_STATUS"])
    pheno_col = first_existing(work, ["PHENOTYPE_SCORE"])
    omim_col = first_existing(work, ["OMIM", "NEEDLR_OMIM"])
    gencc_col = first_existing(work, ["GENCC", "NEEDLR_GENCC"])
    ann_col = first_existing(work, ["ANNotsv_Classification", "AnnotSV_Classification"])

    if callers_col:
        caller_text = work[callers_col]
        work["Sniffles2"] = caller_present(caller_text, "Sniffles2").astype(int)
        work["cuteSV"] = caller_present(caller_text, "cuteSV").astype(int)
        work["Manta"] = caller_present(caller_text, "Manta").astype(int)
        work["Delly"] = caller_present(caller_text, "Delly").astype(int)
    else:
        work[["Sniffles2", "cuteSV", "Manta", "Delly"]] = np.nan

    caller_count = numeric(work[caller_count_col]) if caller_count_col else pd.Series(np.nan, index=work.index)
    work["Multi-caller"] = (caller_count >= 2).astype(float).where(caller_count.notna())

    if af_col:
        af = numeric(work[af_col])
        work["needLR evaluable"] = af.notna().astype(int)
        work[f"needLR AF ≤ {args.rare_af:g}"] = ((af <= args.rare_af) & af.notna()).astype(int)
        work["needLR AF = 0"] = (af.eq(0) & af.notna()).astype(int)
    else:
        af = pd.Series(np.nan, index=work.index)
        work["needLR evaluable"] = 0
        work[f"needLR AF ≤ {args.rare_af:g}"] = 0
        work["needLR AF = 0"] = 0

    if panel_col:
        ptxt = work[panel_col].fillna("").astype(str).str.upper()
        work["Panel gene"] = ptxt.str.strip().isin(["PANEL_GENE", "YES"]).astype(int)
    else:
        work["Panel gene"] = 0

    work["OMIM evidence"] = text_present(work[omim_col]).astype(int) if omim_col else 0
    work["GenCC evidence"] = text_present(work[gencc_col]).astype(int) if gencc_col else 0
    work["AnnotSV class"] = text_present(work[ann_col]).astype(int) if ann_col else 0
    pheno = numeric(work[pheno_col]).fillna(0) if pheno_col else pd.Series(0, index=work.index)
    work["Phenotype overlap"] = (pheno > 0).astype(int)

    optional_flags = []
    for source_col, display in [
        ("STRAGLR_MATCH", "Straglr match"),
        ("TLDR_MATCH", "TLDR match"),
        ("LONGPHASE_MATCH", "LongPhase match"),
    ]:
        col = first_existing(work, [source_col])
        if col:
            work[display] = yes_flag(work[col])
            optional_flags.append(display)

    whatshap_count_col = first_existing(work, ["WHATSHAP_PHASED_HET_COUNT"])
    if whatshap_count_col:
        work["Nearby WhatsHap phase"] = (numeric(work[whatshap_count_col]).fillna(0) > 0).astype(int)
        optional_flags.append("Nearby WhatsHap phase")

    methylation_context_col = first_existing(work, ["METHYLATION_CONTEXT"])
    if methylation_context_col:
        work["Methylation data available"] = (
            work[methylation_context_col].fillna("").astype(str).str.upper().eq("EVALUATED")
        ).astype(int)
        optional_flags.append("Methylation data available")

    mitocarta_encoding_col = first_existing(work, ["MITOCARTA_ENCODING"])
    if mitocarta_encoding_col:
        work["Nuclear mitochondrial gene"] = (
            work[mitocarta_encoding_col]
            .fillna("")
            .astype(str)
            .str.upper()
            .eq("NUCLEAR_MITOCHONDRIAL_GENE")
        ).astype(int)
        optional_flags.append("Nuclear mitochondrial gene")

    mito_on_col = first_existing(work, ["MITO_ON_CONTEXT"])
    if mito_on_col:
        work["Mito + ON context"] = (
            work[mito_on_col].fillna("").astype(str).str.upper().eq("YES")
        ).astype(int)
        optional_flags.append("Mito + ON context")

    score_col = first_existing(
        work,
        [
            "EVENT_GENE_RELEVANCE_SCORE",
            "INTEGRATED_DISCOVERY_SCORE",
            "integrated_discovery_score",
            "PHENOTYPE_SCORE",
        ],
    )
    # Match the documented gene-discovery ranking. Database presence, panel
    # membership and optional analyses must not silently create another score.
    work["PLOT_ORDER_BASIS"] = score_col or (pheno_col or "UNRANKED")
    work["PLOT_ORDER_SCORE"] = numeric(work[score_col]) if score_col else (
        numeric(work[pheno_col]) if pheno_col else np.nan
    )
    work["_caller_count"] = caller_count
    work["_pheno"] = pheno
    work["_af"] = af
    work = (
        work.sort_values(["PLOT_ORDER_SCORE", "_pheno", "_caller_count", id_col, gene_col],
                         ascending=[False, False, False, True, True], na_position="last")
        .drop_duplicates(subset=[id_col, gene_col], keep="first")
        .head(args.top_n)
        .copy()
    )
    # Keep the full SV_ID in the exported TSV; the figure uses a shorter label
    # so candidate rows remain readable at thesis scale.

    evidence_cols = [
        "Sniffles2", "cuteSV", "Manta", "Delly", "Multi-caller",
        "needLR evaluable", f"needLR AF ≤ {args.rare_af:g}", "needLR AF = 0",
        "Panel gene", "AnnotSV class", "OMIM evidence", "GenCC evidence", "Phenotype overlap",
    ] + optional_flags

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work[[id_col, gene_col, "PLOT_ORDER_BASIS", "PLOT_ORDER_SCORE"] + evidence_cols].to_csv(
        prefix.with_name(prefix.name + "_matrix.tsv"), sep="\t", index=False
    )

    matrix = work[evidence_cols].astype(float).to_numpy()
    n = len(work)
    fig_h = max(7.2, 0.42 * n + 2.6)
    fig, ax = plt.subplots(figsize=(17.0, fig_h))
    cmap = ListedColormap(["#F3F4F4", "#0B6E69"])
    cmap.set_bad("#AAB2BA")
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(evidence_cols)))
    wrapped_evidence = ["\n".join(textwrap.wrap(x, width=14)) for x in evidence_cols]
    ax.set_xticklabels(wrapped_evidence, rotation=35, ha="right", fontsize=8.5)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(work["_label"], fontsize=9)
    ax.set_xlabel("Evidence layer")
    ax.set_ylabel("SV | overlapping gene")
    ax.legend(handles=[Patch(facecolor="#0B6E69", label="Reported"),
                       Patch(facecolor="#F3F4F4", edgecolor="#BBBBBB", label="Not reported"),
                       Patch(facecolor="#AAB2BA", label="Unavailable/unknown")],
              loc="upper left", bbox_to_anchor=(0, 1.08), ncol=3, fontsize=8)
    ax.set_xticks(np.arange(-0.5, len(evidence_cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    add_panel_label(ax, "A")

    for i, (_, row) in enumerate(work.iterrows()):
        af_txt = "AF n/e" if pd.isna(row["_af"]) else f"AF={row['_af']:.3g}"
        ax.text(
            len(evidence_cols) - 0.15,
            i,
            f"  {af_txt}; phenotype={row['_pheno']:.2f}",
            va="center",
            ha="left",
            fontsize=8,
            clip_on=False,
        )
    ax.set_xlim(-0.5, len(evidence_cols) + 3.4)

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.008,
        "Order follows the available gene-relevance score (fallback: phenotype score). Straglr, TLDR and phasing are complementary same-dataset analyses; methylation indicates local data availability. Context presence is not independent SV confirmation or a pathogenicity classification.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
