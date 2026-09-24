#!/usr/bin/env python3
"""Create an integrated technical/biological evidence matrix for top SV-gene candidates."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd

from plot_utils import add_panel_label, first_existing, numeric, read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description="Plot integrated evidence for prioritized SV-gene pairs.")
    p.add_argument("--input", required=True, help="Integrated or extended SV-gene analysis TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--rare-af", type=float, default=0.01)
    p.add_argument("--title", default="Integrated evidence for prioritized SV-gene candidates")
    return p.parse_args()


def text_present(series: pd.Series) -> pd.Series:
    s = series.fillna("").astype(str).str.strip()
    return ~s.isin(["", ".", "NA", "N/A", "None", "nan", "NaN"])


def caller_present(series: pd.Series, caller: str) -> pd.Series:
    return series.fillna("").astype(str).str.contains(re.escape(caller), case=False, regex=True)


def yes_flag(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.upper().isin(["YES", "TRUE", "1"])


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
        work["_label"] = (
            work[chrom_col].astype(str)
            + ":"
            + work[start_col].astype(str)
            + " "
            + work[type_col].astype(str)
            + " | "
            + work["_gene"]
        )
    else:
        work["_label"] = work[id_col].fillna(".").astype(str) + " | " + work["_gene"]

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
        work["Delly"] = caller_present(caller_text, "Delly").astype(int)
    else:
        work[["Sniffles2", "cuteSV", "Delly"]] = 0

    caller_count = numeric(work[caller_count_col]) if caller_count_col else pd.Series(np.nan, index=work.index)
    work["Multi-caller"] = (caller_count >= 2).astype(int)

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
        work["Panel gene"] = ptxt.str.contains("PANEL_GENE|^YES$", regex=True).astype(int)
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
            work[display] = yes_flag(work[col]).astype(int)
            optional_flags.append(display)

    whatshap_count_col = first_existing(work, ["WHATSHAP_PHASED_HET_COUNT"])
    if whatshap_count_col:
        work["Nearby WhatsHap phase"] = (numeric(work[whatshap_count_col]).fillna(0) > 0).astype(int)
        optional_flags.append("Nearby WhatsHap phase")

    methylation_context_col = first_existing(work, ["METHYLATION_CONTEXT"])
    if methylation_context_col:
        work["Methylation context"] = (
            work[methylation_context_col].fillna("").astype(str).str.upper().eq("EVALUATED")
        ).astype(int)
        optional_flags.append("Methylation context")

    work["_plot_priority"] = (
        work["Multi-caller"] * 2
        + work[f"needLR AF ≤ {args.rare_af:g}"] * 2
        + work["Panel gene"] * 2
        + work["OMIM evidence"]
        + work["GenCC evidence"]
        + work["Phenotype overlap"] * 2
        + sum(work[x] for x in optional_flags)
        + pheno.rank(pct=True).fillna(0)
    )
    work["_pheno"] = pheno
    work["_af"] = af
    work = (
        work.sort_values(["_plot_priority", "_pheno"], ascending=False)
        .drop_duplicates(subset=["_label"], keep="first")
        .head(args.top_n)
        .copy()
    )

    evidence_cols = [
        "Sniffles2", "cuteSV", "Delly", "Multi-caller",
        "needLR evaluable", f"needLR AF ≤ {args.rare_af:g}", "needLR AF = 0",
        "Panel gene", "AnnotSV class", "OMIM evidence", "GenCC evidence", "Phenotype overlap",
    ] + optional_flags

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    work[[id_col, gene_col] + evidence_cols].to_csv(
        prefix.with_name(prefix.name + "_matrix.tsv"), sep="\t", index=False
    )

    matrix = work[evidence_cols].astype(float).to_numpy()
    n = len(work)
    fig_h = max(7.2, 0.42 * n + 2.6)
    fig, ax = plt.subplots(figsize=(15.5, fig_h))
    cmap = ListedColormap(["#F3F4F4", "#0B6E69"])
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(evidence_cols)))
    ax.set_xticklabels(evidence_cols, rotation=42, ha="right", fontsize=9)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(work["_label"], fontsize=9)
    ax.set_xlabel("Evidence layer")
    ax.set_ylabel("SV | overlapping gene")
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
        "Dark teal indicates evidence/context presence. WhatsHap and methylation are local context layers, not SV confirmations. Ordering is a review aid, not a pathogenicity classification.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
