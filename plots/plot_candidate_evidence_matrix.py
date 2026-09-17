#!/usr/bin/env python3
"""Create an integrated technical/biological evidence matrix for top SV-gene candidates."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import add_panel_label, first_existing, numeric, read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description="Plot integrated evidence for prioritized SV-gene pairs.")
    p.add_argument("--input", required=True, help="*_integrated_SV_gene_analysis.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--top-n", type=int, default=30)
    p.add_argument("--rare-af", type=float, default=0.01)
    p.add_argument("--title", default="Integrated evidence for prioritized SV-gene candidates")
    return p.parse_args()


def text_present(series: pd.Series) -> pd.Series:
    s = series.fillna("").astype(str).str.strip()
    return ~s.isin(["", ".", "NA", "N/A", "None", "nan", "NaN"])


def caller_present(series: pd.Series, caller: str) -> pd.Series:
    pattern = re.escape(caller)
    return series.fillna("").astype(str).str.contains(pattern, case=False, regex=True)


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    if df.empty:
        raise ValueError("Integrated SV-gene table is empty.")

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "Gene", "GENE", "ANNotsv_Gene"])
    if id_col is None or gene_col is None:
        raise ValueError("Input needs an SV ID and gene column.")

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    work["_label"] = work[id_col].fillna(".").astype(str) + " | " + work["_gene"]

    callers_col = first_existing(work, ["CALLERS"])
    caller_count_col = first_existing(work, ["CALLER_COUNT", "SUPP"])
    af_col = first_existing(work, ["NEEDLR_AF"])
    panel_col = first_existing(work, ["PANEL_STATUS"])
    pheno_col = first_existing(work, ["PHENOTYPE_SCORE"])
    omim_col = first_existing(work, ["OMIM", "NEEDLR_OMIM"])
    gencc_col = first_existing(work, ["GENCC", "NEEDLR_GENCC"])
    ann_col = first_existing(work, ["ANNotsv_Classification", "AnnotSV_Classification"])
    candidate_col = first_existing(work, ["CANDIDATE_CLASS"])

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
        work["needLR evaluable"] = 0
        work[f"needLR AF ≤ {args.rare_af:g}"] = 0
        work["needLR AF = 0"] = 0
        af = pd.Series(np.nan, index=work.index)

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

    # Ranking is deliberately transparent and is not a pathogenicity score.
    work["_plot_priority"] = (
        work["Multi-caller"] * 2
        + work[f"needLR AF ≤ {args.rare_af:g}"] * 2
        + work["Panel gene"] * 2
        + work["OMIM evidence"]
        + work["GenCC evidence"]
        + work["Phenotype overlap"] * 2
        + pheno.rank(pct=True).fillna(0)
    )
    work["_pheno"] = pheno
    work["_af"] = af
    work = work.sort_values(["_plot_priority", "_pheno"], ascending=False).head(args.top_n).copy()
    work = work.drop_duplicates(subset=["_label"], keep="first")

    evidence_cols = [
        "Sniffles2",
        "cuteSV",
        "Delly",
        "Multi-caller",
        "needLR evaluable",
        f"needLR AF ≤ {args.rare_af:g}",
        "needLR AF = 0",
        "Panel gene",
        "AnnotSV class",
        "OMIM evidence",
        "GenCC evidence",
        "Phenotype overlap",
    ]

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out_cols = [id_col, gene_col] + [c for c in [caller_count_col, af_col, panel_col, pheno_col, ann_col, candidate_col] if c] + evidence_cols
    work[out_cols].to_csv(prefix.with_name(prefix.name + "_matrix.tsv"), sep="\t", index=False)

    matrix = work[evidence_cols].astype(float).to_numpy()
    n = len(work)
    fig_h = max(6.0, 0.30 * n + 2.0)
    fig, ax = plt.subplots(figsize=(13.4, fig_h))
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(evidence_cols)))
    ax.set_xticklabels(evidence_cols, rotation=42, ha="right")
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(work["_label"], fontsize=7.5)
    ax.set_xlabel("Evidence layer")
    ax.set_ylabel("Master SV | overlapping gene")
    ax.set_xticks(np.arange(-0.5, len(evidence_cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    add_panel_label(ax, "A")

    # Add compact quantitative annotations to the right without converting them to binary evidence.
    for i, (_, row) in enumerate(work.iterrows()):
        af_txt = "AF n/e" if pd.isna(row["_af"]) else f"AF={row['_af']:.3g}"
        pheno_txt = f"phenotype={row['_pheno']:.2f}"
        ax.text(len(evidence_cols) - 0.15, i, f"  {af_txt}; {pheno_txt}", va="center", ha="left", fontsize=6.8, clip_on=False)
    ax.set_xlim(-0.5, len(evidence_cols) + 3.0)

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.text(0.5, 0.01, "Filled cells indicate evidence presence only. The display ordering is a review aid and must not be interpreted as a pathogenicity classification.", ha="center", fontsize=8.2)
    fig.tight_layout(rect=[0, 0.025, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
