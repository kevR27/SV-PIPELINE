#!/usr/bin/env python3
"""Compare completed LRS samples using the same master SV-gene evidence schema."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import SVTYPE_COLORS, first_existing, normalize_svtype, numeric, read_tsv, save_figure, set_thesis_style, style_axis, unique_master_svs


def parse_args():
    p = argparse.ArgumentParser(description="Compare master SV landscapes across completed samples.")
    p.add_argument("--sample-input", action="append", required=True, help="SAMPLE=/path/to/integrated.tsv; repeat per sample")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--rare-af", type=float, default=0.01)
    p.add_argument("--title", default="Cross-sample structural-variant comparison")
    return p.parse_args()


def parse_mapping(items):
    out = []
    for item in items:
        if "=" not in item:
            raise ValueError("--sample-input must be SAMPLE=/path/to/file")
        sample, path = item.split("=", 1)
        out.append((sample, path))
    return out


def main():
    args = parse_args()
    set_thesis_style()

    summaries = []
    svtype_rows = []
    for sample, path in parse_mapping(args.sample_input):
        raw = read_tsv(path)
        master = unique_master_svs(raw)
        id_col = first_existing(raw, ["SV_ID", "ID"])
        gene_col = first_existing(raw, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
        type_col = first_existing(master, ["SVTYPE", "SV_type"])
        caller_col = first_existing(master, ["CALLER_COUNT", "SUPP"])
        af_col = first_existing(master, ["NEEDLR_AF"])
        panel_col = first_existing(raw, ["PANEL_STATUS", "panel_gene"])
        pheno_col = first_existing(raw, ["PHENOTYPE_SCORE", "phenotype_score"])
        disease_col = first_existing(raw, ["GENE_DISEASE_EVIDENCE_SCORE", "gene_disease_evidence_score"])

        svtypes = normalize_svtype(master[type_col]) if type_col else pd.Series("OTHER", index=master.index)
        for svtype, count in svtypes.value_counts().items():
            svtype_rows.append({"sample": sample, "SVTYPE": svtype, "count": int(count)})

        callers = numeric(master[caller_col]) if caller_col else pd.Series(np.nan, index=master.index)
        af = numeric(master[af_col]) if af_col else pd.Series(np.nan, index=master.index)

        pairs = raw.copy()
        if gene_col:
            pairs["_gene"] = pairs[gene_col].fillna(".").astype(str)
        else:
            pairs["_gene"] = "."
        pairs = pairs.drop_duplicates([id_col, gene_col] if gene_col else [id_col])

        if panel_col:
            ptxt = pairs[panel_col].fillna("").astype(str).str.upper().str.strip()
            panel = ptxt.isin(["PANEL_GENE", "YES", "TRUE", "1"])
        else:
            panel = pd.Series(False, index=pairs.index)
        pheno = numeric(pairs[pheno_col]).fillna(0) if pheno_col else pd.Series(0, index=pairs.index)
        disease = numeric(pairs[disease_col]).fillna(0) if disease_col else pd.Series(0, index=pairs.index)
        gene_valid = ~pairs["_gene"].isin(["", ".", "NA", "N/A", "nan", "None"])

        summaries.append(
            {
                "sample": sample,
                "master_SVs": len(master),
                "multicaller_SVs": int((callers >= 2).sum()) if caller_col else 0,
                "needLR_evaluable_SVs": int(af.notna().sum()),
                "needLR_rare_SVs": int(((af <= args.rare_af) & af.notna()).sum()),
                "gene_overlapping_SVs": int(pairs.loc[gene_valid, id_col].nunique()),
                "panel_gene_SVs": int(pairs.loc[panel & gene_valid, id_col].nunique()),
                "nonpanel_HPO_or_disease_SVs": int(pairs.loc[(~panel) & gene_valid & ((pheno > 0) | (disease > 0)), id_col].nunique()),
                "unique_overlapping_genes": int(pairs.loc[gene_valid, "_gene"].nunique()),
            }
        )

    summary = pd.DataFrame(summaries)
    svtypes = pd.DataFrame(svtype_rows)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)
    svtypes.to_csv(prefix.with_name(prefix.name + "_svtypes.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14.8, 7.2))
    ax1, ax2 = axes

    samples = summary["sample"].tolist()
    pivot = svtypes.pivot(index="sample", columns="SVTYPE", values="count").fillna(0).reindex(samples)
    bottom = np.zeros(len(samples))
    order = [x for x in ["DEL", "INS", "DUP", "INV", "BND", "CNV", "OTHER"] if x in pivot.columns]
    for svtype in order:
        vals = pivot[svtype].to_numpy()
        ax1.bar(samples, vals, bottom=bottom, label=svtype, color=SVTYPE_COLORS.get(svtype, "#999999"))
        bottom += vals
    ax1.set_ylabel("Unique master SVs")
    ax1.set_title("SV-type composition")
    ax1.legend(frameon=False, ncol=3, fontsize=8)
    style_axis(ax1, "y")

    metrics = [
        ("multicaller_SVs", "≥2 callers"),
        ("needLR_rare_SVs", f"needLR AF ≤ {args.rare_af:g}"),
        ("gene_overlapping_SVs", "Gene-overlapping"),
        ("panel_gene_SVs", "Panel-gene"),
        ("nonpanel_HPO_or_disease_SVs", "Non-panel HPO/disease"),
    ]
    x = np.arange(len(metrics))
    width = 0.8 / max(len(samples), 1)
    for i, sample in enumerate(samples):
        row = summary[summary["sample"].eq(sample)].iloc[0]
        ax2.bar(x + (i - (len(samples)-1)/2) * width, [row[m] for m, _ in metrics], width=width, label=sample)
    ax2.set_xticks(x)
    ax2.set_xticklabels([label for _, label in metrics], rotation=28, ha="right")
    ax2.set_ylabel("Unique SVs")
    ax2.set_title("Evidence-defined SV subsets")
    ax2.legend(frameon=False)
    style_axis(ax2, "y")

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.012,
        "All counts are deduplicated by master SV_ID. Evidence-defined subsets are descriptive and are not sequential filtering stages.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.055, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
