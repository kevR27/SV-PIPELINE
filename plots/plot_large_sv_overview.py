#!/usr/bin/env python3
"""Summarize very large structural variants without discarding them.

One row represents one master SV. Gene counts are descriptive interval-overlap
context. For inversions and BNDs, interval overlap must not be interpreted as
direct disruption of every gene inside the event.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import SVTYPE_COLORS, first_existing, normalize_svtype, numeric, read_tsv, save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Integrated SV-gene TSV")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-size", type=int, default=10_000_000)
    return p.parse_args()


def mechanism_scope(svtype: str) -> str:
    svtype = str(svtype).upper()
    if svtype == "DEL":
        return "INTERVAL_DOSAGE_LOSS_CONTEXT"
    if svtype in {"DUP", "CNV"}:
        return "INTERVAL_COPY_GAIN_CONTEXT"
    if svtype == "INV":
        return "BREAKPOINT_PRIMARY_INTERVAL_CONTEXT_ONLY"
    if svtype in {"BND", "TRA"}:
        return "BREAKPOINT_PRIMARY"
    if svtype == "INS":
        return "INSERTION_SITE_PRIMARY"
    return "REVIEW_REQUIRED"


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    chrom_col = first_existing(df, ["CHROM", "chrom"])
    start_col = first_existing(df, ["START", "POS"])
    end_col = first_existing(df, ["END"])
    type_col = first_existing(df, ["SVTYPE", "SV_type"])
    len_col = first_existing(df, ["SVLEN", "SV_Length"])
    caller_col = first_existing(df, ["CALLER_COUNT", "SUPP"])
    panel_col = first_existing(df, ["PANEL_STATUS", "panel_gene"])
    af_col = first_existing(df, ["NEEDLR_AF"])

    if None in (id_col, chrom_col, start_col, type_col):
        raise ValueError("Input needs SV_ID, chromosome, start and SVTYPE columns.")

    work = df.copy()
    work["_start"] = numeric(work[start_col])
    work["_end"] = numeric(work[end_col]) if end_col else work["_start"]
    if len_col:
        work["_size"] = numeric(work[len_col]).abs()
    else:
        work["_size"] = (work["_end"] - work["_start"]).abs()
    fallback = (work["_end"] - work["_start"]).abs()
    work["_size"] = work["_size"].fillna(fallback)
    work["_svtype"] = normalize_svtype(work[type_col])
    work["_gene"] = work[gene_col].fillna(".").astype(str) if gene_col else "."
    work["_caller_count"] = numeric(work[caller_col]).fillna(0) if caller_col else 0
    work["_af"] = numeric(work[af_col]) if af_col else np.nan

    if panel_col:
        ptxt = work[panel_col].fillna("").astype(str).str.upper().str.strip()
        work["_panel"] = ptxt.isin(["PANEL_GENE", "YES", "TRUE", "1"])
    else:
        work["_panel"] = False

    work = work[work["_size"] >= args.min_size].copy()

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    if work.empty:
        pd.DataFrame(columns=["SV_ID", "SVTYPE", "size_bp"]).to_csv(
            prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False
        )
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.axis("off")
        ax.text(0.5, 0.55, f"No master SVs ≥ {args.min_size/1e6:.0f} Mb", ha="center", fontsize=14)
        outputs = save_figure(fig, prefix)
        plt.close(fig)
        print("[OK] no very-large SVs")
        print("[OK]", *outputs, sep="\n")
        return

    rows = []
    for sv_id, group in work.groupby(id_col, sort=False):
        first = group.iloc[0]
        genes = sorted({
            str(x) for x in group["_gene"]
            if str(x) not in {"", ".", "NA", "N/A", "nan", "None"}
        })
        panel_genes = sorted({
            str(row["_gene"]) for _, row in group.iterrows()
            if bool(row["_panel"]) and str(row["_gene"]) not in {"", ".", "nan", "None"}
        })
        af_values = group["_af"].dropna()
        rows.append({
            "SV_ID": sv_id,
            "CHROM": first[chrom_col],
            "START": int(first["_start"]) if pd.notna(first["_start"]) else ".",
            "END": int(first["_end"]) if pd.notna(first["_end"]) else ".",
            "SVTYPE": first["_svtype"],
            "size_bp": float(first["_size"]),
            "size_Mb": float(first["_size"]) / 1e6,
            "caller_count": int(first["_caller_count"]),
            "overlapping_gene_count": len(genes),
            "panel_gene_count": len(panel_genes),
            "panel_genes": ";".join(panel_genes) if panel_genes else ".",
            "needLR_AF_min": float(af_values.min()) if not af_values.empty else np.nan,
            "needLR_evaluable": "YES" if not af_values.empty else "NO",
            "interpretation_scope": mechanism_scope(first["_svtype"]),
        })

    summary = pd.DataFrame(rows).sort_values("size_Mb", ascending=False).reset_index(drop=True)
    summary.to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)

    plot = summary.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(14.5, max(7.0, 0.48 * len(plot) + 2.8)))
    y = np.arange(len(plot))
    colors = [SVTYPE_COLORS.get(x, "#999999") for x in plot["SVTYPE"]]
    bars = ax.barh(y, plot["size_Mb"], color=colors)

    labels = [
        f"{row.SVTYPE} {row.CHROM}:{row.START}-{row.END}"
        for row in plot.itertuples()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("SV span (Mb)")
    ax.set_ylabel("Very-large master SV")
    ax.set_title(f"Master structural variants ≥ {args.min_size/1e6:.0f} Mb")
    style_axis(ax, "x")

    xmax = max(plot["size_Mb"].max(), 1)
    for bar, (_, row) in zip(bars, plot.iterrows()):
        panel_txt = f"; panel={int(row['panel_gene_count'])}" if row["panel_gene_count"] else ""
        ax.text(
            bar.get_width() + 0.012 * xmax,
            bar.get_y() + bar.get_height() / 2,
            f"{row['caller_count']} caller(s); genes={int(row['overlapping_gene_count'])}{panel_txt}",
            va="center",
            fontsize=8,
        )

    ax.set_xlim(0, xmax * 1.38)
    fig.text(
        0.5, 0.012,
        "Gene counts describe interval overlap. DEL/DUP intervals can alter dosage if real; for INV/BND, direct disruption is primarily a breakpoint question and internal gene overlap alone is not equivalent to gene disruption.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print(f"[OK] very_large_master_SVs={len(summary)}")
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
