#!/usr/bin/env python3
"""Summarize binned read depth around large deletions and duplications."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def choose_plot_bin(size_bp: int) -> int:
    """Choose a display bin from SV size.

    Coverage calculations still use the original 10 kb mosdepth windows.
    These larger bins are used only for a cleaner plot.
    """
    if size_bp >= 50_000_000:
        return 500_000
    if size_bp >= 10_000_000:
        return 250_000
    if size_bp >= 1_000_000:
        return 100_000
    return 50_000


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--coverage", required=True, help="mosdepth regions BED.gz")
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--bins-output", required=True)
    parser.add_argument("--min-size", type=int, default=100_000)
    parser.add_argument("--max-flank", type=int, default=5_000_000)
    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates, sep="\t", dtype=str, low_memory=False)
    coverage = pd.read_csv(
        args.coverage,
        sep="\t",
        header=None,
        names=["CHROM", "START", "END", "DEPTH"],
        compression="infer",
    )
    coverage["START"] = pd.to_numeric(coverage["START"], errors="coerce")
    coverage["END"] = pd.to_numeric(coverage["END"], errors="coerce")
    coverage["DEPTH"] = pd.to_numeric(coverage["DEPTH"], errors="coerce")

    candidates["SVLEN_NUM"] = pd.to_numeric(candidates["SVLEN"], errors="coerce").abs()
    candidates["START_NUM"] = pd.to_numeric(candidates["START"], errors="coerce")
    candidates["END_NUM"] = pd.to_numeric(candidates["END"], errors="coerce")

    large = candidates[
        candidates["SVTYPE"].isin(["DEL", "DUP"])
        & candidates["SVLEN_NUM"].ge(args.min_size)
    ].copy()

    summary_rows = []
    bin_rows = []

    for sv_id, event_rows in large.groupby("SV_ID", sort=False):
        event_rows = event_rows.copy()
        score = pd.to_numeric(
            event_rows.get("GENE_RELEVANCE_SCORE", pd.Series(index=event_rows.index, dtype=float)),
            errors="coerce",
        ).fillna(0)
        panel = (
            event_rows.get("PANEL_STATUS", pd.Series("", index=event_rows.index))
            .fillna("")
            .astype(str)
            .eq("PANEL_GENE")
        )
        mito = (
            event_rows.get("MITOCARTA", pd.Series("", index=event_rows.index))
            .fillna("")
            .astype(str)
            .isin(["NUCLEAR_MITOCHONDRIAL_GENE", "MTDNA_ENCODED_GENE"])
        )
        event_rows["_display_priority"] = panel.astype(int) * 100 + mito.astype(int) * 10 + score
        event_rows = event_rows.sort_values(
            ["_display_priority", "GENE"],
            ascending=[False, True],
        )
        sv = event_rows.iloc[0]
        genes = sorted({
            str(x)
            for x in event_rows["GENE"].dropna().astype(str)
            if str(x) not in {"", ".", "nan", "None"}
        })
        panel_genes = sorted(set(event_rows.loc[panel, "GENE"].dropna().astype(str)))
        mito_genes = sorted(set(event_rows.loc[mito, "GENE"].dropna().astype(str)))
        top_genes = event_rows["GENE"].dropna().astype(str).drop_duplicates().head(12).tolist()

        chrom = sv["CHROM"]
        start = int(sv["START_NUM"])
        end = int(sv["END_NUM"])
        if end < start:
            start, end = end, start

        size = end - start
        flank = min(args.max_flank, max(500_000, int(size * 0.25)))
        region_start = max(0, start - flank)
        region_end = end + flank
        plot_bin = choose_plot_bin(size)

        local = coverage[
            (coverage["CHROM"] == chrom)
            & (coverage["END"] > region_start)
            & (coverage["START"] < region_end)
        ].copy()
        if local.empty:
            continue

        local["REGION"] = np.where(
            local["END"] <= start,
            "LEFT_FLANK",
            np.where(local["START"] >= end, "RIGHT_FLANK", "SV"),
        )
        flank_depth = local.loc[local["REGION"] != "SV", "DEPTH"].median()
        inside_depth = local.loc[local["REGION"] == "SV", "DEPTH"].median()
        depth_ratio = inside_depth / flank_depth if pd.notna(flank_depth) and flank_depth > 0 else np.nan

        if sv["SVTYPE"] == "DEL":
            pattern = "CONSISTENT_WITH_LOSS" if pd.notna(depth_ratio) and depth_ratio < 0.75 else "NOT_CLEAR"
        else:
            pattern = "CONSISTENT_WITH_GAIN" if pd.notna(depth_ratio) and depth_ratio > 1.25 else "NOT_CLEAR"

        summary_rows.append({
            "SV_ID": sv_id,
            "GENE": sv.get("GENE", "."),
            "LEAD_GENE": sv.get("GENE", "."),
            "GENE_COUNT": len(genes),
            "GENES_IN_SV": ";".join(genes) if genes else ".",
            "PANEL_GENES_IN_SV": ";".join(panel_genes) if panel_genes else ".",
            "MITOCARTA_GENES_IN_SV": ";".join(mito_genes) if mito_genes else ".",
            "TOP_RELEVANT_GENES": ";".join(top_genes) if top_genes else ".",
            "CHROM": chrom,
            "START": start,
            "END": end,
            "SVTYPE": sv["SVTYPE"],
            "SV_SIZE_BP": size,
            "MEDIAN_DEPTH_INSIDE": round(float(inside_depth), 3) if pd.notna(inside_depth) else ".",
            "MEDIAN_DEPTH_FLANKS": round(float(flank_depth), 3) if pd.notna(flank_depth) else ".",
            "DEPTH_RATIO": round(float(depth_ratio), 4) if pd.notna(depth_ratio) else ".",
            "DEPTH_PATTERN": pattern,
            "PLOT_BIN_BP": plot_bin,
        })

        local["PLOT_START"] = (local["START"] // plot_bin) * plot_bin
        grouped = (
            local.groupby(["CHROM", "PLOT_START"], as_index=False)
            .agg(DEPTH=("DEPTH", "median"))
        )
        grouped["PLOT_END"] = grouped["PLOT_START"] + plot_bin
        grouped["SV_ID"] = sv_id
        grouped["GENE"] = sv.get("GENE", ".")
        grouped["LEAD_GENE"] = sv.get("GENE", ".")
        grouped["SVTYPE"] = sv["SVTYPE"]
        grouped["SV_START"] = start
        grouped["SV_END"] = end
        grouped["NORMALIZED_DEPTH"] = (
            grouped["DEPTH"] / flank_depth
            if pd.notna(flank_depth) and flank_depth > 0
            else np.nan
        )
        bin_rows.extend(grouped.to_dict("records"))

    summary = pd.DataFrame(
        summary_rows,
        columns=[
            "SV_ID", "GENE", "LEAD_GENE", "GENE_COUNT", "GENES_IN_SV",
            "PANEL_GENES_IN_SV", "MITOCARTA_GENES_IN_SV", "TOP_RELEVANT_GENES",
            "CHROM", "START", "END", "SVTYPE",
            "SV_SIZE_BP", "MEDIAN_DEPTH_INSIDE", "MEDIAN_DEPTH_FLANKS",
            "DEPTH_RATIO", "DEPTH_PATTERN", "PLOT_BIN_BP",
        ],
    )
    bins = pd.DataFrame(
        bin_rows,
        columns=[
            "CHROM", "PLOT_START", "DEPTH", "PLOT_END", "SV_ID", "GENE", "LEAD_GENE",
            "SVTYPE", "SV_START", "SV_END", "NORMALIZED_DEPTH",
        ],
    )

    summary_path = Path(args.summary_output)
    bins_path = Path(args.bins_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    bins_path.parent.mkdir(parents=True, exist_ok=True)

    summary.to_csv(summary_path, sep="\t", index=False)
    bins.to_csv(bins_path, sep="\t", index=False)
    print(f"[OK] large_sv_depth={len(summary)} output={summary_path}")
    print(f"[OK] plotted_bins={len(bins)} output={bins_path}")


if __name__ == "__main__":
    main()
