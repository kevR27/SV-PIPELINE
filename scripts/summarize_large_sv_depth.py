#!/usr/bin/env python3
"""Summarize binned read depth around large deletions and duplications."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def choose_plot_bin(size_bp: int) -> int:
    if size_bp >= 10_000_000:
        return 100_000
    if size_bp >= 1_000_000:
        return 50_000
    return 10_000


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
    ].drop_duplicates("SV_ID")

    summary_rows = []
    bin_rows = []

    for _, sv in large.iterrows():
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
            "SV_ID": sv["SV_ID"],
            "GENE": sv.get("GENE", "."),
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
        grouped["SV_ID"] = sv["SV_ID"]
        grouped["GENE"] = sv.get("GENE", ".")
        grouped["SVTYPE"] = sv["SVTYPE"]
        grouped["SV_START"] = start
        grouped["SV_END"] = end
        grouped["NORMALIZED_DEPTH"] = (
            grouped["DEPTH"] / flank_depth
            if pd.notna(flank_depth) and flank_depth > 0
            else np.nan
        )
        bin_rows.extend(grouped.to_dict("records"))

    summary = pd.DataFrame(summary_rows)
    bins = pd.DataFrame(bin_rows)

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
