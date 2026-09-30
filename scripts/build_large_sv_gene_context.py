#!/usr/bin/env python3
"""Build an event-aligned gene-position table for large DEL/DUP depth plots."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def normalize_chrom(value):
    text = str(value)
    return text if text.startswith("chr") else "chr" + text


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--depth-summary", required=True)
    p.add_argument("--gene-bed", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    cand = pd.read_csv(args.candidates, sep="\t", dtype=str, low_memory=False)
    depth = pd.read_csv(args.depth_summary, sep="\t", dtype=str, low_memory=False)

    columns = [
        "SV_ID", "CHROM", "SV_START", "SV_END", "SVTYPE", "GENE",
        "GENE_START", "GENE_END", "GENE_MID", "GENE_RELEVANCE_SCORE",
        "PANEL_STATUS", "MITOCARTA", "SV_GENE_EFFECT", "DEPTH_PATTERN",
    ]

    if cand.empty or depth.empty:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=columns).to_csv(args.output, sep="\t", index=False)
        return

    bed = pd.read_csv(args.gene_bed, sep="\t", comment="#", header=None, dtype=str, low_memory=False)
    if bed.shape[1] < 4:
        raise ValueError("Gene BED must contain chromosome, start, end and gene columns.")
    bed = bed.iloc[:, :4].copy()
    bed.columns = ["CHROM", "GENE_START", "GENE_END", "GENE"]
    bed["CHROM"] = bed["CHROM"].map(normalize_chrom)
    bed["GENE_START"] = pd.to_numeric(bed["GENE_START"], errors="coerce")
    bed["GENE_END"] = pd.to_numeric(bed["GENE_END"], errors="coerce")
    bed = bed.dropna(subset=["GENE_START", "GENE_END", "GENE"]).drop_duplicates(["CHROM", "GENE"])

    keep_ids = set(depth["SV_ID"].astype(str))
    work = cand[cand["SV_ID"].astype(str).isin(keep_ids)].copy()
    work = work.merge(
        bed,
        on=["CHROM", "GENE"],
        how="left",
    )
    work = work.merge(
        depth[["SV_ID", "DEPTH_PATTERN"]].drop_duplicates("SV_ID"),
        on="SV_ID",
        how="left",
    )

    work["SV_START"] = pd.to_numeric(work["START"], errors="coerce")
    work["SV_END"] = pd.to_numeric(work["END"], errors="coerce")
    work["GENE_MID"] = (work["GENE_START"] + work["GENE_END"]) / 2

    out = work.rename(columns={"SV_GENE_EFFECT": "SV_GENE_EFFECT"})[
        [
            "SV_ID", "CHROM", "SV_START", "SV_END", "SVTYPE", "GENE",
            "GENE_START", "GENE_END", "GENE_MID", "GENE_RELEVANCE_SCORE",
            "PANEL_STATUS", "MITOCARTA", "SV_GENE_EFFECT", "DEPTH_PATTERN",
        ]
    ].drop_duplicates(["SV_ID", "GENE"])

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)
    print(f"[OK] large_sv_gene_rows={len(out)} output={output}")


if __name__ == "__main__":
    main()
