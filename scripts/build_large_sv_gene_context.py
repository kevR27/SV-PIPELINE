#!/usr/bin/env python3
"""Build a complete gene-position table for large DEL/DUP depth plots.

All genes are obtained directly from the configured gene BED for each large SV
interval. Candidate-table annotations are then added when that SV-gene pair is
available. This avoids mistaking the candidate subset for the complete set of
genes physically overlapped by a large event.
"""
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
        "GENE_START", "GENE_END", "GENE_MID", "GENE_RELATIONSHIP_TO_EVENT",
        "GENE_RELEVANCE_SCORE", "PANEL_STATUS", "MITOCARTA",
        "SV_GENE_EFFECT", "DEPTH_PATTERN",
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    if depth.empty:
        pd.DataFrame(columns=columns).to_csv(output, sep="\t", index=False)
        print(f"[OK] large_sv_gene_rows=0 output={output}")
        return

    bed = pd.read_csv(
        args.gene_bed,
        sep="\t",
        comment="#",
        header=None,
        dtype=str,
        low_memory=False,
    )
    if bed.shape[1] < 4:
        raise ValueError(
            "Gene BED must contain chromosome, start, end and gene columns."
        )

    bed = bed.iloc[:, :4].copy()
    bed.columns = ["CHROM", "GENE_START", "GENE_END", "GENE"]
    bed["CHROM"] = bed["CHROM"].map(normalize_chrom)
    bed["GENE_START"] = pd.to_numeric(bed["GENE_START"], errors="coerce")
    bed["GENE_END"] = pd.to_numeric(bed["GENE_END"], errors="coerce")
    bed = (
        bed.dropna(subset=["GENE_START", "GENE_END", "GENE"])
        .sort_values(["CHROM", "GENE_START", "GENE_END", "GENE"])
        .drop_duplicates(["CHROM", "GENE"])
    )

    candidate_cols = [
        "SV_ID", "GENE", "GENE_RELEVANCE_SCORE", "PANEL_STATUS",
        "MITOCARTA", "SV_GENE_EFFECT",
    ]
    candidate_cols = [col for col in candidate_cols if col in cand.columns]
    candidate_meta = (
        cand[candidate_cols].drop_duplicates(["SV_ID", "GENE"])
        if {"SV_ID", "GENE"}.issubset(candidate_cols)
        else pd.DataFrame(columns=["SV_ID", "GENE"])
    )

    rows = []
    for _, event in depth.iterrows():
        sv_id = str(event["SV_ID"])
        chrom = normalize_chrom(event["CHROM"])
        start = int(float(event["START"]))
        end = int(float(event["END"]))
        if end < start:
            start, end = end, start

        genes = bed[
            bed["CHROM"].eq(chrom)
            & bed["GENE_END"].gt(start)
            & bed["GENE_START"].lt(end)
        ].copy()

        if genes.empty:
            continue

        genes["SV_ID"] = sv_id
        genes["SV_START"] = start
        genes["SV_END"] = end
        genes["SVTYPE"] = str(event["SVTYPE"])
        genes["GENE_MID"] = (genes["GENE_START"] + genes["GENE_END"]) / 2
        genes["DEPTH_PATTERN"] = str(event.get("DEPTH_PATTERN", "."))

        genes["GENE_RELATIONSHIP_TO_EVENT"] = "INSIDE_EVENT"
        genes.loc[
            (genes["GENE_START"] < start) | (genes["GENE_END"] > end),
            "GENE_RELATIONSHIP_TO_EVENT",
        ] = "BOUNDARY_OVERLAP"

        genes = genes.merge(
            candidate_meta,
            on=["SV_ID", "GENE"],
            how="left",
        )

        for col in [
            "GENE_RELEVANCE_SCORE",
            "PANEL_STATUS",
            "MITOCARTA",
            "SV_GENE_EFFECT",
        ]:
            if col not in genes:
                genes[col] = "."
            genes[col] = genes[col].fillna(".")

        rows.append(genes[columns])

    out = (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(columns=columns)
    )
    out.to_csv(output, sep="\t", index=False)

    print(
        f"[OK] large_sv_gene_rows={len(out)} "
        f"events={out['SV_ID'].nunique() if not out.empty else 0} "
        f"unique_genes={out['GENE'].nunique() if not out.empty else 0} "
        f"output={output}"
    )


if __name__ == "__main__":
    main()
