#!/usr/bin/env python3
"""Attach orthogonal Straglr/TLDR evidence to the integrated Jasmine SV-gene table.

The master SV universe is not changed. This script only adds evidence flags and
matched locus identifiers. Straglr and TLDR are not treated as extra Jasmine callers.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from plot_utils import first_existing, normalize_svtype, read_tsv


def parse_args():
    p = argparse.ArgumentParser(description="Intersect master SVs with orthogonal repeat/MEI evidence.")
    p.add_argument("--integrated", required=True, help="*_integrated_SV_gene_analysis.tsv")
    p.add_argument("--straglr", default=None, help="Optional *_straglr.annotated.tsv")
    p.add_argument("--tldr", default=None, help="Optional TLDR *.table.txt")
    p.add_argument("--straglr-breakpoint-tol", type=int, default=500)
    p.add_argument("--tldr-breakpoint-tol", type=int, default=500)
    p.add_argument("--output", required=True)
    return p.parse_args()


def to_int(value):
    try:
        return int(float(value))
    except Exception:
        return None


def interval_overlap(a_start, a_end, b_start, b_end) -> int:
    if None in (a_start, a_end, b_start, b_end):
        return 0
    alo, ahi = sorted((a_start, a_end))
    blo, bhi = sorted((b_start, b_end))
    return max(0, min(ahi, bhi) - max(alo, blo) + 1)


def normalize_chrom(x) -> str:
    x = str(x)
    return x if x.startswith("chr") else "chr" + x


def load_straglr(path: str) -> list[dict]:
    df = read_tsv(path)
    chrom = first_existing(df, ["chrom", "CHROM", "Chr"])
    start = first_existing(df, ["start", "START", "POS"])
    end = first_existing(df, ["end", "END"])
    locus = first_existing(df, ["locus"])
    gene = first_existing(df, ["overlapping_genes"])
    copy_number = first_existing(df, ["copy_number"])
    support = first_existing(df, ["supporting_reads"])

    if chrom is None or start is None or end is None:
        raise ValueError("Straglr table needs chromosome/start/end columns.")

    rows = []
    for i, row in df.iterrows():
        rows.append(
            {
                "chrom": normalize_chrom(row[chrom]),
                "start": to_int(row[start]),
                "end": to_int(row[end]),
                "label": str(row[locus]) if locus else f"STRAGLR_{i+1}",
                "genes": str(row[gene]) if gene else ".",
                "copy_number": str(row[copy_number]) if copy_number else ".",
                "supporting_reads": str(row[support]) if support else ".",
            }
        )
    return rows


def load_tldr(path: str) -> list[dict]:
    df = read_tsv(path)
    chrom = first_existing(df, ["Chrom", "chrom", "CHROM", "chr", "Chr"])
    pos = first_existing(df, ["Start", "start", "START", "pos", "POS", "position", "POSITION"])
    end = first_existing(df, ["End", "end", "END"])
    family = first_existing(df, ["Family", "family", "FAMILY", "repeat_family", "element", "ELEMENT"])
    subfamily = first_existing(df, ["Subfamily", "subfamily", "SUBFAMILY", "repeat_name", "REPEAT_NAME"])
    used_reads = first_existing(df, ["UsedReads", "used_reads", "support", "SUPPORT"])
    span_reads = first_existing(df, ["SpanReads", "span_reads", "spanning_reads"])
    uuid = first_existing(df, ["UUID", "uuid"])

    if chrom is None or pos is None:
        raise ValueError("TLDR table needs chromosome and insertion-position columns.")

    rows = []
    for i, row in df.iterrows():
        label_parts = [
            str(row[x])
            for x in [family, subfamily]
            if x and str(row[x]) not in {"", ".", "nan"}
        ]
        rows.append(
            {
                "chrom": normalize_chrom(row[chrom]),
                "start": to_int(row[pos]),
                "end": to_int(row[end]) if end else to_int(row[pos]),
                "label": ":".join(label_parts) if label_parts else f"TLDR_{i+1}",
                "uuid": str(row[uuid]) if uuid else ".",
                "used_reads": str(row[used_reads]) if used_reads else ".",
                "span_reads": str(row[span_reads]) if span_reads else ".",
            }
        )
    return rows


def main():
    args = parse_args()
    df = read_tsv(args.integrated)

    chrom_col = first_existing(df, ["CHROM", "Chr", "chrom"])
    start_col = first_existing(df, ["START", "POS", "SV_start"])
    end_col = first_existing(df, ["END", "SV_end"])
    type_col = first_existing(df, ["SVTYPE", "SV_type", "Type"])
    if None in (chrom_col, start_col, type_col):
        raise ValueError("Integrated table needs CHROM, START/POS and SVTYPE columns.")

    straglr = load_straglr(args.straglr) if args.straglr else []
    tldr = load_tldr(args.tldr) if args.tldr else []

    out = df.copy()
    out["STRAGLR_MATCH"] = "NO"
    out["STRAGLR_LOCI"] = "."
    out["STRAGLR_GENES"] = "."
    out["STRAGLR_COPY_NUMBER"] = "."
    out["STRAGLR_SUPPORTING_READS"] = "."
    out["TLDR_MATCH"] = "NO"
    out["TLDR_INSERTIONS"] = "."
    out["TLDR_UUID"] = "."
    out["TLDR_USED_READS"] = "."
    out["TLDR_SPAN_READS"] = "."

    svtypes = normalize_svtype(out[type_col])

    for idx, row in out.iterrows():
        chrom = normalize_chrom(row[chrom_col])
        start = to_int(row[start_col])
        end = to_int(row[end_col]) if end_col else start
        svtype = svtypes.loc[idx]

        if start is None:
            continue
        if end is None:
            end = start

        smatches = []
        sgenes = []
        scn = []
        ssupport = []
        for s in straglr:
            if s["chrom"] != chrom or s["start"] is None or s["end"] is None:
                continue

            overlap = interval_overlap(start, end, s["start"], s["end"])
            near_bp = (
                min(
                    abs(start - s["start"]),
                    abs(start - s["end"]),
                    abs(end - s["start"]),
                    abs(end - s["end"]),
                )
                <= args.straglr_breakpoint_tol
            )

            if overlap > 0 or (svtype in {"INS", "BND"} and near_bp):
                smatches.append(s["label"])
                if s["genes"] not in {"", ".", "nan"}:
                    sgenes.append(s["genes"])
                if s["copy_number"] not in {"", ".", "nan"}:
                    scn.append(s["copy_number"])
                if s["supporting_reads"] not in {"", ".", "nan"}:
                    ssupport.append(s["supporting_reads"])

        if smatches:
            out.at[idx, "STRAGLR_MATCH"] = "YES"
            out.at[idx, "STRAGLR_LOCI"] = ";".join(sorted(set(smatches)))
            out.at[idx, "STRAGLR_GENES"] = ";".join(sorted(set(sgenes))) if sgenes else "."
            out.at[idx, "STRAGLR_COPY_NUMBER"] = ";".join(sorted(set(scn))) if scn else "."
            out.at[idx, "STRAGLR_SUPPORTING_READS"] = ";".join(sorted(set(ssupport))) if ssupport else "."

        if svtype in {"INS", "BND"}:
            tmatches = []
            tuuids = []
            tused = []
            tspan = []

            for t in tldr:
                if t["chrom"] != chrom or t["start"] is None:
                    continue

                d1 = abs(start - t["start"])
                d2 = abs(end - t["start"]) if end is not None else d1
                if min(d1, d2) <= args.tldr_breakpoint_tol:
                    tmatches.append(t["label"])
                    if t["uuid"] not in {"", ".", "nan"}:
                        tuuids.append(t["uuid"])
                    if t["used_reads"] not in {"", ".", "nan"}:
                        tused.append(t["used_reads"])
                    if t["span_reads"] not in {"", ".", "nan"}:
                        tspan.append(t["span_reads"])

            if tmatches:
                out.at[idx, "TLDR_MATCH"] = "YES"
                out.at[idx, "TLDR_INSERTIONS"] = ";".join(sorted(set(tmatches)))
                out.at[idx, "TLDR_UUID"] = ";".join(sorted(set(tuuids))) if tuuids else "."
                out.at[idx, "TLDR_USED_READS"] = ";".join(sorted(set(tused))) if tused else "."
                out.at[idx, "TLDR_SPAN_READS"] = ";".join(sorted(set(tspan))) if tspan else "."

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)

    print(
        f"[OK] rows={len(out)} "
        f"straglr_matches={(out['STRAGLR_MATCH'] == 'YES').sum()} "
        f"tldr_matches={(out['TLDR_MATCH'] == 'YES').sum()} "
        f"output={output}"
    )


if __name__ == "__main__":
    main()
