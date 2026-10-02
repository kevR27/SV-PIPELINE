#!/usr/bin/env python3
"""Normalize CNVpytor read-depth CNV calls and BAF genotyping evidence.

CNVpytor calls are retained rather than hard-filtered. Q0, pN and distance to
large reference gaps are exposed as soft evidence flags so candidate review
can distinguish weak genomic context from absence of a CNV.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


MISSING = "."


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_region(text: str):
    match = re.match(r"^([^:]+):(\d+)-(\d+)$", text.replace(",", ""))
    if not match:
        raise ValueError(f"Unrecognized CNVpytor region: {text!r}")
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def read_genotypes(path: str) -> dict[str, list[str]]:
    result = {}
    if not path:
        return result
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2 or ":" not in fields[0]:
                continue
            result[fields[0].replace(",", "")] = fields[1:]
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calls", required=True)
    parser.add_argument("--genotypes", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--max-q0", type=float, default=0.5)
    parser.add_argument("--max-pn", type=float, default=0.5)
    parser.add_argument("--min-gap-distance", type=float, default=100000.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    genotypes = read_genotypes(args.genotypes)
    rows = []
    with open(args.calls, encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 9 or ":" not in fields[1]:
                continue

            cnv_type = fields[0].lower()
            if cnv_type not in {"deletion", "duplication"}:
                continue

            region = fields[1].replace(",", "")
            chrom, start, end = parse_region(region)
            svtype = "DEL" if cnv_type == "deletion" else "DUP"
            length = max(1, end - start + 1)
            signed_length = -length if svtype == "DEL" else length

            # Standard CNVpytor call columns: type, region, size, level,
            # e-val1..4, q0, pN, dG. Keep any additional columns separately.
            padded = fields + [MISSING] * max(0, 12 - len(fields))
            level = padded[3]
            e1, e2, e3, e4 = padded[4:8]
            q0, pn, dg = padded[8:11]

            flags = []
            q0_value = as_float(q0)
            pn_value = as_float(pn)
            dg_value = as_float(dg)
            if q0_value is not None and q0_value > args.max_q0:
                flags.append("HIGH_Q0_FRACTION")
            if pn_value is not None and pn_value > args.max_pn:
                flags.append("HIGH_REFERENCE_N_FRACTION")
            if dg_value is not None and dg_value < args.min_gap_distance:
                flags.append("NEAR_LARGE_REFERENCE_GAP")

            gt = genotypes.get(region, [])
            gt = gt + [MISSING] * max(0, 11 - len(gt))
            # For -genotype BIN -a output, the region is followed by:
            # level,e1,e2,q0,pN,dG,bin proportion,hom,het,BAF shift,BAF p.
            row = {
                "SAMPLE": args.sample,
                "CALLER": "CNVpytor",
                "SV_ID": f"CNVPYTOR_{args.sample}_{index}",
                "CHROM": chrom,
                "START": str(start),
                "END": str(end),
                "SVTYPE": svtype,
                "SVLEN": str(signed_length),
                "CNVPYTOR_SIZE_REPORTED": padded[2],
                "CNVPYTOR_LEVEL": level,
                "CNVPYTOR_EVAL1": e1,
                "CNVPYTOR_EVAL2": e2,
                "CNVPYTOR_EVAL3": e3,
                "CNVPYTOR_EVAL4": e4,
                "CNVPYTOR_Q0": q0,
                "CNVPYTOR_PN": pn,
                "CNVPYTOR_DG": dg,
                "CNVPYTOR_GT_LEVEL": gt[0],
                "CNVPYTOR_GT_EVAL1": gt[1],
                "CNVPYTOR_GT_EVAL2": gt[2],
                "CNVPYTOR_GT_Q0": gt[3],
                "CNVPYTOR_GT_PN": gt[4],
                "CNVPYTOR_GT_DG": gt[5],
                "CNVPYTOR_BIN_PROPORTION": gt[6],
                "CNVPYTOR_HOM_COUNT": gt[7],
                "CNVPYTOR_HET_COUNT": gt[8],
                "CNVPYTOR_BAF_SHIFT": gt[9],
                "CNVPYTOR_BAF_PVALUE": gt[10],
                "EVIDENCE_STATUS": "RETAINED",
                "EVIDENCE_FLAGS": ";".join(flags) if flags else ".",
            }
            rows.append(row)

    columns = [
        "SAMPLE","CALLER","SV_ID","CHROM","START","END","SVTYPE","SVLEN",
        "CNVPYTOR_SIZE_REPORTED","CNVPYTOR_LEVEL","CNVPYTOR_EVAL1","CNVPYTOR_EVAL2",
        "CNVPYTOR_EVAL3","CNVPYTOR_EVAL4","CNVPYTOR_Q0","CNVPYTOR_PN","CNVPYTOR_DG",
        "CNVPYTOR_GT_LEVEL","CNVPYTOR_GT_EVAL1","CNVPYTOR_GT_EVAL2","CNVPYTOR_GT_Q0",
        "CNVPYTOR_GT_PN","CNVPYTOR_GT_DG","CNVPYTOR_BIN_PROPORTION",
        "CNVPYTOR_HOM_COUNT","CNVPYTOR_HET_COUNT","CNVPYTOR_BAF_SHIFT",
        "CNVPYTOR_BAF_PVALUE","EVIDENCE_STATUS","EVIDENCE_FLAGS",
    ]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] CNVpytor calls={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
