#!/usr/bin/env python3
"""Convert normalized CNVpytor DEL/DUP calls to a site-only VCF for AnnotSV.

No diploid genotype is invented. CNVpytor's normalized read-depth level and
BAF metrics remain INFO evidence.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import pysam


MISSING = "."


def finite(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def info_number(key: str, value: str):
    x = finite(value)
    return f"{key}={x:g}" if x is not None else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    fasta = pysam.FastaFile(args.reference)
    with open(args.tsv, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=CNVpytor_to_AnnotSV\n")
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Signed SV length">\n')
        out.write('##INFO=<ID=CNVPYTOR_LEVEL,Number=1,Type=Float,Description="CNVpytor read-depth level normalized to 1">\n')
        out.write('##INFO=<ID=CNVPYTOR_EVAL1,Number=1,Type=Float,Description="CNVpytor e-value 1">\n')
        out.write('##INFO=<ID=CNVPYTOR_Q0,Number=1,Type=Float,Description="CNVpytor q0 read fraction">\n')
        out.write('##INFO=<ID=CNVPYTOR_PN,Number=1,Type=Float,Description="CNVpytor reference-N fraction">\n')
        out.write('##INFO=<ID=CNVPYTOR_DG,Number=1,Type=Float,Description="Distance to nearest large reference gap">\n')
        out.write('##INFO=<ID=CNVPYTOR_BAFSHIFT,Number=1,Type=Float,Description="CNVpytor BAF shift from 0.5">\n')
        out.write('##INFO=<ID=CNVPYTOR_BAFP,Number=1,Type=Float,Description="CNVpytor BAF p-value">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for row in rows:
            chrom = row["CHROM"]
            start = int(row["START"])
            end = int(row["END"])
            svtype = row["SVTYPE"]
            svlen = int(row["SVLEN"])
            try:
                ref = fasta.fetch(chrom, start - 1, start).upper() or "N"
            except (KeyError, ValueError):
                ref = "N"

            info = [f"END={end}", f"SVTYPE={svtype}", f"SVLEN={svlen}"]
            for key, column in (
                ("CNVPYTOR_LEVEL","CNVPYTOR_LEVEL"),
                ("CNVPYTOR_EVAL1","CNVPYTOR_EVAL1"),
                ("CNVPYTOR_Q0","CNVPYTOR_Q0"),
                ("CNVPYTOR_PN","CNVPYTOR_PN"),
                ("CNVPYTOR_DG","CNVPYTOR_DG"),
                ("CNVPYTOR_BAFSHIFT","CNVPYTOR_BAF_SHIFT"),
                ("CNVPYTOR_BAFP","CNVPYTOR_BAF_PVALUE"),
            ):
                value = info_number(key, row.get(column, MISSING))
                if value:
                    info.append(value)

            out.write(
                f"{chrom}\t{start}\t{row['SV_ID']}\t{ref}\t<{svtype}>\t.\t.\t"
                + ";".join(info) + "\n"
            )

    fasta.close()
    print(f"[OK] records={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
