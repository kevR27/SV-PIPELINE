#!/usr/bin/env python3
"""Read-length QC directly from the BAM used by the LRS workflow.

Only primary alignments are counted, so supplementary/secondary records do not
inflate read counts. Metrics are reported for all primary reads, the configured
long-read subset, and the complementary short-read tail.
"""
from __future__ import annotations

import argparse
import statistics
from pathlib import Path

import pysam


def n50(lengths):
    if not lengths:
        return 0
    values = sorted(lengths, reverse=True)
    half = sum(values) / 2
    cumulative = 0
    for length in values:
        cumulative += length
        if cumulative >= half:
            return int(length)
    return 0


def mean_value(values):
    return round(sum(values) / len(values), 2) if values else 0


def median_value(values):
    return round(float(statistics.median(values)), 2) if values else 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--min-length", type=int, default=1000)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    lengths = []
    seen = set()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue

            name = read.query_name
            if name and name in seen:
                continue
            if name:
                seen.add(name)

            length = read.infer_read_length()
            if length is None:
                length = read.query_length
            if length is None or length < 0:
                continue
            lengths.append(int(length))

    long_reads = [x for x in lengths if x >= args.min_length]
    short_reads = [x for x in lengths if x < args.min_length]

    groups = {
        "all_primary_reads": lengths,
        f"long_reads_ge_{args.min_length}": long_reads,
        f"short_reads_lt_{args.min_length}": short_reads,
    }

    total_reads = len(lengths)
    total_bases = sum(lengths)

    metrics = {
        "n_reads": {name: len(vals) for name, vals in groups.items()},
        "total_bases": {name: sum(vals) for name, vals in groups.items()},
        "mean_read_length": {name: mean_value(vals) for name, vals in groups.items()},
        "median_read_length": {name: median_value(vals) for name, vals in groups.items()},
        "N50": {name: n50(vals) for name, vals in groups.items()},
        "read_fraction_percent": {
            name: round(100 * len(vals) / total_reads, 3) if total_reads else 0
            for name, vals in groups.items()
        },
        "base_fraction_percent": {
            name: round(100 * sum(vals) / total_bases, 3) if total_bases else 0
            for name, vals in groups.items()
        },
    }

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as out:
        names = list(groups)
        out.write("metric\t" + "\t".join(names) + "\n")
        for metric, values in metrics.items():
            out.write(
                metric + "\t" + "\t".join(str(values[name]) for name in names) + "\n"
            )

    print(
        f"[OK] primary_reads={total_reads} long_reads={len(long_reads)} "
        f"short_reads={len(short_reads)} output={output}"
    )


if __name__ == "__main__":
    main()
