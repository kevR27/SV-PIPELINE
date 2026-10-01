#!/usr/bin/env python3
"""Read-length QC from a Dorado sequencing-summary TSV.

Reports N50 and basic length statistics for all reads, a configurable long-read
subset, and the complementary short-read tail. The threshold separates the
subsets; N50 itself is a distribution summary and is not a read filter.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys


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
    p.add_argument("--summary", required=True, help="Dorado summary TSV")
    p.add_argument(
        "--min-length",
        type=int,
        default=1000,
        help="Minimum read length in bp for the long-read QC subset",
    )
    p.add_argument("--out", required=True)
    a = p.parse_args()

    length_col = "sequence_length_template"
    seen = set()
    lengths = []

    with open(a.summary, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if length_col not in (reader.fieldnames or []):
            sys.exit(
                f"Column '{length_col}' not found in {a.summary}. "
                "Use a Dorado sequencing-summary file produced from compatible "
                "Dorado BAM/basecalling data."
            )

        for row in reader:
            read_id = str(row.get("read_id", "")).strip()
            if read_id:
                if read_id in seen:
                    continue
                seen.add(read_id)

            raw = str(row.get(length_col, "")).strip()
            if not raw:
                continue
            try:
                length = int(float(raw))
            except ValueError:
                continue
            if length >= 0:
                lengths.append(length)

    long_reads = [x for x in lengths if x >= a.min_length]
    short_reads = [x for x in lengths if x < a.min_length]

    groups = {
        "all_reads": lengths,
        f"long_reads_ge_{a.min_length}": long_reads,
        f"short_reads_lt_{a.min_length}": short_reads,
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

    group_names = list(groups)
    with open(a.out, "w", encoding="utf-8") as out:
        out.write("metric\t" + "\t".join(group_names) + "\n")
        for metric, values in metrics.items():
            out.write(
                metric
                + "\t"
                + "\t".join(str(values[name]) for name in group_names)
                + "\n"
            )


if __name__ == "__main__":
    main()
