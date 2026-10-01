#!/usr/bin/env python3
"""N50 (all reads and long reads only) from a `dorado summary` TSV."""

import argparse
import csv
import sys


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    cum = 0
    for length in lengths:
        cum += length
        if cum >= half:
            return length
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", required=True, help="dorado summary TSV")
    p.add_argument("--min-length", type=int, default=1000,
                   help="minimum read length (bp) for the long-read N50")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    col = "sequence_length_template"
    seen, lengths = set(), []
    with open(a.summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if col not in (reader.fieldnames or []):
            sys.exit(f"Column '{col}' not found in {a.summary}")
        for row in reader:
            rid = row.get("read_id")
            if rid is not None:          # one row per read, even if the BAM
                if rid in seen:          # has secondary/supplementary records
                    continue
                seen.add(rid)
            lengths.append(int(row[col]))

    long_reads = [x for x in lengths if x >= a.min_length]

    with open(a.out, "w") as out:
        out.write(f"metric\tall_reads\tlong_reads_ge_{a.min_length}\n")
        out.write(f"n_reads\t{len(lengths)}\t{len(long_reads)}\n")
        out.write(f"total_bases\t{sum(lengths)}\t{sum(long_reads)}\n")
        out.write(f"N50\t{n50(lengths)}\t{n50(long_reads)}\n")


if __name__ == "__main__":
    main()
