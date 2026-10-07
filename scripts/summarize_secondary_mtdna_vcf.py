#!/usr/bin/env python3
"""Create a readable table for secondary mtDNA candidate review.

This script reports variants already present in a caller VCF. It does not
reclassify them and does not claim that an allele fraction is a validated
heteroplasmy measurement.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_vcf(path: str):
    """Open either a plain-text or bgzip-compressed VCF."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def split_key_values(text: str, item_separator: str) -> dict[str, str]:
    """Read KEY=VALUE fields such as VCF INFO entries."""
    values: dict[str, str] = {}
    for item in text.split(item_separator):
        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value
    return values


def sample_values(format_field: str, sample_field: str) -> dict[str, str]:
    """Match VCF FORMAT labels to values from the first sample column."""
    if not format_field or not sample_field:
        return {}
    return dict(zip(format_field.split(":"), sample_field.split(":")))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    with open_vcf(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            info = split_key_values(fields[7], ";")
            sample = sample_values(
                fields[8] if len(fields) > 8 else "",
                fields[9] if len(fields) > 9 else "",
            )

            rows.append(
                {
                    "CHROM": fields[0],
                    "POS": fields[1],
                    "ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "QUAL": fields[5],
                    "FILTER": fields[6],
                    "GT": sample.get("GT", "."),
                    "CANDIDATE_ALLELE_FRACTION": sample.get(
                        "AF", info.get("AF", ".")
                    ),
                    "DP": sample.get("DP", info.get("DP", ".")),
                    "SOURCE": args.source,
                    "ANALYSIS_ROLE": "SECONDARY_MTDNA_REVIEW",
                    "INTERPRETATION_LIMIT": (
                        "Long-read small-variant candidate; not a validated "
                        "heteroplasmy or complete mtDNA analysis"
                    ),
                }
            )

    columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "GT",
        "CANDIDATE_ALLELE_FRACTION", "DP", "SOURCE", "ANALYSIS_ROLE",
        "INTERPRETATION_LIMIT",
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] secondary_mtDNA_candidates={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
