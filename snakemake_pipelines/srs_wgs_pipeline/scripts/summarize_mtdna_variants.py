#!/usr/bin/env python3
"""
Convert the Mutserve2 mtDNA VCF into a small, readable TSV.

The table keeps the key fields needed for thesis review:
position, reference/alternate allele, genotype, heteroplasmy fraction and depth.

Important:
Mutserve2 VCF output is used here for mtDNA SNV/heteroplasmy analysis.
It is not treated as a complete mtDNA indel workflow.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(text: str) -> dict[str, str]:
    result: dict[str, str] = {}

    for item in str(text).split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "True"

    return result


def parse_sample(format_text: str, sample_text: str) -> dict[str, str]:
    if not format_text or not sample_text:
        return {}

    return dict(
        zip(
            format_text.split(":"),
            sample_text.split(":"),
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize Mutserve2 mtDNA variants."
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []

    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            info = parse_info(fields[7])

            sample = {}
            if len(fields) > 9:
                sample = parse_sample(
                    fields[8],
                    fields[9],
                )

            rows.append(
                {
                    "CHROM": fields[0],
                    "POS": fields[1],
                    "ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "FILTER": fields[6],
                    "GT": sample.get("GT", "."),
                    "HETEROPLASMY_AF": sample.get(
                        "AF",
                        info.get("AF", "."),
                    ),
                    "DP": sample.get(
                        "DP",
                        info.get("DP", "."),
                    ),
                    "ANALYSIS_SCOPE": (
                        "mtDNA SNV/heteroplasmy; "
                        "not a complete mtDNA indel analysis"
                    ),
                }
            )

    columns = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "FILTER",
        "GT",
        "HETEROPLASMY_AF",
        "DP",
        "ANALYSIS_SCOPE",
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[OK] mtDNA_variants={len(rows)} output={output}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
