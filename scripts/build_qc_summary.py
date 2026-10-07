#!/usr/bin/env python3
"""Summarize the main alignment and coverage checks for one SRS sample."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def read_flagstat(path: str) -> tuple[int | None, int | None]:
    """Return total and mapped read counts from samtools flagstat."""
    total_reads = None
    mapped_reads = None

    lines = Path(path).read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    for line in lines:
        match = re.match(r"(\d+) \+ \d+ (.+)", line)
        if not match:
            continue

        value = int(match.group(1))
        label = match.group(2)

        if label.startswith("in total"):
            total_reads = value
        elif label.startswith("mapped ("):
            mapped_reads = value

    return total_reads, mapped_reads


def read_mean_coverage(path: str) -> float | None:
    """Read whole-genome mean coverage from a mosdepth summary table."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        return None

    total_row = next(
        (row for row in rows if row.get("chrom") == "total"),
        rows[-1],
    )

    try:
        return float(total_row["mean"])
    except (KeyError, TypeError, ValueError):
        return None


def qc_flags(
    mean_coverage: float | None,
    mapped_percent: float | None,
    minimum_coverage: float,
    minimum_mapped_percent: float,
) -> list[str]:
    """Describe missing or below-threshold QC measurements."""
    flags = []

    if mean_coverage is None:
        flags.append("MISSING_MEAN_COVERAGE")
    elif mean_coverage < minimum_coverage:
        flags.append("LOW_MEAN_COVERAGE")

    if mapped_percent is None:
        flags.append("MISSING_MAPPED_PERCENT")
    elif mapped_percent < minimum_mapped_percent:
        flags.append("LOW_MAPPED_PERCENT")

    return flags


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--flagstat", required=True)
    parser.add_argument("--mosdepth-summary", required=True)
    parser.add_argument("--min-mean-coverage", type=float, default=25)
    parser.add_argument("--min-mapped-percent", type=float, default=95)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    total_reads, mapped_reads = read_flagstat(args.flagstat)
    mean_coverage = read_mean_coverage(args.mosdepth_summary)

    if total_reads and mapped_reads is not None:
        mapped_percent = 100.0 * mapped_reads / total_reads
    else:
        mapped_percent = None

    flags = qc_flags(
        mean_coverage,
        mapped_percent,
        args.min_mean_coverage,
        args.min_mapped_percent,
    )

    summary = {
        "sample": args.sample,
        "qc_status": "PASS" if not flags else "REVIEW",
        "qc_flags": ";".join(flags) if flags else ".",
        "mean_coverage": (
            "." if mean_coverage is None else f"{mean_coverage:.4g}"
        ),
        "minimum_mean_coverage": args.min_mean_coverage,
        "mapped_percent": (
            "." if mapped_percent is None else f"{mapped_percent:.4g}"
        ),
        "minimum_mapped_percent": args.min_mapped_percent,
        "total_reads": total_reads if total_reads is not None else ".",
        "mapped_reads": mapped_reads if mapped_reads is not None else ".",
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(summary)

    print(f"[OK] sample={args.sample} qc_status={summary['qc_status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
