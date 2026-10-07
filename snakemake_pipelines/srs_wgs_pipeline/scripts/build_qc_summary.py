#!/usr/bin/env python3
"""Create a compact WGS QC gate from mosdepth and samtools outputs."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def flagstat(path: str) -> tuple[int | None, int | None]:
    total = mapped = None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"(\d+) \+ \d+ (.+)", line)
        if not m:
            continue
        value, label = int(m.group(1)), m.group(2)
        if label.startswith("in total"):
            total = value
        elif label.startswith("mapped ("):
            mapped = value
    return total, mapped


def mean_coverage(path: str) -> float | None:
    with open(path, encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row = next((r for r in rows if r.get("chrom") == "total"), rows[-1] if rows else None)
    try:
        return float(row["mean"]) if row else None
    except (KeyError, TypeError, ValueError):
        return None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True); p.add_argument("--flagstat", required=True)
    p.add_argument("--mosdepth-summary", required=True); p.add_argument("--min-mean-coverage", type=float, default=25)
    p.add_argument("--min-mapped-percent", type=float, default=95); p.add_argument("--output", required=True)
    a = p.parse_args()
    total, mapped = flagstat(a.flagstat); coverage = mean_coverage(a.mosdepth_summary)
    mapped_pct = (100.0 * mapped / total) if total and mapped is not None else None
    flags = []
    if coverage is None: flags.append("MISSING_MEAN_COVERAGE")
    elif coverage < a.min_mean_coverage: flags.append("LOW_MEAN_COVERAGE")
    if mapped_pct is None: flags.append("MISSING_MAPPED_PERCENT")
    elif mapped_pct < a.min_mapped_percent: flags.append("LOW_MAPPED_PERCENT")
    row = {"sample": a.sample, "qc_status": "PASS" if not flags else "REVIEW", "qc_flags": ";".join(flags) if flags else ".",
           "mean_coverage": "." if coverage is None else f"{coverage:.4g}", "minimum_mean_coverage": a.min_mean_coverage,
           "mapped_percent": "." if mapped_pct is None else f"{mapped_pct:.4g}", "minimum_mapped_percent": a.min_mapped_percent,
           "total_reads": total if total is not None else ".", "mapped_reads": mapped if mapped is not None else "."}
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
