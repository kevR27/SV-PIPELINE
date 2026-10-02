#!/usr/bin/env python3
"""Create a compact Illumina WGS QC table from samtools and mosdepth outputs.

This script is descriptive only. It deliberately does not convert QC metrics
into a PASS/FAIL clinical label because acceptable thresholds depend on the
sequencing protocol and the validated laboratory workflow.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


MISSING = "."


def first_number(text: str) -> str:
    match = re.match(r"^\s*([0-9]+(?:\.[0-9]+)?)", text)
    return match.group(1) if match else MISSING


def parse_flagstat(path: str) -> dict[str, str]:
    out = {
        "TOTAL_READS": MISSING,
        "MAPPED_READS": MISSING,
        "MAPPED_PERCENT": MISSING,
        "PROPERLY_PAIRED_READS": MISSING,
        "PROPERLY_PAIRED_PERCENT": MISSING,
        "DUPLICATE_READS": MISSING,
    }
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if " in total (QC-passed reads + QC-failed reads)" in line:
                out["TOTAL_READS"] = first_number(line)
            elif " mapped (" in line and "primary mapped" not in line:
                out["MAPPED_READS"] = first_number(line)
                m = re.search(r"\(([-0-9.]+)%", line)
                if m:
                    out["MAPPED_PERCENT"] = m.group(1)
            elif " properly paired (" in line:
                out["PROPERLY_PAIRED_READS"] = first_number(line)
                m = re.search(r"\(([-0-9.]+)%", line)
                if m:
                    out["PROPERLY_PAIRED_PERCENT"] = m.group(1)
            elif " duplicates" in line:
                out["DUPLICATE_READS"] = first_number(line)
    return out


def parse_stats(path: str) -> dict[str, str]:
    wanted = {
        "average length": "AVERAGE_READ_LENGTH",
        "insert size average": "MEAN_INSERT_SIZE",
        "insert size standard deviation": "SD_INSERT_SIZE",
        "error rate": "ERROR_RATE",
    }
    out = {column: MISSING for column in wanted.values()}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("SN\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            key = fields[1].rstrip(":")
            if key in wanted:
                out[wanted[key]] = fields[2]
    return out


def parse_mosdepth_global_dist(path: str) -> dict[str, str]:
    """Read mosdepth cumulative depth distribution for the total genome."""
    depths = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0].lower() != "total":
                continue
            try:
                depth = int(float(fields[1]))
                fraction = float(fields[2])
            except ValueError:
                continue
            depths[depth] = fraction

    out = {
        "MEDIAN_COVERAGE_FROM_GLOBAL_DIST": MISSING,
        "FRACTION_BASES_GE_10X": MISSING,
        "FRACTION_BASES_GE_20X": MISSING,
        "FRACTION_BASES_GE_30X": MISSING,
        "PERCENT_BASES_GE_10X": MISSING,
        "PERCENT_BASES_GE_20X": MISSING,
        "PERCENT_BASES_GE_30X": MISSING,
    }
    if not depths:
        return out

    # mosdepth global.dist is cumulative: fraction of bases at or above depth.
    median_candidates = [depth for depth, fraction in depths.items() if fraction >= 0.5]
    if median_candidates:
        out["MEDIAN_COVERAGE_FROM_GLOBAL_DIST"] = str(max(median_candidates))

    for threshold in (10, 20, 30):
        fraction = depths.get(threshold)
        if fraction is None:
            continue
        out[f"FRACTION_BASES_GE_{threshold}X"] = f"{fraction:.8g}"
        out[f"PERCENT_BASES_GE_{threshold}X"] = f"{100.0 * fraction:.6f}"
    return out


def parse_mosdepth_summary(path: str) -> dict[str, str]:
    out = {
        "AUTOSOMAL_OR_TOTAL_MEAN_COVERAGE": MISSING,
        "MOSDEPTH_TOTAL_LENGTH": MISSING,
        "MOSDEPTH_TOTAL_BASES": MISSING,
    }
    with open(path, encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        return out

    total = next((r for r in rows if str(r.get("chrom", "")).lower() == "total"), None)
    if total is None:
        total = rows[-1]

    out["AUTOSOMAL_OR_TOTAL_MEAN_COVERAGE"] = str(total.get("mean", MISSING))
    out["MOSDEPTH_TOTAL_LENGTH"] = str(total.get("length", MISSING))
    out["MOSDEPTH_TOTAL_BASES"] = str(total.get("bases", MISSING))
    return out


def duplicate_percent(duplicates: str, total: str) -> str:
    try:
        denominator = float(total)
        if denominator <= 0:
            return MISSING
        return f"{100.0 * float(duplicates) / denominator:.6f}"
    except (TypeError, ValueError):
        return MISSING


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--flagstat", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--mosdepth-summary", required=True)
    parser.add_argument("--mosdepth-global-dist", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    row = {"SAMPLE": args.sample}
    row.update(parse_flagstat(args.flagstat))
    row.update(parse_stats(args.stats))
    row.update(parse_mosdepth_summary(args.mosdepth_summary))
    row.update(parse_mosdepth_global_dist(args.mosdepth_global_dist))
    row["DUPLICATE_PERCENT_OF_TOTAL"] = duplicate_percent(
        row["DUPLICATE_READS"], row["TOTAL_READS"]
    )

    columns = list(row)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)

    print(f"[OK] sample={args.sample} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
