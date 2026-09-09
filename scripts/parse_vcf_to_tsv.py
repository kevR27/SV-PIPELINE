#!/usr/bin/env python3
"""Expand a VCF/VCF.GZ into a TSV without losing INFO or FORMAT data."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import TextIO


MISSING = "."


# ============================================================
# FILE HANDLING
# ============================================================

def open_text(path: str) -> TextIO:
    """Open a plain-text or gzipped text file."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


# ============================================================
# FIELD PARSING
# ============================================================

def clean_key(key: str) -> str:
    """Make a VCF key safe for use as a TSV column name."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", key).strip("_") or "UNKNOWN"


def parse_info(raw: str) -> OrderedDict[str, str]:
    """Parse the VCF INFO column into an ordered dictionary."""
    out = OrderedDict()
    if not raw or raw == MISSING:
        return out

    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = "True"

    return out


def parse_format(raw: str, samples: list[str], values: list[str]) -> OrderedDict[str, str]:
    """Expand FORMAT fields into sample-specific TSV columns."""
    out = OrderedDict()
    if not raw or raw == MISSING:
        return out

    keys = raw.split(":")
    for i, sample in enumerate(samples):
        if i >= len(values):
            break

        vals = values[i].split(":")
        sample_key = clean_key(sample)

        for j, key in enumerate(keys):
            column = f"FORMAT_{sample_key}_{clean_key(key)}"
            out[column] = vals[j] if j < len(vals) else MISSING

    return out


def first(info: dict[str, str], *keys: str) -> str:
    """Return the first non-missing value among the requested keys."""
    for key in keys:
        value = info.get(key, MISSING)
        if value not in ("", MISSING):
            return value
    return MISSING


# ============================================================
# SV INFORMATION INFERENCE
# ============================================================

def infer_svtype(info: dict[str, str], alt: str) -> str:
    """Infer SVTYPE from INFO or ALT when necessary."""
    value = first(info, "SVTYPE")
    if value != MISSING:
        return value
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1]
    return MISSING


def infer_end(pos: int, ref: str, info: dict[str, str]) -> str:
    """Infer END from INFO/END or the reference allele length."""
    value = first(info, "END")
    if value != MISSING:
        return value
    return str(pos + max(len(ref), 1) - 1) if ref not in ("", MISSING) else MISSING


# ============================================================
# VCF READING
# ============================================================

def read_vcf_header(vcf_path: str) -> tuple[OrderedDict[str, None], list[str]]:
    """Read INFO keys and sample names from the VCF header."""
    info_keys = OrderedDict()
    samples = []

    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##INFO=<"):
                match = re.search(r"ID=([^,>]+)", line)
                if match:
                    info_keys.setdefault(match.group(1), None)
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break

    return info_keys, samples


def read_vcf_records(
    vcf_path: str,
    samples: list[str],
    info_keys: OrderedDict[str, None],
) -> tuple[list[OrderedDict[str, str]], OrderedDict[str, None]]:
    """Read VCF records and convert them into ordered dictionaries."""
    rows = []
    format_keys = OrderedDict()

    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line.rstrip()}")

            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            format_raw = fields[8] if len(fields) > 8 else MISSING
            sample_values = fields[9:] if len(fields) > 9 else []

            try:
                pos_int = int(pos)
            except ValueError as error:
                raise ValueError(f"Invalid POS {pos!r}") from error

            info = parse_info(info_raw)
            for key in info:
                info_keys.setdefault(key, None)

            fmt = parse_format(format_raw, samples, sample_values)
            for key in fmt:
                format_keys.setdefault(key, None)

            row = OrderedDict({
                "SV_ID": sv_id,
                "CHROM": chrom,
                "START": pos,
                "END": infer_end(pos_int, ref, info),
                "SVTYPE": infer_svtype(info, alt),
                "SVLEN": first(info, "SVLEN"),
                "QUAL": qual,
                "FILTER": filt,
                "REF": ref,
                "ALT": alt,
                "INFO_RAW": info_raw,
                "FORMAT_RAW": format_raw,
            })

            for key in info_keys:
                row[f"INFO_{clean_key(key)}"] = info.get(key, MISSING)
            for key in format_keys:
                row[key] = fmt.get(key, MISSING)

            rows.append(row)

    return rows, format_keys


# ============================================================
# TSV WRITING
# ============================================================

def build_output_columns(
    info_keys: OrderedDict[str, None],
    format_keys: OrderedDict[str, None],
    caller: str,
) -> list[str]:
    """Build the final TSV column order."""
    columns = ["CALLER"] if caller else []
    columns += [
        "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN",
        "QUAL", "FILTER", "REF", "ALT", "INFO_RAW", "FORMAT_RAW",
    ]
    columns += [f"INFO_{clean_key(key)}" for key in info_keys]
    columns += list(format_keys)
    return columns


def write_tsv(output_path: str, rows: list[OrderedDict[str, str]], columns: list[str]) -> None:
    """Write the parsed records to a TSV file."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(column, MISSING) for column in columns) + "\n")


# ============================================================
# MAIN
# ============================================================

def main() -> int:
    """Run the VCF-to-TSV conversion."""
    parser = argparse.ArgumentParser(description="Expand every VCF INFO/FORMAT key into TSV columns.")
    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ file.")
    parser.add_argument("--output", required=True, help="Output TSV file.")
    parser.add_argument("--caller", default="", help="Optional SV caller name.")
    args = parser.parse_args()

    info_keys, samples = read_vcf_header(args.vcf)
    rows, format_keys = read_vcf_records(args.vcf, samples, info_keys)

    if args.caller:
        for row in rows:
            row["CALLER"] = args.caller

    columns = build_output_columns(info_keys, format_keys, args.caller)
    write_tsv(args.output, rows, columns)

    print(f"[OK] records={len(rows)} info_fields={len(info_keys)} format_fields={len(format_keys)}")
    print(f"[OK] output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
