#!/usr/bin/env python3
"""
Create a caller-aware SV evidence table from a VCF/VCF.GZ.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import TextIO


MISSING = "."


# ============================================================
# FILE HANDLING
# ============================================================

def open_text(path: str) -> TextIO:
    """Open plain-text or gzip-compressed text files."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


# ============================================================
# VCF FIELD PARSING
# ============================================================

def parse_info(raw: str) -> dict[str, str]:
    """
    Parse a VCF INFO field into a dictionary.

    Example:
        SVTYPE=DEL;END=1000;SUPPORT=8
    becomes:
        {"SVTYPE": "DEL", "END": "1000", "SUPPORT": "8"}

    Flag-style INFO fields are assigned the value "True".
    """
    info = {}

    if not raw or raw == MISSING:
        return info

    for item in raw.split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "True"

    return info


def first(info: dict[str, str], *keys: str) -> str:
    """Return the first non-empty value found among the supplied keys."""
    for key in keys:
        value = info.get(key, MISSING)
        if value not in {"", MISSING}:
            return value

    return MISSING


def infer_svtype(info: dict[str, str], alt: str) -> str:
    """Infer SVTYPE from INFO or the ALT representation."""
    value = first(info, "SVTYPE")

    if value != MISSING:
        return value
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1]

    return MISSING


# ============================================================
# CALLER-SPECIFIC EVIDENCE
# ============================================================

def evidence_sniffles2(info: dict[str, str]) -> dict[str, str]:
    """Extract Sniffles2-specific evidence."""
    return {
        "CALLER_SUPPORT": first(info, "SUPPORT", "RE", "SUPP"),
        "CALLER_RNAMES": first(info, "RNAMES"),
        "CALLER_STRANDS": first(info, "STRANDS"),
        "CALLER_IMPRECISE": first(info, "IMPRECISE"),
        "CALLER_MOSAIC": first(info, "MOSAIC"),
    }


def evidence_cutesv(info: dict[str, str]) -> dict[str, str]:
    """Extract cuteSV-specific evidence."""
    return {
        "CALLER_SUPPORT": first(info, "RE", "SUPPORT", "SUPP"),
        "CALLER_RNAMES": first(info, "RNAMES"),
        "CALLER_STRANDS": first(info, "STRANDS"),
    }


def evidence_delly(info: dict[str, str]) -> dict[str, str]:
    """Extract Delly-specific evidence."""
    return {
        "CALLER_SUPPORT": first(info, "SU", "SUPPORT"),
        "CALLER_PE": first(info, "PE"),
        "CALLER_SR": first(info, "SR"),
        "CALLER_PRECISE": first(info, "PRECISE"),
    }


def evidence_manta(info: dict[str, str]) -> dict[str, str]:
    """Extract Manta-specific evidence."""
    return {
        "CALLER_SUPPORT": first(info, "SU", "PR", "SR", "SUPPORT"),
        "CALLER_PR": first(info, "PR"),
        "CALLER_SR": first(info, "SR"),
    }


def evidence_jasmine(info: dict[str, str]) -> dict[str, str]:
    """Extract Jasmine/SURVIVOR merged-call evidence."""
    return {
        "CALLER_SUPPORT": first(info, "SUPP", "SUPPORT"),
        "CALLER_SUPPORT_VECTOR": first(info, "SUPP_VEC"),
    }


def evidence_generic(info: dict[str, str]) -> dict[str, str]:
    """Extract generic support information."""
    return {
        "CALLER_SUPPORT": first(info, "SUPPORT", "SUPP", "RE", "SU"),
    }


def evidence(caller: str, info: dict[str, str]) -> dict[str, str]:
    """Return caller-specific evidence fields."""
    caller_name = caller.lower()

    if caller_name == "sniffles2":
        return evidence_sniffles2(info)
    if caller_name == "cutesv":
        return evidence_cutesv(info)
    if caller_name == "delly":
        return evidence_delly(info)
    if caller_name == "manta":
        return evidence_manta(info)
    if caller_name in {"jasmine", "survivor"}:
        return evidence_jasmine(info)

    return evidence_generic(info)


# ============================================================
# MAIN
# ============================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse caller-specific SV evidence into a TSV."
    )

    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ")
    parser.add_argument("--caller", required=True, help="Caller name, e.g. sniffles2, cutesv, manta")
    parser.add_argument("--output", required=True, help="Output TSV")

    args = parser.parse_args()
    rows = []

    # --------------------------------------------------------
    # Read VCF records
    # --------------------------------------------------------

    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line.rstrip()}")

            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            info = parse_info(info_raw)

            row = {
                "CALLER": args.caller,
                "SV_ID": sv_id,
                "CHROM": chrom,
                "START": pos,
                "END": first(info, "END"),
                "SVTYPE": infer_svtype(info, alt),
                "SVLEN": first(info, "SVLEN"),
                "QUAL": qual,
                "FILTER": filt,
                "REF": ref,
                "ALT": alt,
                "INFO_RAW": info_raw,
            }

            row.update(evidence(args.caller, info))
            rows.append(row)

    # --------------------------------------------------------
    # Define output columns
    # --------------------------------------------------------

    columns = [
        "CALLER",
        "SV_ID",
        "CHROM",
        "START",
        "END",
        "SVTYPE",
        "SVLEN",
        "QUAL",
        "FILTER",
        "REF",
        "ALT",
        "INFO_RAW",
        "CALLER_SUPPORT",
        "CALLER_RNAMES",
        "CALLER_STRANDS",
        "CALLER_IMPRECISE",
        "CALLER_MOSAIC",
        "CALLER_PE",
        "CALLER_SR",
        "CALLER_PR",
        "CALLER_PRECISE",
        "CALLER_SUPPORT_VECTOR",
    ]

    # --------------------------------------------------------
    # Write TSV
    # --------------------------------------------------------

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] caller={args.caller} records={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
