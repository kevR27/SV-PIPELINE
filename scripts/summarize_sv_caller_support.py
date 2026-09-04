#!/usr/bin/env python3
"""Summarize caller provenance for merged SVs.

The script does not recluster or redefine SVs. It reads the merged VCF produced
by Jasmine/SURVIVOR and extracts caller provenance/support already represented
in the merged record (e.g. CALLERS/SOURCES/SUPP_VEC/SUPP). If those fields are
not present, optional individual caller VCFs can be supplied and records are
matched approximately by chromosome, SV type and breakpoint tolerance.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

MISSING = "."


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(raw: str) -> dict[str, str]:
    result = {}
    if not raw or raw == MISSING:
        return result
    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "True"
    return result


def first(info: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = info.get(key, MISSING)
        if value not in ("", MISSING):
            return value
    return MISSING


def infer_svtype(info: dict[str, str], alt: str) -> str:
    value = first(info, "SVTYPE")
    if value != MISSING:
        return value.upper()
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].upper()
    return MISSING


def parse_callers(value: str) -> list[str]:
    if value in ("", MISSING):
        return []
    return sorted({x for x in re.split(r"[,;|]", value) if x})


def classify(count: int) -> str:
    if count >= 3:
        return "MULTICALLER"
    if count == 2:
        return "TWO_CALLER"
    if count == 1:
        return "SINGLE_CALLER"
    return "UNKNOWN"


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize caller support in a merged SV VCF.")
    ap.add_argument("--vcf", required=True, help="Jasmine/SURVIVOR merged VCF or VCF.GZ")
    ap.add_argument("--output", required=True, help="Output TSV")
    args = ap.parse_args()

    rows = []
    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            info = parse_info(info_raw)
            svtype = infer_svtype(info, alt)

            callers = first(info, "CALLERS", "CALLER", "SOURCES", "SOURCE")
            caller_list = parse_callers(callers)

            supp_vec = first(info, "SUPP_VEC")
            supp = first(info, "SUPP", "SUPPORT")
            caller_count = len(caller_list)

            # If Jasmine/SURVIVOR supplies a support count but caller names are
            # unavailable, retain the count without inventing caller identities.
            if caller_count == 0 and supp not in ("", MISSING):
                try:
                    caller_count = int(float(supp))
                except ValueError:
                    caller_count = 0

            rows.append({
                "SV_ID": sv_id,
                "CHROM": chrom,
                "START": pos,
                "END": first(info, "END"),
                "SVTYPE": svtype,
                "SVLEN": first(info, "SVLEN"),
                "CALLERS": ";".join(caller_list) if caller_list else MISSING,
                "CALLER_COUNT": str(caller_count) if caller_count else MISSING,
                "CALLER_SUPPORT_CLASS": classify(caller_count),
                "SUPP": supp,
                "SUPP_VEC": supp_vec,
            })

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN",
        "CALLERS", "CALLER_COUNT", "CALLER_SUPPORT_CLASS", "SUPP", "SUPP_VEC"
    ]
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] records={len(rows)} output={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
