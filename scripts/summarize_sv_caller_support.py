#!/usr/bin/env python3
"""Summarize caller provenance from a per-patient Jasmine/SURVIVOR merge.

The merge input order is fixed and biologically meaningful. SUPP_VEC therefore
provides a stable caller-presence vector even when the merged record does not
carry caller names explicitly (LRS: Sniffles2,cuteSV,Delly; SRS: Manta,Delly).

This script derives CALLERS and CALLER_COUNT from SUPP_VEC first, falling back
to CALLERS/SOURCES or SUPP only when necessary.  The optional high-confidence
VCF is a *companion evidence set*; it never replaces the complete master VCF.
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
    result: dict[str, str] = {}
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
    return sorted({x.strip() for x in re.split(r"[,;|]", value) if x.strip()})


def callers_from_supp_vec(supp_vec: str, caller_order: list[str]) -> list[str]:
    if supp_vec in ("", MISSING):
        return []
    vec = str(supp_vec).strip()
    if len(vec) != len(caller_order) or any(bit not in "01" for bit in vec):
        return []
    return [caller for caller, bit in zip(caller_order, vec) if bit == "1"]


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
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--caller-order",
        default="Sniffles2,cuteSV,Delly",
        help="Comma-separated caller order used to create the merge file list.",
    )
    ap.add_argument(
        "--min-callers",
        type=int,
        default=None,
        help="Also create a companion VCF with at least this many callers.",
    )
    ap.add_argument("--high-confidence-vcf", default=None)
    args = ap.parse_args()

    if args.min_callers is not None and not args.high_confidence_vcf:
        ap.error("--high-confidence-vcf is required when --min-callers is set")

    caller_order = [x.strip() for x in args.caller_order.split(",") if x.strip()]
    if not caller_order:
        ap.error("--caller-order must contain at least one caller")

    rows: list[dict[str, str]] = []
    header_lines: list[str] = []
    records_for_hc: list[tuple[str, int]] = []

    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                if args.min_callers is not None:
                    header_lines.append(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            info = parse_info(info_raw)
            svtype = infer_svtype(info, alt)

            supp_vec = first(info, "SUPP_VEC")
            supp = first(info, "SUPP", "SUPPORT")

            caller_list = callers_from_supp_vec(supp_vec, caller_order)
            provenance_source = "SUPP_VEC" if caller_list else MISSING

            if not caller_list:
                explicit = first(info, "CALLERS", "CALLER", "SOURCES", "SOURCE")
                caller_list = parse_callers(explicit)
                if caller_list:
                    provenance_source = "CALLERS_FIELD"

            caller_count = len(caller_list)
            if caller_count == 0 and supp_vec not in ("", MISSING):
                # Even if the vector length no longer matches the configured
                # caller order, its number of 1s still gives a support count.
                if all(bit in "01" for bit in supp_vec):
                    caller_count = supp_vec.count("1")
                    provenance_source = "SUPP_VEC_COUNT_ONLY"

            if caller_count == 0 and supp not in ("", MISSING):
                try:
                    caller_count = int(float(supp))
                    provenance_source = "SUPP_COUNT_ONLY"
                except ValueError:
                    caller_count = 0

            rows.append(
                {
                    "SV_ID": sv_id,
                    "CHROM": chrom,
                    "START": pos,
                    "END": first(info, "END"),
                    "SVTYPE": svtype,
                    "SVLEN": first(info, "SVLEN"),
                    "CALLERS": ";".join(caller_list) if caller_list else MISSING,
                    "CALLER_COUNT": str(caller_count) if caller_count else MISSING,
                    "CALLER_SUPPORT_CLASS": classify(caller_count),
                    "CALLER_PROVENANCE_SOURCE": provenance_source,
                    "SUPP": supp,
                    "SUPP_VEC": supp_vec,
                }
            )

            if args.min_callers is not None:
                records_for_hc.append((line, caller_count))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "SV_ID",
        "CHROM",
        "START",
        "END",
        "SVTYPE",
        "SVLEN",
        "CALLERS",
        "CALLER_COUNT",
        "CALLER_SUPPORT_CLASS",
        "CALLER_PROVENANCE_SOURCE",
        "SUPP",
        "SUPP_VEC",
    ]
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] records={len(rows)} output={out}")

    if args.min_callers is not None:
        hc_out = Path(args.high_confidence_vcf)
        hc_out.parent.mkdir(parents=True, exist_ok=True)
        n_kept = 0
        with hc_out.open("w", encoding="utf-8", newline="") as fh:
            fh.writelines(header_lines)
            for line, caller_count in records_for_hc:
                if caller_count >= args.min_callers:
                    fh.write(line)
                    n_kept += 1
        print(
            f"[OK] high_confidence min_callers={args.min_callers} "
            f"kept={n_kept}/{len(rows)} output={hc_out}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
