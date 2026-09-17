#!/usr/bin/env python3
"""Summarize caller provenance in a Jasmine/SURVIVOR merged SV VCF.

This script never reclusters variants. It reads the caller-support metadata
already present in the merged record. For Jasmine, SUPP_VEC is decoded using
an explicit caller order matching the VCF list supplied to Jasmine.

The optional high-confidence VCF is a companion view only: the full merged VCF
remains the source of truth for annotation and clinical review.
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
    out: dict[str, str] = {}
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


def first(info: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = info.get(key, MISSING)
        if value not in ("", MISSING, None):
            return str(value)
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


def parse_named_callers(value: str) -> list[str]:
    if value in ("", MISSING):
        return []
    return sorted({x.strip() for x in re.split(r"[,;|]", value) if x.strip()})


def decode_supp_vec(value: str, caller_order: list[str]) -> list[str]:
    """Decode Jasmine SUPP_VEC against the exact file-list caller order."""
    if value in ("", MISSING):
        return []
    bits = value.strip()
    if not caller_order:
        return []
    if not re.fullmatch(r"[01]+", bits):
        return []
    if len(bits) != len(caller_order):
        raise ValueError(
            f"SUPP_VEC length ({len(bits)}) does not match --caller-order "
            f"length ({len(caller_order)}): {bits}"
        )
    return [caller for caller, bit in zip(caller_order, bits) if bit == "1"]


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
        default="",
        help="Comma-separated caller order used to build Jasmine file_list, e.g. Sniffles2,cuteSV,delly",
    )
    ap.add_argument(
        "--min-callers", type=int, default=None,
        help="Write a companion VCF with at least this many distinct callers.",
    )
    ap.add_argument("--high-confidence-vcf", default=None)
    args = ap.parse_args()

    if args.min_callers is not None and not args.high_confidence_vcf:
        ap.error("--high-confidence-vcf is required when --min-callers is set")

    caller_order = [x.strip() for x in args.caller_order.split(",") if x.strip()]
    rows: list[dict[str, str]] = []
    header_lines: list[str] = []
    kept_lines: list[str] = []

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

            named = parse_named_callers(first(info, "CALLERS", "CALLER", "SOURCES", "SOURCE"))
            supp_vec = first(info, "SUPP_VEC")
            vec_callers = decode_supp_vec(supp_vec, caller_order) if caller_order else []
            callers = sorted(set(named) | set(vec_callers))

            supp = first(info, "SUPP", "SUPPORT")
            if callers:
                caller_count = len(callers)
            else:
                try:
                    caller_count = int(float(supp)) if supp != MISSING else 0
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
                    "CALLERS": ";".join(callers) if callers else MISSING,
                    "CALLER_COUNT": str(caller_count) if caller_count else MISSING,
                    "CALLER_SUPPORT_CLASS": classify(caller_count),
                    "SUPP": supp,
                    "SUPP_VEC": supp_vec,
                }
            )

            if args.min_callers is not None and caller_count >= args.min_callers:
                # Preserve every qualifying record in input order. Do not key by
                # SV_ID: paired BND representations can legitimately require more
                # than one record and must never be silently overwritten.
                kept_lines.append(line)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN",
        "CALLERS", "CALLER_COUNT", "CALLER_SUPPORT_CLASS", "SUPP", "SUPP_VEC",
    ]
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] records={len(rows)} output={out}")

    if args.min_callers is not None:
        hc = Path(args.high_confidence_vcf)
        hc.parent.mkdir(parents=True, exist_ok=True)
        with hc.open("w", encoding="utf-8", newline="") as fh:
            fh.writelines(header_lines)
            fh.writelines(kept_lines)
        print(
            f"[OK] high_confidence min_callers={args.min_callers} "
            f"kept={len(kept_lines)}/{len(rows)} output={hc}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
