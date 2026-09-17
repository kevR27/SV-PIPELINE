#!/usr/bin/env python3
"""Summarize caller provenance for Jasmine-merged SVs.

For this LRS workflow the Jasmine input list is ordered as
Sniffles2,cuteSV,delly. Jasmine writes SUPP_VEC and IDLIST in that same input
order, so caller names can be recovered without inventing them.

The optional high-confidence VCF remains a companion view only; the complete
Jasmine VCF stays the master callset for downstream annotation.
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


def classify(count: int) -> str:
    if count >= 3:
        return "MULTICALLER"
    if count == 2:
        return "TWO_CALLER"
    if count == 1:
        return "SINGLE_CALLER"
    return "UNKNOWN"


def decode_provenance(info: dict[str, str], caller_order: list[str]):
    supp_vec = first(info, "SUPP_VEC")
    id_list = first(info, "IDLIST")
    if supp_vec in ("", MISSING):
        return [], {}, "NO_SUPP_VEC"

    callers = [
        caller_order[i]
        for i, bit in enumerate(supp_vec)
        if bit == "1" and i < len(caller_order)
    ]

    ids = [] if id_list in ("", MISSING) else [x for x in id_list.split(",") if x]
    supporting_idx = [i for i, bit in enumerate(supp_vec) if bit == "1"]
    provenance: dict[str, str] = {}
    status = "DECODED"

    if ids and len(ids) == len(supporting_idx):
        for idx, sv_id in zip(supporting_idx, ids):
            if idx < len(caller_order):
                provenance[caller_order[idx]] = sv_id
    elif ids:
        status = "IDLIST_LENGTH_MISMATCH"

    return callers, provenance, status


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize caller support in a Jasmine merged VCF.")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--caller-order",
        default="Sniffles2,cuteSV,delly",
        help="Comma-separated Jasmine file-list order used to decode SUPP_VEC/IDLIST",
    )
    ap.add_argument(
        "--min-callers",
        type=int,
        default=None,
        help="Write a companion VCF containing records supported by at least this many caller inputs",
    )
    ap.add_argument("--high-confidence-vcf", default=None)
    args = ap.parse_args()

    if args.min_callers is not None and not args.high_confidence_vcf:
        ap.error("--high-confidence-vcf is required when --min-callers is set")

    caller_order = [x.strip() for x in args.caller_order.split(",") if x.strip()]
    rows = []
    header_lines: list[str] = []
    kept_lines: list[tuple[str, int]] = []

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
            callers, provenance, decode_status = decode_provenance(info, caller_order)

            supp = first(info, "SUPP", "SUPPORT")
            caller_count = len(callers)
            if caller_count == 0 and supp not in ("", MISSING):
                try:
                    caller_count = int(float(supp))
                except ValueError:
                    caller_count = 0

            row = {
                "SV_ID": sv_id,
                "CHROM": chrom,
                "START": pos,
                "END": first(info, "END"),
                "SVTYPE": infer_svtype(info, alt),
                "SVLEN": first(info, "SVLEN"),
                "CALLERS": ";".join(callers) if callers else MISSING,
                "CALLER_COUNT": str(caller_count) if caller_count else MISSING,
                "CALLER_SUPPORT_CLASS": classify(caller_count),
                "SUPP": supp,
                "SUPP_VEC": first(info, "SUPP_VEC"),
                "IDLIST": first(info, "IDLIST"),
                "PROVENANCE_STATUS": decode_status,
            }
            for caller in caller_order:
                row[f"{caller.upper()}_ID"] = provenance.get(caller, MISSING)
            rows.append(row)

            if args.min_callers is not None:
                kept_lines.append((line, caller_count))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN",
        "CALLERS", "CALLER_COUNT", "CALLER_SUPPORT_CLASS", "SUPP", "SUPP_VEC",
        "IDLIST", "PROVENANCE_STATUS",
    ] + [f"{caller.upper()}_ID" for caller in caller_order]

    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] records={len(rows)} output={out}")

    if args.min_callers is not None:
        hc = Path(args.high_confidence_vcf)
        hc.parent.mkdir(parents=True, exist_ok=True)
        n_kept = 0
        with hc.open("w", encoding="utf-8", newline="") as fh:
            fh.writelines(header_lines)
            for line, caller_count in kept_lines:
                if caller_count >= args.min_callers:
                    fh.write(line)
                    n_kept += 1
        print(
            f"[OK] high_confidence min_callers={args.min_callers} "
            f"kept={n_kept}/{len(rows)} output={hc}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
