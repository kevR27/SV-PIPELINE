#!/usr/bin/env python3
"""Create a caller-aware SV evidence table from VCF/VCF.GZ.

The output normalizes evidence used by the pre-Jasmine QC filter while retaining
caller-specific fields needed for later review.  SV IDs are required to be
present and unique because the filtered PASS-ID list is used to select records
from the original caller VCF and Jasmine later records those IDs in IDLIST.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

MISSING = "."


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(raw: str) -> dict[str, str]:
    info: dict[str, str] = {}
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
    for key in keys:
        value = info.get(key, MISSING)
        if value not in ("", MISSING):
            return value
    return MISSING


def infer_svtype(info: dict[str, str], alt: str) -> str:
    svtype = first(info, "SVTYPE")
    if svtype != MISSING:
        return svtype
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1]
    return MISSING


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def summed_support(*values: str) -> str:
    parsed = [number(v) for v in values if v not in ("", MISSING, None)]
    parsed = [v for v in parsed if v is not None]
    if not parsed:
        return MISSING
    total = sum(parsed)
    return str(int(total)) if total.is_integer() else str(total)


def evidence(caller: str, info: dict[str, str]) -> dict[str, str]:
    caller_lc = caller.lower()

    if caller_lc == "sniffles2":
        return {
            "CALLER_SUPPORT": first(info, "SUPPORT", "RE", "SUPP"),
            "CALLER_RNAMES": first(info, "RNAMES"),
            "CALLER_STRANDS": first(info, "STRANDS"),
            "CALLER_IMPRECISE": first(info, "IMPRECISE"),
            "CALLER_MOSAIC": first(info, "MOSAIC"),
        }

    if caller_lc == "cutesv":
        return {
            "CALLER_SUPPORT": first(info, "RE", "SUPPORT", "SUPP"),
            "CALLER_RNAMES": first(info, "RNAMES"),
            "CALLER_STRANDS": first(info, "STRANDS"),
        }

    if caller_lc == "delly":
        pe = first(info, "PE")
        sr = first(info, "SR")
        support = first(info, "SU", "SUPPORT")
        if support == MISSING:
            # DELLY long-read VCFs commonly expose PE/SR rather than SU.
            # Treat their sum as the record-level supporting evidence used by
            # the caller-agnostic pre-merge filter.  PE and SR are retained as
            # separate columns so the derivation stays transparent.
            support = summed_support(pe, sr)
        return {
            "CALLER_SUPPORT": support,
            "CALLER_PE": pe,
            "CALLER_SR": sr,
            "CALLER_PRECISE": first(info, "PRECISE"),
        }

    if caller_lc == "manta":
        return {
            "CALLER_SUPPORT": first(info, "SU", "PR", "SR", "SUPPORT"),
            "CALLER_PR": first(info, "PR"),
            "CALLER_SR": first(info, "SR"),
        }

    if caller_lc in {"jasmine", "survivor"}:
        return {
            "CALLER_SUPPORT": first(info, "SUPP", "SUPPORT"),
            "CALLER_SUPPORT_VECTOR": first(info, "SUPP_VEC"),
        }

    return {"CALLER_SUPPORT": first(info, "SUPPORT", "SUPP", "RE", "SU")}


def parse_format(fmt_raw: str, sample_raw: str) -> dict[str, str]:
    out = {
        "CALLER_GT": MISSING,
        "CALLER_GQ": MISSING,
        "CALLER_DP": MISSING,
        "CALLER_DV": MISSING,
        "CALLER_RV": MISSING,
    }

    if not fmt_raw or fmt_raw == MISSING or not sample_raw or sample_raw == MISSING:
        return out

    keys = fmt_raw.split(":")
    values = sample_raw.split(":")
    fields = dict(zip(keys, values))

    for key, out_key in (
        ("GT", "CALLER_GT"),
        ("GQ", "CALLER_GQ"),
        ("DV", "CALLER_DV"),
        ("RV", "CALLER_RV"),
    ):
        value = fields.get(key, MISSING)
        if value not in ("", MISSING):
            out[out_key] = value

    for depth_key in ("DP", "DR"):
        value = fields.get(depth_key, MISSING)
        if value not in ("", MISSING):
            out["CALLER_DP"] = value
            break

    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Parse caller-specific SV evidence into a uniform TSV."
    )
    ap.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ")
    ap.add_argument("--caller", required=True, help="Caller name")
    ap.add_argument("--output", required=True, help="Output TSV")
    args = ap.parse_args()

    rows: list[dict[str, str]] = []
    seen_ids: set[str] = set()

    with open_text(args.vcf) as fh:
        for lineno, line in enumerate(fh, start=1):
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at line {lineno}: {line.rstrip()}")

            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]

            if sv_id in {"", MISSING}:
                raise ValueError(
                    f"{args.caller} VCF contains a record without an ID at {chrom}:{pos}. "
                    "Stable unique IDs are required for PASS-ID filtering and Jasmine provenance."
                )
            if sv_id in seen_ids:
                raise ValueError(
                    f"{args.caller} VCF contains duplicate ID {sv_id!r} at {chrom}:{pos}. "
                    "Normalize caller IDs before running the evidence filter/Jasmine."
                )
            seen_ids.add(sv_id)

            info = parse_info(info_raw)
            fmt_raw = fields[8] if len(fields) > 8 else MISSING
            sample_raw = fields[9] if len(fields) > 9 else MISSING

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
            row.update(parse_format(fmt_raw, sample_raw))
            rows.append(row)

    columns = [
        "CALLER", "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN",
        "QUAL", "FILTER", "REF", "ALT", "INFO_RAW",
        "CALLER_SUPPORT", "CALLER_RNAMES", "CALLER_STRANDS",
        "CALLER_IMPRECISE", "CALLER_MOSAIC", "CALLER_PE", "CALLER_SR",
        "CALLER_PR", "CALLER_PRECISE", "CALLER_SUPPORT_VECTOR",
        "CALLER_GT", "CALLER_GQ", "CALLER_DP", "CALLER_DV", "CALLER_RV",
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] caller={args.caller} records={len(rows)} unique_ids={len(seen_ids)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
