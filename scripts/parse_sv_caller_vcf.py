#!/usr/bin/env python3
"""Create a caller-aware structural-variant evidence table from a VCF.

The output is intentionally caller-agnostic for the fields used by the
pre-Jasmine evidence filter, while retaining useful caller-specific evidence.
A stable, unique VCF ID is required because filter_sv_evidence.py writes PASS
IDs which are subsequently used by bcftools to subset the original VCF.
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


def first(mapping: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = mapping.get(key, MISSING)
        if value not in ("", MISSING, None):
            return str(value)
    return MISSING


def parse_format(fmt_raw: str, sample_raw: str) -> dict[str, str]:
    if not fmt_raw or fmt_raw == MISSING or not sample_raw or sample_raw == MISSING:
        return {}
    keys = fmt_raw.split(":")
    values = sample_raw.split(":")
    return dict(zip(keys, values))


def infer_svtype(info: dict[str, str], alt: str) -> str:
    svtype = first(info, "SVTYPE")
    if svtype != MISSING:
        return svtype.upper()
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].upper()
    return MISSING


def caller_evidence(caller: str, info: dict[str, str], fmt: dict[str, str]) -> dict[str, str]:
    """Return normalized evidence plus caller-specific fields.

    In particular, DELLY commonly reports alternate-read evidence in FORMAT/DV
    (and sometimes RV) rather than INFO/SU. DR is reference-supporting depth and
    must not be used as alternate support.
    """
    caller_lc = caller.lower()

    base = {
        "CALLER_SUPPORT": MISSING,
        "CALLER_RNAMES": MISSING,
        "CALLER_STRANDS": MISSING,
        "CALLER_IMPRECISE": MISSING,
        "CALLER_MOSAIC": MISSING,
        "CALLER_PE": MISSING,
        "CALLER_SR": MISSING,
        "CALLER_PR": MISSING,
        "CALLER_PRECISE": MISSING,
        "CALLER_SUPPORT_VECTOR": MISSING,
        "CALLER_DV": first(fmt, "DV"),
        "CALLER_DR": first(fmt, "DR"),
        "CALLER_RV": first(fmt, "RV"),
        "CALLER_RR": first(fmt, "RR"),
    }

    if caller_lc == "sniffles2":
        base.update(
            CALLER_SUPPORT=first(info, "SUPPORT", "RE", "SUPP"),
            CALLER_RNAMES=first(info, "RNAMES"),
            CALLER_STRANDS=first(info, "STRANDS"),
            CALLER_IMPRECISE=first(info, "IMPRECISE"),
            CALLER_MOSAIC=first(info, "MOSAIC"),
        )
    elif caller_lc == "cutesv":
        base.update(
            CALLER_SUPPORT=first(info, "RE", "SUPPORT", "SUPP"),
            CALLER_RNAMES=first(info, "RNAMES"),
            CALLER_STRANDS=first(info, "STRANDS"),
        )
    elif caller_lc == "delly":
        # Prefer an explicit total support INFO value when it exists. Otherwise
        # use the alternate-read FORMAT fields emitted by DELLY long-read mode.
        base.update(
            CALLER_SUPPORT=first(info, "SU", "SUPPORT")
            if first(info, "SU", "SUPPORT") != MISSING
            else first(fmt, "DV", "RV"),
            CALLER_PE=first(info, "PE"),
            CALLER_SR=first(info, "SR"),
            CALLER_PRECISE=first(info, "PRECISE"),
        )
    elif caller_lc == "manta":
        base.update(
            CALLER_SUPPORT=first(info, "SU", "SUPPORT"),
            CALLER_PR=first(info, "PR"),
            CALLER_SR=first(info, "SR"),
        )
    elif caller_lc in {"jasmine", "survivor"}:
        base.update(
            CALLER_SUPPORT=first(info, "SUPP", "SUPPORT"),
            CALLER_SUPPORT_VECTOR=first(info, "SUPP_VEC"),
        )
    else:
        base["CALLER_SUPPORT"] = first(info, "SUPPORT", "SUPP", "RE", "SU")

    return base


def main() -> int:
    ap = argparse.ArgumentParser(description="Parse caller-specific SV evidence into a uniform TSV.")
    ap.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ")
    ap.add_argument("--caller", required=True, help="Caller name")
    ap.add_argument("--output", required=True, help="Output TSV")
    args = ap.parse_args()

    rows: list[dict[str, str]] = []
    seen_ids: set[str] = set()

    with open_text(args.vcf) as fh:
        for line_number, line in enumerate(fh, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at line {line_number}: {line.rstrip()}")

            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            if sv_id in ("", MISSING):
                raise ValueError(
                    f"{args.caller}: missing VCF ID at {chrom}:{pos}. "
                    "Stable unique IDs are required for PASS-ID filtering before Jasmine."
                )
            if sv_id in seen_ids:
                raise ValueError(
                    f"{args.caller}: duplicate VCF ID '{sv_id}'. "
                    "IDs must be unique before PASS-ID filtering."
                )
            seen_ids.add(sv_id)

            info = parse_info(info_raw)
            fmt_raw = fields[8] if len(fields) > 8 else MISSING
            sample_raw = fields[9] if len(fields) > 9 else MISSING
            fmt = parse_format(fmt_raw, sample_raw)

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
                "CALLER_GT": first(fmt, "GT"),
                "CALLER_GQ": first(fmt, "GQ"),
                "CALLER_DP": first(fmt, "DP"),
            }
            row.update(caller_evidence(args.caller, info, fmt))
            rows.append(row)

    columns = [
        "CALLER", "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN",
        "QUAL", "FILTER", "REF", "ALT", "INFO_RAW",
        "CALLER_SUPPORT", "CALLER_RNAMES", "CALLER_STRANDS",
        "CALLER_IMPRECISE", "CALLER_MOSAIC", "CALLER_PE", "CALLER_SR",
        "CALLER_PR", "CALLER_PRECISE", "CALLER_SUPPORT_VECTOR",
        "CALLER_GT", "CALLER_GQ", "CALLER_DP",
        "CALLER_DV", "CALLER_DR", "CALLER_RV", "CALLER_RR",
    ]

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t",
            extrasaction="ignore", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] caller={args.caller} records={len(rows)} unique_ids={len(seen_ids)} output={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
