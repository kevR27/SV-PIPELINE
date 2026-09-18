#!/usr/bin/env python3
"""Normalize caller-specific SV evidence into one TSV schema.

The output is used by the transparent pre-Jasmine evidence filter.  Read
support is therefore normalized conservatively but caller-aware:

* Sniffles2: INFO/SUPPORT (fallback RE/SUPP)
* cuteSV:    INFO/RE (fallback SUPPORT/SUPP)
* Delly:     INFO/PE + INFO/SR; if absent, FORMAT/DV + FORMAT/RV

The script deliberately fails on missing or duplicated VCF IDs because the
next pipeline step selects PASS records from the original VCF by ID.  Silent
ID ambiguity would otherwise remove or retain the wrong records.
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


def numeric(value: str):
    if value in (None, "", MISSING):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def alt_count(value: str) -> float | None:
    """Return alternate-read support from a scalar or ref,alt count pair.

    Manta stores paired-read and split-read counts in sample FORMAT/PR and
    FORMAT/SR as ``ref,alt``.  Other callers may expose a scalar.  For a pair,
    only the alternate count is variant support; summing both values would
    incorrectly reward reference evidence.
    """
    if value in (None, "", MISSING):
        return None
    parts = str(value).split(",")
    target = parts[-1]
    return numeric(target)


def sum_fields(mapping: dict[str, str], *keys: str) -> str:
    values = [numeric(mapping.get(key, MISSING)) for key in keys]
    values = [v for v in values if v is not None]
    if not values:
        return MISSING
    total = sum(values)
    return str(int(total)) if float(total).is_integer() else str(total)


def infer_svtype(info: dict[str, str], alt: str) -> str:
    svtype = first(info, "SVTYPE")
    if svtype != MISSING:
        return svtype.upper()
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].upper()
    return MISSING


def parse_format(fmt_raw: str, sample_raw: str) -> dict[str, str]:
    out = {
        "CALLER_GT": MISSING,
        "CALLER_GQ": MISSING,
        "CALLER_DP": MISSING,
        "CALLER_DR": MISSING,
        "CALLER_DV": MISSING,
        "CALLER_RR": MISSING,
        "CALLER_RV": MISSING,
        "CALLER_PR": MISSING,
        "CALLER_SR": MISSING,
    }
    if not fmt_raw or fmt_raw == MISSING or not sample_raw or sample_raw == MISSING:
        return out

    keys = fmt_raw.split(":")
    values = sample_raw.split(":")
    fields = dict(zip(keys, values))

    for source, target in (
        ("GT", "CALLER_GT"),
        ("GQ", "CALLER_GQ"),
        ("DP", "CALLER_DP"),
        ("DR", "CALLER_DR"),
        ("DV", "CALLER_DV"),
        ("RR", "CALLER_RR"),
        ("RV", "CALLER_RV"),
        ("PR", "CALLER_PR"),
        ("SR", "CALLER_SR"),
    ):
        value = fields.get(source, MISSING)
        if value not in ("", MISSING):
            out[target] = value

    if out["CALLER_DP"] == MISSING:
        # Delly often provides reference/variant support instead of DP.
        parts = [
            numeric(out["CALLER_DR"]),
            numeric(out["CALLER_DV"]),
            numeric(out["CALLER_RR"]),
            numeric(out["CALLER_RV"]),
        ]
        parts = [v for v in parts if v is not None]
        if parts:
            total = sum(parts)
            out["CALLER_DP"] = str(int(total)) if float(total).is_integer() else str(total)

    return out


def caller_evidence(caller: str, info: dict[str, str], fmt: dict[str, str]) -> dict[str, str]:
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
            "CALLER_IMPRECISE": first(info, "IMPRECISE"),
        }

    if caller_lc == "delly":
        # DELLY declares PE and SR as INFO support fields.  The former parser
        # looked for INFO/SU, which DELLY does not use, causing genuine DELLY
        # calls to be labelled LOW_SUPPORT before Jasmine.
        support = sum_fields(info, "PE", "SR")
        if support == MISSING:
            support = sum_fields(
                {
                    "DV": fmt.get("CALLER_DV", MISSING),
                    "RV": fmt.get("CALLER_RV", MISSING),
                },
                "DV",
                "RV",
            )
        if support == MISSING:
            support = first(info, "SUPP", "SUPPORT")
        return {
            "CALLER_SUPPORT": support,
            "CALLER_PE": first(info, "PE"),
            "CALLER_SR": first(info, "SR"),
            "CALLER_PRECISE": first(info, "PRECISE"),
            "CALLER_IMPRECISE": first(info, "IMPRECISE"),
            "CALLER_STRANDS": first(info, "CT"),
        }

    if caller_lc == "manta":
        pr_alt = alt_count(fmt.get("CALLER_PR", MISSING))
        sr_alt = alt_count(fmt.get("CALLER_SR", MISSING))
        components = [value for value in (pr_alt, sr_alt) if value is not None]
        if components:
            support_value = sum(components)
            support = (
                str(int(support_value))
                if float(support_value).is_integer()
                else str(support_value)
            )
        else:
            # Retain INFO fallbacks for converted or non-standard Manta VCFs.
            support = first(info, "SU", "SUPPORT")
        return {
            "CALLER_SUPPORT": support,
            "CALLER_PR": first(fmt, "CALLER_PR"),
            "CALLER_SR": first(fmt, "CALLER_SR"),
        }

    if caller_lc in {"jasmine", "survivor"}:
        return {
            "CALLER_SUPPORT": first(info, "SUPP", "SUPPORT"),
            "CALLER_SUPPORT_VECTOR": first(info, "SUPP_VEC"),
        }

    return {"CALLER_SUPPORT": first(info, "SUPPORT", "SUPP", "RE", "SU")}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse caller-specific SV evidence into a uniform TSV."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ")
    parser.add_argument("--caller", required=True, help="Caller name")
    parser.add_argument("--output", required=True, help="Output TSV")
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    seen_ids: set[str] = set()

    with open_text(args.vcf) as fh:
        for line_no, line in enumerate(fh, start=1):
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at line {line_no}: {line.rstrip()}")

            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            if sv_id in ("", MISSING):
                raise ValueError(
                    f"{args.caller} VCF contains a record without a usable ID at "
                    f"{chrom}:{pos}. PASS-ID filtering would be unsafe."
                )
            if sv_id in seen_ids:
                raise ValueError(
                    f"{args.caller} VCF contains duplicate ID {sv_id!r}. "
                    "PASS-ID filtering requires unique IDs."
                )
            seen_ids.add(sv_id)

            info = parse_info(info_raw)
            fmt_raw = fields[8] if len(fields) > 8 else MISSING
            sample_raw = fields[9] if len(fields) > 9 else MISSING
            fmt = parse_format(fmt_raw, sample_raw)

            row: dict[str, str] = {
                "CALLER": args.caller,
                "SV_ID": sv_id,
                "CHROM": chrom,
                "START": pos,
                "END": first(info, "END"),
                "CHR2": first(info, "CHR2"),
                "POS2": first(info, "POS2"),
                "SVTYPE": infer_svtype(info, alt),
                "SVLEN": first(info, "SVLEN"),
                "QUAL": qual,
                "FILTER": filt,
                "REF": ref,
                "ALT": alt,
                "INFO_RAW": info_raw,
            }
            row.update(fmt)
            row.update(caller_evidence(args.caller, info, fmt))
            rows.append(row)

    columns = [
        "CALLER", "SV_ID", "CHROM", "START", "END", "CHR2", "POS2",
        "SVTYPE", "SVLEN", "QUAL", "FILTER", "REF", "ALT", "INFO_RAW",
        "CALLER_SUPPORT", "CALLER_RNAMES", "CALLER_STRANDS",
        "CALLER_IMPRECISE", "CALLER_MOSAIC", "CALLER_PE", "CALLER_SR",
        "CALLER_PR", "CALLER_PRECISE", "CALLER_SUPPORT_VECTOR",
        "CALLER_GT", "CALLER_GQ", "CALLER_DP", "CALLER_DR", "CALLER_DV",
        "CALLER_RR", "CALLER_RV",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] caller={args.caller} records={len(rows)} output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
