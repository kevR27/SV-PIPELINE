#!/usr/bin/env python3
"""
parse_sv_caller_vcf.py

Create a compact, caller-aware SV evidence table from a caller's VCF
(or bgzipped VCF.GZ) file.

Different SV callers use different INFO/FORMAT field names for the same
underlying evidence (e.g. supporting read count is SUPPORT in Sniffles2,
RE in cuteSV, SU in Delly). This script normalizes those into a single
set of CALLER_* columns so downstream steps (filtering, merging, review)
can work with one consistent schema regardless of which caller produced
the record.

Output: one row per VCF record, written as a tab-separated file.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

# Placeholder used everywhere a value is absent, to keep the output
# table free of blank cells.
MISSING = "."


def open_text(path: str):
    """Open a plain-text or gzip-compressed VCF for reading, as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(raw: str) -> dict:
    """
    Parse a VCF INFO field string (e.g. "SVTYPE=DEL;SVLEN=-500;IMPRECISE")
    into a dict. Flag-style fields with no "=" are recorded as 'True'.
    """
    info = {} # empty dict
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

    return info #paste the found character to dict


def first(info: dict, *keys: str) -> str:
    """
    Return the value of the first key (in order) that exists in `info`
    and is non-empty; otherwise MISSING. Lets us handle the fact that
    different callers name the same evidence field differently.
    """
    for key in keys:
        value = info.get(key, MISSING)
        if value not in ("", MISSING):
            return value

    return MISSING


def infer_svtype(info: dict, alt: str) -> str:
    """
    Determine the SV type. Prefer the INFO/SVTYPE field; fall back to
    parsing it out of the ALT allele for breakends (BND) or symbolic alleles
    that lack an explicit SVTYPE
    """
    svtype = first(info, "SVTYPE")
    if svtype != MISSING:
        return svtype

    if "[" in alt or "]" in alt:
        return "BND"

    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1]

    return MISSING


def evidence(caller: str, info: dict) -> dict:
    """
    Extract caller-specific supporting evidence from the INFO dict and
    normalize it into CALLER_* keys. Each caller reports read support
    (and other evidence) under different field names, so this is a
    per-caller lookup table.
    """
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
        return {
            "CALLER_SUPPORT": first(info, "SU", "SUPPORT"),
            "CALLER_PE": first(info, "PE"),
            "CALLER_SR": first(info, "SR"),
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

    # Unknown caller: try the common field names as a best-effort fallback.
    return {"CALLER_SUPPORT": first(info, "SUPPORT", "SUPP", "RE", "SU")}


def parse_format(fmt_raw: str, sample_raw: str) -> dict:
    """
    Pull GT/GQ/DP out of the first sample column, given the FORMAT and
    sample value strings (e.g. FORMAT="GT:GQ:DP", sample="0/1:40:20").
    Returns MISSING for any field not present in this record's FORMAT.
    """
    out = {"CALLER_GT": MISSING, "CALLER_GQ": MISSING, "CALLER_DP": MISSING}

    if not fmt_raw or not sample_raw:
        return out

    keys = fmt_raw.split(":")
    values = sample_raw.split(":")
    fields = dict(zip(keys, values))

    if "GT" in fields:
        out["CALLER_GT"] = fields["GT"]

    if "GQ" in fields and fields["GQ"] not in ("", MISSING):
        out["CALLER_GQ"] = fields["GQ"]

    # Depth is called DP by most callers, DR (reference-supporting depth)
    # by some; take whichever is present, preferring DP.
    for depth_key in ("DP", "DR"):
        if depth_key in fields and fields[depth_key] not in ("", MISSING):
            out["CALLER_DP"] = fields[depth_key]
            break

    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse caller-specific SV evidence into a uniform TSV."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ file")
    parser.add_argument("--caller", required=True, help="Caller name (e.g. Sniffles2, cuteSV, delly)")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    rows = []

    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line.rstrip()}")

            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
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
        "CALLER_GT", "CALLER_GQ", "CALLER_DP",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t",
            extrasaction="ignore", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] caller={args.caller} records={len(rows)} output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
