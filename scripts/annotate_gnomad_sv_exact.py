#!/usr/bin/env python3
"""Annotate SV-gene events with conservative gnomAD-SV v4.1 site matches.

The script reports gnomAD AF/AC/AN only for exact coordinate/type matches.
For DEL/DUP/INV, POS and END must match exactly. For INS, POS and SV length
must match when both lengths are available. Sequence identity is not inferred
for symbolic insertions. BND/CTX exact matching is intentionally not attempted.

This is gnomAD-specific site-frequency evidence and is kept separate from
AnnotSV benign-region overlap AF.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import math
import pandas as pd
import pysam


MISSING = {None, "", ".", "NA", "N/A", "nan", "None"}


def number(value):
    try:
        if value in MISSING:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def first_number(value):
    if value in MISSING:
        return None
    if isinstance(value, (tuple, list)):
        if not value:
            return None
        return number(value[0])
    return number(value)


def normalize_chrom(value):
    text = str(value)
    return text if text.startswith("chr") else "chr" + text


def record_svtype(record):
    value = record.info.get("SVTYPE")
    if value not in MISSING:
        return str(value).upper()
    alt = record.alts[0] if record.alts else ""
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].split(":")[0].upper()
    return "."


def record_svlen(record):
    value = record.info.get("SVLEN")
    if isinstance(value, (tuple, list)):
        value = value[0] if value else None
    value = number(value)
    if value is not None:
        return abs(int(value))
    if record_svtype(record) in {"DEL", "DUP", "INV"}:
        return abs(int(record.stop) - int(record.pos))
    return None


def info_value(record, key):
    try:
        return record.info.get(key)
    except Exception:
        return None


def af_value(record):
    # gnomAD SV releases provide site-level frequency annotations. AF is the
    # preferred global field; fallback to AC/AN only when both are available.
    af = first_number(info_value(record, "AF"))
    ac = first_number(info_value(record, "AC"))
    an = first_number(info_value(record, "AN"))
    if af is None and ac is not None and an not in (None, 0):
        af = ac / an
    return af, ac, an


def filters(record):
    values = list(record.filter.keys())
    return ";".join(values) if values else "PASS"


def exact_match(query, record):
    qtype = str(query["SVTYPE"]).upper()
    rtype = record_svtype(record)
    if qtype != rtype:
        return False, "TYPE_MISMATCH"

    qstart = int(float(query["START"]))
    if int(record.pos) != qstart:
        return False, "POS_MISMATCH"

    if qtype in {"DEL", "DUP", "INV"}:
        qend = int(float(query["END"]))
        if int(record.stop) != qend:
            return False, "END_MISMATCH"
        return True, "EXACT_POS_END_TYPE"

    if qtype == "INS":
        qlen = number(query.get("SV_SPAN_BP"))
        if qlen is None:
            qlen = number(query.get("SVLEN"))
        if qlen is not None:
            qlen = abs(int(qlen))
        rlen = record_svlen(record)
        if qlen is not None and rlen is not None and qlen != rlen:
            return False, "SVLEN_MISMATCH"
        return True, (
            "EXACT_POS_TYPE_SIZE_NOT_SEQUENCE_CONFIRMED"
            if qlen is not None and rlen is not None
            else "EXACT_POS_TYPE_SIZE_UNAVAILABLE"
        )

    return False, "SVTYPE_NOT_SUPPORTED_FOR_EXACT_MATCH"


def fetch_records(vcf, chrom, pos):
    names = [chrom]
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    if alt not in names:
        names.append(alt)

    records = []
    for name in names:
        try:
            records.extend(list(vcf.fetch(name, max(0, pos - 1), pos)))
        except (ValueError, KeyError):
            continue
    return records


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--gnomad-vcf", default=None)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False)

    added = [
        "GNOMAD_SV_EXACT_MATCH",
        "GNOMAD_SV_ID",
        "GNOMAD_SV_AF",
        "GNOMAD_SV_AC",
        "GNOMAD_SV_AN",
        "GNOMAD_SV_FILTER",
        "GNOMAD_SV_MATCH_SCOPE",
        "GNOMAD_SV_RESOURCE",
    ]
    for col in added:
        df[col] = "."

    if not args.gnomad_vcf:
        df["GNOMAD_SV_EXACT_MATCH"] = "RESOURCE_NOT_CONFIGURED"
        df["GNOMAD_SV_MATCH_SCOPE"] = (
            "NO_GNOMAD_VCF_CONFIGURED;ANNOTSV_OVERLAP_IS_SEPARATE"
        )
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, sep="\t", index=False)
        print(f"[OK] gnomad_resource=not_configured rows={len(df)} output={args.output}")
        return

    vcf_path = Path(args.gnomad_vcf)
    if not vcf_path.exists():
        raise FileNotFoundError(f"Configured gnomAD-SV VCF does not exist: {vcf_path}")

    vcf = pysam.VariantFile(str(vcf_path))
    cache = {}
    results_by_id = {}
    unique_svs = int(df["SV_ID"].nunique())

    for sv_number, (sv_id, group) in enumerate(df.groupby("SV_ID", sort=False), start=1):
        if sv_number == 1 or sv_number % 1000 == 0 or sv_number == unique_svs:
            print(
                f"[INFO] gnomAD-SV exact matching {sv_number}/{unique_svs} unique SVs",
                flush=True,
            )
        row = group.iloc[0]
        qtype = str(row.get("SVTYPE", ".")).upper()

        result = {
            "GNOMAD_SV_EXACT_MATCH": "NO",
            "GNOMAD_SV_ID": ".",
            "GNOMAD_SV_AF": ".",
            "GNOMAD_SV_AC": ".",
            "GNOMAD_SV_AN": ".",
            "GNOMAD_SV_FILTER": ".",
            "GNOMAD_SV_MATCH_SCOPE": (
                "BND_CTX_EXACT_MATCH_NOT_ATTEMPTED"
                if qtype in {"BND", "TRA", "CTX"}
                else "NO_EXACT_COORDINATE_TYPE_MATCH"
            ),
            "GNOMAD_SV_RESOURCE": vcf_path.name,
        }

        if qtype not in {"BND", "TRA", "CTX"}:
            chrom = normalize_chrom(row["CHROM"])
            pos = int(float(row["START"]))
            key = (chrom, pos)
            if key not in cache:
                cache[key] = fetch_records(vcf, chrom, pos)

            matches = []
            for record in cache[key]:
                ok, scope = exact_match(row, record)
                if ok:
                    matches.append((record, scope))

            if matches:
                # Exact coordinate/type duplicates are unusual. If multiple
                # records exist, retain the record with the highest AF but
                # explicitly report the multiplicity.
                ranked = []
                for record, scope in matches:
                    af, ac, an = af_value(record)
                    ranked.append((af if af is not None else -1.0, record, scope, ac, an))
                ranked.sort(key=lambda x: x[0], reverse=True)
                af, record, scope, ac, an = ranked[0]
                if af == -1.0:
                    af = None
                result.update({
                    "GNOMAD_SV_EXACT_MATCH": "YES",
                    "GNOMAD_SV_ID": record.id or ".",
                    "GNOMAD_SV_AF": "." if af is None else f"{af:.8g}",
                    "GNOMAD_SV_AC": "." if ac is None else f"{ac:.8g}",
                    "GNOMAD_SV_AN": "." if an is None else f"{an:.8g}",
                    "GNOMAD_SV_FILTER": filters(record),
                    "GNOMAD_SV_MATCH_SCOPE": (
                        scope
                        if len(matches) == 1
                        else scope + f";MULTIPLE_EXACT_RECORDS={len(matches)}"
                    ),
                })

        results_by_id[str(sv_id)] = result

    vcf.close()

    # Apply one annotation result per master SV in a single vectorized pass.
    # The input may contain multiple SV-gene rows for the same SV; all such
    # rows receive the same gnomAD site annotation without repeatedly scanning
    # the full dataframe for every unique SV.
    sv_ids = df["SV_ID"].astype(str)
    for col in added:
        df[col] = sv_ids.map(
            lambda sv_id: results_by_id.get(sv_id, {}).get(col, ".")
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False)

    exact = int(df.drop_duplicates("SV_ID")["GNOMAD_SV_EXACT_MATCH"].eq("YES").sum())
    print(f"[OK] exact_gnomad_matches={exact} svs={df['SV_ID'].nunique()} output={output}")


if __name__ == "__main__":
    main()
