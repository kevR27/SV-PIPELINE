#!/usr/bin/env python3
"""Convert a single-sample normalized mity VCF to a transparent mtDNA TSV."""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


MISSING = "."


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else open(path, encoding="utf-8", errors="replace")


def parse_pairs(raw: str, sep=";"):
    out = {}
    for item in str(raw).split(sep):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        elif item:
            out[item] = "True"
    return out


def pick_alt_value(value: str, index: int):
    if value in (None, "", MISSING):
        return MISSING
    values = str(value).split(",")
    return values[index] if index < len(values) else values[-1]


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt_raw, qual, filt, info_raw = fields[:8]
            info = parse_pairs(info_raw)
            fmt = {}
            if len(fields) >= 10:
                fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))

            alts = alt_raw.split(",")
            for alt_index, alt in enumerate(alts):
                ao = pick_alt_value(fmt.get("AO", info.get("AO", MISSING)), alt_index)
                ro = fmt.get("RO", info.get("RO", MISSING))
                dp = fmt.get("DP", info.get("DP", MISSING))

                # mity normalise explicitly writes FORMAT/VAF = AO / DP.
                # Use that value first; only derive a fallback when VAF is
                # absent so this parser does not redefine mity's heteroplasmy.
                vaf = pick_alt_value(fmt.get("VAF", MISSING), alt_index)
                if number(vaf) is None:
                    ao_n, dp_n = number(ao), number(dp)
                    if ao_n is not None and dp_n is not None and dp_n > 0:
                        vaf = f"{ao_n / dp_n:.8g}"
                    else:
                        vaf = MISSING

                rows.append({
                    "SAMPLE": args.sample,
                    "CHROM": chrom,
                    "POS": pos,
                    "ID": vid,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": filt,
                    "TOTAL_DEPTH": dp,
                    "REF_DEPTH": ro,
                    "ALT_DEPTH": ao,
                    "HETEROPLASMY_VAF": vaf,
                    "MITY_TIER": pick_alt_value(fmt.get("tier", MISSING), alt_index),
                    "MITY_Q": pick_alt_value(fmt.get("q", MISSING), alt_index),
                    "MITY_POS_FILTER": pick_alt_value(fmt.get("POS_filter", MISSING), alt_index),
                    "MITY_SBR_FILTER": pick_alt_value(fmt.get("SBR_filter", MISSING), alt_index),
                    "MITY_SBA_FILTER": pick_alt_value(fmt.get("SBA_filter", MISSING), alt_index),
                    "MITY_MQMR_FILTER": pick_alt_value(fmt.get("MQMR_filter", MISSING), alt_index),
                    "MITY_AQR_FILTER": pick_alt_value(fmt.get("AQR_filter", MISSING), alt_index),
                    "MITY_INFO_RAW": info_raw,
                })

    columns = [
        "SAMPLE","CHROM","POS","ID","REF","ALT","QUAL","FILTER","TOTAL_DEPTH",
        "REF_DEPTH","ALT_DEPTH","HETEROPLASMY_VAF","MITY_TIER","MITY_Q",
        "MITY_POS_FILTER","MITY_SBR_FILTER","MITY_SBA_FILTER",
        "MITY_MQMR_FILTER","MITY_AQR_FILTER","MITY_INFO_RAW",
    ]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] mtDNA alleles={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
