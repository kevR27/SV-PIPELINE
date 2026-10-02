#!/usr/bin/env python3
"""Collapse reciprocal GRIDSS breakend records into one descriptive event table.

GRIDSS is retained as breakpoint/local-assembly evidence. This script does not
reinterpret its BND output as DEL/DUP/INV calls and does not apply a new
pathogenicity or confidence classifier.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


MISSING = "."


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else open(path, encoding="utf-8", errors="replace")


def parse_info(raw: str):
    out = {}
    for item in str(raw).split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = "True"
    return out


def parse_format(fmt_raw: str, sample_raw: str):
    if not fmt_raw or fmt_raw == MISSING or not sample_raw or sample_raw == MISSING:
        return {}
    return dict(zip(fmt_raw.split(":"), sample_raw.split(":")))


def first(info, fmt, *keys):
    for key in keys:
        for source in (fmt, info):
            value = source.get(key, MISSING)
            if value not in ("", MISSING, None):
                return str(value)
    return MISSING


def remote_breakend(alt: str):
    match = re.search(r"([\[\]])([^\[\]]+):(\d+)[\[\]]", alt)
    if not match:
        return MISSING, MISSING, MISSING
    chrom2, pos2 = match.group(2), match.group(3)
    local = "-" if alt.startswith(("[", "]")) else "+"
    remote = "+" if match.group(1) == "]" else "-"
    return chrom2, pos2, local + remote


def numeric_max(values):
    parsed = []
    for value in values:
        try:
            parsed.append(float(value))
        except (TypeError, ValueError):
            pass
    if not parsed:
        return MISSING
    value = max(parsed)
    return str(int(value)) if value.is_integer() else str(value)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    records = {}
    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, record_id, ref, alt, qual, filt, info_raw = fields[:8]
            info = parse_info(info_raw)
            fmt = parse_format(fields[8], fields[9]) if len(fields) >= 10 else {}
            chrom2, pos2, orientation = remote_breakend(alt)
            if chrom2 == MISSING:
                chrom2 = info.get("CHR2", MISSING)
                pos2 = info.get("END", info.get("POS2", MISSING))
            records[record_id] = {
                "id": record_id,
                "mate": info.get("MATEID", MISSING).split(",")[0],
                "event": info.get("EVENT", MISSING),
                "chrom1": chrom,
                "pos1": pos,
                "chrom2": chrom2,
                "pos2": pos2,
                "orientation": orientation,
                "qual": qual,
                "filter": filt,
                "AS": first(info, fmt, "AS"),
                "RAS": first(info, fmt, "RAS"),
                "CAS": first(info, fmt, "CAS"),
                "ASSR": first(info, fmt, "ASSR"),
                "ASRP": first(info, fmt, "ASRP"),
                "SR": first(info, fmt, "SR"),
                "RP": first(info, fmt, "RP"),
                "VF": first(info, fmt, "VF"),
            }

    used = set()
    rows = []
    for record_id, record in records.items():
        if record_id in used:
            continue
        mate = records.get(record["mate"])
        group = [record]
        if mate is not None and mate["id"] not in used:
            group.append(mate)
        used.update(x["id"] for x in group)

        # Each GRIDSS breakend record fully specifies both endpoints. Prefer the
        # current record; the reciprocal record is used to aggregate evidence.
        event_id = record["event"] if record["event"] != MISSING else "|".join(sorted(x["id"] for x in group))
        filters = sorted({x["filter"] for x in group if x["filter"] not in ("", MISSING)})
        row = {
            "GRIDSS_EVENT_ID": event_id,
            "GRIDSS_RECORD_IDS": ";".join(sorted(x["id"] for x in group)),
            "CHROM1": record["chrom1"],
            "POS1": record["pos1"],
            "CHROM2": record["chrom2"],
            "POS2": record["pos2"],
            "ORIENTATION": record["orientation"],
            "QUAL_MAX": numeric_max([x["qual"] for x in group]),
            "FILTERS": ";".join(filters) if filters else MISSING,
        }
        for key in ("AS","RAS","CAS","ASSR","ASRP","SR","RP","VF"):
            row[key] = numeric_max([x[key] for x in group])
        rows.append(row)

    columns = [
        "GRIDSS_EVENT_ID","GRIDSS_RECORD_IDS","CHROM1","POS1","CHROM2","POS2",
        "ORIENTATION","QUAL_MAX","FILTERS","AS","RAS","CAS","ASSR","ASRP","SR","RP","VF",
    ]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] GRIDSS events={len(rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
