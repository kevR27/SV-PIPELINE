#!/usr/bin/env python3
"""Attach MELT and ExpansionHunter locus evidence to the SRS master table.

The row-defining Manta/DELLY/SURVIVOR SV universe is never changed. MELT is
matched only to insertion/BND-like master calls by breakpoint proximity.
ExpansionHunter records are reported as locus overlaps, not as extra SV-caller
support, because its catalog-driven repeat genotypes are a distinct variant
class.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

MISSING = "."


def open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(raw: str) -> dict[str, str]:
    result: dict[str, str] = {}
    if raw in ("", MISSING):
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


def to_int(value):
    if value in (None, "", MISSING):
        return None
    try:
        return int(float(str(value).replace(",", "")))
    except (TypeError, ValueError):
        return None


def normalize_chrom(value: str) -> str:
    text = str(value)
    return text if text.startswith("chr") else "chr" + text


def normalize_svtype(value: str) -> str:
    text = str(value or MISSING).upper()
    if "BND" in text or text == "TRA":
        return "BND"
    for svtype in ("DEL", "DUP", "INS", "INV", "CNV"):
        if svtype in text:
            return svtype
    return text


def parse_sample(format_raw: str, sample_raw: str) -> dict[str, str]:
    if format_raw in ("", MISSING) or sample_raw in ("", MISSING):
        return {}
    return dict(zip(format_raw.split(":"), sample_raw.split(":")))


def vcf_records(path: str):
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            sample = parse_sample(
                fields[8] if len(fields) > 8 else MISSING,
                fields[9] if len(fields) > 9 else MISSING,
            )
            yield {
                "chrom": normalize_chrom(fields[0]),
                "pos": to_int(fields[1]),
                "id": fields[2],
                "alt": fields[4],
                "filter": fields[6],
                "info": info,
                "sample": sample,
            }


def load_melt(path: str | None) -> dict[str, list[dict]]:
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    if not path:
        return by_chrom
    for index, record in enumerate(vcf_records(path), start=1):
        info = record["info"]
        family = info.get("SVTYPE") or info.get("MEINFO") or "MEI"
        family = str(family).split(",", 1)[0]
        ident = record["id"]
        if ident in ("", MISSING):
            ident = f"{family}:{record['chrom']}:{record['pos']}:{index}"
        record["label"] = ident
        record["family"] = family
        by_chrom[record["chrom"]].append(record)
    return by_chrom


def load_expansionhunter(path: str) -> dict[str, list[dict]]:
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for index, record in enumerate(vcf_records(path), start=1):
        info = record["info"]
        sample = record["sample"]
        end = to_int(info.get("END")) or record["pos"]
        repeat_id = (
            info.get("REPID")
            or info.get("VARID")
            or record["id"]
            or f"EH_{index}"
        )
        if repeat_id in ("", MISSING):
            repeat_id = f"EH_{index}"
        details = [str(repeat_id)]
        if sample.get("REPCN") not in (None, "", MISSING):
            details.append(f"REPCN={sample['REPCN']}")
        if sample.get("REPCI") not in (None, "", MISSING):
            details.append(f"REPCI={sample['REPCI']}")
        record["end"] = end
        record["label"] = "|".join(details)
        by_chrom[record["chrom"]].append(record)
    return by_chrom


def interval_overlap(a_start, a_end, b_start, b_end) -> int:
    if None in (a_start, a_end, b_start, b_end):
        return 0
    alo, ahi = sorted((a_start, a_end))
    blo, bhi = sorted((b_start, b_end))
    return max(0, min(ahi, bhi) - max(alo, blo) + 1)


def find_column(columns: list[str], *names: str) -> str | None:
    direct = {name: name for name in columns}
    lower = {name.lower(): name for name in columns}
    for name in names:
        if name in direct:
            return direct[name]
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Attach orthogonal SRS evidence")
    parser.add_argument("--integrated", required=True)
    parser.add_argument("--expansionhunter", required=True)
    parser.add_argument("--melt")
    parser.add_argument("--breakpoint-tolerance", type=int, default=500)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open_text(args.integrated) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        input_columns = list(reader.fieldnames or [])

    chrom_col = find_column(input_columns, "CHROM", "Chr", "chrom")
    start_col = find_column(input_columns, "START", "POS", "SV_start")
    end_col = find_column(input_columns, "END", "SV_end")
    type_col = find_column(input_columns, "SVTYPE", "SV_type", "Type")
    if not all((chrom_col, start_col, type_col)):
        raise ValueError("Integrated table requires CHROM, START/POS and SVTYPE")

    melt = load_melt(args.melt)
    eh = load_expansionhunter(args.expansionhunter)

    for row in rows:
        chrom = normalize_chrom(row[chrom_col])
        start = to_int(row[start_col])
        end = to_int(row[end_col]) if end_col else start
        if end is None:
            end = start
        svtype = normalize_svtype(row[type_col])

        melt_matches: list[str] = []
        if args.melt and start is not None and svtype in {"INS", "BND"}:
            for record in melt.get(chrom, []):
                if record["pos"] is None:
                    continue
                if abs(start - record["pos"]) <= args.breakpoint_tolerance:
                    melt_matches.append(
                        f"{record['label']}|{record['family']}|{record['filter']}"
                    )
        row["MELT_MATCH"] = (
            "YES" if melt_matches else ("NO" if args.melt else "NOT_RUN")
        )
        row["MELT_INSERTIONS"] = (
            ";".join(sorted(set(melt_matches))) if melt_matches else MISSING
        )

        eh_matches: list[str] = []
        if start is not None:
            for record in eh.get(chrom, []):
                pos = record["pos"]
                locus_end = record["end"]
                if pos is None or locus_end is None:
                    continue
                overlap = interval_overlap(start, end, pos, locus_end)
                near_breakpoint = min(abs(start - pos), abs(start - locus_end))
                if overlap > 0 or (
                    svtype in {"INS", "BND"}
                    and near_breakpoint <= args.breakpoint_tolerance
                ):
                    eh_matches.append(record["label"])
        row["EXPANSIONHUNTER_LOCUS_OVERLAP"] = "YES" if eh_matches else "NO"
        row["EXPANSIONHUNTER_LOCI"] = (
            ";".join(sorted(set(eh_matches))) if eh_matches else MISSING
        )

    output_columns = input_columns + [
        "MELT_MATCH",
        "MELT_INSERTIONS",
        "EXPANSIONHUNTER_LOCUS_OVERLAP",
        "EXPANSIONHUNTER_LOCI",
    ]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    melt_count = sum(row["MELT_MATCH"] == "YES" for row in rows)
    eh_count = sum(
        row["EXPANSIONHUNTER_LOCUS_OVERLAP"] == "YES" for row in rows
    )
    print(
        f"[OK] rows={len(rows)} melt_matches={melt_count} "
        f"expansionhunter_overlaps={eh_count} output={output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
