#!/usr/bin/env python3
"""
Add independent supporting evidence to the final SRS SV/gene table.

This script does not change which SVs exist in the integrated callset.
It only adds evidence columns showing whether a candidate is also supported by:

- CNVpytor read depth;
- GRIDSS breakpoint/assembly evidence;
- MELT mobile-element insertion calls;
- ExpansionHunter repeat-expansion loci.

The final TECHNICAL_SUPPORT column is descriptive only. It is not an ACMG
classification and it does not imply pathogenicity.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

MISSING = "."


def open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(text: str) -> dict[str, str]:
    result: dict[str, str] = {}

    for item in str(text).split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "True"

    return result


def as_int(value):
    try:
        return int(float(str(value).replace(",", "")))
    except (TypeError, ValueError):
        return None


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def normalize_chromosome(value: str) -> str:
    value = str(value)

    if value == "MT":
        return "chrM"

    if value.startswith("chr"):
        return value

    return "chr" + value


def normalize_svtype(value: str) -> str:
    value = str(value or MISSING).upper()

    if "BND" in value or value == "TRA":
        return "BND"

    for svtype in ("DEL", "DUP", "INS", "INV", "CNV"):
        if svtype in value:
            return svtype

    return value


def read_vcf(path: str) -> list[dict]:
    records: list[dict] = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            info = parse_info(fields[7])

            sample = {}
            if len(fields) > 9:
                sample = dict(
                    zip(
                        fields[8].split(":"),
                        fields[9].split(":"),
                    )
                )

            records.append(
                {
                    "chrom": normalize_chromosome(fields[0]),
                    "pos": as_int(fields[1]),
                    "id": fields[2],
                    "alt": fields[4],
                    "qual": fields[5],
                    "filter": fields[6],
                    "info": info,
                    "sample": sample,
                }
            )

    return records


def parse_breakend(record: dict):
    """Extract the two GRIDSS/BND breakpoints when possible."""
    match = re.search(
        r"([\[\]])([^\[\]]+):(\d+)[\[\]]",
        record["alt"],
    )

    if match:
        return (
            record["chrom"],
            record["pos"],
            normalize_chromosome(match.group(2)),
            int(match.group(3)),
        )

    chromosome_2 = record["info"].get("CHR2", MISSING)
    if chromosome_2 != MISSING:
        chromosome_2 = normalize_chromosome(chromosome_2)

    position_2 = as_int(
        record["info"].get(
            "POS2",
            record["info"].get("END"),
        )
    )

    return (
        record["chrom"],
        record["pos"],
        chromosome_2,
        position_2,
    )


def interval_overlap(
    start_a,
    end_a,
    start_b,
    end_b,
) -> int:
    if None in (start_a, end_a, start_b, end_b):
        return 0

    return max(
        0,
        min(max(start_a, end_a), max(start_b, end_b))
        - max(min(start_a, end_a), min(start_b, end_b))
        + 1,
    )


def find_value(row: dict, *names: str) -> str:
    """Find a column while tolerating capitalization differences."""
    lower = {key.lower(): key for key in row}

    for name in names:
        key = name if name in row else lower.get(name.lower())

        if key is None:
            continue

        value = str(row.get(key, "")).strip()
        if value not in {"", MISSING}:
            return value

    return MISSING


def load_integrated_vcf(path: str) -> dict[str, dict]:
    """Index integrated SV records by SV_ID."""
    return {
        record["id"]: record
        for record in read_vcf(path)
    }


def load_gridss(path: str | None, min_qual: float | None = None):
    """Index GRIDSS breakends by first chromosome."""
    index = defaultdict(list)

    if not path:
        return index

    for record in read_vcf(path):
        if record["filter"] not in {"PASS", "."}:
            continue
        qual = as_float(record.get("qual"))
        if min_qual is not None and (qual is None or qual < min_qual):
            continue
        chromosome_1, position_1, chromosome_2, position_2 = parse_breakend(record)

        if position_1 is not None:
            index[chromosome_1].append(
                (
                    position_1,
                    chromosome_2,
                    position_2,
                    record,
                )
            )

    return index


def load_melt(path: str | None):
    """Index MELT insertions by chromosome."""
    index = defaultdict(list)

    if not path:
        return index

    for record in read_vcf(path):
        index[record["chrom"]].append(record)

    return index


def load_expansionhunter(path: str):
    """Index ExpansionHunter loci by chromosome."""
    index = defaultdict(list)

    for record in read_vcf(path):
        record["end"] = (
            as_int(record["info"].get("END"))
            or record["pos"]
        )
        index[record["chrom"]].append(record)

    return index


def gridss_support(
    chromosome: str,
    start: int | None,
    end: int | None,
    svtype: str,
    chromosome_2: str,
    position_2: int | None,
    gridss_index,
    tolerance: int,
):
    """
    Compare candidate breakpoints with GRIDSS breakends.

    For DEL/DUP/INV/BND events with two defined breakpoints, strong GRIDSS
    support requires both breakpoints to match the same GRIDSS record within
    the configured distance. A one-breakpoint match is reported separately
    and is not treated as equivalent to full event support.

    For an event represented by only one local breakpoint, one matching
    GRIDSS breakend is sufficient to report single-breakpoint support.
    """
    if start is None:
        return "NO", []

    if svtype == "BND" and position_2 is not None:
        candidate_1 = (chromosome, start)
        candidate_2 = (chromosome_2, position_2)
    elif end is not None and end != start:
        candidate_1 = (chromosome, start)
        candidate_2 = (chromosome, end)
    else:
        candidate_1 = (chromosome, start)
        candidate_2 = None

    full_matches: list[str] = []
    partial_matches: list[str] = []

    # GRIDSS records are indexed by the chromosome of their first breakend.
    # Search records starting near the first candidate breakpoint.
    for gridss_position_1, gridss_chromosome_2, gridss_position_2, record in gridss_index.get(
        candidate_1[0],
        [],
    ):
        if gridss_position_1 is None:
            continue

        first_matches = (
            abs(gridss_position_1 - candidate_1[1]) <= tolerance
        )

        # A paired GRIDSS breakend can also be represented in the opposite
        # direction, so test both direct and swapped breakpoint order.
        direct_second_matches = False
        swapped_matches = False

        if candidate_2 is not None and gridss_position_2 is not None:
            direct_second_matches = (
                gridss_chromosome_2 == candidate_2[0]
                and abs(gridss_position_2 - candidate_2[1]) <= tolerance
            )

            swapped_matches = (
                gridss_chromosome_2 == candidate_1[0]
                and abs(gridss_position_2 - candidate_1[1]) <= tolerance
                and record["chrom"] == candidate_2[0]
                and record["pos"] is not None
                and abs(record["pos"] - candidate_2[1]) <= tolerance
            )

        record_id = record["id"]
        if record_id in {"", MISSING}:
            record_id = f"{record['chrom']}:{record['pos']}"

        if candidate_2 is None:
            if first_matches:
                full_matches.append(record_id)
            continue

        if (first_matches and direct_second_matches) or swapped_matches:
            full_matches.append(record_id)
        elif first_matches or direct_second_matches:
            partial_matches.append(record_id)

    # Also search the second candidate chromosome because GRIDSS can store the
    # pair with either breakend first.
    if candidate_2 is not None:
        for gridss_position_1, gridss_chromosome_2, gridss_position_2, record in gridss_index.get(
            candidate_2[0],
            [],
        ):
            if gridss_position_1 is None:
                continue

            first_matches_second = (
                abs(gridss_position_1 - candidate_2[1]) <= tolerance
            )
            remote_matches_first = (
                gridss_position_2 is not None
                and gridss_chromosome_2 == candidate_1[0]
                and abs(gridss_position_2 - candidate_1[1]) <= tolerance
            )

            record_id = record["id"]
            if record_id in {"", MISSING}:
                record_id = f"{record['chrom']}:{record['pos']}"

            if first_matches_second and remote_matches_first:
                full_matches.append(record_id)
            elif first_matches_second or remote_matches_first:
                partial_matches.append(record_id)

    full_matches = sorted(set(full_matches))
    partial_matches = sorted(set(partial_matches) - set(full_matches))

    if full_matches:
        if candidate_2 is None:
            return "YES_SINGLE_BREAKPOINT", full_matches
        return "YES_BOTH_BREAKPOINTS", full_matches

    if partial_matches:
        return "PARTIAL_ONE_BREAKPOINT", partial_matches

    return "NO", []

def melt_matches(
    chromosome: str,
    start: int | None,
    svtype: str,
    melt_index,
    tolerance: int,
) -> list[str]:
    """Find MELT insertions close to INS/BND candidates."""
    if start is None or svtype not in {"INS", "BND"}:
        return []

    matches: list[str] = []

    for record in melt_index.get(chromosome, []):
        if record["pos"] is None:
            continue

        if abs(record["pos"] - start) <= tolerance:
            matches.append(record["id"])

    return sorted(set(matches))


def expansionhunter_matches(
    chromosome: str,
    start: int | None,
    end: int | None,
    svtype: str,
    expansionhunter_index,
) -> list[str]:
    """
    Report repeat loci overlapping the candidate interval.

    This is context only. ExpansionHunter is not counted as another general
    SV caller.
    """
    if start is None:
        return []

    local_end = start if end is None or svtype == "BND" else end
    matches: list[str] = []

    for record in expansionhunter_index.get(chromosome, []):
        if interval_overlap(
            start,
            local_end,
            record["pos"],
            record["end"],
        ) == 0:
            continue

        repeat_id = record["info"].get(
            "REPID",
            record["id"],
        )
        matches.append(repeat_id)

    return sorted(set(matches))


def technical_support_label(
    caller_count: float,
    cnvpytor_supported: bool,
    gridss_supported: bool,
    cnvpytor_only: bool,
) -> str:
    """Return a readable description of the technical support."""
    if cnvpytor_only:
        return "READ_DEPTH_ONLY_CNV"

    if caller_count >= 2 and cnvpytor_supported:
        return "MULTI_CALLER_PLUS_READ_DEPTH"

    if caller_count >= 2 and gridss_supported:
        return "MULTI_CALLER_PLUS_GRIDSS"

    if caller_count >= 2:
        return "MULTI_CALLER"

    if caller_count >= 1 and (
        cnvpytor_supported or gridss_supported
    ):
        return "SINGLE_CALLER_PLUS_SUPPORTING_EVIDENCE"

    if caller_count >= 1:
        return "SINGLE_CALLER"

    return "REVIEW_REQUIRED"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Add supporting SRS evidence to the final SV/gene table."
    )
    parser.add_argument("--integrated", required=True)
    parser.add_argument("--integrated-vcf", required=True)
    parser.add_argument("--gridss")
    parser.add_argument("--gridss-min-qual", type=float)
    parser.add_argument("--melt")
    parser.add_argument("--expansionhunter", required=True)
    parser.add_argument(
        "--breakpoint-tolerance",
        type=int,
        default=500,
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open_text(args.integrated) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        input_columns = list(reader.fieldnames or [])

    integrated_vcf = load_integrated_vcf(args.integrated_vcf)
    gridss_index = load_gridss(args.gridss, args.gridss_min_qual)
    melt_index = load_melt(args.melt)
    expansionhunter_index = load_expansionhunter(
        args.expansionhunter
    )

    for row in rows:
        sv_id = find_value(row, "SV_ID", "ID")
        vcf_record = integrated_vcf.get(sv_id)
        vcf_info = vcf_record["info"] if vcf_record else {}

        if "CNVPYTOR_ONLY" in vcf_info:
            cnvpytor_match = "DEPTH_ONLY"
        else:
            cnvpytor_match = vcf_info.get(
                "CNVPYTOR_MATCH",
                "NO",
            )

        row["CNVPYTOR_READ_DEPTH_MATCH"] = cnvpytor_match
        row["CNVPYTOR_STATUS"] = vcf_info.get(
            "CNVPYTOR_STATUS",
            MISSING,
        )
        row["CNVPYTOR_LEVEL"] = vcf_info.get(
            "CNVPYTOR_LEVEL",
            MISSING,
        )
        row["CNVPYTOR_EVAL1"] = vcf_info.get(
            "CNVPYTOR_EVAL1",
            MISSING,
        )
        row["CNVPYTOR_Q0"] = vcf_info.get(
            "CNVPYTOR_Q0",
            MISSING,
        )
        row["CNVPYTOR_PN"] = vcf_info.get(
            "CNVPYTOR_PN",
            MISSING,
        )
        row["CNVPYTOR_FLAGS"] = vcf_info.get(
            "CNVPYTOR_FLAGS",
            MISSING,
        )

        chromosome = normalize_chromosome(
            find_value(row, "CHROM", "Chr")
        )
        start = as_int(
            find_value(row, "START", "POS", "SV_start")
        )
        end = as_int(
            find_value(row, "END", "SV_end")
        )
        svtype = normalize_svtype(
            find_value(row, "SVTYPE", "SV_type")
        )

        chromosome_2 = chromosome
        if find_value(row, "CHR2") != MISSING:
            chromosome_2 = normalize_chromosome(
                find_value(row, "CHR2")
            )

        position_2 = as_int(
            find_value(row, "POS2")
        )

        gridss_status, gridss_records = gridss_support(
            chromosome,
            start,
            end,
            svtype,
            chromosome_2,
            position_2,
            gridss_index,
            args.breakpoint_tolerance,
        )

        row["GRIDSS_SUPPORT"] = (
            gridss_status
            if args.gridss
            else "NOT_RUN"
        )
        row["GRIDSS_RECORDS"] = (
            ";".join(gridss_records)
            if gridss_records
            else MISSING
        )

        melt_records = melt_matches(
            chromosome,
            start,
            svtype,
            melt_index,
            args.breakpoint_tolerance,
        )

        row["MELT_SUPPORT"] = (
            "YES"
            if melt_records
            else ("NO" if args.melt else "NOT_RUN")
        )
        row["MELT_RECORDS"] = (
            ";".join(melt_records)
            if melt_records
            else MISSING
        )

        repeat_loci = expansionhunter_matches(
            chromosome,
            start,
            end,
            svtype,
            expansionhunter_index,
        )

        row["EXPANSIONHUNTER_LOCUS_OVERLAP"] = (
            "YES" if repeat_loci else "NO"
        )
        row["EXPANSIONHUNTER_LOCI"] = (
            ";".join(repeat_loci)
            if repeat_loci
            else MISSING
        )

        caller_count = (
            as_float(
                find_value(
                    row,
                    "CALLER_COUNT",
                    "SUPP",
                )
            )
            or 0
        )

        cnvpytor_supported = (
            row["CNVPYTOR_STATUS"] == "PASS"
            or row["CNVPYTOR_READ_DEPTH_MATCH"] == "DEPTH_ONLY"
        )
        gridss_supported = row["GRIDSS_SUPPORT"] in {
            "YES_BOTH_BREAKPOINTS",
            "YES_SINGLE_BREAKPOINT",
        }
        cnvpytor_only = "CNVPYTOR_ONLY" in vcf_info

        row["TECHNICAL_SUPPORT"] = technical_support_label(
            caller_count,
            cnvpytor_supported,
            gridss_supported,
            cnvpytor_only,
        )

        row["TECHNICAL_SUPPORT_NOTE"] = (
            "Technical evidence only; not an ACMG or pathogenicity classification."
        )

    added_columns = [
        "CNVPYTOR_READ_DEPTH_MATCH",
        "CNVPYTOR_STATUS",
        "CNVPYTOR_LEVEL",
        "CNVPYTOR_EVAL1",
        "CNVPYTOR_Q0",
        "CNVPYTOR_PN",
        "CNVPYTOR_FLAGS",
        "GRIDSS_SUPPORT",
        "GRIDSS_RECORDS",
        "MELT_SUPPORT",
        "MELT_RECORDS",
        "EXPANSIONHUNTER_LOCUS_OVERLAP",
        "EXPANSIONHUNTER_LOCI",
        "TECHNICAL_SUPPORT",
        "TECHNICAL_SUPPORT_NOTE",
    ]

    output_columns = input_columns + [
        column
        for column in added_columns
        if column not in input_columns
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[OK] rows={len(rows)} output={output}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
