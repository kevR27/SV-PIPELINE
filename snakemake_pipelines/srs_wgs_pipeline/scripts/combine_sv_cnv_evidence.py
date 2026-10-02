#!/usr/bin/env python3
"""
Combine breakpoint-based SV calls with CNVpytor read-depth CNV calls.

Why this script exists
----------------------
Manta and DELLY are good at detecting SV breakpoints from paired-end and
split-read evidence. CNVpytor provides a different type of evidence: changes
in read depth.

The aim is therefore not to make CNVpytor behave like another breakpoint
caller. Instead:

1. Keep the Manta+DELLY SURVIVOR callset as the breakpoint-based SV callset.
2. For deletions and duplications, check whether CNVpytor detects an
   overlapping read-depth change.
3. Add the CNVpytor information to matching SVs.
4. Retain CNVpytor-only deletions/duplications when they pass the configured
   quality filters, so a real depth-only CNV is not lost.

The output is the integrated SRS structural-variant VCF used by AnnotSV/VEP.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
from pathlib import Path

MISSING = "."


def open_text(path: str):
    """Open a normal or gzipped text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(info_text: str) -> dict[str, str]:
    """Convert a VCF INFO field into a dictionary."""
    result: dict[str, str] = {}

    for item in str(info_text).split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "True"

    return result


def as_number(value):
    """Return a finite float or None."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(number):
        return None

    return number


def normalize_chromosome(chromosome: str) -> str:
    """Use UCSC-style chromosome names to match hg38 BAM/VCF files."""
    chromosome = str(chromosome)

    if chromosome == "MT":
        return "chrM"

    if chromosome.startswith("chr"):
        return chromosome

    return "chr" + chromosome


def normalize_svtype(value: str) -> str:
    """Reduce CNV/SV labels to the types needed for this integration."""
    value = str(value).upper()

    if "DEL" in value:
        return "DEL"

    if "DUP" in value or "GAIN" in value:
        return "DUP"

    return value


def reciprocal_overlap(
    start_a: int,
    end_a: int,
    start_b: int,
    end_b: int,
) -> float:
    """
    Return the smaller reciprocal overlap between two genomic intervals.

    A value of 0.5 means that at least 50% of both events overlap.
    """
    overlap_start = max(min(start_a, end_a), min(start_b, end_b))
    overlap_end = min(max(start_a, end_a), max(start_b, end_b))
    overlap_bp = max(0, overlap_end - overlap_start + 1)

    if overlap_bp == 0:
        return 0.0

    length_a = max(1, abs(end_a - start_a) + 1)
    length_b = max(1, abs(end_b - start_b) + 1)

    return min(overlap_bp / length_a, overlap_bp / length_b)


def read_cnvpytor_calls(
    path: str,
    max_q0: float,
    max_pn: float,
    max_eval1: float,
) -> list[dict]:
    """
    Read the standard CNVpytor '-call' TSV output.

    CNVpytor columns used here:
      1  CNV type
      2  genomic region
      3  size
      4  normalized read-depth level
      5  e-val1
      9  q0
      10 pN
      11 distance from large reference gap

    The quality thresholds are used only to decide whether an unmatched
    CNVpytor call is strong enough to remain as a depth-only CNV.
    """
    calls: list[dict] = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split()
            if len(fields) < 4:
                continue

            svtype = normalize_svtype(fields[0])
            if svtype not in {"DEL", "DUP"}:
                continue

            region = re.match(
                r"([^:]+):(\d+)-(\d+)",
                fields[1].replace(",", ""),
            )
            if region is None:
                continue

            chromosome = normalize_chromosome(region.group(1))
            start = int(region.group(2))
            end = int(region.group(3))

            size = as_number(fields[2])
            level = as_number(fields[3])
            eval1 = as_number(fields[4]) if len(fields) > 4 else None
            eval2 = as_number(fields[5]) if len(fields) > 5 else None
            eval3 = as_number(fields[6]) if len(fields) > 6 else None
            eval4 = as_number(fields[7]) if len(fields) > 7 else None
            q0 = as_number(fields[8]) if len(fields) > 8 else None
            pN = as_number(fields[9]) if len(fields) > 9 else None
            distance_to_gap = as_number(fields[10]) if len(fields) > 10 else None

            flags: list[str] = []

            if q0 is not None and q0 > max_q0:
                flags.append("HIGH_Q0")

            if pN is not None and pN > max_pn:
                flags.append("HIGH_REFERENCE_N_CONTENT")

            if eval1 is not None and eval1 > max_eval1:
                flags.append("WEAK_READ_DEPTH_EVALUE")

            status = "PASS" if not flags else "REVIEW"

            calls.append(
                {
                    "index": line_number,
                    "chrom": chromosome,
                    "start": start,
                    "end": end,
                    "svtype": svtype,
                    "size": size,
                    "level": level,
                    "eval1": eval1,
                    "eval2": eval2,
                    "eval3": eval3,
                    "eval4": eval4,
                    "q0": q0,
                    "pN": pN,
                    "distance_to_gap": distance_to_gap,
                    "status": status,
                    "flags": ";".join(flags) if flags else MISSING,
                }
            )

    return calls


def read_vcf(path: str):
    """Read VCF metadata, header and records without changing the genotypes."""
    metadata: list[str] = []
    header = None
    records: list[dict] = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                metadata.append(line)
                continue

            if line.startswith("#CHROM"):
                header = line
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            info = parse_info(fields[7])
            svtype = normalize_svtype(
                info.get("SVTYPE", fields[4].strip("<>"))
            )

            end = as_number(info.get("END"))
            start = int(fields[1])

            records.append(
                {
                    "fields": fields,
                    "info": info,
                    "chrom": normalize_chromosome(fields[0]),
                    "start": start,
                    "end": int(end) if end is not None else start,
                    "svtype": svtype,
                }
            )

    if header is None:
        raise ValueError("Input VCF does not contain a #CHROM header.")

    return metadata, header, records


def format_value(value) -> str:
    """Write numeric evidence compactly in the VCF INFO field."""
    if value is None:
        return MISSING

    if isinstance(value, float):
        return f"{value:.8g}"

    return str(value)


def add_cnvpytor_evidence_to_breakpoint_calls(
    breakpoint_records: list[dict],
    cnv_calls: list[dict],
    minimum_overlap: float,
) -> set[int]:
    """Attach the best overlapping same-type CNVpytor call to DEL/DUP events."""
    matched_cnvpytor: set[int] = set()

    for record in breakpoint_records:
        if record["svtype"] not in {"DEL", "DUP"}:
            continue

        best_index = None
        best_overlap = 0.0

        for index, cnv in enumerate(cnv_calls):
            if cnv["chrom"] != record["chrom"]:
                continue

            if cnv["svtype"] != record["svtype"]:
                continue

            overlap = reciprocal_overlap(
                record["start"],
                record["end"],
                cnv["start"],
                cnv["end"],
            )

            if overlap >= minimum_overlap and overlap > best_overlap:
                best_index = index
                best_overlap = overlap

        if best_index is None:
            continue

        cnv = cnv_calls[best_index]
        matched_cnvpytor.add(best_index)

        cnvpytor_info = {
            "CNVPYTOR_MATCH": "YES",
            "CNVPYTOR_STATUS": cnv["status"],
            "CNVPYTOR_RECIPROCAL_OVERLAP": format_value(best_overlap),
            "CNVPYTOR_LEVEL": format_value(cnv["level"]),
            "CNVPYTOR_EVAL1": format_value(cnv["eval1"]),
            "CNVPYTOR_Q0": format_value(cnv["q0"]),
            "CNVPYTOR_PN": format_value(cnv["pN"]),
            "CNVPYTOR_FLAGS": cnv["flags"],
        }

        record["info"].update(cnvpytor_info)
        record["fields"][7] = ";".join(
            key if value == "True" else f"{key}={value}"
            for key, value in record["info"].items()
        )

    return matched_cnvpytor


def add_depth_only_cnvs(
    records: list[dict],
    cnv_calls: list[dict],
    matched_cnvpytor: set[int],
    vcf_header: str,
) -> int:
    """Add PASS CNVpytor calls that did not match a Manta/DELLY event."""
    added = 0
    sample_columns = max(
        0,
        len(vcf_header.rstrip("\n").split("\t")) - 9,
    )

    for index, cnv in enumerate(cnv_calls):
        if index in matched_cnvpytor:
            continue

        if cnv["status"] != "PASS":
            continue

        length = cnv["end"] - cnv["start"] + 1
        svlen = -length if cnv["svtype"] == "DEL" else length

        info = {
            "SVTYPE": cnv["svtype"],
            "END": str(cnv["end"]),
            "SVLEN": str(svlen),
            "SRS_SOURCE": "CNVPytor",
            "CNVPYTOR_ONLY": "True",
            "CNVPYTOR_STATUS": cnv["status"],
            "CNVPYTOR_LEVEL": format_value(cnv["level"]),
            "CNVPYTOR_EVAL1": format_value(cnv["eval1"]),
            "CNVPYTOR_Q0": format_value(cnv["q0"]),
            "CNVPYTOR_PN": format_value(cnv["pN"]),
            "CNVPYTOR_FLAGS": cnv["flags"],
        }

        fields = [
            cnv["chrom"],
            str(cnv["start"]),
            f"CNVPYTOR_ONLY_{cnv['index']}",
            "N",
            f"<{cnv['svtype']}>",
            ".",
            "PASS",
            ";".join(
                key if value == "True" else f"{key}={value}"
                for key, value in info.items()
            ),
        ]

        if sample_columns:
            fields += ["GT"] + ["./."] * sample_columns

        records.append(
            {
                "fields": fields,
                "info": info,
                "chrom": cnv["chrom"],
                "start": cnv["start"],
                "end": cnv["end"],
                "svtype": cnv["svtype"],
            }
        )
        added += 1

    return added


def add_vcf_info_headers(metadata: list[str]) -> list[str]:
    """Add definitions for the new CNVpytor INFO fields."""
    definitions = [
        '##INFO=<ID=SRS_SOURCE,Number=1,Type=String,Description="Additional SRS evidence source">\n',
        '##INFO=<ID=CNVPYTOR_ONLY,Number=0,Type=Flag,Description="Depth-only CNV added from CNVpytor">\n',
        '##INFO=<ID=CNVPYTOR_MATCH,Number=1,Type=String,Description="CNVpytor read-depth call overlaps breakpoint event">\n',
        '##INFO=<ID=CNVPYTOR_STATUS,Number=1,Type=String,Description="CNVpytor read-depth evidence QC status">\n',
        '##INFO=<ID=CNVPYTOR_RECIPROCAL_OVERLAP,Number=1,Type=Float,Description="Minimum reciprocal overlap with CNVpytor call">\n',
        '##INFO=<ID=CNVPYTOR_LEVEL,Number=1,Type=Float,Description="CNVpytor normalized read-depth level">\n',
        '##INFO=<ID=CNVPYTOR_EVAL1,Number=1,Type=Float,Description="CNVpytor primary read-depth e-value">\n',
        '##INFO=<ID=CNVPYTOR_Q0,Number=1,Type=Float,Description="CNVpytor q0 fraction">\n',
        '##INFO=<ID=CNVPYTOR_PN,Number=1,Type=Float,Description="CNVpytor reference N fraction">\n',
        '##INFO=<ID=CNVPYTOR_FLAGS,Number=1,Type=String,Description="CNVpytor review flags">\n',
    ]

    existing = "".join(metadata)

    for definition in definitions:
        info_id = definition.split("ID=", 1)[1].split(",", 1)[0]
        if info_id not in existing:
            metadata.append(definition)

    return metadata


def chromosome_order(reference_fai: str) -> dict[str, int]:
    """Return chromosome order from the same FASTA index used for alignment."""
    order: dict[str, int] = {}

    with open(reference_fai, "r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            chromosome = line.split("\t", 1)[0]
            order[chromosome] = index

    return order


def write_vcf(
    path: str,
    metadata: list[str],
    header: str,
    records: list[dict],
    reference_fai: str,
):
    """Sort records by reference order and write the integrated VCF."""
    order = chromosome_order(reference_fai)

    records.sort(
        key=lambda record: (
            order.get(record["chrom"], 10**9),
            record["start"],
            record["end"],
            record["fields"][2],
        )
    )

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as handle:
        handle.writelines(metadata)
        handle.write(header if header.endswith("\n") else header + "\n")

        for record in records:
            handle.write("\t".join(record["fields"]) + "\n")


def write_integration_table(
    path: str,
    cnv_calls: list[dict],
    matched_cnvpytor: set[int],
):
    """Write an audit table showing what happened to every CNVpytor call."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as handle:
        handle.write(
            "CNVPYTOR_INDEX\tCHROM\tSTART\tEND\tSVTYPE\t"
            "STATUS\tFLAGS\tMATCHED_TO_BREAKPOINT_CALL\n"
        )

        for index, cnv in enumerate(cnv_calls):
            matched = "YES" if index in matched_cnvpytor else "NO"

            handle.write(
                f"{cnv['index']}\t{cnv['chrom']}\t{cnv['start']}\t"
                f"{cnv['end']}\t{cnv['svtype']}\t{cnv['status']}\t"
                f"{cnv['flags']}\t{matched}\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Combine Manta/DELLY SVs with CNVpytor read-depth CNVs."
    )
    parser.add_argument("--breakpoint-vcf", required=True)
    parser.add_argument("--cnvpytor-calls", required=True)
    parser.add_argument("--reference-fai", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit-output", required=True)
    parser.add_argument("--min-reciprocal-overlap", type=float, default=0.5)
    parser.add_argument("--max-q0", type=float, default=0.5)
    parser.add_argument("--max-pn", type=float, default=0.5)
    parser.add_argument("--max-eval1", type=float, default=1e-4)
    args = parser.parse_args()

    metadata, header, breakpoint_records = read_vcf(args.breakpoint_vcf)

    cnv_calls = read_cnvpytor_calls(
        args.cnvpytor_calls,
        max_q0=args.max_q0,
        max_pn=args.max_pn,
        max_eval1=args.max_eval1,
    )

    matched_cnvpytor = add_cnvpytor_evidence_to_breakpoint_calls(
        breakpoint_records,
        cnv_calls,
        minimum_overlap=args.min_reciprocal_overlap,
    )

    depth_only_added = add_depth_only_cnvs(
        breakpoint_records,
        cnv_calls,
        matched_cnvpytor,
        header,
    )

    metadata = add_vcf_info_headers(metadata)

    write_vcf(
        args.output,
        metadata,
        header,
        breakpoint_records,
        args.reference_fai,
    )

    write_integration_table(
        args.audit_output,
        cnv_calls,
        matched_cnvpytor,
    )

    print(
        "[OK] "
        f"CNVpytor_calls={len(cnv_calls)} "
        f"matched_to_breakpoint_calls={len(matched_cnvpytor)} "
        f"depth_only_CNVS_added={depth_only_added} "
        f"output={args.output}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
