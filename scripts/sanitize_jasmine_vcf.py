#!/usr/bin/env python3
"""Repair Jasmine VCF metadata without changing variant records.

Jasmine can propagate INFO/FORMAT/FILTER tags from caller records without
always carrying every corresponding header declaration.  Downstream tools such
as bcftools, VEP and AnnotSV are stricter.  This script scans the complete VCF,
adds only missing metadata definitions, adds missing contig declarations from
the reference FAI, and writes a plain-text VCF for subsequent bcftools sorting.

Important: variant coordinates, alleles, IDs, INFO values and sample columns are
preserved exactly.  In particular, BND records are not rewritten and no caller
annotations are discarded.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


def open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def read_fai(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                try:
                    lengths[fields[0]] = int(fields[1])
                except ValueError:
                    pass
    return lengths


def header_ids(lines: list[str], kind: str) -> set[str]:
    out: set[str] = set()
    prefix = f"##{kind}=<ID="
    for line in lines:
        if line.startswith(prefix):
            match = re.match(rf"##{re.escape(kind)}=<ID=([^,>]+)", line)
            if match:
                out.add(match.group(1))
    return out


def infer_scalar_type(value: str) -> str:
    if re.fullmatch(r"-?\d+", value or ""):
        return "Integer"
    if re.fullmatch(r"-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", value or ""):
        return "Float"
    return "String"


def info_definition(key: str, values: list[str], is_flag: bool) -> str:
    known = {
        "END": ('Number=1,Type=Integer', "End position of the structural variant"),
        "SVLEN": ('Number=1,Type=Integer', "Structural variant length"),
        "SVTYPE": ('Number=1,Type=String', "Structural variant type"),
        "CHR2": ('Number=1,Type=String', "Chromosome of second breakpoint"),
        "POS2": ('Number=1,Type=Integer', "Position of second breakpoint"),
        "SUPP": ('Number=1,Type=Integer', "Number of input callsets supporting the merged SV"),
        "SUPP_VEC": ('Number=1,Type=String', "Input-callset support vector in Jasmine file-list order"),
        "IDLIST": ('Number=.,Type=String', "IDs of input records represented by the merged SV"),
        "CIPOS": ('Number=2,Type=Integer', "Confidence interval around POS"),
        "CIEND": ('Number=2,Type=Integer', "Confidence interval around END"),
        "CILEN": ('Number=2,Type=Integer', "Confidence interval around SV length"),
        "PRECISE": ('Number=0,Type=Flag', "Precise structural variant"),
        "IMPRECISE": ('Number=0,Type=Flag', "Imprecise structural variant"),
    }
    if key in known:
        spec, description = known[key]
        return f'##INFO=<ID={key},{spec},Description="{description}">'

    if is_flag:
        return (
            f'##INFO=<ID={key},Number=0,Type=Flag,'
            'Description="Automatically restored metadata for a flag present in a Jasmine record">'
        )

    sample_value = next((v for v in values if v not in {"", "."}), "")
    if "," in sample_value:
        number = "."
        value_type = "String"
    else:
        number = "1"
        value_type = infer_scalar_type(sample_value)
    return (
        f'##INFO=<ID={key},Number={number},Type={value_type},'
        'Description="Automatically restored metadata for a field present in a Jasmine record">'
    )


def format_definition(key: str, values: list[str]) -> str:
    known = {
        "GT": ('Number=1,Type=String', "Genotype"),
        "GQ": ('Number=1,Type=Integer', "Genotype quality"),
        "DP": ('Number=1,Type=Integer', "Read depth"),
        "DR": ('Number=1,Type=Integer', "Reference-supporting reads"),
        "DV": ('Number=1,Type=Integer', "Variant-supporting reads"),
        "RR": ('Number=1,Type=Integer', "Reference junction reads"),
        "RV": ('Number=1,Type=Integer', "Variant junction reads"),
        "PS": ('Number=1,Type=Integer', "Phase set"),
    }
    if key in known:
        spec, description = known[key]
        return f'##FORMAT=<ID={key},{spec},Description="{description}">'

    sample_value = next((v for v in values if v not in {"", "."}), "")
    if "," in sample_value:
        number = "."
        value_type = "String"
    else:
        number = "1"
        value_type = infer_scalar_type(sample_value)
    return (
        f'##FORMAT=<ID={key},Number={number},Type={value_type},'
        'Description="Automatically restored metadata for a FORMAT field present in a Jasmine record">'
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--reference-fai", required=True)
    args = parser.parse_args()

    metadata: list[str] = []
    chrom_header: str | None = None
    records: list[str] = []

    used_info: dict[str, list[str]] = {}
    flag_info: set[str] = set()
    used_format: dict[str, list[str]] = {}
    used_filter: set[str] = set()
    used_contigs: set[str] = set()

    with open_text(args.input) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                metadata.append(line)
                continue
            if line.startswith("#CHROM"):
                chrom_header = line
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed Jasmine VCF record: {line}")
            records.append(line)
            used_contigs.add(fields[0])

            filt = fields[6]
            if filt not in {"", ".", "PASS"}:
                used_filter.update(x for x in filt.split(";") if x)

            info = fields[7]
            if info not in {"", "."}:
                for item in info.split(";"):
                    if not item:
                        continue
                    if "=" in item:
                        key, value = item.split("=", 1)
                        used_info.setdefault(key, [])
                        if len(used_info[key]) < 20:
                            used_info[key].append(value)
                    else:
                        flag_info.add(item)
                        used_info.setdefault(item, [])

            if len(fields) >= 9 and fields[8] not in {"", "."}:
                fmt_keys = fields[8].split(":")
                for idx, key in enumerate(fmt_keys):
                    used_format.setdefault(key, [])
                    for sample_field in fields[9:]:
                        values = sample_field.split(":")
                        if idx < len(values) and len(used_format[key]) < 20:
                            used_format[key].append(values[idx])

    if chrom_header is None:
        raise ValueError("Input VCF has no #CHROM header line")

    declared_info = header_ids(metadata, "INFO")
    declared_format = header_ids(metadata, "FORMAT")
    declared_filter = header_ids(metadata, "FILTER")
    declared_contig = header_ids(metadata, "contig")

    additions: list[str] = []
    for key in sorted(used_info):
        if key not in declared_info:
            additions.append(info_definition(key, used_info[key], key in flag_info))
    for key in sorted(used_format):
        if key not in declared_format:
            additions.append(format_definition(key, used_format[key]))
    for key in sorted(used_filter):
        if key not in declared_filter:
            additions.append(
                f'##FILTER=<ID={key},Description="Automatically restored FILTER metadata">'
            )

    ref_lengths = read_fai(args.reference_fai)
    for contig in sorted(used_contigs):
        if contig in declared_contig:
            continue
        if contig not in ref_lengths:
            raise ValueError(
                f"Contig {contig!r} occurs in Jasmine VCF but is absent from {args.reference_fai}"
            )
        additions.append(f"##contig=<ID={contig},length={ref_lengths[contig]}>")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as out:
        for line in metadata:
            out.write(line + "\n")
        for line in additions:
            out.write(line + "\n")
        out.write(chrom_header + "\n")
        for line in records:
            out.write(line + "\n")

    print(
        f"[OK] records={len(records)} added_metadata={len(additions)} "
        f"output={output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
