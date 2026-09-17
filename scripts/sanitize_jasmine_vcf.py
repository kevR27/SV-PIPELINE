#!/usr/bin/env python3
"""Repair Jasmine VCF headers without changing variant records.

Jasmine writes its merged header from the first VCF in the input list. When
multiple callers are merged, INFO/FILTER/FORMAT tags used by a consensus record
can therefore be absent from that header even though the record itself is
valid. bcftools/AnnotSV may then reject the VCF.

This script performs a two-pass, text-based repair:
  * preserve every record exactly as Jasmine wrote it;
  * preserve every existing header declaration;
  * add declarations only for metadata that is used but missing;
  * add missing contig declarations from the reference .fai when supplied.

No INFO values are removed. In particular BND partner information is retained.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_declared(header_lines: list[str], prefix: str) -> set[str]:
    out: set[str] = set()
    rx = re.compile(rf"##{re.escape(prefix)}=<ID=([^,>]+)")
    for line in header_lines:
        m = rx.match(line)
        if m:
            out.add(m.group(1))
    return out


def read_fai(path: str | None) -> dict[str, int]:
    if not path:
        return {}
    lengths: dict[str, int] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                try:
                    lengths[fields[0]] = int(fields[1])
                except ValueError:
                    pass
    return lengths


KNOWN_INFO = {
    "SVTYPE": ('1', 'String'),
    "END": ('1', 'Integer'),
    "SVLEN": ('1', 'String'),
    "SUPP": ('1', 'Integer'),
    "SUPP_EXT": ('1', 'Integer'),
    "SUPP_VEC": ('1', 'String'),
    "SUPP_VEC_EXT": ('1', 'String'),
    "IDLIST": ('.', 'String'),
    "IDLIST_EXT": ('.', 'String'),
    "INTRASAMPLE_IDLIST": ('.', 'String'),
    "ALLVARS_EXT": ('.', 'String'),
    "VARCALLS": ('1', 'Integer'),
    "SVMETHOD": ('1', 'String'),
    "STARTVARIANCE": ('1', 'String'),
    "ENDVARIANCE": ('1', 'String'),
    "AVG_LEN": ('1', 'String'),
    "AVG_START": ('1', 'String'),
    "AVG_END": ('1', 'String'),
    "CHR2": ('1', 'String'),
    "CIPOS": ('2', 'Integer'),
    "CIEND": ('2', 'Integer'),
    "CILEN": ('2', 'Integer'),
    "RE": ('1', 'Integer'),
    "SUPPORT": ('1', 'Integer'),
    "AF": ('A', 'Float'),
    "STRANDS": ('1', 'String'),
}

KNOWN_FORMAT = {
    "GT": ('1', 'String'),
    "GQ": ('1', 'Integer'),
    "DP": ('1', 'Integer'),
    "DR": ('1', 'Integer'),
    "DV": ('1', 'Integer'),
    "RE": ('1', 'Integer'),
    "AD": ('R', 'Integer'),
    "PL": ('G', 'Integer'),
}


def definition(kind: str, key: str, flag: bool = False) -> str:
    if kind == "INFO":
        if flag:
            number, typ = "0", "Flag"
        else:
            number, typ = KNOWN_INFO.get(key, (".", "String"))
    else:
        number, typ = KNOWN_FORMAT.get(key, (".", "String"))
    return (
        f'##{kind}=<ID={key},Number={number},Type={typ},'
        f'Description="Added by sanitize_jasmine_vcf.py because Jasmine output uses this field">'
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--reference-fai",
        default=None,
        help="Optional FASTA .fai used to add missing contig declarations",
    )
    args = ap.parse_args()

    header_meta: list[str] = []
    chrom_header: str | None = None

    with open_text(args.input) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                header_meta.append(line)
            elif line.startswith("#CHROM"):
                chrom_header = line
                break

    if chrom_header is None:
        raise SystemExit("ERROR: input VCF has no #CHROM header")

    declared_info = parse_declared(header_meta, "INFO")
    declared_format = parse_declared(header_meta, "FORMAT")
    declared_filter = parse_declared(header_meta, "FILTER")
    declared_contig = parse_declared(header_meta, "contig")

    used_info: set[str] = set()
    info_flags: set[str] = set()
    used_format: set[str] = set()
    used_filter: set[str] = set()
    used_contig: set[str] = set()

    with open_text(args.input) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise SystemExit(f"ERROR: malformed Jasmine VCF record: {line.rstrip()}")

            used_contig.add(fields[0])

            filt = fields[6]
            if filt not in {"", ".", "PASS"}:
                used_filter.update(x for x in filt.split(";") if x)

            info = fields[7]
            if info not in {"", "."}:
                for item in info.split(";"):
                    if not item:
                        continue
                    if "=" in item:
                        key = item.split("=", 1)[0]
                        used_info.add(key)
                    else:
                        used_info.add(item)
                        info_flags.add(item)

            if len(fields) >= 9 and fields[8] not in {"", "."}:
                used_format.update(x for x in fields[8].split(":") if x)

    repaired = list(header_meta)

    for key in sorted(used_info - declared_info):
        repaired.append(definition("INFO", key, flag=key in info_flags))

    for key in sorted(used_format - declared_format):
        repaired.append(definition("FORMAT", key))

    for key in sorted(used_filter - declared_filter):
        repaired.append(
            f'##FILTER=<ID={key},Description="Added by sanitize_jasmine_vcf.py because Jasmine output uses this FILTER">'
        )

    contig_lengths = read_fai(args.reference_fai)
    for contig in sorted(used_contig - declared_contig):
        if contig_lengths:
            if contig not in contig_lengths:
                raise SystemExit(
                    f"ERROR: VCF uses contig {contig!r}, but it is absent from {args.reference_fai}"
                )
            repaired.append(f"##contig=<ID={contig},length={contig_lengths[contig]}>")
        else:
            repaired.append(f"##contig=<ID={contig}>")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="") as out:
        for line in repaired:
            out.write(line + "\n")
        out.write(chrom_header + "\n")

        with open_text(args.input) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                out.write(line if line.endswith("\n") else line + "\n")

    print(
        "[OK] Jasmine header repaired: "
        f"INFO+{len(used_info - declared_info)} "
        f"FORMAT+{len(used_format - declared_format)} "
        f"FILTER+{len(used_filter - declared_filter)} "
        f"contig+{len(used_contig - declared_contig)}"
    )
    print(f"[OK] output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
