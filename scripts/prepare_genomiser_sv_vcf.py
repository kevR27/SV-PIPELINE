#!/usr/bin/env python3
"""Prepare a Genomiser-only indexed SV VCF with symbolic large deletions.

Exomiser/Jannovar can exhaust Java heap while constructing an HGVS string for
a multi-megabase sequence-encoded REF deletion. Convert only unambiguous
single-allele, left-anchored deletions; never change the source caller VCF.

Original coordinates, allele lengths, and existing END/SVLEN are recorded in
a separate audit TSV. Unrelated SVs, breakends, FORMAT, and genotypes remain
unchanged. This is VCF representation normalization, not variant filtering.
"""
import argparse
import csv
import gzip
import subprocess
from pathlib import Path


def _run(command):
    subprocess.run([str(value) for value in command], check=True)


def _info_dict(value):
    fields = [] if value in ("", ".") else value.split(";")
    return fields


def _get_info(fields, name):
    for field in fields:
        if field.startswith(name + "="):
            return field.split("=", 1)[1]
    return "."


def _set_info(fields, name, value):
    return [field for field in fields if not field.startswith(name + "=")] + [f"{name}={value}"]


def prepare_sv_vcf(source, output, audit, min_ref_length=1000, maximum_unsupported_length=1000000):
    source = Path(source)
    output = Path(output)
    audit = Path(audit)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.parent.mkdir(parents=True, exist_ok=True)
    plain = output.with_suffix("")  # local intermediate .vcf, NOT the caller source
    if plain == source:
        raise ValueError("Refusing to overwrite the source SV VCF")
    if output.resolve() == source.resolve():
        raise ValueError("Genomiser output and input SV VCF must be distinct")

    converted = 0
    inspected = 0
    headings = set()
    fieldnames = [
        "CHROM", "POS", "ID", "OLD_REF_LENGTH", "OLD_ALT_LENGTH",
        "NEW_END", "NEW_SVLEN", "OLD_END", "OLD_SVLEN",
        "END_LENGTH_DISAGREEMENT", "SVLEN_DISAGREEMENT",
    ]
    try:
        with gzip.open(source, "rt", encoding="utf-8") as infile, \\
             plain.open("w", encoding="utf-8") as outfile, \\
             audit.open("w", encoding="utf-8", newline="") as log:
            writer = csv.DictWriter(log, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for line in infile:
                if line.startswith("##"):
                    if line.startswith("##ALT=<ID=DEL,"):
                        headings.add("ALT")
                    if line.startswith("##INFO=<ID=SVTYPE,"):
                        headings.add("SVTYPE")
                    if line.startswith("##INFO=<ID=SVLEN,"):
                        headings.add("SVLEN")
                    if line.startswith("##INFO=<ID=END,"):
                        headings.add("END")
                    outfile.write(line)
                    continue
                if line.startswith("#CHROM"):
                    if "ALT" not in headings:
                        outfile.write('##ALT=<ID=DEL,Description="Deletion">\n')
                    if "SVTYPE" not in headings:
                        outfile.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                    if "SVLEN" not in headings:
                        outfile.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length of structural variant">\n')
                    if "END" not in headings:
                        outfile.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
                    outfile.write(line)
                    continue
                if line.startswith("#"):
                    outfile.write(line)
                    continue
                inspected += 1
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    raise ValueError(f"Malformed VCF data line {inspected}: fewer than 8 columns")
                chrom, pos_str, variant_id, ref, alt = parts[:5]
                info = _info_dict(parts[7])
                svtype = _get_info(info, "SVTYPE")
                is_simple_deletion = (
                    len(ref) >= min_ref_length
                    and len(alt) == 1
                    and alt.upper() == ref[0].upper()
                    and svtype in ("DEL", ".")
                )
                if is_simple_deletion:
                    pos = int(pos_str)
                    old_end = _get_info(info, "END")
                    old_svlen = _get_info(info, "SVLEN")
                    end = pos + len(ref) - 1
                    svlen = 1 - len(ref)
                    writer.writerow({
                        "CHROM": chrom, "POS": pos, "ID": variant_id,
                        "OLD_REF_LENGTH": len(ref), "OLD_ALT_LENGTH": len(alt),
                        "NEW_END": end, "NEW_SVLEN": svlen,
                        "OLD_END": old_end, "OLD_SVLEN": old_svlen,
                        "END_LENGTH_DISAGREEMENT": str(old_end not in (".", str(end))),
                        "SVLEN_DISAGREEMENT": str(old_svlen not in (".", str(svlen))),
                    })
                    info = _set_info(info, "SVTYPE", "DEL")
                    info = _set_info(info, "END", end)
                    info = _set_info(info, "SVLEN", svlen)
                    parts[3] = ref[0]
                    parts[4] = "<DEL>"
                    parts[7] = ";".join(info)
                    converted += 1
                    outfile.write("\t".join(parts) + "\n")
                else:
                    if max(len(ref), max(len(a) for a in alt.split(","))) > maximum_unsupported_length:
                        raise ValueError(
                            f"Unconverted explicit allele over {maximum_unsupported_length} bp "
                            f"at {chrom}:{pos_str} ({variant_id}), "
                            "inspect before Genomiser to avoid Java heap exhaustion"
                        )
                    outfile.write(line)

        _run(["bcftools", "view", "-Oz", "-o", output, plain])
        _run(["bcftools", "index", "-f", "-t", output])
        print(f"[Genomiser SV normalization] scanned={inspected} symbolic_deletions={converted} "
              f"audit={audit} output={output}", flush=True)
        return {"scanned": inspected, "symbolic_deletions": converted, "audit": str(audit)}
    finally:
        plain.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input single-sample, bgzipped SV VCF")
    parser.add_argument("--output", required=True, help="Separate BGZF Genomiser-only SV VCF")
    parser.add_argument("--audit", required=True, help="TSV mapping of converted deletion records")
    parser.add_argument("--min-ref-length", type=int, default=1000)
    args = parser.parse_args()
    if args.min_ref_length < 50:
        parser.error("--min-ref-length must be at least 50")
    prepare_sv_vcf(args.input, args.output, args.audit, args.min_ref_length)


if __name__ == "__main__":
    main()
