#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile


REMOVE_INFO = {
    "AF",
    "DR",
    "DV",
}

REMOVE_FORMAT = {
    "DR",
    "DV",
}


def run(cmd):
    print("[CMD]", " ".join(cmd), flush=True)

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        sys.exit(
            f"ERROR: command failed with exit code {result.returncode}"
        )

    if result.stderr:
        print(result.stderr, file=sys.stderr)


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fai(path):
    """
    Read a FASTA .fai index and return:
        {chromosome: length}
    """
    if not os.path.exists(path):
        sys.exit(f"ERROR: FASTA index not found: {path}")

    lengths = {}

    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                continue

            chrom = fields[0]

            try:
                length = int(fields[1])
            except ValueError:
                continue

            lengths[chrom] = length

    return lengths


def infer_info_definition(key, values):
    value = values[0] if values else ""

    if re.fullmatch(r"-?\d+", value or ""):
        vtype = "Integer"
    elif re.fullmatch(r"-?\d+\.\d+", value or ""):
        vtype = "Float"
    else:
        vtype = "String"

    return (
        f'##INFO=<ID={key},Number=1,Type={vtype},'
        f'Description="Automatically added missing INFO field">'
    )


def infer_format_definition(key, values):
    value = values[0] if values else ""

    if re.fullmatch(r"-?\d+", value or ""):
        vtype = "Integer"
    elif re.fullmatch(r"-?\d+\.\d+", value or ""):
        vtype = "Float"
    else:
        vtype = "String"

    return (
        f'##FORMAT=<ID={key},Number=1,Type={vtype},'
        f'Description="Automatically added missing FORMAT field">'
    )


def infer_filter_definition(key):
    return (
        f'##FILTER=<ID={key},'
        f'Description="Automatically added missing FILTER field">'
    )


def collect_missing_metadata(input_vcf, header_lines):

    declared_info = set()
    declared_format = set()
    declared_filter = set()
    declared_contig = set()

    for line in header_lines:

        if line.startswith("##INFO=<ID="):
            m = re.match(r"##INFO=<ID=([^,>]+)", line)
            if m:
                declared_info.add(m.group(1))

        elif line.startswith("##FORMAT=<ID="):
            m = re.match(r"##FORMAT=<ID=([^,>]+)", line)
            if m:
                declared_format.add(m.group(1))

        elif line.startswith("##FILTER=<ID="):
            m = re.match(r"##FILTER=<ID=([^,>]+)", line)
            if m:
                declared_filter.add(m.group(1))

        elif line.startswith("##contig=<ID="):
            m = re.match(r"##contig=<ID=([^,>]+)", line)
            if m:
                declared_contig.add(m.group(1))

    used_info = {}
    used_format = {}
    used_filter = set()
    used_contigs = set()

    with open_vcf(input_vcf) as f:

        for line in f:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom = fields[0]
            filt = fields[6]
            info = fields[7]

            used_contigs.add(chrom)

            if filt not in {".", "PASS"}:
                for value in filt.split(";"):
                    if value:
                        used_filter.add(value)

            if info not in {"", "."}:
                for item in info.split(";"):

                    if "=" in item:
                        key, value = item.split("=", 1)

                        if key not in used_info:
                            used_info[key] = []

                        if len(used_info[key]) < 10:
                            used_info[key].append(value)

                    else:
                        key = item

                        if key not in used_info:
                            used_info[key] = []

            if len(fields) >= 9:

                format_field = fields[8]

                if format_field not in {"", "."}:

                    format_keys = format_field.split(":")

                    for idx, key in enumerate(format_keys):

                        if key not in used_format:
                            used_format[key] = []

                        for sample in fields[9:]:

                            sample_values = sample.split(":")

                            if idx < len(sample_values):

                                value = sample_values[idx]

                                if (
                                    value not in {"", "."}
                                    and len(used_format[key]) < 10
                                ):
                                    used_format[key].append(value)

    missing_info = {
        key: values
        for key, values in used_info.items()
        if key not in declared_info
    }

    missing_format = {
        key: values
        for key, values in used_format.items()
        if key not in declared_format
    }

    missing_filter = {
        key
        for key in used_filter
        if key not in declared_filter
    }

    missing_contigs = {
        chrom
        for chrom in used_contigs
        if chrom not in declared_contig
    }

    return (
        missing_info,
        missing_format,
        missing_filter,
        missing_contigs,
    )


def clean_info(info):
    if info in {"", "."}:
        return "."

    kept = []

    for item in info.split(";"):

        key = item.split("=", 1)[0]

        if key in REMOVE_INFO:
            continue

        kept.append(item)

    return ";".join(kept) if kept else "."


def clean_record(line):

    fields = line.rstrip("\n").split("\t")

    if len(fields) < 8:
        return None

    # ---------------------------------------------------------
    # INFO
    # ---------------------------------------------------------

    fields[7] = clean_info(fields[7])

    # ---------------------------------------------------------
    # FORMAT + samples
    # ---------------------------------------------------------

    if len(fields) >= 9:

        format_keys = fields[8].split(":")

        keep_indices = [
            i
            for i, key in enumerate(format_keys)
            if key not in REMOVE_FORMAT
        ]

        kept_format_keys = [
            format_keys[i]
            for i in keep_indices
        ]

        fields[8] = (
            ":".join(kept_format_keys)
            if kept_format_keys
            else "."
        )

        for sample_idx in range(9, len(fields)):

            sample_values = fields[sample_idx].split(":")

            kept_values = [
                sample_values[i]
                if i < len(sample_values)
                else "."
                for i in keep_indices
            ]

            fields[sample_idx] = (
                ":".join(kept_values)
                if kept_values
                else "."
            )

    return "\t".join(fields)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prepare a Jasmine merged SV VCF for NeedLR."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    parser.add_argument(
        "--reference-fai",
        required=True,
        help="FASTA .fai index used to add contig lengths",
    )

    args = parser.parse_args()

    if shutil.which("bcftools") is None:
        sys.exit("ERROR: bcftools not found in PATH.")

    if shutil.which("tabix") is None:
        sys.exit("ERROR: tabix not found in PATH.")

    input_vcf = os.path.abspath(args.input)
    output_vcf = os.path.abspath(args.output)
    reference_fai = os.path.abspath(args.reference_fai)

    if not os.path.exists(input_vcf):
        sys.exit(
            f"ERROR: Input VCF not found: {input_vcf}"
        )

    contig_lengths = read_fai(reference_fai)

    print(
        f"Loaded {len(contig_lengths)} contig lengths "
        f"from {reference_fai}"
    )

    os.makedirs(
        os.path.dirname(output_vcf),
        exist_ok=True,
    )

    with tempfile.TemporaryDirectory(
        prefix="needlr_"
    ) as tmp:

        repaired_vcf = os.path.join(
            tmp,
            "repaired.vcf",
        )

        normalized_vcf = os.path.join(
            tmp,
            "normalized.vcf.gz",
        )

        sorted_vcf = os.path.join(
            tmp,
            "sorted.vcf.gz",
        )

        # ---------------------------------------------------------
        # Read header
        # ---------------------------------------------------------

        header = []

        with open_vcf(input_vcf) as f:

            for line in f:

                if line.startswith("##"):
                    header.append(
                        line.rstrip("\n")
                    )

                elif line.startswith("#CHROM"):
                    header.append(
                        line.rstrip("\n")
                    )
                    break

        (
            missing_info,
            missing_format,
            missing_filter,
            missing_contigs,
        ) = collect_missing_metadata(
            input_vcf,
            header,
        )

        print(
            f"Missing INFO definitions: "
            f"{len(missing_info)}"
        )

        print(
            f"Missing FORMAT definitions: "
            f"{len(missing_format)}"
        )

        print(
            f"Missing FILTER definitions: "
            f"{len(missing_filter)}"
        )

        print(
            f"Missing contig definitions: "
            f"{len(missing_contigs)}"
        )

        # ---------------------------------------------------------
        # Split #CHROM from metadata
        # ---------------------------------------------------------

        chrom_header = header.pop()

        # ---------------------------------------------------------
        # Add missing INFO definitions
        # ---------------------------------------------------------

        for key, values in sorted(
            missing_info.items()
        ):

            if key not in REMOVE_INFO:

                header.append(
                    infer_info_definition(
                        key,
                        values,
                    )
                )

        # ---------------------------------------------------------
        # Add missing FORMAT definitions
        # ---------------------------------------------------------

        for key, values in sorted(
            missing_format.items()
        ):

            if key not in REMOVE_FORMAT:

                header.append(
                    infer_format_definition(
                        key,
                        values,
                    )
                )

        # ---------------------------------------------------------
        # Add missing FILTER definitions
        # ---------------------------------------------------------

        for key in sorted(missing_filter):

            header.append(
                infer_filter_definition(key)
            )

        # ---------------------------------------------------------
        # Add missing contigs WITH LENGTH
        # ---------------------------------------------------------

        for chrom in sorted(missing_contigs):

            if chrom not in contig_lengths:

                sys.exit(
                    f"ERROR: Contig {chrom} is present in the VCF "
                    f"but not found in reference FAI: "
                    f"{reference_fai}"
                )

            header.append(
                f"##contig=<ID={chrom},"
                f"length={contig_lengths[chrom]}>"
            )

        header.append(chrom_header)

        # ---------------------------------------------------------
        # Remove conflicting original INFO/FORMAT definitions
        # ---------------------------------------------------------

        filtered_header = []

        for line in header:

            remove = False

            if line.startswith("##INFO=<ID="):

                m = re.match(
                    r"##INFO=<ID=([^,>]+)",
                    line,
                )

                if m and m.group(1) in REMOVE_INFO:
                    remove = True

            elif line.startswith("##FORMAT=<ID="):

                m = re.match(
                    r"##FORMAT=<ID=([^,>]+)",
                    line,
                )

                if m and m.group(1) in REMOVE_FORMAT:
                    remove = True

            if not remove:
                filtered_header.append(line)

        header = filtered_header

        # ---------------------------------------------------------
        # Write cleaned VCF
        # ---------------------------------------------------------

        print(
            "Writing repaired and cleaned VCF..."
        )

        with open(
            repaired_vcf,
            "w",
        ) as out:

            for line in header:
                out.write(line + "\n")

            with open_vcf(input_vcf) as f:

                for line in f:

                    if line.startswith("#"):
                        continue

                    cleaned = clean_record(line)

                    if cleaned is not None:
                        out.write(
                            cleaned + "\n"
                        )

        # ---------------------------------------------------------
        # Compress
        # ---------------------------------------------------------

        run([
            "bcftools",
            "view",
            "-Oz",
            "-o",
            normalized_vcf,
            repaired_vcf,
        ])

        # ---------------------------------------------------------
        # Sort
        # ---------------------------------------------------------

        print("Sorting VCF...")

        run([
            "bcftools",
            "sort",
            normalized_vcf,
            "-Oz",
            "-o",
            sorted_vcf,
        ])

        # ---------------------------------------------------------
        # Copy final VCF
        # ---------------------------------------------------------

        shutil.copy2(
            sorted_vcf,
            output_vcf,
        )

        # ---------------------------------------------------------
        # Index
        # ---------------------------------------------------------

        print("Indexing VCF...")

        run([
            "tabix",
            "-f",
            "-p",
            "vcf",
            output_vcf,
        ])

    print()
    print("SUCCESS")
    print(f"Output VCF   : {output_vcf}")
    print(f"Output index : {output_vcf}.tbi")


if __name__ == "__main__":
    main()
