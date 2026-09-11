#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile


# ============================================================
# Fields to remove because they are not needed by needLR or
# may conflict with the merged VCF representation.
# ============================================================

REMOVE_INFO = {
    "AF",
    "DR",
    "DV",
}

REMOVE_FORMAT = {
    "DR",
    "DV",
}


# ============================================================
# Canonical human chromosomes accepted by this pipeline.
#
# This is a safety-net filter:
# if upstream processing already removed non-canonical
# chromosomes, this filter simply retains all remaining records.
# ============================================================

ALLOWED_CHROMOSOMES = {
    *(f"chr{i}" for i in range(1, 23)),
    "chrX",
    "chrY",
    "chrM",
}


# ============================================================
# Utility functions
# ============================================================

def run(cmd):
    """Run an external command and stop on failure."""

    print("[CMD]", " ".join(cmd), flush=True)

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:

        if result.stderr:
            print(
                result.stderr,
                file=sys.stderr,
            )

        sys.exit(
            f"ERROR: command failed with exit code "
            f"{result.returncode}"
        )

    if result.stderr:
        print(
            result.stderr,
            file=sys.stderr,
        )


def open_vcf(path):
    """Open plain-text or gzip-compressed VCF."""

    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "r")


def read_fai(path):
    """
    Read a FASTA .fai index.

    Returns:
        dictionary:
            chromosome -> chromosome length
    """

    if not os.path.exists(path):
        sys.exit(
            f"ERROR: FASTA index not found: {path}"
        )

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


# ============================================================
# Header metadata inference
# ============================================================

def infer_info_definition(key, values):
    """Create a missing INFO header definition."""

    value = values[0] if values else ""

    if re.fullmatch(r"-?\d+", value or ""):

        vtype = "Integer"

    elif re.fullmatch(
        r"-?\d+(?:\.\d+)?",
        value or "",
    ):

        vtype = "Float"

    else:

        vtype = "String"

    return (
        f'##INFO=<ID={key},Number=1,Type={vtype},'
        f'Description="Automatically added missing INFO field">'
    )


def infer_format_definition(key, values):
    """Create a missing FORMAT header definition."""

    value = values[0] if values else ""

    if re.fullmatch(r"-?\d+", value or ""):

        vtype = "Integer"

    elif re.fullmatch(
        r"-?\d+(?:\.\d+)?",
        value or "",
    ):

        vtype = "Float"

    else:

        vtype = "String"

    return (
        f'##FORMAT=<ID={key},Number=1,Type={vtype},'
        f'Description="Automatically added missing FORMAT field">'
    )


def infer_filter_definition(key):
    """Create a missing FILTER header definition."""

    return (
        f'##FILTER=<ID={key},'
        f'Description="Automatically added missing FILTER field">'
    )


# ============================================================
# Scan VCF metadata
# ============================================================

def collect_missing_metadata(
    input_vcf,
    header_lines,
):
    """
    Identify INFO, FORMAT, FILTER and contig definitions used
    by the VCF but missing from the header.

    Only used metadata is collected.
    """

    declared_info = set()
    declared_format = set()
    declared_filter = set()
    declared_contig = set()

    # --------------------------------------------------------
    # Read definitions already present in the header
    # --------------------------------------------------------

    for line in header_lines:

        if line.startswith("##INFO=<ID="):

            match = re.match(
                r"##INFO=<ID=([^,>]+)",
                line,
            )

            if match:
                declared_info.add(match.group(1))

        elif line.startswith("##FORMAT=<ID="):

            match = re.match(
                r"##FORMAT=<ID=([^,>]+)",
                line,
            )

            if match:
                declared_format.add(match.group(1))

        elif line.startswith("##FILTER=<ID="):

            match = re.match(
                r"##FILTER=<ID=([^,>]+)",
                line,
            )

            if match:
                declared_filter.add(match.group(1))

        elif line.startswith("##contig=<ID="):

            match = re.match(
                r"##contig=<ID=([^,>]+)",
                line,
            )

            if match:
                declared_contig.add(match.group(1))

    # --------------------------------------------------------
    # Scan records
    # --------------------------------------------------------

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

            # ------------------------------------------------
            # FILTER
            # ------------------------------------------------

            if filt not in {".", "PASS"}:

                for value in filt.split(";"):

                    if value:
                        used_filter.add(value)

            # ------------------------------------------------
            # INFO
            # ------------------------------------------------

            if info not in {"", "."}:

                for item in info.split(";"):

                    if "=" in item:

                        key, value = item.split(
                            "=",
                            1,
                        )

                        if key not in used_info:
                            used_info[key] = []

                        if (
                            value not in {"", "."}
                            and len(used_info[key]) < 10
                        ):

                            used_info[key].append(
                                value
                            )

                    else:

                        key = item

                        if key not in used_info:
                            used_info[key] = []

            # ------------------------------------------------
            # FORMAT
            # ------------------------------------------------

            if len(fields) >= 9:

                format_field = fields[8]

                if format_field not in {"", "."}:

                    format_keys = format_field.split(":")

                    for idx, key in enumerate(
                        format_keys
                    ):

                        if key not in used_format:
                            used_format[key] = []

                        for sample in fields[9:]:

                            sample_values = sample.split(":")

                            if idx >= len(sample_values):
                                continue

                            value = sample_values[idx]

                            if (
                                value not in {"", "."}
                                and len(used_format[key]) < 10
                            ):

                                used_format[key].append(
                                    value
                                )

    # --------------------------------------------------------
    # Determine missing metadata
    # --------------------------------------------------------

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


# ============================================================
# INFO cleaning
# ============================================================

def clean_info(info):
    """Remove unwanted INFO fields."""

    if info in {"", "."}:
        return "."

    kept = []

    for item in info.split(";"):

        key = item.split("=", 1)[0]

        if key in REMOVE_INFO:
            continue

        kept.append(item)

    return ";".join(kept) if kept else "."


# ============================================================
# Record cleaning
# ============================================================

def clean_record(line):
    """
    Remove unwanted INFO and FORMAT fields from one VCF record.
    """

    fields = line.rstrip("\n").split("\t")

    if len(fields) < 8:
        return None

    # --------------------------------------------------------
    # INFO
    # --------------------------------------------------------

    fields[7] = clean_info(fields[7])

    # --------------------------------------------------------
    # FORMAT + samples
    # --------------------------------------------------------

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

        for sample_idx in range(
            9,
            len(fields),
        ):

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


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prepare a Jasmine merged SV VCF for needLR "
            "by repairing metadata, removing unwanted fields, "
            "filtering chromosomes, sorting and indexing."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input VCF or VCF.GZ",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output VCF.GZ",
    )

    parser.add_argument(
        "--reference-fai",
        required=True,
        help=(
            "FASTA .fai index corresponding to the reference "
            "used by the VCF, e.g. hg38.fa.fai"
        ),
    )

    args = parser.parse_args()

    # --------------------------------------------------------
    # Check dependencies
    # --------------------------------------------------------

    if shutil.which("bcftools") is None:

        sys.exit(
            "ERROR: bcftools not found in PATH."
        )

    if shutil.which("tabix") is None:

        sys.exit(
            "ERROR: tabix not found in PATH."
        )

    # --------------------------------------------------------
    # Resolve paths
    # --------------------------------------------------------

    input_vcf = os.path.abspath(
        args.input
    )

    output_vcf = os.path.abspath(
        args.output
    )

    reference_fai = os.path.abspath(
        args.reference_fai
    )

    if not os.path.exists(input_vcf):

        sys.exit(
            f"ERROR: Input VCF not found: "
            f"{input_vcf}"
        )

    # --------------------------------------------------------
    # Read reference chromosome lengths
    # --------------------------------------------------------

    contig_lengths = read_fai(
        reference_fai
    )

    print(
        f"Loaded {len(contig_lengths)} contig lengths "
        f"from {reference_fai}"
    )

    # --------------------------------------------------------
    # Create output directory
    # --------------------------------------------------------

    os.makedirs(
        os.path.dirname(output_vcf),
        exist_ok=True,
    )

    # --------------------------------------------------------
    # Counters
    # --------------------------------------------------------

    chromosome_filtered_records = 0
    chromosome_retained_records = 0

    # ========================================================
    # Temporary working directory
    # ========================================================

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

        # ====================================================
        # Read header
        # ====================================================

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

        if not header:

            sys.exit(
                "ERROR: VCF header could not be read."
            )

        # ====================================================
        # Collect missing metadata
        # ====================================================

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

        # ====================================================
        # Separate #CHROM from metadata
        # ====================================================

        chrom_header = header.pop()

        # ====================================================
        # Add missing INFO definitions
        # ====================================================

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

        # ====================================================
        # Add missing FORMAT definitions
        # ====================================================

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

        # ====================================================
        # Add missing FILTER definitions
        # ====================================================

        for key in sorted(
            missing_filter
        ):

            header.append(
                infer_filter_definition(
                    key
                )
            )

        # ====================================================
        # Determine chromosomes present in the final VCF
        #
        # Only canonical chromosomes are added to the header.
        # Non-canonical contigs are intentionally omitted.
        # ====================================================

        for chrom in sorted(
            missing_contigs
        ):

            if chrom not in ALLOWED_CHROMOSOMES:
                continue

            if chrom not in contig_lengths:

                sys.exit(
                    f"ERROR: Canonical contig {chrom} "
                    f"is present in the VCF but not found "
                    f"in reference FAI: "
                    f"{reference_fai}"
                )

            header.append(
                f"##contig=<ID={chrom},"
                f"length={contig_lengths[chrom]}>"
            )

        # ====================================================
        # Remove original conflicting INFO/FORMAT definitions
        # ====================================================

        filtered_header = []

        for line in header:

            remove = False

            if line.startswith(
                "##INFO=<ID="
            ):

                match = re.match(
                    r"##INFO=<ID=([^,>]+)",
                    line,
                )

                if (
                    match
                    and match.group(1)
                    in REMOVE_INFO
                ):

                    remove = True

            elif line.startswith(
                "##FORMAT=<ID="
            ):

                match = re.match(
                    r"##FORMAT=<ID=([^,>]+)",
                    line,
                )

                if (
                    match
                    and match.group(1)
                    in REMOVE_FORMAT
                ):

                    remove = True

            # ------------------------------------------------
            # Remove original contig declarations for
            # non-canonical chromosomes.
            # ------------------------------------------------

            elif line.startswith(
                "##contig=<ID="
            ):

                match = re.match(
                    r"##contig=<ID=([^,>]+)",
                    line,
                )

                if (
                    match
                    and match.group(1)
                    not in ALLOWED_CHROMOSOMES
                ):

                    remove = True

            if not remove:

                filtered_header.append(
                    line
                )

        header = filtered_header

        # ====================================================
        # Add #CHROM header as the final header line
        # ====================================================

        header.append(
            chrom_header
        )

        # ====================================================
        # Write repaired and filtered VCF
        # ====================================================

        print(
            "Writing repaired and cleaned VCF..."
        )

        with open(
            repaired_vcf,
            "w",
        ) as out:

            # ------------------------------------------------
            # Write header
            # ------------------------------------------------

            for line in header:

                out.write(
                    line + "\n"
                )

            # ------------------------------------------------
            # Read records
            # ------------------------------------------------

            with open_vcf(input_vcf) as f:

                for line in f:

                    if line.startswith("#"):
                        continue

                    fields = line.rstrip(
                        "\n"
                    ).split("\t")

                    if len(fields) < 8:
                        continue

                    chrom = fields[0]

                    # ========================================
                    # Safety-net chromosome filter
                    #
                    # Keep only:
                    #   chr1 - chr22
                    #   chrX
                    #   chrY
                    #   chrM
                    #
                    # If upstream filtering has already removed
                    # non-canonical chromosomes, this block
                    # removes nothing.
                    # ========================================

                    if chrom not in ALLOWED_CHROMOSOMES:

                        chromosome_filtered_records += 1

                        continue

                    chromosome_retained_records += 1

                    # ----------------------------------------
                    # Clean INFO/FORMAT
                    # ----------------------------------------

                    cleaned = clean_record(
                        line
                    )

                    if cleaned is not None:

                        out.write(
                            cleaned + "\n"
                        )

        # ====================================================
        # Compress
        # ====================================================

        print(
            "Compressing VCF..."
        )

        run([
            "bcftools",
            "view",
            "-Oz",
            "-o",
            normalized_vcf,
            repaired_vcf,
        ])

        # ====================================================
        # Sort
        # ====================================================

        print(
            "Sorting VCF..."
        )

        run([
            "bcftools",
            "sort",
            normalized_vcf,
            "-Oz",
            "-o",
            sorted_vcf,
        ])

        # ====================================================
        # Copy final VCF
        # ====================================================

        shutil.copy2(
            sorted_vcf,
            output_vcf,
        )

        # ====================================================
        # Index final VCF
        # ====================================================

        print(
            "Indexing VCF..."
        )

        run([
            "tabix",
            "-f",
            "-p",
            "vcf",
            output_vcf,
        ])

    # ========================================================
    # Summary
    # ========================================================

    print()
    print(
        "Chromosome filtering"
    )
    print(
        "--------------------"
    )

    print(
        f"Retained canonical records : "
        f"{chromosome_retained_records:,}"
    )

    print(
        f"Removed non-canonical records : "
        f"{chromosome_filtered_records:,}"
    )

    print()
    print(
        "SUCCESS"
    )

    print(
        f"Output VCF   : "
        f"{output_vcf}"
    )

    print(
        f"Output index : "
        f"{output_vcf}.tbi"
    )


if __name__ == "__main__":
    main()
