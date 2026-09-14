#!/usr/bin/env python3

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def open_text(path):
    """Open plain or gzipped text file."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def run_command(cmd):
    """Run command and stop on error."""
    print("\n[CMD]", " ".join(cmd))
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print("\nERROR: command failed.")
        sys.exit(result.returncode)


def get_header_info(vcf):
    """Read VCF header and determine number of declared samples."""
    samples = []

    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n\r").split("\t")

                if len(fields) < 8:
                    raise ValueError(
                        "Invalid VCF header: fewer than 8 columns."
                    )

                # FORMAT is column 9, samples start at column 10
                if len(fields) >= 10:
                    samples = fields[9:]

                return fields, samples

            if not line.startswith("#"):
                break

    raise ValueError("Could not find #CHROM header line.")


def inspect_vcf(vcf, expected_columns):
    """
    Check all records.

    Returns:
        total_records
        malformed_records
        column_distribution
    """

    total = 0
    malformed = 0
    distribution = {}

    with open_text(vcf) as fh:
        for lineno, line in enumerate(fh, start=1):

            if line.startswith("#"):
                continue

            total += 1

            fields = line.rstrip("\n\r").split("\t")
            nf = len(fields)

            distribution[nf] = distribution.get(nf, 0) + 1

            if nf != expected_columns:
                malformed += 1

    return total, malformed, distribution


def repair_vcf(input_vcf, output_vcf, expected_columns, n_samples):
    """
    Repair records with missing sample columns.

    Missing samples are filled with ./.:. etc., depending on FORMAT.
    Existing columns are never modified.
    """

    print("\n[REPAIR] Creating corrected VCF:")
    print(output_vcf)

    with open_text(input_vcf) as fin, gzip.open(output_vcf, "wt") as fout:

        for line in fin:

            if line.startswith("#"):
                fout.write(line)
                continue

            line = line.rstrip("\n\r")
            fields = line.split("\t")

            # Already correct
            if len(fields) == expected_columns:
                fout.write(line + "\n")
                continue

            # We only repair records with missing sample columns.
            if len(fields) < 9:
                raise ValueError(
                    f"Invalid VCF record with fewer than 9 columns:\n{line}"
                )

            format_field = fields[8]

            if n_samples == 0:
                raise ValueError(
                    "Header contains no samples, but records require repair."
                )

            # Number of FORMAT subfields
            n_format_fields = len(format_field.split(":"))

            # Missing sample representation
            missing_sample = ":".join(
                ["."] * n_format_fields
            )

            current_samples = len(fields) - 9

            missing_samples = n_samples - current_samples

            if missing_samples < 0:
                raise ValueError(
                    "Record contains MORE samples than declared in header:\n"
                    + line
                )

            fields.extend([missing_sample] * missing_samples)

            if len(fields) != expected_columns:
                raise ValueError(
                    f"Could not repair record at {fields[0]}:{fields[1]} "
                    f"(resulting columns={len(fields)}, "
                    f"expected={expected_columns})"
                )

            fout.write("\t".join(fields) + "\n")


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Validate and, if necessary, repair a VCF by adding missing "
            "sample columns. The original VCF is never modified."
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input VCF or VCF.GZ"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output corrected VCF.GZ"
    )

    parser.add_argument(
        "--bcftools",
        default="bcftools",
        help="Path to bcftools (default: bcftools)"
    )

    args = parser.parse_args()

    input_vcf = Path(args.input)
    output_vcf = Path(args.output)

    if not input_vcf.exists():
        print(f"ERROR: input VCF does not exist: {input_vcf}")
        sys.exit(1)

    if output_vcf.exists():
        print(f"ERROR: output already exists: {output_vcf}")
        print("Choose another output name or remove the existing file.")
        sys.exit(1)

    # ------------------------------------------------------------
    # STEP 1: Read header
    # ------------------------------------------------------------

    print("=" * 70)
    print("VCF VALIDATION / REPAIR FOR AnnotSV")
    print("=" * 70)

    print(f"\nInput: {input_vcf}")

    header, samples = get_header_info(input_vcf)

    n_samples = len(samples)

    print("\n[HEADER]")
    print(f"Number of samples declared: {n_samples}")

    if samples:
        for i, sample in enumerate(samples, start=1):
            print(f"  Sample {i}: {sample}")
    else:
        print("  No samples declared.")

    # VCF mandatory columns = 8
    # FORMAT = column 9
    # Each sample adds one column
    if n_samples > 0:
        expected_columns = 9 + n_samples
    else:
        expected_columns = 8

    print(f"\nExpected columns per record: {expected_columns}")

    # ------------------------------------------------------------
    # STEP 2: Validate records
    # ------------------------------------------------------------

    print("\n[CHECKING RECORDS]")

    total, malformed, distribution = inspect_vcf(
        input_vcf,
        expected_columns
    )

    print(f"Total variant records: {total}")

    print("\nColumn distribution:")

    for nf in sorted(distribution):
        print(f"  NF={nf}: {distribution[nf]}")

    # ------------------------------------------------------------
    # STEP 3: Already valid
    # ------------------------------------------------------------

    if malformed == 0:

        print("\nSTATUS: VCF structure is already valid.")
        print("No repair is necessary.")

        # Copy to requested output.
        if str(input_vcf).endswith(".gz"):
            shutil.copy2(input_vcf, output_vcf)
        else:
            with open(input_vcf, "rb") as fin:
                with gzip.open(output_vcf, "wb") as fout:
                    shutil.copyfileobj(fin, fout)

    # ------------------------------------------------------------
    # STEP 4: Repair
    # ------------------------------------------------------------

    else:

        print("\nSTATUS: VCF contains malformed records.")
        print(f"Records requiring repair: {malformed}")

        if n_samples == 0:
            print(
                "\nERROR: Cannot automatically repair because the VCF "
                "header declares no samples."
            )
            sys.exit(1)

        print(
            "\nRepair strategy:"
            "\n  - Existing columns will NOT be modified."
            "\n  - Missing sample columns will be added."
            "\n  - Missing FORMAT fields will be represented by '.'."
        )

        repair_vcf(
            input_vcf,
            output_vcf,
            expected_columns,
            n_samples
        )

    # ------------------------------------------------------------
    # STEP 5: Validate resulting VCF with bcftools
    # ------------------------------------------------------------

    print("\n" + "=" * 70)
    print("VALIDATING CORRECTED VCF")
    print("=" * 70)

    run_command([
        args.bcftools,
        "view",
        "-h",
        str(output_vcf)
    ])

    # Validate entire file
    print("\n[bcftools validation]")

    result = subprocess.run(
        [
            args.bcftools,
            "view",
            "-Ov",
            "-o",
            "/dev/null",
            str(output_vcf)
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        print("\nERROR: Corrected VCF still fails bcftools validation.")
        print(result.stderr)
        sys.exit(1)

    print("bcftools validation: OK")

    # ------------------------------------------------------------
    # STEP 6: Index
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # STEP 6: Sort and index
    # ------------------------------------------------------------

    print("\n[ SORTING VCF ]")

    sorted_vcf = output_vcf.with_name(
        output_vcf.stem + ".sorted.vcf.gz"
    )

    run_command([
        args.bcftools,
        "sort",
        "-Oz",
        "-o",
        str(sorted_vcf),
        str(output_vcf)
    ])

    print("\n[SORTING] Sorted VCF:")
    print(sorted_vcf)

    print("\n[INDEXING]")

    run_command([
        args.bcftools,
        "index",
        str(sorted_vcf)
    ])

    print("\n" + "=" * 70)
    print("SUCCESS")
    print("=" * 70)

    print("\nCorrected and sorted VCF:")
    print(sorted_vcf)

    print("\nIndex:")
    print(str(sorted_vcf) + ".csi")

    print("\nThe corrected VCF is ready for AnnotSV.")
  

if __name__ == "__main__":
    main()
