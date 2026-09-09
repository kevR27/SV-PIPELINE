```python
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

# needLR-compatible canonical human chromosomes for this pipeline.
ALLOWED_CHROMOSOMES = {
    *(f"chr{i}" for i in range(1, 23)),
    "chrX",
    "chrY",
}

MIN_SV_LENGTH = 50
EXCLUDED_SVTYPES = {"BND"}


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
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_fai(path):
    """Read a FASTA .fai index and return {chromosome: length}."""
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


def infer_value_type(value):
    if re.fullmatch(r"-?\d+", value or ""):
        return "Integer"

    if re.fullmatch(r"-?\d+(?:\.\d+)?", value or ""):
        return "Float"

    return "String"


def infer_info_definition(key, values):
    value = values[0] if values else ""
    value_type = infer_value_type(value)

    return (
        f'##INFO=<ID={key},Number=1,Type={value_type},'
        f'Description="Automatically added missing INFO field">'
    )


def infer_format_definition(key, values):
    value = values[0] if values else ""
    value_type = infer_value_type(value)

    return (
        f'##FORMAT=<ID={key},Number=1,Type={value_type},'
        f'Description="Automatically added missing FORMAT field">'
    )


def infer_filter_definition(key):
    return (
        f'##FILTER=<ID={key},'
        f'Description="Automatically added missing FILTER field">'
    )


def parse_info(info):
    """Return INFO entries as a dictionary."""
    parsed = {}

    if info in {"", "."}:
        return parsed

    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = True

    return parsed


def get_svtype(info_dict):
    value = info_dict.get("SVTYPE")

    if value is True or value in {None, ""}:
        return None

    return str(value).upper()


def get_sv_length(pos, info_dict):
    """
    Return absolute SV length.

    Priority:
      1. SVLEN
      2. END - POS

    Returns None if neither can be determined.
    """
    svlen = info_dict.get("SVLEN")

    if svlen not in {None, "", "."}:
        try:
            value = str(svlen).split(",")[0]
            return abs(int(float(value)))
        except (ValueError, TypeError):
            pass

    end = info_dict.get("END")

    if end not in {None, "", "."}:
        try:
            return abs(int(end) - int(pos))
        except (ValueError, TypeError):
            pass

    return None


def passes_needlr_filter(fields):
    """
    Apply the pre-needLR filters:

      - canonical chromosomes only
      - FILTER = PASS
      - SVTYPE != BND
      - SV length >= 50 bp
      - No MAX SV-length limit is applied
    
    Returns:
        (True, reason) if retained
        (False, reason) if filtered
    """
    chrom = fields[0]
    pos = fields[1]
    filt = fields[6]
    info = fields[7]

    if chrom not in ALLOWED_CHROMOSOMES:
        return False, "noncanonical_chromosome"

    if filt != "PASS":
        return False, "non_PASS"

    info_dict = parse_info(info)

    svtype = get_svtype(info_dict)

    if svtype in EXCLUDED_SVTYPES: #da vedere se è una buona cosa da fare 
        return False, "BND"

    sv_length = get_sv_length(pos, info_dict)

    if sv_length is None:
        return False, "missing_SV_length"

    if sv_length < MIN_SV_LENGTH:
        return False, "SV_below_50bp"


    return True, "PASS"


def collect_metadata_from_record(
    fields,
    used_info,
    used_format,
    used_filter,
    used_contigs,
):
    """Collect metadata only from retained records."""
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

                if (
                    value not in {"", "."}
                    and len(used_info[key]) < 10
                ):
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

                    if idx >= len(sample_values):
                        continue

                    value = sample_values[idx]

                    if (
                        value not in {"", "."}
                        and len(used_format[key]) < 10
                    ):
                        used_format[key].append(value)


def collect_missing_metadata(header, used_info, used_format, used_filter):
    declared_info = set()
    declared_format = set()
    declared_filter = set()
    declared_contig = set()

    for line in header:
        if line.startswith("##INFO=<ID="):
            match = re.match(r"##INFO=<ID=([^,>]+)", line)

            if match:
                declared_info.add(match.group(1))

        elif line.startswith("##FORMAT=<ID="):
            match = re.match(r"##FORMAT=<ID=([^,>]+)", line)

            if match:
                declared_format.add(match.group(1))

        elif line.startswith("##FILTER=<ID="):
            match = re.match(r"##FILTER=<ID=([^,>]+)", line)

            if match:
                declared_filter.add(match.group(1))

        elif line.startswith("##contig=<ID="):
            match = re.match(r"##contig=<ID=([^,>]+)", line)

            if match:
                declared_contig.add(match.group(1))

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

    return (
        missing_info,
        missing_format,
        missing_filter,
        declared_contig,
    )


def clean_info(info):
    if info in {"", "."}:
        return "."

    kept = []

    for item in info.split(";"):
        key = item.split("=", 1)[0]

        if key not in REMOVE_INFO:
            kept.append(item)

    return ";".join(kept) if kept else "."


def clean_record(line):
    fields = line.rstrip("\n").split("\t")

    if len(fields) < 8:
        return None

    fields[7] = clean_info(fields[7])

    if len(fields) >= 9:
        format_keys = fields[8].split(":")

        keep_indices = [
            idx
            for idx, key in enumerate(format_keys)
            if key not in REMOVE_FORMAT
        ]

        kept_format_keys = [
            format_keys[idx]
            for idx in keep_indices
        ]

        fields[8] = ":".join(kept_format_keys) if kept_format_keys else "."

        for sample_idx in range(9, len(fields)):
            sample_values = fields[sample_idx].split(":")

            kept_values = [
                sample_values[idx]
                if idx < len(sample_values)
                else "."
                for idx in keep_indices
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
            "Prepare a Jasmine merged SV VCF for needLR by "
            "repairing metadata and applying needLR-compatible "
            "SV/chromosome filters."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input Jasmine merged VCF/VCF.GZ",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output prepared VCF.GZ",
    )

    parser.add_argument(
        "--reference-fai",
        required=True,
        help="FASTA .fai index used to add contig lengths",
    )

    args = parser.parse_args()

    for tool in ("bcftools", "tabix"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: {tool} not found in PATH.")

    input_vcf = os.path.abspath(args.input)
    output_vcf = os.path.abspath(args.output)
    reference_fai = os.path.abspath(args.reference_fai)

    if not os.path.exists(input_vcf):
        sys.exit(f"ERROR: Input VCF not found: {input_vcf}")

    contig_lengths = read_fai(reference_fai)

    print(
        f"Loaded {len(contig_lengths)} contig lengths "
        f"from {reference_fai}"
    )

    print(
        "Applying filters: "
        "PASS; canonical chromosomes; "
        f"SV >= {MIN_SV_LENGTH} bp; "
        "exclude BND; no upper SV-length limit"
    )

    os.makedirs(
        os.path.dirname(output_vcf),
        exist_ok=True,
    )

    with tempfile.TemporaryDirectory(prefix="needlr_") as tmp:
        repaired_vcf = os.path.join(tmp, "repaired.vcf")
        normalized_vcf = os.path.join(tmp, "normalized.vcf.gz")
        sorted_vcf = os.path.join(tmp, "sorted.vcf.gz")

        header = []
        chrom_header = None

        with open_vcf(input_vcf) as f:
            for line in f:
                if line.startswith("##"):
                    header.append(line.rstrip("\n"))
                elif line.startswith("#CHROM"):
                    chrom_header = line.rstrip("\n")
                    break

        if chrom_header is None:
            sys.exit("ERROR: #CHROM header line not found.")

        declared_contigs = set()

        for line in header:
            if line.startswith("##contig=<ID="):
                match = re.match(r"##contig=<ID=([^,>]+)", line)

                if match:
                    declared_contigs.add(match.group(1))

        used_info = {}
        used_format = {}
        used_filter = set()
        used_contigs = set()

        filter_counts = {}
        total_records = 0
        retained_records = 0

        print("Filtering and writing temporary VCF...")

        with open(repaired_vcf, "w") as out:
            for line in header:
                out.write(line + "\n")

            with open_vcf(input_vcf) as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    fields = line.rstrip("\n").split("\t")

                    if len(fields) < 8:
                        filter_counts["malformed_record"] = (
                            filter_counts.get("malformed_record", 0) + 1
                        )
                        continue

                    total_records += 1

                    keep, reason = passes_needlr_filter(fields)

                    if not keep:
                        filter_counts[reason] = (
                            filter_counts.get(reason, 0) + 1
                        )
                        continue

                    cleaned = clean_record(line)

                    if cleaned is None:
                        filter_counts["malformed_record"] = (
                            filter_counts.get("malformed_record", 0) + 1
                        )
                        continue

                    cleaned_fields = cleaned.split("\t")

                    collect_metadata_from_record(
                        cleaned_fields,
                        used_info,
                        used_format,
                        used_filter,
                        used_contigs,
                    )

                    out.write(cleaned + "\n")
                    retained_records += 1

        (
            missing_info,
            missing_format,
            missing_filter,
            declared_contigs,
        ) = collect_missing_metadata(
            header,
            used_info,
            used_format,
            used_filter,
        )

        # Read the already-written records and rebuild the header
        # with any required missing metadata definitions.
        with open(repaired_vcf) as f:
            body = f.read().splitlines()

        body_records = [
            line
            for line in body
            if not line.startswith("#")
        ]

        new_header = []

        for line in header:
            remove = False

            if line.startswith("##INFO=<ID="):
                match = re.match(r"##INFO=<ID=([^,>]+)", line)

                if match and match.group(1) in REMOVE_INFO:
                    remove = True

            elif line.startswith("##FORMAT=<ID="):
                match = re.match(r"##FORMAT=<ID=([^,>]+)", line)

                if match and match.group(1) in REMOVE_FORMAT:
                    remove = True

            elif line.startswith("##contig=<ID="):
                match = re.match(r"##contig=<ID=([^,>]+)", line)

                if match and match.group(1) not in used_contigs:
                    remove = True

            if not remove:
                new_header.append(line)

        for key, values in sorted(missing_info.items()):
            if key not in REMOVE_INFO:
                new_header.append(
                    infer_info_definition(key, values)
                )

        for key, values in sorted(missing_format.items()):
            if key not in REMOVE_FORMAT:
                new_header.append(
                    infer_format_definition(key, values)
                )

        for key in sorted(missing_filter):
            new_header.append(infer_filter_definition(key))

        for chrom in sorted(used_contigs):
            if chrom in declared_contigs:
                continue

            if chrom not in contig_lengths:
                sys.exit(
                    f"ERROR: Contig {chrom} is present in the filtered "
                    f"VCF but not found in reference FAI: "
                    f"{reference_fai}"
                )

            new_header.append(
                f"##contig=<ID={chrom},"
                f"length={contig_lengths[chrom]}>"
            )

        new_header.append(chrom_header)

        with open(repaired_vcf, "w") as out:
            for line in new_header:
                out.write(line + "\n")

            for line in body_records:
                out.write(line + "\n")

        print()
        print("Filtering summary")
        print("-----------------")
        print(f"Input SV records      : {total_records:,}")
        print(f"Retained SV records   : {retained_records:,}")
        print(
            f"Removed SV records    : "
            f"{total_records - retained_records:,}"
        )

        if total_records:
            percentage = (retained_records / total_records) * 100
            print(f"Retention             : {percentage:.2f}%")

        print()
        print("Removal reasons:")

        if filter_counts:
            for reason, count in sorted(filter_counts.items()):
                print(f"  {reason:<24} {count:,}")
        else:
            print("  None")

        print()
        print("Compressing VCF...")

        run([
            "bcftools",
            "view",
            "-Oz",
            "-o",
            normalized_vcf,
            repaired_vcf,
        ])

        print("Sorting VCF...")

        run([
            "bcftools",
            "sort",
            normalized_vcf,
            "-Oz",
            "-o",
            sorted_vcf,
        ])

        shutil.copy2(sorted_vcf, output_vcf)

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
    print(f"Prepared VCF : {output_vcf}")
    print(f"Index        : {output_vcf}.tbi")
    print(f"SV records   : {retained_records:,}")


if __name__ == "__main__":
    main()
```
