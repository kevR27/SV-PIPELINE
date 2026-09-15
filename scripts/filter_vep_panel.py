#!/usr/bin/env python3
"""
Filter VEP's default tab-delimited text output to only the rows whose
gene SYMBOL (parsed from the Extra column) is in a supplied gene list,
and write the result as a clean TSV with the Extra column expanded
into its own named fields (SYMBOL, IMPACT, BIOTYPE, ...).

Usage:
    python3 filter_and_clean_vep.py \
        --vep vep_output.txt \
        --genes gene_panel.txt \
        --output filtered_clean.tsv
"""

import argparse


def load_genes(path):
    with open(path) as f:
        return {
            line.strip()
            for line in f
            if line.strip() and not line.startswith("#")
        }


def parse_extra(extra_field):
    """Parse the semicolon-separated key=value Extra column into a dict."""
    extra = {}
    if extra_field == "-" or extra_field == "":
        return extra
    for pair in extra_field.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            extra[key] = value
    return extra


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep", required=True, help="VEP tab-delimited text output")
    parser.add_argument("--genes", required=True, help="Gene list, one symbol per line")
    parser.add_argument("--output", required=True, help="Filtered clean TSV output path")
    args = parser.parse_args()

    genes = load_genes(args.genes)

    base_fields = None
    extra_idx = None
    extra_keys = []
    extra_keys_seen = set()
    matched_rows = []  # list of (fields, extra_dict) for genes of interest

    # single pass over the raw VEP file: keep only rows matching the gene
    # list, while tracking every Extra key seen among the matched rows
    with open(args.vep) as fin:
        for line in fin:
            if line.startswith("##"):
                continue

            if line.startswith("#"):
                base_fields = line.lstrip("#").rstrip("\n").split("\t")
                if "Extra" in base_fields:
                    extra_idx = base_fields.index("Extra")
                continue

            if base_fields is None or extra_idx is None:
                # no header seen yet, or file has no Extra column -> can't filter
                continue

            fields = line.rstrip("\n").split("\t")
            if extra_idx >= len(fields):
                continue

            extra_dict = parse_extra(fields[extra_idx])
            symbol = extra_dict.get("SYMBOL")

            if symbol is None or symbol not in genes:
                continue

            for key in extra_dict:
                if key not in extra_keys_seen:
                    extra_keys_seen.add(key)
                    extra_keys.append(key)

            matched_rows.append((fields, extra_dict))

    if base_fields is None:
        raise SystemExit("No header line found in input file.")

    # base columns minus Extra, since Extra gets expanded into extra_keys
    out_base_fields = [f for i, f in enumerate(base_fields) if i != extra_idx]

    with open(args.output, "w") as fout:
        fout.write("\t".join(out_base_fields + extra_keys) + "\n")

        for fields, extra_dict in matched_rows:
            out_fields = [f for i, f in enumerate(fields) if i != extra_idx]
            out_fields += [extra_dict.get(key, "") for key in extra_keys]
            fout.write("\t".join(out_fields) + "\n")


if __name__ == "__main__":
    main()
