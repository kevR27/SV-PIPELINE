#!/usr/bin/env python3

import argparse
import pysam


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    invcf = pysam.VariantFile(args.input)

    # Add INFO definitions that may be present in caller records
    # but lost from Jasmine's merged header.
    if "CIPOS" not in invcf.header.info:
        invcf.header.info.add(
            "CIPOS",
            number=2,
            type="Integer",
            description="Confidence interval around POS"
        )

    if "CILEN" not in invcf.header.info:
        invcf.header.info.add(
            "CILEN",
            number=2,
            type="Integer",
            description="Confidence interval around SV length"
        )

    if "RE" not in invcf.header.info:
        invcf.header.info.add(
            "RE",
            number=1,
            type="Integer",
            description="Number of supporting reads"
        )

    if "AF" not in invcf.header.info:
        invcf.header.info.add(
            "AF",
            number=1,
            type="Float",
            description="Allele frequency"
        )

    outvcf = pysam.VariantFile(
        args.output,
        "w",
        header=invcf.header
    )

    for record in invcf:

        svtype = record.info.get("SVTYPE")

        # For BNDs, END may represent the partner breakpoint
        # rather than a genomic interval endpoint.
        # The partner position is already encoded in ALT / CHR2.
        if svtype == "BND" and "END" in record.info:
            del record.info["END"]

        outvcf.write(record)

    outvcf.close()
    invcf.close()


if __name__ == "__main__":
    main()
