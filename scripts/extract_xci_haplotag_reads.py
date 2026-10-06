#!/usr/bin/env python3
"""Extract informative chrX read haplotags and modified-base QC from a BAM.

The XCI analysis requires read-level haplotype labels (HP/PS) and native ONT
modified-base tags (MM/ML) on the same reads. WhatsHap haplotagging is used as
the primary read-level haplotype assignment because it writes HP/PS tags while
preserving the original BAM record tags.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd
import pysam


def chromosome_name(bam: pysam.AlignmentFile, requested: str) -> str:
    refs = set(bam.references)
    candidates = [requested]
    if requested.startswith("chr"):
        candidates.append(requested[3:])
    else:
        candidates.append("chr" + requested)
    for candidate in candidates:
        if candidate in refs:
            return candidate
    raise ValueError(
        f"Chromosome {requested!r} is not present in BAM. "
        f"Examples: {', '.join(list(bam.references)[:8])}"
    )


def has_mod_tag(read) -> bool:
    return (
        read.has_tag("MM")
        or read.has_tag("Mm")
        or read.has_tag("ML")
        or read.has_tag("Ml")
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--chromosome", default="chrX")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--output", required=True)
    p.add_argument("--summary", required=True)
    args = p.parse_args()

    bam_path = Path(args.bam)
    if not bam_path.exists():
        raise FileNotFoundError(bam_path)

    # pysam/htslib requires an index for regional fetch. Create it only when
    # absent; this does not rewrite the BAM.
    bai1 = Path(str(bam_path) + ".bai")
    bai2 = bam_path.with_suffix(".bai")
    if not bai1.exists() and not bai2.exists():
        pysam.index(str(bam_path))

    out_path = Path(args.output)
    summary_path = Path(args.summary)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    counts = {
        "CHR_X_PRIMARY_READS": 0,
        "CHR_X_MAPQ_PASS_READS": 0,
        "CHR_X_MOD_TAGGED_READS": 0,
        "CHR_X_HAPLOTAGGED_READS": 0,
        "CHR_X_INFORMATIVE_XCI_READS": 0,
    }

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        chrom = chromosome_name(bam, args.chromosome)

        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            counts["CHR_X_PRIMARY_READS"] += 1
            if read.mapping_quality < args.min_mapq:
                continue

            counts["CHR_X_MAPQ_PASS_READS"] += 1
            mod_tag = has_mod_tag(read)
            if mod_tag:
                counts["CHR_X_MOD_TAGGED_READS"] += 1

            hp = read.get_tag("HP") if read.has_tag("HP") else None
            ps = read.get_tag("PS") if read.has_tag("PS") else None
            haplotagged = hp in {1, 2} and ps not in {None, 0, "."}

            if haplotagged:
                counts["CHR_X_HAPLOTAGGED_READS"] += 1

            if not (mod_tag and haplotagged):
                continue

            counts["CHR_X_INFORMATIVE_XCI_READS"] += 1
            rows.append(
                {
                    "read_name": read.query_name,
                    "HP": int(hp),
                    "PS": str(ps),
                    "CHROM": chrom,
                    "START": int(read.reference_start),
                    "END": int(read.reference_end or read.reference_start),
                    "MAPQ": int(read.mapping_quality),
                    "HAS_MODIFICATION_TAG": "YES",
                }
            )

    columns = [
        "read_name", "HP", "PS", "CHROM", "START", "END", "MAPQ",
        "HAS_MODIFICATION_TAG",
    ]
    hap = pd.DataFrame(rows, columns=columns)
    if not hap.empty:
        hap = hap.sort_values(["PS", "HP", "START", "read_name"])
        hap = hap.drop_duplicates("read_name", keep="first")

    compression = "gzip" if out_path.suffix == ".gz" else None
    hap.to_csv(out_path, sep="\t", index=False, compression=compression)

    primary = counts["CHR_X_PRIMARY_READS"]
    mapq_pass = counts["CHR_X_MAPQ_PASS_READS"]
    mod_fraction = (
        counts["CHR_X_MOD_TAGGED_READS"] / mapq_pass if mapq_pass else 0.0
    )
    hap_fraction = (
        counts["CHR_X_HAPLOTAGGED_READS"] / mapq_pass if mapq_pass else 0.0
    )

    if counts["CHR_X_INFORMATIVE_XCI_READS"] == 0:
        status = "NO_READS_WITH_BOTH_MODIFICATION_AND_HP_PS"
    elif mod_fraction < 0.10:
        status = "LOW_MODIFICATION_TAG_FRACTION_REVIEW"
    elif hap_fraction < 0.05:
        status = "LOW_CHRX_HAPLOTAG_FRACTION_REVIEW"
    else:
        status = "PASS_FOR_XCI_CLUSTERING"

    summary = pd.DataFrame(
        [
            {
                "BAM": str(bam_path),
                "CHROMOSOME": args.chromosome,
                "MIN_MAPQ": args.min_mapq,
                **counts,
                "MOD_TAGGED_FRACTION_OF_MAPQ_PASS": round(mod_fraction, 6),
                "HAPLOTAGGED_FRACTION_OF_MAPQ_PASS": round(hap_fraction, 6),
                "XCI_INPUT_STATUS": status,
                "INTERPRETATION": (
                    "XCI skew requires the same native-DNA reads to retain "
                    "methylation tags and informative HP/PS haplotags."
                ),
            }
        ]
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] XCI haplotags={len(hap)} status={status} "
        f"output={out_path}"
    )


if __name__ == "__main__":
    main()
