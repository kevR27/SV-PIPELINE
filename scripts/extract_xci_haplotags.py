#!/usr/bin/env python3
"""Extract chrX WhatsHap read haplotags for X-inactivation analysis.

The output contains only primary mapped reads with HP=1/2 and a phase-set tag.
These read-level labels are later combined with CpG-island methylation clusters.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd
import pysam


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return path.open("w", encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--chrom", default="chrX")
    p.add_argument("--output", required=True)
    p.add_argument("--summary-output", required=True)
    args = p.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    if args.chrom not in bam.references:
        raise ValueError(
            f"{args.chrom} is not present in BAM reference dictionary."
        )

    records = []
    total_primary = 0
    hp_reads = 0
    ps_reads = 0

    for read in bam.fetch(args.chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        total_primary += 1

        if not read.has_tag("HP"):
            continue
        hp = int(read.get_tag("HP"))
        if hp not in (1, 2):
            continue
        hp_reads += 1

        ps = "."
        if read.has_tag("PS"):
            ps = str(read.get_tag("PS"))
            ps_reads += 1

        records.append(
            {
                "read_name": read.query_name,
                "HP": hp,
                "PS": ps,
                "chrom": args.chrom,
                "alignment_start": int(read.reference_start),
                "alignment_end": int(read.reference_end or read.reference_start),
                "mapq": int(read.mapping_quality),
            }
        )

    bam.close()

    out = pd.DataFrame(
        records,
        columns=[
            "read_name",
            "HP",
            "PS",
            "chrom",
            "alignment_start",
            "alignment_end",
            "mapq",
        ],
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False, compression="infer")

    summary = pd.DataFrame(
        [
            {
                "chrom": args.chrom,
                "primary_mapped_reads": total_primary,
                "HP1_reads": int((out["HP"] == 1).sum()) if not out.empty else 0,
                "HP2_reads": int((out["HP"] == 2).sum()) if not out.empty else 0,
                "haplotagged_reads": hp_reads,
                "reads_with_PS": ps_reads,
                "haplotag_fraction": (
                    round(hp_reads / total_primary, 6)
                    if total_primary
                    else 0.0
                ),
            }
        ]
    )
    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] chr={args.chrom} primary={total_primary} "
        f"haplotagged={hp_reads} output={output}"
    )


if __name__ == "__main__":
    main()
