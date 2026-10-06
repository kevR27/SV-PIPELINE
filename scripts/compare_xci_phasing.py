#!/usr/bin/env python3
"""Compare WhatsHap and LongPhase chrX phasing with block-wise flip tolerance.

Haplotype labels are arbitrary between phase blocks and between phasing tools.
For each overlapping WhatsHap/LongPhase phase-block pair, the script therefore
reports the better of SAME versus FLIPPED orientation rather than requiring
identical 0|1 / 1|0 labels globally.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


def first_sample(vcf: pysam.VariantFile) -> str:
    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError("VCF has no sample column.")
    return samples[0]


def phased_chrX_snvs(path: str, chrom: str):
    vcf = pysam.VariantFile(path)
    sample_name = first_sample(vcf)
    records = {}

    for rec in vcf.fetch(chrom):
        if len(rec.ref) != 1 or len(rec.alts or ()) != 1 or len(rec.alts[0]) != 1:
            continue
        sample = rec.samples[sample_name]
        gt = sample.get("GT")
        if not sample.phased or gt not in ((0, 1), (1, 0)):
            continue
        ps = sample.get("PS")
        if ps is None:
            ps = rec.pos
        alt_haplotype = 2 if gt == (0, 1) else 1
        key = (rec.chrom, int(rec.pos), rec.ref, rec.alts[0])
        records[key] = {
            "PS": str(ps),
            "ALT_HAPLOTYPE": alt_haplotype,
        }

    vcf.close()
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--whatshap-vcf", required=True)
    p.add_argument("--longphase-vcf", required=True)
    p.add_argument("--chrom", default="chrX")
    p.add_argument("--blocks-output", required=True)
    p.add_argument("--summary-output", required=True)
    args = p.parse_args()

    wh = phased_chrX_snvs(args.whatshap_vcf, args.chrom)
    lp = phased_chrX_snvs(args.longphase_vcf, args.chrom)

    shared = sorted(set(wh) & set(lp))
    grouped = defaultdict(list)

    for key in shared:
        wh_rec = wh[key]
        lp_rec = lp[key]
        grouped[(wh_rec["PS"], lp_rec["PS"])].append(
            (
                key,
                wh_rec["ALT_HAPLOTYPE"],
                lp_rec["ALT_HAPLOTYPE"],
            )
        )

    rows = []
    total_same = 0
    total_flipped = 0
    total_oriented = 0

    for (wh_ps, lp_ps), items in grouped.items():
        same = sum(wh_h == lp_h for _, wh_h, lp_h in items)
        flipped = sum(wh_h != lp_h for _, wh_h, lp_h in items)
        orientation = "SAME" if same >= flipped else "FLIPPED"
        concordant = max(same, flipped)
        n = len(items)
        positions = [x[0][1] for x in items]

        rows.append(
            {
                "chrom": args.chrom,
                "WHATSHAP_PS": wh_ps,
                "LONGPHASE_PS": lp_ps,
                "BLOCK_START": min(positions),
                "BLOCK_END": max(positions),
                "SHARED_PHASED_SNVS": n,
                "SAME_ORIENTATION_SNVS": same,
                "FLIPPED_ORIENTATION_SNVS": flipped,
                "BEST_ORIENTATION": orientation,
                "ORIENTATION_CONCORDANCE": round(concordant / n, 6),
            }
        )
        total_same += same
        total_flipped += flipped
        total_oriented += concordant

    blocks = pd.DataFrame(
        rows,
        columns=[
            "chrom",
            "WHATSHAP_PS",
            "LONGPHASE_PS",
            "BLOCK_START",
            "BLOCK_END",
            "SHARED_PHASED_SNVS",
            "SAME_ORIENTATION_SNVS",
            "FLIPPED_ORIENTATION_SNVS",
            "BEST_ORIENTATION",
            "ORIENTATION_CONCORDANCE",
        ],
    )
    if not blocks.empty:
        blocks = blocks.sort_values(
            ["BLOCK_START", "BLOCK_END", "WHATSHAP_PS", "LONGPHASE_PS"]
        )

    blocks_path = Path(args.blocks_output)
    blocks_path.parent.mkdir(parents=True, exist_ok=True)
    blocks.to_csv(blocks_path, sep="\t", index=False)

    overall = (
        total_oriented / len(shared)
        if shared
        else np.nan
    )
    summary = pd.DataFrame(
        [
            {
                "chrom": args.chrom,
                "WHATSHAP_PHASED_SNVS": len(wh),
                "LONGPHASE_PHASED_SNVS": len(lp),
                "SHARED_PHASED_SNVS": len(shared),
                "OVERLAPPING_PHASE_BLOCK_PAIRS": len(blocks),
                "FLIP_TOLERANT_PHASE_CONCORDANCE": (
                    round(overall, 6) if np.isfinite(overall) else "."
                ),
                "INTERPRETATION": (
                    "Concordance is calculated after allowing a SAME/FLIPPED "
                    "orientation independently for each overlapping phase-block pair. "
                    "This is a phasing-consistency QC metric, not an inheritance result."
                ),
            }
        ]
    )
    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] shared_chrX_phased_snvs={len(shared)} "
        f"block_pairs={len(blocks)} concordance="
        f"{overall if np.isfinite(overall) else 'NA'}"
    )


if __name__ == "__main__":
    main()
