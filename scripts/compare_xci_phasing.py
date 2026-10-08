#!/usr/bin/env python3
"""Compare WhatsHap and LongPhase chrX phasing with block-wise flip tolerance.

Haplotype labels are arbitrary between phase blocks and between phasing tools.
For each overlapping WhatsHap/LongPhase phase-block pair, concordance is
therefore calculated after allowing either SAME or FLIPPED orientation.

The output also reports the actual WhatsHap/LongPhase block boundaries,
their genomic intersection, variant counts per overlap, and chrX
heterozygosity counts for XCI QC.
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


def chrX_variants(path: str, chrom: str):
    """Return phased het SNVs plus chrX heterozygosity/phasing counts."""
    vcf = pysam.VariantFile(path)
    sample_name = first_sample(vcf)

    phased = {}
    het_snvs = 0
    phased_het_snvs = 0
    block_positions = defaultdict(list)

    for rec in vcf.fetch(chrom):
        if len(rec.ref) != 1 or len(rec.alts or ()) != 1 or len(rec.alts[0]) != 1:
            continue

        sample = rec.samples[sample_name]
        gt = sample.get("GT")
        if gt is None or len(gt) != 2 or None in gt or gt[0] == gt[1]:
            continue

        het_snvs += 1
        if not sample.phased or gt not in ((0, 1), (1, 0)):
            continue

        phased_het_snvs += 1
        ps = sample.get("PS")
        if ps is None:
            ps = rec.pos
        ps = str(ps)

        alt_haplotype = 2 if gt == (0, 1) else 1
        key = (rec.chrom, int(rec.pos), rec.ref, rec.alts[0])
        phased[key] = {
            "PS": ps,
            "ALT_HAPLOTYPE": alt_haplotype,
        }
        block_positions[ps].append(int(rec.pos))

    vcf.close()

    block_bounds = {
        ps: (min(pos), max(pos), len(pos))
        for ps, pos in block_positions.items()
    }
    return phased, block_bounds, het_snvs, phased_het_snvs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--whatshap-vcf", required=True)
    p.add_argument("--longphase-vcf", required=True)
    p.add_argument("--chrom", default="chrX")
    p.add_argument("--blocks-output", required=True)
    p.add_argument("--summary-output", required=True)
    args = p.parse_args()

    wh, wh_bounds, wh_het, wh_phased_het = chrX_variants(
        args.whatshap_vcf, args.chrom
    )
    lp, lp_bounds, lp_het, lp_phased_het = chrX_variants(
        args.longphase_vcf, args.chrom
    )

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
    total_oriented = 0

    for (wh_ps, lp_ps), items in grouped.items():
        same = sum(wh_h == lp_h for _, wh_h, lp_h in items)
        flipped = sum(wh_h != lp_h for _, wh_h, lp_h in items)
        orientation = "SAME" if same >= flipped else "FLIPPED"
        concordant = max(same, flipped)
        n = len(items)

        wh_start, wh_end, wh_n = wh_bounds[wh_ps]
        lp_start, lp_end, lp_n = lp_bounds[lp_ps]
        intersection_start = max(wh_start, lp_start)
        intersection_end = min(wh_end, lp_end)

        # The pair is defined by shared phased SNVs carrying these two PS IDs,
        # so the comparison is already restricted to their overlapping
        # genomic span. Report that span explicitly for auditability.
        shared_positions = [
            x[0][1]
            for x in items
            if intersection_start <= x[0][1] <= intersection_end
        ]
        n_intersection = len(shared_positions)

        rows.append(
            {
                "chrom": args.chrom,
                "WHATSHAP_PS": wh_ps,
                "LONGPHASE_PS": lp_ps,
                "WHATSHAP_BLOCK_START": wh_start,
                "WHATSHAP_BLOCK_END": wh_end,
                "WHATSHAP_PHASED_SNVS_IN_BLOCK": wh_n,
                "LONGPHASE_BLOCK_START": lp_start,
                "LONGPHASE_BLOCK_END": lp_end,
                "LONGPHASE_PHASED_SNVS_IN_BLOCK": lp_n,
                "INTERSECTION_START": intersection_start,
                "INTERSECTION_END": intersection_end,
                "SHARED_PHASED_SNVS": n_intersection,
                "SAME_ORIENTATION_SNVS": same,
                "FLIPPED_ORIENTATION_SNVS": flipped,
                "BEST_CONCORDANT_SNVS": concordant,
                "BEST_ORIENTATION": orientation,
                "ORIENTATION_CONCORDANCE": (
                    round(concordant / n, 6) if n else np.nan
                ),
            }
        )
        total_oriented += concordant

    blocks = pd.DataFrame(rows)
    if not blocks.empty:
        blocks = blocks.sort_values(
            [
                "INTERSECTION_START",
                "INTERSECTION_END",
                "WHATSHAP_PS",
                "LONGPHASE_PS",
            ]
        )

    blocks_path = Path(args.blocks_output)
    blocks_path.parent.mkdir(parents=True, exist_ok=True)
    blocks.to_csv(blocks_path, sep="\t", index=False)

    overall = total_oriented / len(shared) if shared else np.nan
    summary = pd.DataFrame(
        [
            {
                "chrom": args.chrom,
                "WHATSHAP_CHRX_HET_SNVS": wh_het,
                "WHATSHAP_CHRX_PHASED_HET_SNVS": wh_phased_het,
                "WHATSHAP_CHRX_PHASED_HET_FRACTION": (
                    round(wh_phased_het / wh_het, 6) if wh_het else "."
                ),
                "LONGPHASE_CHRX_HET_SNVS": lp_het,
                "LONGPHASE_CHRX_PHASED_HET_SNVS": lp_phased_het,
                "LONGPHASE_CHRX_PHASED_HET_FRACTION": (
                    round(lp_phased_het / lp_het, 6) if lp_het else "."
                ),
                "SHARED_PHASED_SNVS": len(shared),
                "OVERLAPPING_PHASE_BLOCK_PAIRS": len(blocks),
                "FLIP_TOLERANT_PHASE_CONCORDANCE": (
                    round(overall, 6) if np.isfinite(overall) else "."
                ),
                "INTERPRETATION": (
                    "Concordance uses only shared phased chrX SNVs within "
                    "overlapping WhatsHap/LongPhase block spans and allows a "
                    "SAME/FLIPPED orientation independently for each block pair. "
                    "This is a phase-consistency QC metric, not an XCI direction "
                    "or inheritance result."
                ),
            }
        ]
    )

    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] whatshap_chrX_hets={wh_het} "
        f"whatshap_phased={wh_phased_het} "
        f"shared_chrX_phased_snvs={len(shared)} "
        f"block_pairs={len(blocks)} "
        f"concordance={overall if np.isfinite(overall) else 'NA'}"
    )


if __name__ == "__main__":
    main()
