#!/usr/bin/env python3
"""Compare chrX phasing between WhatsHap and LongPhase.

WhatsHap HP/PS tags are the primary read-level labels for the XCI assay.
LongPhase is used as an independent phase-consistency layer. Phase labels may
flip between blocks, so concordance is calculated after allowing a block-pair
orientation flip.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pysam


def resolve_chrom(vcf: pysam.VariantFile, requested: str) -> str:
    contigs = set(vcf.header.contigs)
    candidates = [requested]
    if requested.startswith("chr"):
        candidates.append(requested[3:])
    else:
        candidates.append("chr" + requested)
    for candidate in candidates:
        if candidate in contigs:
            return candidate
    raise ValueError(f"{requested} not found in VCF contigs")


def phased_hets(path: str, requested_chrom: str) -> pd.DataFrame:
    rows = []
    with pysam.VariantFile(path) as vcf:
        chrom = resolve_chrom(vcf, requested_chrom)
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError(f"No sample genotype in {path}")
        sample = samples[0]

        for rec in vcf.fetch(chrom):
            if len(rec.alts or ()) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue

            call = rec.samples[sample]
            gt = call.get("GT")
            if not call.phased or gt not in {(0, 1), (1, 0)}:
                continue

            ps = call.get("PS")
            if ps in {None, "."}:
                # A phased record without PS cannot define a comparable block.
                continue

            rows.append(
                {
                    "CHROM": chrom,
                    "POS": int(rec.pos),
                    "REF": rec.ref,
                    "ALT": rec.alts[0],
                    "GT": f"{gt[0]}|{gt[1]}",
                    "PS": str(ps),
                }
            )

    return pd.DataFrame(
        rows,
        columns=["CHROM", "POS", "REF", "ALT", "GT", "PS"],
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--whatshap-vcf", required=True)
    p.add_argument("--longphase-vcf", required=True)
    p.add_argument("--chromosome", default="chrX")
    p.add_argument("--min-variants-per-block-pair", type=int, default=2)
    p.add_argument("--detail-output", required=True)
    p.add_argument("--summary-output", required=True)
    args = p.parse_args()

    wh = phased_hets(args.whatshap_vcf, args.chromosome).rename(
        columns={"GT": "WHATSHAP_GT", "PS": "WHATSHAP_PS"}
    )
    lp = phased_hets(args.longphase_vcf, args.chromosome).rename(
        columns={"GT": "LONGPHASE_GT", "PS": "LONGPHASE_PS"}
    )

    keys = ["CHROM", "POS", "REF", "ALT"]
    shared = wh.merge(lp, on=keys, how="inner")
    if not shared.empty:
        shared["RAW_ORIENTATION"] = (
            shared["WHATSHAP_GT"] == shared["LONGPHASE_GT"]
        ).map({True: "SAME", False: "FLIPPED"})

    block_rows = []
    if not shared.empty:
        for (wh_ps, lp_ps), group in shared.groupby(
            ["WHATSHAP_PS", "LONGPHASE_PS"],
            sort=False,
        ):
            same = int((group["RAW_ORIENTATION"] == "SAME").sum())
            flipped = int((group["RAW_ORIENTATION"] == "FLIPPED").sum())
            total = len(group)
            chosen = "SAME" if same >= flipped else "FLIPPED"
            concordant = max(same, flipped)
            block_rows.append(
                {
                    "WHATSHAP_PS": wh_ps,
                    "LONGPHASE_PS": lp_ps,
                    "SHARED_PHASED_HET_SNPS": total,
                    "SAME_ORIENTATION_SNPS": same,
                    "FLIPPED_ORIENTATION_SNPS": flipped,
                    "BEST_BLOCK_ORIENTATION": chosen,
                    "CONCORDANT_AFTER_ALLOWING_FLIP": concordant,
                    "PHASE_CONCORDANCE": (
                        round(concordant / total, 6) if total else 0.0
                    ),
                    "BLOCK_PAIR_INFORMATIVE": (
                        "YES"
                        if total >= args.min_variants_per_block_pair
                        else "NO"
                    ),
                    "BLOCK_START": int(group["POS"].min()),
                    "BLOCK_END": int(group["POS"].max()),
                }
            )

    detail = pd.DataFrame(
        block_rows,
        columns=[
            "WHATSHAP_PS", "LONGPHASE_PS", "SHARED_PHASED_HET_SNPS",
            "SAME_ORIENTATION_SNPS", "FLIPPED_ORIENTATION_SNPS",
            "BEST_BLOCK_ORIENTATION", "CONCORDANT_AFTER_ALLOWING_FLIP",
            "PHASE_CONCORDANCE", "BLOCK_PAIR_INFORMATIVE",
            "BLOCK_START", "BLOCK_END",
        ],
    )

    informative = (
        detail[detail["BLOCK_PAIR_INFORMATIVE"].eq("YES")]
        if not detail.empty
        else detail
    )
    denom = (
        int(informative["SHARED_PHASED_HET_SNPS"].sum())
        if not informative.empty
        else 0
    )
    numer = (
        int(informative["CONCORDANT_AFTER_ALLOWING_FLIP"].sum())
        if not informative.empty
        else 0
    )
    concordance = numer / denom if denom else None

    if shared.empty:
        status = "NO_SHARED_PHASED_HETEROZYGOUS_CHRX_SNPS"
    elif informative.empty:
        status = "NO_INFORMATIVE_SHARED_PHASE_BLOCKS"
    elif concordance is not None and concordance < 0.90:
        status = "PHASING_METHOD_DISCORDANCE_REVIEW"
    else:
        status = "PHASING_METHODS_CONCORDANT"

    summary = pd.DataFrame(
        [
            {
                "CHROMOSOME": args.chromosome,
                "WHATSHAP_PHASED_HET_SNPS": len(wh),
                "LONGPHASE_PHASED_HET_SNPS": len(lp),
                "SHARED_PHASED_HET_SNPS": len(shared),
                "INFORMATIVE_BLOCK_PAIRS": len(informative),
                "WEIGHTED_PHASE_CONCORDANCE": (
                    round(concordance, 6) if concordance is not None else "."
                ),
                "PHASE_CONCORDANCE_STATUS": status,
                "INTERPRETATION": (
                    "WhatsHap provides read HP/PS labels for XCI. LongPhase is "
                    "an independent phase-consistency check; block orientation "
                    "flips are allowed before concordance is calculated."
                ),
            }
        ]
    )

    detail_path = Path(args.detail_output)
    summary_path = Path(args.summary_output)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(detail_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] shared_chrX_phased_snps={len(shared)} "
        f"informative_block_pairs={len(informative)} status={status}"
    )


if __name__ == "__main__":
    main()
