#!/usr/bin/env python3
"""Calculate global folded X-inactivation skew from haplotype-block read counts.

The likelihood follows the folded-binomial model used by Gocuk et al. /
SkewX. A value of 0.5 represents balanced XCI; values closer to 0 represent
increasing skew after folding the arbitrary phase-block orientation.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp


def read_chr_coverage(path: str, chrom: str) -> float | None:
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    lower = {str(c).lower(): c for c in df.columns}
    chrom_col = lower.get("chrom")
    mean_col = lower.get("mean")
    if chrom_col is None or mean_col is None:
        return None
    hit = df[df[chrom_col].astype(str).eq(chrom)]
    if hit.empty:
        return None
    try:
        value = float(hit.iloc[0][mean_col])
        return value if np.isfinite(value) else None
    except Exception:
        return None


def folded_logpmf(x: int, n: int, p: float) -> float:
    """Log PMF of the folded binomial used by SkewX."""
    if n <= 0:
        return 0.0
    eps = 1e-12
    p = min(max(float(p), eps), 0.5)
    q = 1.0 - p

    log_choose = (
        gammaln(n + 1)
        - gammaln(x + 1)
        - gammaln(n - x + 1)
    )
    a = x * math.log(p) + (n - x) * math.log(q)
    b = (n - x) * math.log(p) + x * math.log(q)
    delta = 1 if x == n - x else 0
    return math.log(1.0 - 0.5 * delta) + log_choose + logsumexp([a, b])


def negative_log_likelihood(p: float, xs, ns) -> float:
    return -sum(
        folded_logpmf(int(x), int(n), p)
        for x, n in zip(xs, ns)
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--block-skew", required=True)
    p.add_argument("--mosdepth-summary", required=True)
    p.add_argument("--haplotag-summary", required=True)
    p.add_argument("--phase-summary", required=True)
    p.add_argument("--chrom", default="chrX")
    p.add_argument("--min-chrx-coverage", type=float, default=15.0)
    p.add_argument("--min-block-reads", type=int, default=5)
    p.add_argument("--neutral-threshold", type=float, default=0.40)
    p.add_argument("--high-skew-threshold", type=float, default=0.20)
    p.add_argument("--blocks-output", required=True)
    p.add_argument("--summary-output", required=True)
    args = p.parse_args()

    blocks = pd.read_csv(args.block_skew, sep="\t", dtype=str)
    numeric_cols = [
        "H1_Xa",
        "H1_Xi",
        "H2_Xa",
        "H2_Xi",
        "H1_Xa_SKEW",
        "BLOCK_START",
        "BLOCK_END",
        "INFORMATIVE_CPG_ISLANDS",
        "UNIQUE_READS",
    ]
    for col in numeric_cols:
        if col in blocks.columns:
            blocks[col] = pd.to_numeric(blocks[col], errors="coerce")

    if blocks.empty:
        filtered = blocks.copy()
    else:
        for col in ("H1_Xa", "H1_Xi", "H2_Xa", "H2_Xi"):
            if col not in blocks.columns:
                blocks[col] = 0
        blocks["TRIALS"] = (
            blocks["H1_Xa"].fillna(0)
            + blocks["H1_Xi"].fillna(0)
            + blocks["H2_Xa"].fillna(0)
            + blocks["H2_Xi"].fillna(0)
        ).astype(int)
        blocks["SUCCESS_H1_XA"] = (
            blocks["H1_Xa"].fillna(0)
            + blocks["H2_Xi"].fillna(0)
        ).astype(int)
        blocks["SUCCESS_FOLDED"] = np.minimum(
            blocks["SUCCESS_H1_XA"],
            blocks["TRIALS"] - blocks["SUCCESS_H1_XA"],
        )
        blocks["FOLDED_BLOCK_SKEW"] = np.where(
            blocks["TRIALS"] > 0,
            blocks["SUCCESS_FOLDED"] / blocks["TRIALS"],
            np.nan,
        )
        blocks["BLOCK_XCI_ORIENTATION"] = np.where(
            blocks["H1_Xa_SKEW"] > 0.5,
            "H1_PREFERENTIAL_Xa",
            np.where(
                blocks["H1_Xa_SKEW"] < 0.5,
                "H2_PREFERENTIAL_Xa",
                "BALANCED_OR_UNRESOLVED",
            ),
        )
        filtered = blocks[
            blocks["TRIALS"].ge(args.min_block_reads)
        ].copy()

    coverage = read_chr_coverage(args.mosdepth_summary, args.chrom)

    hap_summary = pd.read_csv(args.haplotag_summary, sep="\t")
    haplotagged_reads = (
        int(pd.to_numeric(hap_summary.get("haplotagged_reads"), errors="coerce").fillna(0).iloc[0])
        if not hap_summary.empty and "haplotagged_reads" in hap_summary.columns
        else 0
    )

    phase_summary = pd.read_csv(args.phase_summary, sep="\t")
    phase_concordance = "."
    shared_phased = 0
    if not phase_summary.empty:
        if "FLIP_TOLERANT_PHASE_CONCORDANCE" in phase_summary.columns:
            phase_concordance = phase_summary.iloc[0][
                "FLIP_TOLERANT_PHASE_CONCORDANCE"
            ]
        if "SHARED_PHASED_SNVS" in phase_summary.columns:
            try:
                shared_phased = int(
                    float(phase_summary.iloc[0]["SHARED_PHASED_SNVS"])
                )
            except Exception:
                shared_phased = 0

    mle = np.nan
    if not filtered.empty:
        xs = filtered["SUCCESS_FOLDED"].astype(int).to_numpy()
        ns = filtered["TRIALS"].astype(int).to_numpy()
        result = minimize_scalar(
            negative_log_likelihood,
            bounds=(1e-6, 0.5),
            method="bounded",
            args=(xs, ns),
            options={"xatol": 1e-8},
        )
        if result.success:
            mle = float(result.x)

    coverage_ok = coverage is not None and coverage >= args.min_chrx_coverage
    blocks_ok = len(filtered) > 0
    estimate_ok = np.isfinite(mle)

    if not coverage_ok:
        eligibility = "LOW_CHRX_COVERAGE_REVIEW"
    elif not blocks_ok:
        eligibility = "NO_INFORMATIVE_PHASED_XCI_BLOCKS"
    elif not estimate_ok:
        eligibility = "MLE_FAILED"
    else:
        eligibility = "PASS"

    if estimate_ok:
        major = 100.0 * (1.0 - mle)
        minor = 100.0 * mle
        if mle <= args.high_skew_threshold:
            status = "HIGH_SKEW_P_LE_0.20"
        elif mle < args.neutral_threshold:
            status = "MODERATE_SKEW_0.20_TO_0.40"
        else:
            status = "BALANCED_OR_LOW_SKEW_P_GE_0.40"
        ratio = f"{major:.1f}:{minor:.1f}"
    else:
        major = np.nan
        minor = np.nan
        status = "NOT_ESTIMATED"
        ratio = "."

    filtered["GLOBAL_FOLDED_SKEW_P"] = (
        round(mle, 6) if estimate_ok else np.nan
    )

    # Per-block orientation evidence, following the likelihood-ratio idea used
    # by SkewX. The published report code defines a signed log10 odds quantity
    # whose sign is opposite to the intuitive H1-Xa direction. We retain that
    # raw value for traceability and also report an explicitly oriented value:
    #
    #   > 0 : H1 is more likely the preferential Xa (H2 preferential Xi)
    #   < 0 : H2 is more likely the preferential Xa (H1 preferential Xi)
    #
    # |log10 odds| = 1 corresponds to 10:1 orientation odds.
    if estimate_ok and not filtered.empty:
        eps = 1e-12
        p_minor = min(max(float(mle), eps), 0.5 - eps)
        successes = filtered["SUCCESS_H1_XA"].astype(float)
        trials = filtered["TRIALS"].astype(float)

        skewx_raw = (
            (2.0 * successes - trials) * math.log10(p_minor)
            + (trials - 2.0 * successes) * math.log10(1.0 - p_minor)
        )
        h1_xa_log10_odds = -skewx_raw

        filtered["SKEWX_LOG10_ORIENTATION_ODDS_RAW"] = skewx_raw
        filtered["LOG10_ODDS_H1_XA_VS_H2_XA"] = h1_xa_log10_odds

        # Convert log10 odds into an easily interpretable posterior-like
        # orientation probability under equal prior odds for the two possible
        # block orientations. This is a likelihood-based orientation
        # plausibility, not a disease probability.
        clipped = np.clip(h1_xa_log10_odds.to_numpy(dtype=float), -300, 300)
        odds = np.power(10.0, clipped)
        prob_h1_xa = odds / (1.0 + odds)

        filtered["P_H1_PREFERENTIAL_XA"] = prob_h1_xa
        filtered["P_H2_PREFERENTIAL_XA"] = 1.0 - prob_h1_xa
        filtered["PREFERRED_XA_HAPLOTYPE"] = np.where(
            prob_h1_xa > 0.5,
            "H1",
            np.where(prob_h1_xa < 0.5, "H2", "UNRESOLVED"),
        )
        filtered["PREFERRED_XI_HAPLOTYPE"] = np.where(
            prob_h1_xa > 0.5,
            "H2",
            np.where(prob_h1_xa < 0.5, "H1", "UNRESOLVED"),
        )
        filtered["ORIENTATION_PLAUSIBILITY"] = np.maximum(
            prob_h1_xa,
            1.0 - prob_h1_xa,
        )

        abs_logodds = np.abs(h1_xa_log10_odds)
        filtered["ORIENTATION_EVIDENCE"] = np.select(
            [
                abs_logodds >= 2.0,
                abs_logodds >= 1.0,
                abs_logodds >= 0.5,
            ],
            [
                "VERY_STRONG_GE_100_TO_1",
                "STRONG_GE_10_TO_1",
                "MODERATE_GE_3_TO_1",
            ],
            default="WEAK_LT_3_TO_1",
        )
    else:
        filtered["SKEWX_LOG10_ORIENTATION_ODDS_RAW"] = np.nan
        filtered["LOG10_ODDS_H1_XA_VS_H2_XA"] = np.nan
        filtered["P_H1_PREFERENTIAL_XA"] = np.nan
        filtered["P_H2_PREFERENTIAL_XA"] = np.nan
        filtered["PREFERRED_XA_HAPLOTYPE"] = "."
        filtered["PREFERRED_XI_HAPLOTYPE"] = "."
        filtered["ORIENTATION_PLAUSIBILITY"] = np.nan
        filtered["ORIENTATION_EVIDENCE"] = "NOT_ESTIMATED"

    out_blocks = Path(args.blocks_output)
    out_blocks.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(out_blocks, sep="\t", index=False)

    informative_cgis = (
        int(filtered["INFORMATIVE_CPG_ISLANDS"].fillna(0).sum())
        if "INFORMATIVE_CPG_ISLANDS" in filtered.columns
        else 0
    )
    informative_reads = (
        int(filtered["TRIALS"].sum())
        if "TRIALS" in filtered.columns
        else 0
    )

    strong_orientation_blocks = 0
    very_strong_orientation_blocks = 0
    median_abs_log10_odds = np.nan
    median_orientation_plausibility = np.nan
    if "LOG10_ODDS_H1_XA_VS_H2_XA" in filtered.columns and not filtered.empty:
        lod = pd.to_numeric(
            filtered["LOG10_ODDS_H1_XA_VS_H2_XA"],
            errors="coerce",
        )
        plaus = pd.to_numeric(
            filtered["ORIENTATION_PLAUSIBILITY"],
            errors="coerce",
        )
        if lod.notna().any():
            abs_lod = lod.abs()
            strong_orientation_blocks = int((abs_lod >= 1.0).sum())
            very_strong_orientation_blocks = int((abs_lod >= 2.0).sum())
            median_abs_log10_odds = float(abs_lod.median())
        if plaus.notna().any():
            median_orientation_plausibility = float(plaus.median())

    summary = pd.DataFrame(
        [
            {
                "chrom": args.chrom,
                "XCI_ANALYSIS_STATUS": eligibility,
                "CHRX_MEAN_COVERAGE": (
                    round(coverage, 3) if coverage is not None else "."
                ),
                "MIN_CHRX_COVERAGE": args.min_chrx_coverage,
                "HAPLOTAGGED_CHRX_READS": haplotagged_reads,
                "SHARED_WHATSHAP_LONGPHASE_PHASED_SNVS": shared_phased,
                "FLIP_TOLERANT_WHATSHAP_LONGPHASE_CONCORDANCE": phase_concordance,
                "INFORMATIVE_PHASE_BLOCKS": len(filtered),
                "INFORMATIVE_CPG_ISLAND_SUM": informative_cgis,
                "INFORMATIVE_READS": informative_reads,
                "MIN_READS_PER_BLOCK": args.min_block_reads,
                "STRONG_ORIENTATION_BLOCKS_LOG10_ODDS_GE_1": (
                    strong_orientation_blocks
                ),
                "VERY_STRONG_ORIENTATION_BLOCKS_LOG10_ODDS_GE_2": (
                    very_strong_orientation_blocks
                ),
                "MEDIAN_ABS_LOG10_ORIENTATION_ODDS": (
                    round(median_abs_log10_odds, 4)
                    if np.isfinite(median_abs_log10_odds)
                    else "."
                ),
                "MEDIAN_ORIENTATION_PLAUSIBILITY": (
                    round(median_orientation_plausibility, 4)
                    if np.isfinite(median_orientation_plausibility)
                    else "."
                ),
                "GLOBAL_FOLDED_SKEW_P": (
                    round(mle, 6) if estimate_ok else "."
                ),
                "MAJOR_X_PERCENT": (
                    round(major, 2) if estimate_ok else "."
                ),
                "MINOR_X_PERCENT": (
                    round(minor, 2) if estimate_ok else "."
                ),
                "XCI_MAJOR_MINOR_RATIO": ratio,
                "XCI_SKEW_STATUS": status,
                "HIGH_SKEW_THRESHOLD_P": args.high_skew_threshold,
                "NEUTRAL_THRESHOLD_P": args.neutral_threshold,
                "METHOD": (
                    "CpG-island Xa/Xi clustering + WhatsHap HP/PS + "
                    "folded-binomial MLE + per-block log10 orientation odds; "
                    "LongPhase used as phase-consistency QC"
                ),
                "INTERPRETATION": (
                    "P=0.5 is balanced after folding; values closer to 0 "
                    "represent stronger skew. Per-block |log10 odds| >=1 means "
                    "at least 10:1 support for one Xa/Xi orientation; >=2 means "
                    "at least 100:1. These are orientation likelihoods, not "
                    "pathogenicity probabilities. Haplotype labels remain local "
                    "unless separately linked to a disease allele."
                ),
            }
        ]
    )
    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] status={eligibility} blocks={len(filtered)} "
        f"P={mle if estimate_ok else 'NA'} ratio={ratio}"
    )


if __name__ == "__main__":
    main()
