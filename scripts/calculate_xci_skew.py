#!/usr/bin/env python3
"""Estimate global X-chromosome inactivation skew from phased methylation blocks.

The primary estimate follows the folded-binomial likelihood used by SkewX:
HP1/HP2 orientation may flip between phase blocks, so the likelihood is
symmetric around 0.5 and estimates the minor-X fraction P in [0, 0.5].

This implementation keeps that estimator unchanged and adds:
  * 95% profile-likelihood confidence interval for P;
  * likelihood-ratio evidence against balanced XCI (P=0.5);
  * minimum-read sensitivity analysis;
  * WhatsHap/LongPhase concordance as a QC/sensitivity filter;
  * explicit per-block Xa-orientation log10 likelihood ratios;
  * chrX depth/heterozygosity QC for XX compatibility review.

The result is a research XCI estimate, not a pathogenicity classification.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.stats import chi2


def read_coverage_qc(path: str, chrom: str):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None, None, None

    lower = {str(c).lower(): c for c in df.columns}
    chrom_col = lower.get("chrom")
    mean_col = lower.get("mean")
    if chrom_col is None or mean_col is None:
        return None, None, None

    work = df[[chrom_col, mean_col]].copy()
    work[mean_col] = pd.to_numeric(work[mean_col], errors="coerce")

    hit = work[work[chrom_col].astype(str).eq(chrom)]
    x_cov = None
    if not hit.empty:
        value = hit.iloc[0][mean_col]
        if pd.notna(value) and np.isfinite(value):
            x_cov = float(value)

    autosomes = work[
        work[chrom_col].astype(str).str.match(
            r"^(chr)?([1-9]|1[0-9]|2[0-2])$"
        )
    ][mean_col].dropna()
    auto_median = float(autosomes.median()) if not autosomes.empty else None
    ratio = (
        x_cov / auto_median
        if x_cov is not None and auto_median not in (None, 0)
        else None
    )
    return x_cov, auto_median, ratio


def folded_logpmf(x: int, n: int, p_minor: float) -> float:
    """Log PMF of the folded binomial used by SkewX."""
    if n <= 0:
        return 0.0

    eps = 1e-12
    p = min(max(float(p_minor), eps), 0.5)
    q = 1.0 - p

    log_choose = (
        gammaln(n + 1)
        - gammaln(x + 1)
        - gammaln(n - x + 1)
    )
    a = x * math.log(p) + (n - x) * math.log(q)
    b = (n - x) * math.log(p) + x * math.log(q)
    delta = 1 if x == n - x else 0

    return (
        math.log(1.0 - 0.5 * delta)
        + log_choose
        + logsumexp([a, b])
    )


def negative_log_likelihood(p_minor: float, xs, ns) -> float:
    return -sum(
        folded_logpmf(int(x), int(n), p_minor)
        for x, n in zip(xs, ns)
    )


def fit_folded_mle(df: pd.DataFrame) -> float:
    if df.empty:
        return np.nan

    xs = df["SUCCESS_FOLDED"].astype(int).to_numpy()
    ns = df["TRIALS"].astype(int).to_numpy()

    result = minimize_scalar(
        negative_log_likelihood,
        bounds=(1e-6, 0.5),
        method="bounded",
        args=(xs, ns),
        options={"xatol": 1e-10},
    )
    return float(result.x) if result.success else np.nan


def profile_likelihood_ci(
    df: pd.DataFrame,
    mle: float,
    level: float = 0.95,
) -> tuple[float, float]:
    """Profile-likelihood CI using 2*DeltaNLL <= chi2_1(level)."""
    if df.empty or not np.isfinite(mle):
        return np.nan, np.nan

    xs = df["SUCCESS_FOLDED"].astype(int).to_numpy()
    ns = df["TRIALS"].astype(int).to_numpy()
    nll_min = negative_log_likelihood(mle, xs, ns)
    target = nll_min + 0.5 * chi2.ppf(level, df=1)
    eps = 1e-8

    def f(value):
        return negative_log_likelihood(value, xs, ns) - target

    if mle <= eps:
        low = eps
    elif f(eps) <= 0:
        low = eps
    else:
        low = brentq(f, eps, mle)

    upper_bound = 0.5
    if mle >= upper_bound - eps:
        high = upper_bound
    elif f(upper_bound) <= 0:
        high = upper_bound
    else:
        high = brentq(f, mle, upper_bound)

    return float(low), float(high)


def balanced_lrt(df: pd.DataFrame, mle: float):
    """LR test of balanced XCI P=0.5 versus P<0.5.

    Because P=0.5 is the boundary of the folded parameter space, the
    asymptotic null is the 50:50 chi-square(0)/chi-square(1) mixture.
    """
    if df.empty or not np.isfinite(mle):
        return np.nan, np.nan

    xs = df["SUCCESS_FOLDED"].astype(int).to_numpy()
    ns = df["TRIALS"].astype(int).to_numpy()
    nll_mle = negative_log_likelihood(mle, xs, ns)
    nll_balanced = negative_log_likelihood(0.5, xs, ns)
    statistic = max(0.0, 2.0 * (nll_balanced - nll_mle))

    if statistic <= 0:
        p_value = 1.0
    else:
        p_value = 0.5 * chi2.sf(statistic, df=1)

    return float(statistic), float(p_value)


def parse_cutoffs(text: str) -> list[int]:
    out = []
    for token in str(text).split(","):
        token = token.strip()
        if token:
            out.append(int(token))
    return sorted(set(out))


def ratio_text(p_minor: float) -> str:
    if not np.isfinite(p_minor):
        return "."
    return f"{100.0 * (1.0 - p_minor):.1f}:{100.0 * p_minor:.1f}"


def p_range_label(
    p_minor: float,
    threshold_70_30: float,
    threshold_80_20: float,
    threshold_90_10: float,
) -> str:
    if not np.isfinite(p_minor):
        return "NOT_ESTIMATED"
    if p_minor >= threshold_70_30:
        return "P_GE_0.30"
    if p_minor >= threshold_80_20:
        return "P_0.20_TO_0.30"
    if p_minor >= threshold_90_10:
        return "P_0.10_TO_0.20"
    return "P_LT_0.10"


def conventional_context(
    p_minor: float,
    threshold_70_30: float,
    threshold_80_20: float,
    threshold_90_10: float,
) -> str:
    if not np.isfinite(p_minor):
        return "NOT_ESTIMATED"
    if p_minor >= threshold_70_30:
        return "BELOW_70_30_SKEW_THRESHOLD"
    if p_minor >= threshold_80_20:
        return "AT_OR_BEYOND_70_30_BUT_BELOW_80_20"
    if p_minor >= threshold_90_10:
        return "AT_OR_BEYOND_80_20_BUT_BELOW_90_10"
    return "AT_OR_BEYOND_90_10"


def load_phase_qc(path: str) -> pd.DataFrame:
    phase = pd.read_csv(path, sep="\t", dtype=str)
    required = {"WHATSHAP_PS", "SHARED_PHASED_SNVS"}
    if phase.empty or not required.issubset(phase.columns):
        return pd.DataFrame(
            columns=[
                "PS",
                "PHASE_QC_SHARED_SNVS",
                "PHASE_QC_CONCORDANT_SNVS",
                "PHASE_QC_CONCORDANCE",
                "PHASE_QC_PAIR_COUNT",
            ]
        )

    phase["SHARED_PHASED_SNVS"] = pd.to_numeric(
        phase["SHARED_PHASED_SNVS"], errors="coerce"
    ).fillna(0)

    if "BEST_CONCORDANT_SNVS" in phase.columns:
        phase["BEST_CONCORDANT_SNVS"] = pd.to_numeric(
            phase["BEST_CONCORDANT_SNVS"], errors="coerce"
        ).fillna(0)
    else:
        phase["BEST_CONCORDANT_SNVS"] = (
            pd.to_numeric(
                phase.get("ORIENTATION_CONCORDANCE"),
                errors="coerce",
            ).fillna(0)
            * phase["SHARED_PHASED_SNVS"]
        )

    grouped = (
        phase.groupby("WHATSHAP_PS", as_index=False)
        .agg(
            PHASE_QC_SHARED_SNVS=("SHARED_PHASED_SNVS", "sum"),
            PHASE_QC_CONCORDANT_SNVS=("BEST_CONCORDANT_SNVS", "sum"),
            PHASE_QC_PAIR_COUNT=("WHATSHAP_PS", "size"),
        )
        .rename(columns={"WHATSHAP_PS": "PS"})
    )
    grouped["PS"] = grouped["PS"].astype(str)
    grouped["PHASE_QC_CONCORDANCE"] = np.where(
        grouped["PHASE_QC_SHARED_SNVS"] > 0,
        grouped["PHASE_QC_CONCORDANT_SNVS"]
        / grouped["PHASE_QC_SHARED_SNVS"],
        np.nan,
    )
    return grouped


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--block-skew", required=True)
    p.add_argument("--mosdepth-summary", required=True)
    p.add_argument("--haplotag-summary", required=True)
    p.add_argument("--phase-summary", required=True)
    p.add_argument("--phase-blocks", required=True)
    p.add_argument("--chrom", default="chrX")
    p.add_argument("--min-chrx-coverage", type=float, default=15.0)
    p.add_argument("--min-block-reads", type=int, default=5)
    p.add_argument("--phase-concordance-threshold", type=float, default=0.90)
    p.add_argument("--phase-min-shared-snvs", type=int, default=3)
    p.add_argument("--min-chrx-het-snvs", type=int, default=50)
    p.add_argument("--threshold-70-30-minor", type=float, default=0.30)
    p.add_argument("--threshold-80-20-minor", type=float, default=0.20)
    p.add_argument("--threshold-90-10-minor", type=float, default=0.10)
    p.add_argument("--orientation-min-log10-odds", type=float, default=1.0)
    p.add_argument("--sensitivity-min-reads", default="5,8,10,15")
    p.add_argument("--blocks-output", required=True)
    p.add_argument("--summary-output", required=True)
    p.add_argument("--sensitivity-output", required=True)
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
        "MULTI_CGI_READS",
        "DISCORDANT_MULTI_CGI_READS",
        "CONSENSUS_READS",
    ]
    for col in numeric_cols:
        if col in blocks.columns:
            blocks[col] = pd.to_numeric(blocks[col], errors="coerce")

    for col in ("H1_Xa", "H1_Xi", "H2_Xa", "H2_Xi"):
        if col not in blocks.columns:
            blocks[col] = 0

    if "PS" not in blocks.columns:
        blocks["PS"] = "."
    blocks["PS"] = blocks["PS"].astype(str)

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
    blocks["H1_XA_PROPORTION_RAW"] = np.where(
        blocks["TRIALS"] > 0,
        blocks["SUCCESS_H1_XA"] / blocks["TRIALS"],
        np.nan,
    )

    phase_qc = load_phase_qc(args.phase_blocks)
    blocks = blocks.merge(phase_qc, on="PS", how="left")

    blocks["PHASE_QC_SHARED_SNVS"] = pd.to_numeric(
        blocks.get("PHASE_QC_SHARED_SNVS"), errors="coerce"
    ).fillna(0)
    blocks["PHASE_QC_CONCORDANCE"] = pd.to_numeric(
        blocks.get("PHASE_QC_CONCORDANCE"), errors="coerce"
    )
    blocks["PHASE_QC_PASS"] = np.where(
        (blocks["PHASE_QC_SHARED_SNVS"] >= args.phase_min_shared_snvs)
        & (
            blocks["PHASE_QC_CONCORDANCE"]
            >= args.phase_concordance_threshold
        ),
        "YES",
        "NO",
    )

    # Primary estimator: preserve the SkewX-style folded-binomial fit across
    # all blocks meeting the read threshold. Phase concordance is QC only.
    primary = blocks[blocks["TRIALS"].ge(args.min_block_reads)].copy()
    phase_filtered = primary[primary["PHASE_QC_PASS"].eq("YES")].copy()

    mle = fit_folded_mle(primary)
    phase_filtered_mle = fit_folded_mle(phase_filtered)
    ci_low, ci_high = profile_likelihood_ci(primary, mle)
    lrt_stat, lrt_p = balanced_lrt(primary, mle)

    x_cov, auto_median, x_auto_ratio = read_coverage_qc(
        args.mosdepth_summary,
        args.chrom,
    )

    hap_summary = pd.read_csv(args.haplotag_summary, sep="\t")
    haplotagged_reads = (
        int(
            pd.to_numeric(
                hap_summary.get("haplotagged_reads"),
                errors="coerce",
            ).fillna(0).iloc[0]
        )
        if not hap_summary.empty
        and "haplotagged_reads" in hap_summary.columns
        else 0
    )

    phase_summary = pd.read_csv(args.phase_summary, sep="\t")
    phase_concordance = "."
    shared_phased = 0
    whatshap_chrX_hets = 0
    whatshap_chrX_phased_hets = 0

    if not phase_summary.empty:
        row = phase_summary.iloc[0]
        if "FLIP_TOLERANT_PHASE_CONCORDANCE" in phase_summary.columns:
            phase_concordance = row[
                "FLIP_TOLERANT_PHASE_CONCORDANCE"
            ]

        def get_int(name):
            if name not in phase_summary.columns:
                return 0
            try:
                return int(float(row[name]))
            except Exception:
                return 0

        shared_phased = get_int("SHARED_PHASED_SNVS")
        whatshap_chrX_hets = get_int("WHATSHAP_CHRX_HET_SNVS")
        whatshap_chrX_phased_hets = get_int(
            "WHATSHAP_CHRX_PHASED_HET_SNVS"
        )

    coverage_ok = (
        x_cov is not None
        and x_cov >= args.min_chrx_coverage
    )
    estimate_ok = np.isfinite(mle)

    if not coverage_ok:
        eligibility = "LOW_CHRX_COVERAGE_REVIEW"
    elif primary.empty:
        eligibility = "NO_INFORMATIVE_PHASED_XCI_BLOCKS"
    elif not estimate_ok:
        eligibility = "MLE_FAILED"
    else:
        eligibility = "PASS"

    if (
        x_auto_ratio is not None
        and 0.75 <= x_auto_ratio <= 1.25
        and whatshap_chrX_hets >= args.min_chrx_het_snvs
    ):
        xx_qc = "XX_COMPATIBLE"
    else:
        xx_qc = "REVIEW_XX_COMPATIBILITY"

    ratio = ratio_text(mle)
    major = 100.0 * (1.0 - mle) if estimate_ok else np.nan
    minor = 100.0 * mle if estimate_ok else np.nan
    range_label = p_range_label(
        mle,
        args.threshold_70_30_minor,
        args.threshold_80_20_minor,
        args.threshold_90_10_minor,
    )
    context_label = conventional_context(
        mle,
        args.threshold_70_30_minor,
        args.threshold_80_20_minor,
        args.threshold_90_10_minor,
    )

    # Per-block Xa orientation likelihood.
    if estimate_ok and not primary.empty:
        eps = 1e-12
        p_minor = min(max(float(mle), eps), 0.5 - eps)
        q_major = 1.0 - p_minor
        successes = primary["SUCCESS_H1_XA"].astype(float)
        trials = primary["TRIALS"].astype(float)

        # log10 LR(H1 major Xa vs H2 major Xa)
        log10_odds = (
            (2.0 * successes - trials)
            * math.log10(q_major / p_minor)
        )

        primary["LOG10_ODDS_H1_XA_VS_H2_XA"] = log10_odds
        primary["PREFERRED_XA_HAPLOTYPE"] = np.where(
            log10_odds > 0,
            "H1",
            np.where(log10_odds < 0, "H2", "UNRESOLVED"),
        )
        primary["PREFERRED_XI_HAPLOTYPE"] = np.where(
            log10_odds > 0,
            "H2",
            np.where(log10_odds < 0, "H1", "UNRESOLVED"),
        )
        abs_lod = np.abs(log10_odds)
        primary["ORIENTATION_EVIDENCE"] = np.select(
            [
                abs_lod >= 2.0,
                abs_lod >= 1.0,
                abs_lod >= 0.5,
            ],
            [
                "VERY_STRONG_GE_100_TO_1",
                "STRONG_GE_10_TO_1",
                "MODERATE_GE_3_TO_1",
            ],
            default="WEAK_LT_3_TO_1",
        )
        primary["ORIENTATION_USABLE"] = np.where(
            (abs_lod >= args.orientation_min_log10_odds)
            & primary["PHASE_QC_PASS"].eq("YES"),
            "YES",
            "NO",
        )
    else:
        primary["LOG10_ODDS_H1_XA_VS_H2_XA"] = np.nan
        primary["PREFERRED_XA_HAPLOTYPE"] = "."
        primary["PREFERRED_XI_HAPLOTYPE"] = "."
        primary["ORIENTATION_EVIDENCE"] = "NOT_ESTIMATED"
        primary["ORIENTATION_USABLE"] = "NO"

    primary["GLOBAL_FOLDED_SKEW_P"] = (
        round(mle, 6) if estimate_ok else np.nan
    )

    out_blocks = Path(args.blocks_output)
    out_blocks.parent.mkdir(parents=True, exist_ok=True)
    primary.to_csv(out_blocks, sep="\t", index=False)

    # Sensitivity: do not replace the primary estimate. Show how P changes
    # with stricter read cutoffs, and separately after phase-QC filtering.
    sensitivity_rows = []
    for cutoff in parse_cutoffs(args.sensitivity_min_reads):
        subset = blocks[blocks["TRIALS"].ge(cutoff)].copy()
        subset_phase = subset[
            subset["PHASE_QC_PASS"].eq("YES")
        ].copy()

        all_est = fit_folded_mle(subset)
        phase_est = fit_folded_mle(subset_phase)

        sensitivity_rows.append(
            {
                "MIN_READS_PER_BLOCK": cutoff,
                "ALL_BLOCKS_N": len(subset),
                "ALL_BLOCKS_FOLDED_MINOR_X_P": all_est,
                "ALL_BLOCKS_RATIO": ratio_text(all_est),
                "PHASE_QC_PASS_BLOCKS_N": len(subset_phase),
                "PHASE_QC_FOLDED_MINOR_X_P": phase_est,
                "PHASE_QC_RATIO": ratio_text(phase_est),
            }
        )

    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity_path = Path(args.sensitivity_output)
    sensitivity_path.parent.mkdir(parents=True, exist_ok=True)
    sensitivity.to_csv(sensitivity_path, sep="\t", index=False)

    informative_cgis = (
        int(primary["INFORMATIVE_CPG_ISLANDS"].fillna(0).sum())
        if "INFORMATIVE_CPG_ISLANDS" in primary.columns
        else 0
    )
    informative_reads = (
        int(primary["TRIALS"].sum())
        if "TRIALS" in primary.columns
        else 0
    )
    discordant_reads = (
        int(
            primary["DISCORDANT_MULTI_CGI_READS"]
            .fillna(0)
            .sum()
        )
        if "DISCORDANT_MULTI_CGI_READS" in primary.columns
        else 0
    )

    lod = pd.to_numeric(
        primary.get("LOG10_ODDS_H1_XA_VS_H2_XA"),
        errors="coerce",
    )
    strong_orientation_blocks = (
        int((lod.abs() >= 1.0).sum())
        if lod.notna().any()
        else 0
    )
    very_strong_orientation_blocks = (
        int((lod.abs() >= 2.0).sum())
        if lod.notna().any()
        else 0
    )

    summary = pd.DataFrame(
        [
            {
                "chrom": args.chrom,
                "XCI_ANALYSIS_STATUS": eligibility,
                "CHRX_MEAN_COVERAGE": (
                    round(x_cov, 3)
                    if x_cov is not None
                    else "."
                ),
                "AUTOSOME_MEDIAN_MEAN_COVERAGE": (
                    round(auto_median, 3)
                    if auto_median is not None
                    else "."
                ),
                "CHRX_AUTOSOME_COVERAGE_RATIO": (
                    round(x_auto_ratio, 4)
                    if x_auto_ratio is not None
                    else "."
                ),
                "WHATSHAP_CHRX_HET_SNVS": whatshap_chrX_hets,
                "WHATSHAP_CHRX_PHASED_HET_SNVS": (
                    whatshap_chrX_phased_hets
                ),
                "XX_COMPATIBILITY_QC": xx_qc,
                "HAPLOTAGGED_CHRX_READS": haplotagged_reads,
                "SHARED_WHATSHAP_LONGPHASE_PHASED_SNVS": (
                    shared_phased
                ),
                "FLIP_TOLERANT_WHATSHAP_LONGPHASE_CONCORDANCE": (
                    phase_concordance
                ),
                "INFORMATIVE_PHASE_BLOCKS": len(primary),
                "PHASE_QC_PASS_BLOCKS": len(phase_filtered),
                "PHASE_CONCORDANCE_THRESHOLD": (
                    args.phase_concordance_threshold
                ),
                "PHASE_MIN_SHARED_SNVS": (
                    args.phase_min_shared_snvs
                ),
                "INFORMATIVE_CPG_ISLAND_SUM": informative_cgis,
                "INFORMATIVE_READS": informative_reads,
                "MIN_READS_PER_BLOCK": args.min_block_reads,
                "DISCORDANT_MULTI_CGI_READS_DROPPED": (
                    discordant_reads
                ),
                "STRONG_ORIENTATION_BLOCKS_LOG10_ODDS_GE_1": (
                    strong_orientation_blocks
                ),
                "VERY_STRONG_ORIENTATION_BLOCKS_LOG10_ODDS_GE_2": (
                    very_strong_orientation_blocks
                ),
                "GLOBAL_FOLDED_SKEW_P": (
                    round(mle, 6)
                    if estimate_ok
                    else "."
                ),
                "PROFILE_LIKELIHOOD_CI95_P_LOW": (
                    round(ci_low, 6)
                    if np.isfinite(ci_low)
                    else "."
                ),
                "PROFILE_LIKELIHOOD_CI95_P_HIGH": (
                    round(ci_high, 6)
                    if np.isfinite(ci_high)
                    else "."
                ),
                "BALANCED_XCI_LRT_STATISTIC": (
                    round(lrt_stat, 6)
                    if np.isfinite(lrt_stat)
                    else "."
                ),
                "BALANCED_XCI_LRT_BOUNDARY_P": (
                    round(lrt_p, 8)
                    if np.isfinite(lrt_p)
                    else "."
                ),
                "PHASE_QC_FILTERED_FOLDED_SKEW_P": (
                    round(phase_filtered_mle, 6)
                    if np.isfinite(phase_filtered_mle)
                    else "."
                ),
                "MAJOR_X_PERCENT": (
                    round(major, 2)
                    if estimate_ok
                    else "."
                ),
                "MINOR_X_PERCENT": (
                    round(minor, 2)
                    if estimate_ok
                    else "."
                ),
                "XCI_MAJOR_MINOR_RATIO": ratio,
                "XCI_P_RANGE": range_label,
                "COMMON_THRESHOLD_CONTEXT": context_label,
                "THRESHOLD_70_30_MINOR_P": (
                    args.threshold_70_30_minor
                ),
                "THRESHOLD_80_20_MINOR_P": (
                    args.threshold_80_20_minor
                ),
                "THRESHOLD_90_10_MINOR_P": (
                    args.threshold_90_10_minor
                ),
                "SKEW_METRIC_DEFINITION": (
                    "GLOBAL_FOLDED_SKEW_P is the folded-binomial "
                    "minor-X fraction: P=0.50 is balanced and values "
                    "closer to 0 indicate stronger skew."
                ),
                "LOG10_ODDS_MODEL": (
                    "Per block log10 LR(H1 major Xa vs H2 major Xa) "
                    "= (2*S-N)*log10((1-P)/P), where "
                    "S=H1_Xa+H2_Xi and N=TRIALS."
                ),
                "LRT_NOTE": (
                    "Balanced P=0.5 lies at the boundary of the folded "
                    "parameter space; the reported asymptotic p-value "
                    "uses the 50:50 chi-square(0)/chi-square(1) mixture."
                ),
                "METHOD": (
                    "CpG-island Xa/Xi read clustering + WhatsHap HP/PS "
                    "+ folded-binomial MLE; profile-likelihood CI; "
                    "balanced-XCI likelihood-ratio test; read-threshold "
                    "and WhatsHap/LongPhase phase-QC sensitivity."
                ),
                "INTERPRETATION": (
                    "Report the major:minor ratio and CI as the primary "
                    "result. XCI cutoffs vary among studies; 70:30, "
                    "75:25, 80:20 and 90:10 conventions are all used. "
                    "Phase concordance is a QC/sensitivity layer and "
                    "does not replace the primary folded-binomial fit."
                ),
            }
        ]
    )

    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] status={eligibility} blocks={len(primary)} "
        f"P={mle if estimate_ok else 'NA'} ratio={ratio} "
        f"profile_CI=({ci_low},{ci_high}) "
        f"LRT_p={lrt_p}"
    )


if __name__ == "__main__":
    main()
