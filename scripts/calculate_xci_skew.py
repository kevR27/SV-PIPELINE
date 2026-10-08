#!/usr/bin/env python3
"""Estimate X-chromosome inactivation skew with uncertainty and phase QC.

The primary estimator remains the folded-binomial maximum-likelihood model
used by Gocuk et al. / SkewX. Folding is necessary because HP1/HP2 labels can
flip between independent phase blocks. This implementation adds:

* WhatsHap-LongPhase block concordance filtering;
* block-bootstrap 95% confidence intervals;
* sensitivity estimates across minimum-read cutoffs;
* explicit major:minor classification thresholds;
* chrX/autosomal coverage and chrX heterozygosity QC;
* an explicit per-block log10 likelihood ratio for Xa orientation.

The output is a research XCI estimate, not a pathogenicity classification.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp


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
        work[chrom_col].astype(str).str.match(r"^(chr)?([1-9]|1[0-9]|2[0-2])$")
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
    return math.log(1.0 - 0.5 * delta) + log_choose + logsumexp([a, b])


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
        options={"xatol": 1e-8},
    )
    return float(result.x) if result.success else np.nan


def bootstrap_ci(
    df: pd.DataFrame,
    replicates: int,
    seed: int,
) -> tuple[float, float, int]:
    if df.empty or replicates <= 0:
        return np.nan, np.nan, 0

    rng = np.random.default_rng(seed)
    values = []
    n = len(df)

    for _ in range(replicates):
        sampled = df.iloc[rng.integers(0, n, size=n)]
        value = fit_folded_mle(sampled)
        if np.isfinite(value):
            values.append(value)

    if not values:
        return np.nan, np.nan, 0

    low, high = np.quantile(values, [0.025, 0.975])
    return float(low), float(high), len(values)


def parse_cutoffs(text: str) -> list[int]:
    values = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        values.append(int(token))
    return sorted(set(values))


def ratio_text(p_minor: float) -> str:
    if not np.isfinite(p_minor):
        return "."
    return f"{100.0 * (1.0 - p_minor):.1f}:{100.0 * p_minor:.1f}"


def classify_xci(
    p_minor: float,
    random_minor: float,
    high_minor: float,
    extreme_minor: float,
) -> str:
    if not np.isfinite(p_minor):
        return "NOT_ESTIMATED"
    if p_minor >= random_minor:
        return "RANDOM_RANGE_MAJOR_LE_70_PERCENT"
    if p_minor >= high_minor:
        return "MODERATE_IMBALANCE_70_TO_80_PERCENT"
    if p_minor >= extreme_minor:
        return "HIGH_SKEW_80_TO_90_PERCENT"
    return "EXTREME_SKEW_MAJOR_GT_90_PERCENT"


def load_phase_qc(path: str) -> pd.DataFrame:
    phase = pd.read_csv(path, sep="\t", dtype=str)
    if phase.empty or "WHATSHAP_PS" not in phase.columns:
        return pd.DataFrame(
            columns=[
                "PS",
                "PHASE_QC_SHARED_SNVS",
                "PHASE_QC_CONCORDANT_SNVS",
                "PHASE_QC_CONCORDANCE",
            ]
        )

    for col in ["SHARED_PHASED_SNVS", "BEST_CONCORDANT_SNVS"]:
        if col in phase.columns:
            phase[col] = pd.to_numeric(phase[col], errors="coerce").fillna(0)
        else:
            phase[col] = 0

    grouped = (
        phase.groupby("WHATSHAP_PS", as_index=False)
        .agg(
            PHASE_QC_SHARED_SNVS=("SHARED_PHASED_SNVS", "sum"),
            PHASE_QC_CONCORDANT_SNVS=("BEST_CONCORDANT_SNVS", "sum"),
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
    p.add_argument("--min-primary-blocks", type=int, default=3)
    p.add_argument("--min-chrx-het-snvs", type=int, default=50)
    p.add_argument("--random-minor-threshold", type=float, default=0.30)
    p.add_argument("--high-skew-minor-threshold", type=float, default=0.20)
    p.add_argument("--extreme-skew-minor-threshold", type=float, default=0.10)
    p.add_argument("--orientation-min-log10-odds", type=float, default=1.0)
    p.add_argument("--bootstrap-replicates", type=int, default=2000)
    p.add_argument("--bootstrap-seed", type=int, default=27)
    p.add_argument("--sensitivity-min-reads", default="5,8,10,15")
    p.add_argument("--blocks-output", required=True)
    p.add_argument("--summary-output", required=True)
    p.add_argument("--sensitivity-output", required=True)
    args = p.parse_args()

    if not (
        0 < args.extreme_skew_minor_threshold
        <= args.high_skew_minor_threshold
        <= args.random_minor_threshold
        <= 0.5
    ):
        raise ValueError(
            "Require 0 < extreme <= high <= random <= 0.5 for minor-X thresholds."
        )

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

    # "Success" means support for the local orientation H1=major Xa.
    blocks["SUCCESS_H1_XA"] = (
        blocks["H1_Xa"].fillna(0)
        + blocks["H2_Xi"].fillna(0)
    ).astype(int)
    blocks["SUCCESS_FOLDED"] = np.minimum(
        blocks["SUCCESS_H1_XA"],
        blocks["TRIALS"] - blocks["SUCCESS_H1_XA"],
    )
    blocks["FOLDED_MINOR_PROPORTION"] = np.where(
        blocks["TRIALS"] > 0,
        blocks["SUCCESS_FOLDED"] / blocks["TRIALS"],
        np.nan,
    )
    # Backward-compatible column name for older plotting/report code.
    blocks["FOLDED_BLOCK_SKEW"] = blocks["FOLDED_MINOR_PROPORTION"]
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
        & (blocks["PHASE_QC_CONCORDANCE"] >= args.phase_concordance_threshold),
        "YES",
        "NO",
    )

    read_filtered = blocks[blocks["TRIALS"].ge(args.min_block_reads)].copy()
    phase_filtered = read_filtered[read_filtered["PHASE_QC_PASS"].eq("YES")].copy()

    all_mle = fit_folded_mle(read_filtered)
    phase_mle = fit_folded_mle(phase_filtered)

    if len(phase_filtered) >= args.min_primary_blocks and np.isfinite(phase_mle):
        primary = phase_filtered.copy()
        mle = phase_mle
        estimate_basis = "PHASE_QC_FILTERED_BLOCKS"
    else:
        primary = read_filtered.copy()
        mle = all_mle
        estimate_basis = "READ_FILTERED_FALLBACK_PHASE_QC_INSUFFICIENT"

    ci_low, ci_high, bootstrap_success = bootstrap_ci(
        primary,
        args.bootstrap_replicates,
        args.bootstrap_seed,
    )

    x_cov, auto_median, x_auto_ratio = read_coverage_qc(
        args.mosdepth_summary,
        args.chrom,
    )

    hap_summary = pd.read_csv(args.haplotag_summary, sep="\t")
    haplotagged_reads = (
        int(
            pd.to_numeric(
                hap_summary.get("haplotagged_reads"), errors="coerce"
            ).fillna(0).iloc[0]
        )
        if not hap_summary.empty and "haplotagged_reads" in hap_summary.columns
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
            phase_concordance = row["FLIP_TOLERANT_PHASE_CONCORDANCE"]
        for key, target in [
            ("SHARED_PHASED_SNVS", "shared"),
            ("WHATSHAP_CHRX_HET_SNVS", "hets"),
            ("WHATSHAP_CHRX_PHASED_HET_SNVS", "phased"),
        ]:
            if key in phase_summary.columns:
                try:
                    value = int(float(row[key]))
                except Exception:
                    value = 0
                if target == "shared":
                    shared_phased = value
                elif target == "hets":
                    whatshap_chrX_hets = value
                else:
                    whatshap_chrX_phased_hets = value

    coverage_ok = x_cov is not None and x_cov >= args.min_chrx_coverage
    blocks_ok = len(primary) > 0
    estimate_ok = np.isfinite(mle)

    if not coverage_ok:
        eligibility = "LOW_CHRX_COVERAGE_REVIEW"
    elif not blocks_ok:
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

    status = classify_xci(
        mle,
        args.random_minor_threshold,
        args.high_skew_minor_threshold,
        args.extreme_skew_minor_threshold,
    )
    ratio = ratio_text(mle)
    major = 100.0 * (1.0 - mle) if estimate_ok else np.nan
    minor = 100.0 * mle if estimate_ok else np.nan

    # Explicit orientation likelihood ratio:
    # H1-major-Xa hypothesis: successes ~ q=1-p, failures ~ p
    # H2-major-Xa hypothesis: successes ~ p, failures ~ q
    # log10 LR = (2*successes - n) * log10(q/p).
    if estimate_ok and not read_filtered.empty:
        eps = 1e-12
        p_minor = min(max(float(mle), eps), 0.5 - eps)
        q_major = 1.0 - p_minor

        successes = read_filtered["SUCCESS_H1_XA"].astype(float)
        trials = read_filtered["TRIALS"].astype(float)
        log10_odds = (
            (2.0 * successes - trials)
            * math.log10(q_major / p_minor)
        )

        read_filtered["LOG10_ODDS_H1_XA_VS_H2_XA"] = log10_odds
        read_filtered["PREFERRED_XA_HAPLOTYPE"] = np.where(
            log10_odds > 0,
            "H1",
            np.where(log10_odds < 0, "H2", "UNRESOLVED"),
        )
        read_filtered["PREFERRED_XI_HAPLOTYPE"] = np.where(
            log10_odds > 0,
            "H2",
            np.where(log10_odds < 0, "H1", "UNRESOLVED"),
        )
        read_filtered["ORIENTED_MAJOR_XA_PROPORTION"] = np.maximum(
            read_filtered["H1_XA_PROPORTION_RAW"],
            1.0 - read_filtered["H1_XA_PROPORTION_RAW"],
        )
        abs_lod = np.abs(log10_odds)
        read_filtered["ORIENTATION_EVIDENCE"] = np.select(
            [abs_lod >= 2.0, abs_lod >= 1.0, abs_lod >= 0.5],
            [
                "VERY_STRONG_GE_100_TO_1",
                "STRONG_GE_10_TO_1",
                "MODERATE_GE_3_TO_1",
            ],
            default="WEAK_LT_3_TO_1",
        )
        read_filtered["ORIENTATION_USABLE"] = np.where(
            (abs_lod >= args.orientation_min_log10_odds)
            & read_filtered["PHASE_QC_PASS"].eq("YES"),
            "YES",
            "NO",
        )
    else:
        read_filtered["LOG10_ODDS_H1_XA_VS_H2_XA"] = np.nan
        read_filtered["PREFERRED_XA_HAPLOTYPE"] = "."
        read_filtered["PREFERRED_XI_HAPLOTYPE"] = "."
        read_filtered["ORIENTED_MAJOR_XA_PROPORTION"] = np.nan
        read_filtered["ORIENTATION_EVIDENCE"] = "NOT_ESTIMATED"
        read_filtered["ORIENTATION_USABLE"] = "NO"

    read_filtered["GLOBAL_FOLDED_MINOR_X_P"] = (
        round(mle, 6) if estimate_ok else np.nan
    )
    read_filtered["GLOBAL_ESTIMATE_BASIS"] = estimate_basis

    out_blocks = Path(args.blocks_output)
    out_blocks.parent.mkdir(parents=True, exist_ok=True)
    read_filtered.to_csv(out_blocks, sep="\t", index=False)

    sensitivity_rows = []
    for i, cutoff in enumerate(parse_cutoffs(args.sensitivity_min_reads)):
        subset_all = blocks[blocks["TRIALS"].ge(cutoff)].copy()
        subset_phase = subset_all[subset_all["PHASE_QC_PASS"].eq("YES")].copy()

        if len(subset_phase) >= args.min_primary_blocks:
            chosen = subset_phase
            basis = "PHASE_QC_FILTERED_BLOCKS"
        else:
            chosen = subset_all
            basis = "READ_FILTERED_FALLBACK_PHASE_QC_INSUFFICIENT"

        estimate = fit_folded_mle(chosen)
        low, high, boot_n = bootstrap_ci(
            chosen,
            args.bootstrap_replicates,
            args.bootstrap_seed + i + 1,
        )
        sensitivity_rows.append(
            {
                "MIN_READS_PER_BLOCK": cutoff,
                "ALL_READ_QUALIFIED_BLOCKS": len(subset_all),
                "PHASE_QC_PASS_BLOCKS": len(subset_phase),
                "ESTIMATE_BASIS": basis,
                "FOLDED_MINOR_X_P": estimate,
                "MAJOR_X_PERCENT": (
                    100.0 * (1.0 - estimate)
                    if np.isfinite(estimate)
                    else np.nan
                ),
                "MINOR_X_PERCENT": (
                    100.0 * estimate
                    if np.isfinite(estimate)
                    else np.nan
                ),
                "XCI_MAJOR_MINOR_RATIO": ratio_text(estimate),
                "BOOTSTRAP_CI95_MINOR_LOW": low,
                "BOOTSTRAP_CI95_MINOR_HIGH": high,
                "BOOTSTRAP_SUCCESSFUL_REPLICATES": boot_n,
                "XCI_SKEW_STATUS": classify_xci(
                    estimate,
                    args.random_minor_threshold,
                    args.high_skew_minor_threshold,
                    args.extreme_skew_minor_threshold,
                ),
            }
        )

    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity_path = Path(args.sensitivity_output)
    sensitivity_path.parent.mkdir(parents=True, exist_ok=True)
    sensitivity.to_csv(sensitivity_path, sep="\t", index=False)

    informative_cgis = (
        int(read_filtered["INFORMATIVE_CPG_ISLANDS"].fillna(0).sum())
        if "INFORMATIVE_CPG_ISLANDS" in read_filtered.columns
        else 0
    )
    informative_reads = (
        int(read_filtered["TRIALS"].sum())
        if "TRIALS" in read_filtered.columns
        else 0
    )

    lod = pd.to_numeric(
        read_filtered.get("LOG10_ODDS_H1_XA_VS_H2_XA"),
        errors="coerce",
    )
    strong_orientation_blocks = int((lod.abs() >= 1.0).sum()) if lod.notna().any() else 0
    very_strong_orientation_blocks = int((lod.abs() >= 2.0).sum()) if lod.notna().any() else 0

    summary = pd.DataFrame(
        [
            {
                "chrom": args.chrom,
                "XCI_ANALYSIS_STATUS": eligibility,
                "PRIMARY_ESTIMATE_BASIS": estimate_basis,
                "CHRX_MEAN_COVERAGE": round(x_cov, 3) if x_cov is not None else ".",
                "AUTOSOME_MEDIAN_MEAN_COVERAGE": (
                    round(auto_median, 3) if auto_median is not None else "."
                ),
                "CHRX_AUTOSOME_COVERAGE_RATIO": (
                    round(x_auto_ratio, 4) if x_auto_ratio is not None else "."
                ),
                "WHATSHAP_CHRX_HET_SNVS": whatshap_chrX_hets,
                "WHATSHAP_CHRX_PHASED_HET_SNVS": whatshap_chrX_phased_hets,
                "XX_COMPATIBILITY_QC": xx_qc,
                "MIN_CHRX_HET_SNVS_FOR_QC": args.min_chrx_het_snvs,
                "HAPLOTAGGED_CHRX_READS": haplotagged_reads,
                "SHARED_WHATSHAP_LONGPHASE_PHASED_SNVS": shared_phased,
                "FLIP_TOLERANT_WHATSHAP_LONGPHASE_CONCORDANCE": phase_concordance,
                "READ_QUALIFIED_PHASE_BLOCKS": len(read_filtered),
                "PHASE_QC_PASS_BLOCKS": len(phase_filtered),
                "PHASE_CONCORDANCE_THRESHOLD": args.phase_concordance_threshold,
                "PHASE_MIN_SHARED_SNVS": args.phase_min_shared_snvs,
                "INFORMATIVE_CPG_ISLAND_SUM": informative_cgis,
                "INFORMATIVE_READS": informative_reads,
                "MIN_READS_PER_BLOCK": args.min_block_reads,
                "STRONG_ORIENTATION_BLOCKS_LOG10_ODDS_GE_1": strong_orientation_blocks,
                "VERY_STRONG_ORIENTATION_BLOCKS_LOG10_ODDS_GE_2": very_strong_orientation_blocks,
                "ORIENTATION_MIN_LOG10_ODDS_FOR_XA_XI_PLOT": (
                    args.orientation_min_log10_odds
                ),
                "ALL_READ_FILTERED_FOLDED_MINOR_X_P": (
                    round(all_mle, 6) if np.isfinite(all_mle) else "."
                ),
                "PHASE_QC_FOLDED_MINOR_X_P": (
                    round(phase_mle, 6) if np.isfinite(phase_mle) else "."
                ),
                "GLOBAL_FOLDED_SKEW_P": (
                    round(mle, 6) if estimate_ok else "."
                ),
                "FOLDED_MINOR_X_P_CI95_LOW": (
                    round(ci_low, 6) if np.isfinite(ci_low) else "."
                ),
                "FOLDED_MINOR_X_P_CI95_HIGH": (
                    round(ci_high, 6) if np.isfinite(ci_high) else "."
                ),
                "BOOTSTRAP_REPLICATES": args.bootstrap_replicates,
                "BOOTSTRAP_SUCCESSFUL_REPLICATES": bootstrap_success,
                "MAJOR_X_PERCENT": round(major, 2) if estimate_ok else ".",
                "MINOR_X_PERCENT": round(minor, 2) if estimate_ok else ".",
                "XCI_MAJOR_MINOR_RATIO": ratio,
                "XCI_SKEW_STATUS": status,
                "RANDOM_RANGE_MINOR_THRESHOLD": args.random_minor_threshold,
                "HIGH_SKEW_MINOR_THRESHOLD": args.high_skew_minor_threshold,
                "EXTREME_SKEW_MINOR_THRESHOLD": args.extreme_skew_minor_threshold,
                "SKEW_METRIC_DEFINITION": (
                    "GLOBAL_FOLDED_SKEW_P is the minor-X fraction after "
                    "folded-binomial maximum-likelihood estimation: 0.50 is "
                    "balanced and values toward 0 indicate stronger skew."
                ),
                "LOG10_ODDS_MODEL": (
                    "Per block: log10 LR(H1 major Xa vs H2 major Xa) = "
                    "(2*S-N)*log10((1-P)/P), where S=H1_Xa+H2_Xi, "
                    "N=TRIALS, and P is the primary folded minor-X MLE."
                ),
                "METHOD": (
                    "CpG-island Xa/Xi read clustering + WhatsHap HP/PS + "
                    "folded-binomial MLE; block-bootstrap CI; minimum-read "
                    "sensitivity; WhatsHap-LongPhase phase-concordance QC; "
                    "per-block orientation likelihood ratios."
                ),
                "INTERPRETATION": (
                    "Classification is descriptive and uses the configured "
                    "70:30/80:20/90:10 major:minor convention. Published XCI "
                    "studies use different cutoffs, so this is not a universal "
                    "clinical threshold. XX compatibility is a sequencing QC "
                    "based on chrX/autosomal depth and chrX heterozygosity, not "
                    "a karyotype diagnosis."
                ),
            }
        ]
    )

    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] status={eligibility} basis={estimate_basis} "
        f"blocks={len(primary)} P={mle if estimate_ok else 'NA'} "
        f"ratio={ratio} CI95=({ci_low},{ci_high})"
    )


if __name__ == "__main__":
    main()
