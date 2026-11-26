#!/usr/bin/env python3
"""Sensitivity of a final SV-gene shortlist to evidence eligibility thresholds.

This does not rerun or alter clinical scoring. Baseline event ranks are kept
fixed; the top-N list is reselected among eligible events per scenario.
Missing caller counts remain UNKNOWN; missing AF remains UNKNOWN, not rare.
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def extract_af(frame):
    # gnomAD-SV values are only used when exact match status permits them.
    af = pd.Series(np.nan, index=frame.index, dtype=float)
    if "GNOMAD_SV_AF" in frame:
        values = numeric(frame["GNOMAD_SV_AF"])
        if "GNOMAD_SV_EXACT_MATCH" in frame:
            values = values.where(
                frame["GNOMAD_SV_EXACT_MATCH"].fillna("").astype(str).str.upper().eq("YES")
            )
        af = af.combine_first(values)
    if "NEEDLR_AF" in frame:
        values = numeric(frame["NEEDLR_AF"])
        if "NEEDLR_STATUS" in frame:
            bad = frame["NEEDLR_STATUS"].fillna("").astype(str).str.upper().str.contains(
                "NO_MATCH|NOT_EVALUABLE|NOT_AVAILABLE|NO_POPULATION_MATCH", regex=True
            )
            values = values.mask(bad)
        af = af.combine_first(values)
    return af.where(af.between(0, 1, inclusive="both"))


def unique_gene_ranks(frame):
    return {
        str(gene): i + 1
        for i, gene in enumerate(frame["GENE"].astype(str).drop_duplicates())
        if gene not in {"", ".", "nan", "NA"}
    }


def analyze(path, top_n=20):
    raw = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    required = {"SV_ID", "GENE", "FINAL_EVENT_RANK_GLOBAL"}
    if not required.issubset(raw.columns):
        raise ValueError(f"Final SV-gene table needs {sorted(required)}; found {list(raw)}")
    df = raw.drop_duplicates(["SV_ID", "GENE"]).copy()
    df["_rank"] = numeric(df["FINAL_EVENT_RANK_GLOBAL"])
    df = df.dropna(subset=["_rank"])
    df = df[df["GENE"].notna() & ~df["GENE"].isin([".", "", "NA"])].copy()
    df = df.sort_values(["_rank", "SV_ID", "GENE"])
    callers_col = next((c for c in ["CALLER_COUNT", "SUPP"] if c in df), None)
    if callers_col:
        df["_caller_count"] = numeric(df[callers_col])
    else:
        df["_caller_count"] = np.nan
    df["_af"] = extract_af(df)
    scenarios = [
        ("baseline", None, None),
        ("af_le_0.01_unknown_review", None, 0.01),
        ("af_le_0.001_unknown_review", None, 0.001),
        ("at_least_2_callers", 2, None),
        ("two_callers_af_le_0.01", 2, 0.01),
        ("two_callers_af_le_0.001", 2, 0.001),
    ]
    baseline_ranks = unique_gene_ranks(df)
    baseline = {g for g, rank in baseline_ranks.items() if rank <= top_n}
    summaries, details = [], []
    for label, min_callers, max_af in scenarios:
        # A caller-count comparison is unevaluable when the column is absent.
        if min_callers is not None and not df["_caller_count"].notna().any():
            summaries.append({
                "scenario": label, "status": "NOT_EVALUABLE_NO_CALLER_COUNTS",
                "genes_selected": ".", "shared_top_n": ".", "jaccard": ".",
                "baseline_top_n_retained": ".", "unknown_af_events_retained": ".",
                "unknown_callers_events": ".",
            })
            continue
        selected = df.copy()
        if min_callers is not None:
            selected = selected[selected["_caller_count"] >= min_callers]
        if max_af is not None:
            # Unknown frequency is retained for review, not assumed AF=0.
            selected = selected[selected["_af"].isna() | (selected["_af"] <= max_af)]
        new_ranks = unique_gene_ranks(selected)
        new_top = {g for g, rank in new_ranks.items() if rank <= top_n}
        shared = baseline & new_top
        union = baseline | new_top
        summaries.append({
            "scenario": label, "status": "EVALUATED",
            "genes_selected": len(new_top), "shared_top_n": len(shared),
            "jaccard": len(shared) / len(union) if union else np.nan,
            "baseline_top_n_retained": len(shared) / len(baseline) if baseline else np.nan,
            "unknown_af_events_retained": int(selected["_af"].isna().sum()),
            "unknown_callers_events": int(selected["_caller_count"].isna().sum()),
        })
        for gene in sorted(baseline | new_top):
            details.append({
                "scenario": label, "gene": gene,
                "baseline_gene_rank": baseline_ranks.get(gene, "."),
                "scenario_gene_rank": new_ranks.get(gene, "."),
                "baseline_top_n": gene in baseline, "scenario_top_n": gene in new_top,
            })
    return pd.DataFrame(summaries), pd.DataFrame(details), len(df), len(baseline_ranks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()
    if args.top_n < 1:
        parser.error("--top-n must be positive")
    summary, detail, n_events, n_genes = analyze(args.input, args.top_n)
    folder = Path(args.outdir)
    folder.mkdir(parents=True, exist_ok=True)
    summary.to_csv(folder / "stability_summary.tsv", sep="\t", index=False)
    detail.to_csv(folder / "stability_genes.tsv", sep="\t", index=False)
    manifest = {
        "source": str(Path(args.input).resolve()),
        "top_n": args.top_n,
        "evaluated_sv_gene_rows": n_events,
        "evaluated_unique_genes": n_genes,
        "analysis": "eligibility_filter_sensitivity_not_reranking_or_truth_benchmark",
        "missing_af_policy": "keep unknown for review; do not assume rarity",
        "missing_caller_policy": "unknown excluded only from multi-caller scenarios",
        "note": "Ranks are re-indexed after filtering; original clinical priority is unchanged.",
    }
    (folder / "stability_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    plot = summary[summary["status"] == "EVALUATED"].copy()
    fig, ax = plt.subplots(figsize=(11, 5))
    y = np.arange(len(plot))
    bars = ax.barh(y, pd.to_numeric(plot["baseline_top_n_retained"], errors="coerce"))
    ax.set_yticks(y, labels=plot["scenario"].str.replace("_", " ", regex=False))
    ax.set_xlim(0, 1)
    ax.set_xlabel(f"Fraction of baseline top-{args.top_n} genes retained")
    ax.invert_yaxis()
    for bar, shared in zip(bars, plot["shared_top_n"]):
        ax.text(min(bar.get_width()+0.015, 0.96), bar.get_y()+bar.get_height()/2,
                str(shared), va="center", fontsize=9)
    ax.set_title("SV-gene shortlist sensitivity to evidence eligibility")
    fig.text(0.5, 0.01,
             "Descriptive filter sensitivity; not a new pathogenicity score, functional validation, or diagnostic benchmark.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(folder / "stability_summary.pdf", bbox_inches="tight")
    fig.savefig(folder / "stability_summary.png", bbox_inches="tight", dpi=250)
    plt.close(fig)
    print(f"[OK] {folder} scenarios={len(summary)}")


if __name__ == "__main__":
    main()
