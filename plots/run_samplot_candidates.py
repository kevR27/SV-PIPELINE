#!/usr/bin/env python3
"""Generate read-level Samplot images for prioritized long-read SV candidates."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd


def first_existing(df, names):
    lower = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        hit = lower.get(name.lower())
        if hit is not None:
            return hit
    return None


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "candidate"


def number(value):
    try:
        return float(value)
    except Exception:
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Integrated/multimodal SV-gene TSV")
    p.add_argument("--bam", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--top-n", type=int, default=12)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--long-read-min", type=int, default=1000)
    p.add_argument("--zoom", type=int, default=20000)
    p.add_argument("--window", type=int, default=5000, help="Local context on each side for small/breakpoint events")
    p.add_argument("--large-sv-threshold", type=int, default=1000000)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False)
    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    chrom_col = first_existing(df, ["CHROM", "chrom"])
    start_col = first_existing(df, ["START", "POS"])
    end_col = first_existing(df, ["END"])
    chr2_col = first_existing(df, ["CHR2"])
    pos2_col = first_existing(df, ["POS2"])
    type_col = first_existing(df, ["SVTYPE", "SV_type"])
    score_col = first_existing(df, ["EVENT_GENE_RELEVANCE_SCORE", "INTEGRATED_DISCOVERY_SCORE", "integrated_discovery_score", "PHENOTYPE_SCORE", "ALLELE_RESEARCH_SCORE"])
    caller_col = first_existing(df, ["CALLER_COUNT", "SUPP"])
    af_col = first_existing(df, ["NEEDLR_AF"])
    if None in (id_col, chrom_col, start_col, type_col):
        raise ValueError("Input needs SV_ID, chromosome, start and SVTYPE columns.")

    work = df.copy()
    work["_score"] = pd.to_numeric(work[score_col], errors="coerce").fillna(0) if score_col else 0
    work["_callers"] = pd.to_numeric(work[caller_col], errors="coerce").fillna(0) if caller_col else 0
    work["_af"] = pd.to_numeric(work[af_col], errors="coerce") if af_col else pd.Series(float("nan"), index=work.index)
    work["_population_rank"] = 0
    work.loc[work["_af"].isna(), "_population_rank"] = 1
    work.loc[work["_af"].notna() & work["_af"].le(0.01), "_population_rank"] = 2
    work["_gene"] = work[gene_col].fillna(".").astype(str) if gene_col else "."
    work["_start"] = pd.to_numeric(work[start_col], errors="coerce")
    work["_end"] = pd.to_numeric(work[end_col], errors="coerce") if end_col else work["_start"]
    work["_end"] = work["_end"].fillna(work["_start"])

    bucket_col = first_existing(work, ["EVENT_REVIEW_BUCKET"])
    bucket_rank_col = first_existing(work, ["EVENT_RANK_WITHIN_BUCKET"])

    collapsed = []
    for sv_id, group in work.groupby(id_col, sort=False):
        sort_cols = ["_score", "_population_rank", "_callers"]
        ascending = [False, False, False]

        if bucket_rank_col:
            group = group.copy()
            group["_event_bucket_rank"] = pd.to_numeric(
                group[bucket_rank_col],
                errors="coerce",
            ).fillna(float("inf"))
            sort_cols = ["_event_bucket_rank"] + sort_cols
            ascending = [True] + ascending

        row = group.sort_values(
            sort_cols,
            ascending=ascending,
        ).iloc[0].copy()
        row["_genes"] = ";".join(
            sorted({
                g
                for g in group["_gene"]
                if g not in {"", ".", "nan", "None"}
            })
        ) or "."
        collapsed.append(row)

    cand = pd.DataFrame(collapsed)

    bucket_order = [
        "BREAKPOINT_GENE_CANDIDATE",
        "LARGE_CNV_GENE_CANDIDATE",
        "VERY_LARGE_CNV_GENE_CONTEXT",
        "INSERTION_GENE_CANDIDATE",
        "SMALL_MEDIUM_CNV_GENE_CANDIDATE",
        "LARGE_COMPLEX_INTERVAL_CONTEXT",
    ]

    if bucket_col and bucket_col in cand.columns:
        per_bucket = max(1, args.top_n // len(bucket_order))
        selected_indices = []

        for bucket in bucket_order:
            sub = cand[cand[bucket_col].eq(bucket)].sort_values(
                ["_score", "_population_rank", "_callers"],
                ascending=[False, False, False],
            )
            selected_indices.extend(
                sub.head(per_bucket).index.tolist()
            )

        selected_indices = list(dict.fromkeys(selected_indices))

        if len(selected_indices) < args.top_n:
            remainder = cand.loc[
                ~cand.index.isin(selected_indices)
            ].sort_values(
                ["_score", "_population_rank", "_callers"],
                ascending=[False, False, False],
            )
            selected_indices.extend(
                remainder.head(
                    args.top_n - len(selected_indices)
                ).index.tolist()
            )

        cand = cand.loc[selected_indices].head(args.top_n)
    else:
        cand = cand.sort_values(
            ["_score", "_population_rank", "_callers"],
            ascending=[False, False, False],
        ).head(args.top_n)

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = []

    def run_one(rank, sv_id, genes, chrom, start, end, svtype, suffix="", context_only=False):
        output = outdir / f"{rank:02d}_{safe_name(genes)}_{safe_name(sv_id)}{suffix}.png"
        command = [
            "samplot", "plot",
            "-n", args.sample,
            "-b", args.bam,
            "-r", args.reference,
            "-o", str(output),
            "-c", str(chrom),
            "-s", str(int(start)),
            "-e", str(int(end)),
            "-q", str(args.min_mapq),
            "--long_read", str(args.long_read_min),
        ]
        if not context_only and svtype in {"DEL", "DUP", "INV"}:
            command += ["-t", svtype]
        span = abs(int(end) - int(start))
        if not context_only and span >= args.large_sv_threshold:
            command += ["--zoom", str(args.zoom)]
        elif not context_only:
            command += ["--window", str(args.window)]
        subprocess.run(command, check=True)
        return output, "BREAKPOINT_CONTEXT_ONLY" if context_only else "SV_SIGNAL_PLOT"

    for rank, (_, row) in enumerate(cand.iterrows(), 1):
        sv_id = str(row[id_col])
        genes = str(row["_genes"])
        chrom = str(row[chrom_col])
        start = number(row["_start"])
        end = number(row["_end"])
        svtype = str(row[type_col]).upper()
        if start is None:
            continue
        if end is None:
            end = start + 1
        if end <= start:
            end = start + 1

        outputs = []
        if svtype == "INS":
            region_start = max(1, start - args.window)
            region_end = start + args.window
            out, _scope = run_one(
                rank, sv_id, genes, chrom, region_start, region_end,
                svtype, "_insertion_context", True,
            )
            outputs.append((str(out), "INSERTION_BREAKPOINT_CONTEXT"))
        elif svtype in {"BND", "TRA"}:
            region_start = max(1, start - args.window)
            region_end = start + args.window
            out1, _scope1 = run_one(
                rank, sv_id, genes, chrom, region_start, region_end,
                svtype, "_bp1", True,
            )
            outputs.append((str(out1), "BREAKPOINT_CONTEXT_ONLY"))
            if chr2_col and pos2_col:
                chr2 = str(row.get(chr2_col, "."))
                pos2 = number(row.get(pos2_col))
                if chr2 not in {"", ".", "nan", "None"} and pos2 is not None:
                    out2, _scope2 = run_one(
                        rank, sv_id, genes, chr2,
                        max(1, pos2 - args.window), pos2 + args.window,
                        svtype, "_bp2", True,
                    )
                    outputs.append((str(out2), "BREAKPOINT_CONTEXT_ONLY"))
        else:
            out, scope = run_one(rank, sv_id, genes, chrom, start, end, svtype)
            outputs.append((str(out), scope))

        for output, scope in outputs:
            manifest.append({
                "rank": rank,
                "SV_ID": sv_id,
                "genes": genes,
                "SVTYPE": svtype,
                "score": row["_score"],
                "score_field": score_col or ".",
                "needLR_AF": row["_af"] if pd.notna(row["_af"]) else ".",
                "caller_count": row["_callers"],
                "event_bucket": row.get("EVENT_REVIEW_BUCKET", "."),
                "gene_relationship": row.get("SV_GENE_RELATIONSHIP", "."),
                "sv_gene_count": row.get("SV_GENE_COUNT", "."),
                "population_tier": row.get("EVENT_POPULATION_TIER", "."),
                "plot_scope": scope,
                "output": output,
            })

    manifest_path = outdir / "samplot_manifest.tsv"
    pd.DataFrame(manifest).to_csv(manifest_path, sep="\t", index=False)
    print(f"[OK] samplot_images={len(manifest)} manifest={manifest_path}")


if __name__ == "__main__":
    main()
