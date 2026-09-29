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
    p.add_argument("--top-n", type=int, default=8)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--long-read-min", type=int, default=1000)
    p.add_argument("--zoom", type=int, default=20000)
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
    score_col = first_existing(df, ["ALLELE_RESEARCH_SCORE", "INTEGRATED_DISCOVERY_SCORE", "PHENOTYPE_SCORE"])
    caller_col = first_existing(df, ["CALLER_COUNT", "SUPP"])
    if None in (id_col, chrom_col, start_col, type_col):
        raise ValueError("Input needs SV_ID, chromosome, start and SVTYPE columns.")

    work = df.copy()
    work["_score"] = pd.to_numeric(work[score_col], errors="coerce").fillna(0) if score_col else 0
    work["_callers"] = pd.to_numeric(work[caller_col], errors="coerce").fillna(0) if caller_col else 0
    work["_gene"] = work[gene_col].fillna(".").astype(str) if gene_col else "."
    work["_start"] = pd.to_numeric(work[start_col], errors="coerce")
    work["_end"] = pd.to_numeric(work[end_col], errors="coerce") if end_col else work["_start"]
    work["_end"] = work["_end"].fillna(work["_start"])

    collapsed = []
    for sv_id, group in work.groupby(id_col, sort=False):
        row = group.sort_values(["_score", "_callers"], ascending=False).iloc[0].copy()
        row["_genes"] = ";".join(sorted({g for g in group["_gene"] if g not in {"", ".", "nan", "None"}})) or "."
        collapsed.append(row)
    cand = pd.DataFrame(collapsed)
    cand = cand.sort_values(["_score", "_callers"], ascending=False).head(args.top_n)

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = []

    def run_one(rank, sv_id, genes, chrom, start, end, svtype, suffix="", context_only=False):
        output = outdir / f"{rank:02d}_{safe_name(genes)}_{safe_name(sv_id)}{suffix}.png"
        command = [
            "samplot", "plot",
            "-n", args.sample,
            "-b", args.bam,
            "-o", str(output),
            "-c", str(chrom),
            "-s", str(int(start)),
            "-e", str(int(end)),
            "-q", str(args.min_mapq),
            "--long_read", str(args.long_read_min),
        ]
        if not context_only and svtype in {"DEL", "DUP", "INV", "INS"}:
            command += ["-t", svtype]
        span = abs(int(end) - int(start))
        if not context_only and span >= args.large_sv_threshold:
            command += ["--zoom", str(args.zoom)]
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
        if svtype in {"BND", "TRA"}:
            out1, scope1 = run_one(rank, sv_id, genes, chrom, start, start + 1, svtype, "_bp1", True)
            outputs.append((str(out1), scope1))
            if chr2_col and pos2_col:
                chr2 = str(row.get(chr2_col, "."))
                pos2 = number(row.get(pos2_col))
                if chr2 not in {"", ".", "nan", "None"} and pos2 is not None:
                    out2, scope2 = run_one(rank, sv_id, genes, chr2, pos2, pos2 + 1, svtype, "_bp2", True)
                    outputs.append((str(out2), scope2))
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
                "caller_count": row["_callers"],
                "plot_scope": scope,
                "output": output,
            })

    manifest_path = outdir / "samplot_manifest.tsv"
    pd.DataFrame(manifest).to_csv(manifest_path, sep="\t", index=False)
    print(f"[OK] samplot_images={len(manifest)} manifest={manifest_path}")


if __name__ == "__main__":
    main()
