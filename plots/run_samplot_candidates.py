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


def safe_name(value, max_length=60):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "candidate"
    return name[:max_length].rstrip("_.-") or "candidate"


def short_gene_label(primary_gene, all_genes):
    genes = [
        gene
        for gene in str(all_genes).split(";")
        if gene not in {"", ".", "nan", "None"}
    ]
    primary = safe_name(primary_gene, 35) if primary_gene not in {"", ".", "nan", "None"} else "gene"
    if len(genes) <= 1:
        return primary
    return f"{primary}_plus{len(genes) - 1}genes"


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
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--long-read-min", type=int, default=1000)
    p.add_argument("--zoom", type=int, default=20000)
    p.add_argument("--window", type=int, default=5000, help="Local context on each side for small/breakpoint events")
    p.add_argument("--large-sv-threshold", type=int, default=1000000)
    p.add_argument(
        "--gene-annotation",
        default=None,
        help="Optional tabix-indexed BED/GFF gene annotation track for Samplot.",
    )
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
        row["_has_panel_event"] = (
            group["PANEL_STATUS"]
            .fillna("")
            .astype(str)
            .eq("PANEL_GENE")
            .any()
            if "PANEL_STATUS" in group.columns
            else False
        )
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

            if "_has_panel_event" in sub.columns and per_bucket >= 2:
                panel_sub = sub[sub["_has_panel_event"].eq(True)]
                nonpanel_sub = sub[~sub["_has_panel_event"].eq(True)]

                chosen = []
                if not panel_sub.empty:
                    chosen.extend(panel_sub.head(1).index.tolist())
                if not nonpanel_sub.empty:
                    chosen.extend(
                        nonpanel_sub.head(
                            per_bucket - len(chosen)
                        ).index.tolist()
                    )

                if len(chosen) < per_bucket:
                    chosen.extend(
                        sub.loc[~sub.index.isin(chosen)]
                        .head(per_bucket - len(chosen))
                        .index.tolist()
                    )

                selected_indices.extend(chosen)
            else:
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

    def run_one(
        rank,
        sv_id,
        primary_gene,
        genes,
        chrom,
        start,
        end,
        svtype,
        suffix="",
        context_only=False,
        plot_title=None,
    ):
        gene_label = short_gene_label(primary_gene, genes)
        output_name = (
            f"{rank:02d}_{safe_name(svtype, 12)}_"
            f"{safe_name(chrom, 16)}_{int(start)}_{int(end)}_"
            f"{safe_name(sv_id, 45)}_{gene_label}{suffix}.png"
        )
        output = outdir / output_name
        command = [
            "samplot", "plot",
            "-n", plot_title or args.sample,
            "-b", args.bam,
            "-r", args.reference,
            "-o", str(output),
            "-c", str(chrom),
            "-s", str(int(start)),
            "-e", str(int(end)),
            "-q", str(args.min_mapq),
            "--long_read", str(args.long_read_min),
        ]
        if args.gene_annotation:
            command += [
                "-A", args.gene_annotation,
                "--annotation_filenames", "Genes",
                "--annotation_fontsize", "7",
            ]
        if not context_only and svtype in {"DEL", "DUP", "INV"}:
            command += ["-t", svtype]
        span = abs(int(end) - int(start))
        if not context_only and span >= args.large_sv_threshold:
            command += ["--zoom", str(args.zoom)]
        elif not context_only:
            command += ["--window", str(args.window)]
        try:
            completed = subprocess.run(
                command,
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
        except Exception as exc:
            return None, f"SAMPLOT_EXECUTION_ERROR: {exc}"

        valid_output = output.exists() and output.stat().st_size > 0
        if completed.returncode != 0 or not valid_output:
            if output.exists() and output.stat().st_size == 0:
                output.unlink()
            message = (completed.stdout or "").strip().splitlines()
            last_message = message[-1] if message else "Samplot did not create an image"
            return None, f"SAMPLOT_FAILED: {last_message}"

        return (
            output,
            "BREAKPOINT_CONTEXT_ONLY" if context_only else "SV_SIGNAL_PLOT",
        )

    for rank, (_, row) in enumerate(cand.iterrows(), 1):
        sv_id = str(row[id_col])
        genes = str(row["_genes"])
        primary_gene = str(row.get("_gene", "."))
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

        gnomad_af = str(row.get("GNOMAD_SV_AF", "."))
        gnomad_match = str(row.get("GNOMAD_SV_EXACT_MATCH", "."))
        gnomad_id = str(row.get("GNOMAD_SV_ID", "."))
        needlr_af = str(row.get("NEEDLR_AF", "."))
        title_parts = [args.sample]
        if needlr_af not in {"", ".", "nan", "None"}:
            title_parts.append(f"needLR AF={needlr_af}")
        if gnomad_match == "YES" and gnomad_af not in {"", ".", "nan", "None"}:
            title_parts.append(f"gnomAD-SV AF={gnomad_af}")
        plot_title = " | ".join(title_parts)

        outputs = []
        failures = []

        def add_plot(result, requested_scope):
            out, scope_or_error = result
            if out is not None:
                outputs.append((str(out), requested_scope or scope_or_error))
                return True
            failures.append(scope_or_error)
            return False

        if svtype == "INS":
            region_start = max(1, start - args.window)
            region_end = start + args.window
            success = add_plot(
                run_one(
                    rank, sv_id, primary_gene, genes, chrom,
                    region_start, region_end, svtype,
                    "_insertion_context", True, plot_title,
                ),
                "INSERTION_BREAKPOINT_CONTEXT",
            )
            if not success:
                add_plot(
                    run_one(
                        rank, sv_id, primary_gene, genes, chrom,
                        max(1, start - 2 * args.window),
                        start + 2 * args.window,
                        svtype,
                        "_insertion_context_retry",
                        True,
                        plot_title,
                    ),
                    "INSERTION_BREAKPOINT_CONTEXT_RETRY",
                )

        elif svtype in {"BND", "TRA"}:
            bp1_success = add_plot(
                run_one(
                    rank, sv_id, primary_gene, genes, chrom,
                    max(1, start - args.window),
                    start + args.window,
                    svtype, "_bp1", True, plot_title,
                ),
                "BREAKPOINT_1_CONTEXT",
            )
            if not bp1_success:
                add_plot(
                    run_one(
                        rank, sv_id, primary_gene, genes, chrom,
                        max(1, start - 2 * args.window),
                        start + 2 * args.window,
                        svtype, "_bp1_retry", True, plot_title,
                    ),
                    "BREAKPOINT_1_CONTEXT_RETRY",
                )

            if chr2_col and pos2_col:
                chr2 = str(row.get(chr2_col, "."))
                pos2 = number(row.get(pos2_col))
                if chr2 not in {"", ".", "nan", "None"} and pos2 is not None:
                    bp2_success = add_plot(
                        run_one(
                            rank, sv_id, primary_gene, genes, chr2,
                            max(1, pos2 - args.window),
                            pos2 + args.window,
                            svtype, "_bp2", True, plot_title,
                        ),
                        "BREAKPOINT_2_CONTEXT",
                    )
                    if not bp2_success:
                        add_plot(
                            run_one(
                                rank, sv_id, primary_gene, genes, chr2,
                                max(1, pos2 - 2 * args.window),
                                pos2 + 2 * args.window,
                                svtype, "_bp2_retry", True, plot_title,
                            ),
                            "BREAKPOINT_2_CONTEXT_RETRY",
                        )

        else:
            span = abs(end - start)

            # Whole-event Samplot views become visually compressed and can
            # produce unreadable axes for multi-megabase events. For large
            # DEL/DUP/INV calls, show both breakpoints locally by default.
            if svtype in {"DEL", "DUP", "INV"} and span >= args.large_sv_threshold:
                bp1_success = add_plot(
                    run_one(
                        rank, sv_id, primary_gene, genes, chrom,
                        max(1, start - args.zoom),
                        start + args.zoom,
                        svtype, "_bp1_large", True, plot_title,
                    ),
                    "BREAKPOINT_1_LARGE_SV",
                )
                bp2_success = add_plot(
                    run_one(
                        rank, sv_id, primary_gene, genes, chrom,
                        max(1, end - args.zoom),
                        end + args.zoom,
                        svtype, "_bp2_large", True, plot_title,
                    ),
                    "BREAKPOINT_2_LARGE_SV",
                )
                if not bp1_success and not bp2_success:
                    raise RuntimeError(
                        f"Samplot could not create local breakpoint views for "
                        f"large event {sv_id}: " + " | ".join(failures)
                    )
            else:
                main_success = add_plot(
                    run_one(
                        rank, sv_id, primary_gene, genes, chrom,
                        start, end, svtype, plot_title=plot_title,
                    ),
                    "SV_SIGNAL_PLOT",
                )

                # If the whole-event plot cannot be drawn, do not discard the SV.
                if not main_success and svtype in {"DEL", "DUP", "INV"}:
                    bp1_success = add_plot(
                        run_one(
                            rank, sv_id, primary_gene, genes, chrom,
                            max(1, start - args.window),
                            start + args.window,
                            svtype, "_bp1_fallback", True, plot_title,
                        ),
                        "BREAKPOINT_1_FALLBACK",
                    )
                    bp2_success = add_plot(
                        run_one(
                            rank, sv_id, primary_gene, genes, chrom,
                            max(1, end - args.window),
                            end + args.window,
                            svtype, "_bp2_fallback", True, plot_title,
                        ),
                        "BREAKPOINT_2_FALLBACK",
                    )

                    if not bp1_success and not bp2_success:
                        raise RuntimeError(
                            f"Samplot could not create a whole-event or breakpoint "
                            f"plot for {sv_id}: " + " | ".join(failures)
                        )

        if not outputs:
            raise RuntimeError(
                f"No Samplot figure was created for {sv_id}: "
                + " | ".join(failures)
            )

        for output, scope in outputs:
            manifest.append({
                "rank": rank,
                "SV_ID": sv_id,
                "primary_gene": primary_gene,
                "genes": genes,
                "SVTYPE": svtype,
                "requested_top_n": args.top_n,
                "selected_sv_rank": rank,
                "score": row["_score"],
                "score_field": score_col or ".",
                "needLR_AF": row["_af"] if pd.notna(row["_af"]) else ".",
                "caller_count": row["_callers"],
                "event_bucket": row.get("EVENT_REVIEW_BUCKET", "."),
                "gene_relationship": row.get("SV_GENE_RELATIONSHIP", "."),
                "sv_gene_count": row.get("SV_GENE_COUNT", "."),
                "population_tier": row.get("EVENT_POPULATION_TIER", "."),
                "contains_panel_gene": "YES" if bool(row.get("_has_panel_event", False)) else "NO",
                "GNOMAD_SV_EXACT_MATCH": gnomad_match,
                "GNOMAD_SV_ID": gnomad_id,
                "GNOMAD_SV_AF": gnomad_af,
                "NEEDLR_AF": needlr_af,
                "plot_scope": scope,
                "plot_status": (
                    "FALLBACK"
                    if "FALLBACK" in scope or "RETRY" in scope
                    else "PRIMARY"
                ),
                "failed_attempts": " | ".join(failures) if failures else ".",
                "output": output,
            })

    manifest_path = outdir / "samplot_manifest.tsv"
    pd.DataFrame(manifest).to_csv(manifest_path, sep="\t", index=False)
    print(
        f"[OK] selected_SVs={len(cand)} requested_top_n={args.top_n} "
        f"samplot_images={len(manifest)} manifest={manifest_path}"
    )


if __name__ == "__main__":
    main()
