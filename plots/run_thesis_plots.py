#!/usr/bin/env python3
"""Run all applicable thesis plots for one completed LRS or SRS sample.

This launcher does not alter biological pipeline outputs. It only discovers files
under the existing per-sample output structure and calls the plotting scripts.
Optional analyses are skipped when their inputs are absent.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def parse_args():
    p = argparse.ArgumentParser(description="Generate all available thesis plots for one WGS sample.")
    p.add_argument("--root", required=True, help="Pipeline output root from the LRS/SRS config path")
    p.add_argument("--sample", required=True)
    p.add_argument("--platform", choices=["lrs", "srs"], default="lrs")
    p.add_argument("--out-dir", default=None, help="Default: <root>/<sample>/plots")
    p.add_argument("--methylation-region", default=None, help="Optional chr:start-end for methylation plot")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def run(command: list[str], required_paths: list[Path], dry_run: bool):
    missing = [p for p in required_paths if not p.exists()]
    if missing:
        print("[SKIP] missing:", ", ".join(str(p) for p in missing))
        return
    print("[RUN]", " ".join(command))
    if not dry_run:
        subprocess.run(command, check=True)


def main():
    args = parse_args()
    root = Path(args.root).expanduser().resolve()
    sample_root = root / args.sample
    out = Path(args.out_dir).expanduser().resolve() if args.out_dir else sample_root / "plots"
    out.mkdir(parents=True, exist_ok=True)
    py = sys.executable
    s = args.sample

    coverage_summary = sample_root / "coverage" / f"{s}.mosdepth.summary.txt"
    coverage_dist = sample_root / "coverage" / f"{s}.mosdepth.global.dist.txt"
    caller_summary = sample_root / "sv" / "merged" / f"{s}_caller_support_summary.tsv"
    integrated_base = sample_root / "gene_discovery" / f"{s}_integrated_SV_gene_analysis.tsv"
    integrated_orthogonal = sample_root / "gene_discovery" / f"{s}_integrated_SV_gene_with_orthogonal_evidence.tsv"
    integrated = (
        integrated_orthogonal
        if args.platform == "srs" and integrated_orthogonal.exists()
        else integrated_base
    )
    needlr = sample_root / "sv" / "needlr" / f"{s}_needLR_RESULTS.tsv"
    ranked = sample_root / "gene_discovery" / f"{s}_ranked_candidates.tsv"
    phenotypes = sample_root / "gene_discovery" / f"{s}_human_gene_phenotypes.tsv"
    straglr = sample_root / "sv" / "straglr" / f"{s}_straglr.annotated.tsv"
    tldr = sample_root / "mei" / "tldr" / f"{s}.tldr.table.txt"
    whatshap = sample_root / "phasing" / f"{s}.phased.vcf.gz"
    longphase = sample_root / "phasing_longphase" / f"{s}.longphase.vcf.gz"
    methylation = sample_root / "methylation" / f"{s}.modkit.tsv"

    caller_order = "Manta,Delly" if args.platform == "srs" else "Sniffles2,cuteSV,Delly"

    jobs = [
        ([py, str(HERE / "plot_lrs_qc.py"), "--summary", str(coverage_summary), "--global-dist", str(coverage_dist), "--out-prefix", str(out / f"{s}_coverage_qc")], [coverage_summary, coverage_dist]),
        ([py, str(HERE / "plot_caller_concordance.py"), "--input", str(caller_summary), "--caller-order", caller_order, "--out-prefix", str(out / f"{s}_caller_concordance")], [caller_summary]),
        ([py, str(HERE / "plot_sv_landscape.py"), "--input", str(integrated), "--out-prefix", str(out / f"{s}_sv_landscape")], [integrated]),
        ([py, str(HERE / "plot_candidate_evidence_matrix.py"), "--input", str(integrated), "--out-prefix", str(out / f"{s}_candidate_evidence")], [integrated]),
        ([py, str(HERE / "plot_candidate_genes.py"), "--input", str(ranked), "--out-prefix", str(out / f"{s}_candidate_genes")], [ranked]),
        ([py, str(HERE / "plot_gene_hpo_heatmap.py"), "--input", str(phenotypes), "--ranking", str(ranked), "--out-prefix", str(out / f"{s}_gene_hpo")], [phenotypes, ranked]),
    ]

    if args.platform == "lrs":
        jobs.extend([
            ([py, str(HERE / "plot_needlr_population.py"), "--input", str(needlr), "--out-prefix", str(out / f"{s}_needlr_population")], [needlr]),
            ([py, str(HERE / "plot_straglr.py"), "--input", str(straglr), "--out-prefix", str(out / f"{s}_straglr")], [straglr]),
            ([py, str(HERE / "plot_mei.py"), "--input", str(tldr), "--out-prefix", str(out / f"{s}_mei")], [tldr]),
        ])

    phase_cmd = [py, str(HERE / "plot_phasing_qc.py"), "--out-prefix", str(out / f"{s}_phasing_qc")]
    phase_required = []
    if whatshap.exists():
        phase_cmd += ["--whatshap-vcf", str(whatshap)]
        phase_required.append(whatshap)
    if args.platform == "lrs" and longphase.exists():
        phase_cmd += ["--longphase-vcf", str(longphase)]
        phase_required.append(longphase)
    if phase_required:
        jobs.append((phase_cmd, phase_required))
    else:
        print("[SKIP] no phased VCF found")

    if args.platform == "lrs" and args.methylation_region:
        jobs.append((
            [py, str(HERE / "plot_methylation.py"), "--input", str(methylation), "--region", args.methylation_region, "--out-prefix", str(out / f"{s}_methylation")],
            [methylation],
        ))
    elif args.platform == "lrs":
        print("[SKIP] methylation plot requires --methylation-region so that a genome-wide read-level file is not plotted indiscriminately")

    for command, required in jobs:
        run(command, required, args.dry_run)

    print(f"[OK] plot directory: {out}")


if __name__ == "__main__":
    main()
