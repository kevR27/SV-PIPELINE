#!/usr/bin/env python3
"""Run all applicable thesis plots for one completed LRS or SRS sample.

The launcher does not alter biological pipeline outputs. Optional analyses are
skipped when their inputs are absent.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def parse_args():
    p = argparse.ArgumentParser(description="Generate all available thesis plots for one WGS sample.")
    p.add_argument("--root", required=True, help="Pipeline output root from config")
    p.add_argument("--sample", required=True)
    p.add_argument("--platform", choices=["lrs", "srs"], default="lrs")
    p.add_argument("--out-dir", default=None, help="Default: <root>/<sample>/plots")
    p.add_argument("--methylation-region", default=None, help="Optional chr:start-end for methylation plot")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def first_existing_path(candidates: list[Path]) -> Path | None:
    for path in candidates:
        if path.exists():
            return path
    return None


def first_glob(root: Path, patterns: list[str]) -> Path | None:
    for pattern in patterns:
        hits = sorted(root.glob(pattern))
        if hits:
            return hits[0]
    return None


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

    coverage_summary = first_existing_path(
        [
            sample_root / "coverage" / f"{s}.mosdepth.summary.txt",
            sample_root / "coverage" / f"{s}.summary.txt",
        ]
    )
    coverage_dist = first_existing_path(
        [
            sample_root / "coverage" / f"{s}.mosdepth.global.dist.txt",
            sample_root / "coverage" / f"{s}.global.dist.txt",
        ]
    )

    caller_summary = sample_root / "sv" / "merged" / f"{s}_caller_support_summary.tsv"
    integrated = sample_root / "gene_discovery" / f"{s}_integrated_SV_gene_analysis.tsv"
    needlr = sample_root / "sv" / "needlr" / f"{s}_needLR_RESULTS.tsv"
    ranked = sample_root / "gene_discovery" / f"{s}_ranked_candidates.tsv"
    phenotypes = sample_root / "gene_discovery" / f"{s}_human_gene_phenotypes.tsv"
    straglr = sample_root / "sv" / "straglr" / f"{s}_straglr.annotated.tsv"
    tldr = sample_root / "mei" / "tldr" / f"{s}.tldr.table.txt"

    clair3_phased = first_glob(
        sample_root,
        [
            "clair3/**/*phased*.vcf.gz",
            "small_variants/**/*phased*.vcf.gz",
            "phasing/**/*clair3*.vcf.gz",
            "phasing/**/*phased*.vcf.gz",
        ],
    )
    longphase = first_glob(
        sample_root,
        [
            "phasing_longphase/*.longphase.vcf.gz",
            "phasing/**/*longphase*.vcf.gz",
            "**/*.longphase.vcf.gz",
        ],
    )

    methylation = first_glob(
        sample_root,
        [
            "methylation/*.cpg.bedmethyl.gz",
            "methylation/*.bedmethyl.gz",
            "**/*.cpg.bedmethyl.gz",
            "**/*.bedmethyl.gz",
        ],
    )

    caller_order = "Manta,Delly" if args.platform == "srs" else "Sniffles2,cuteSV,Delly"

    jobs = []

    if coverage_summary:
        cmd = [
            py, str(HERE / "plot_lrs_qc.py"),
            "--summary", str(coverage_summary),
            "--out-prefix", str(out / f"{s}_coverage_qc"),
        ]
        required = [coverage_summary]
        if coverage_dist:
            cmd += ["--global-dist", str(coverage_dist)]
            required.append(coverage_dist)
        jobs.append((cmd, required))
    else:
        print("[SKIP] no mosdepth summary found")

    jobs.extend(
        [
            (
                [
                    py, str(HERE / "plot_caller_concordance.py"),
                    "--input", str(caller_summary),
                    "--caller-order", caller_order,
                    "--out-prefix", str(out / f"{s}_caller_concordance"),
                ],
                [caller_summary],
            ),
            (
                [
                    py, str(HERE / "plot_sv_landscape.py"),
                    "--input", str(integrated),
                    "--out-prefix", str(out / f"{s}_sv_landscape"),
                ],
                [integrated],
            ),
            (
                [
                    py, str(HERE / "plot_candidate_evidence_matrix.py"),
                    "--input", str(integrated),
                    "--out-prefix", str(out / f"{s}_candidate_evidence"),
                ],
                [integrated],
            ),
            (
                [
                    py, str(HERE / "plot_candidate_genes.py"),
                    "--input", str(ranked),
                    "--out-prefix", str(out / f"{s}_candidate_genes"),
                ],
                [ranked],
            ),
            (
                [
                    py, str(HERE / "plot_gene_hpo_heatmap.py"),
                    "--input", str(phenotypes),
                    "--ranking", str(ranked),
                    "--out-prefix", str(out / f"{s}_gene_hpo"),
                ],
                [phenotypes, ranked],
            ),
        ]
    )

    if args.platform == "lrs":
        jobs.extend(
            [
                (
                    [
                        py, str(HERE / "plot_needlr_population.py"),
                        "--input", str(needlr),
                        "--out-prefix", str(out / f"{s}_needlr_population"),
                    ],
                    [needlr],
                ),
                (
                    [
                        py, str(HERE / "plot_straglr.py"),
                        "--input", str(straglr),
                        "--out-prefix", str(out / f"{s}_straglr"),
                    ],
                    [straglr],
                ),
                (
                    [
                        py, str(HERE / "plot_mei.py"),
                        "--input", str(tldr),
                        "--out-prefix", str(out / f"{s}_mei"),
                    ],
                    [tldr],
                ),
            ]
        )

    phase_cmd = [
        py,
        str(HERE / "plot_phasing_qc.py"),
        "--out-prefix",
        str(out / f"{s}_phasing_qc"),
    ]
    phase_required = []
    if clair3_phased:
        phase_cmd += ["--clair3-vcf", str(clair3_phased)]
        phase_required.append(clair3_phased)
    if args.platform == "lrs" and longphase:
        phase_cmd += ["--longphase-vcf", str(longphase)]
        phase_required.append(longphase)

    if phase_required:
        jobs.append((phase_cmd, phase_required))
    else:
        print("[SKIP] no phased Clair3/LongPhase VCF found")

    if args.platform == "lrs" and args.methylation_region and methylation:
        jobs.append(
            (
                [
                    py,
                    str(HERE / "plot_methylation.py"),
                    "--input",
                    str(methylation),
                    "--format",
                    "bedmethyl",
                    "--region",
                    args.methylation_region,
                    "--out-prefix",
                    str(out / f"{s}_methylation"),
                ],
                [methylation],
            )
        )
    elif args.platform == "lrs" and args.methylation_region:
        print("[SKIP] methylation region requested but no bedMethyl file was found")
    elif args.platform == "lrs":
        print("[SKIP] methylation plot requires --methylation-region")

    for command, required in jobs:
        run(command, required, args.dry_run)

    print(f"[OK] plot directory: {out}")


if __name__ == "__main__":
    main()
