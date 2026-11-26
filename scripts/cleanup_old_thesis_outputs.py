#!/usr/bin/env python3
"""Preview/remove ONLY explicitly named old thesis outputs on any RESULTS disk.

Designed for server RESULTS or an independently mounted backup RESULTS folder.
Run BEFORE copying fresh plots. Default is dry-run. No upstream FASTQ/BAM/VCF,
phasing, coverage, methylation, mtDNA or QC source is ever selected.
"""
import argparse
import csv
import shutil
import sys
from datetime import datetime
from pathlib import Path

DEFAULT_SAMPLES = ("11446", "13048", "19545", "14996")
GENERATED_GENE_FILES = (
    "{s}_allele_assessment.tsv",
    "{s}_allele_disease_hypotheses.tsv",
    "{s}_allele_assessment.manifest.json",
    "{s}_integrated_SV_gene_analysis.final.tsv.gz",
    "{s}_independent_complementary_findings.tsv",
    "{s}_mitocarta_annotation_summary.tsv",
    "{s}_gene_multimodal_evidence_summary.tsv",
    "{s}_mitochondrial_gene_ranking.tsv",
    "{s}_integrated_SV_gene_with_complementary_evidence.tsv",
    "{s}_integrated_SV_gene_with_multimodal_context.tsv",
    "{s}_integrated_SV_gene_with_mitocarta.tsv",
    "{s}_exomiser_comparison.tsv",
)


def candidates(root, samples, scope):
    paths = []
    for s in samples:
        parent = root / s
        if not parent.is_dir() or parent.is_symlink():
            print(f"[SKIP] sample directory missing or symlink: {parent}",file=sys.stderr)
            continue
        figures = [
            parent / "plots",
            parent / "plots_r",
            parent / "plots_reanalysis",
            parent / "downstream_filtering" / "plots",
            parent / "diagnostic_review" / "samplot",
            parent / "diagnostic_review" / f"{s}.samplot_manifest.tsv",
        ]
        if scope in ("plots", "all"):
            paths.extend(figures)
        if scope in ("generated", "all"):
            paths.extend([
                parent / "gene_discovery" / "final",
                parent / "gene_discovery" / "further_filtering",
                parent / "gene_discovery" / ".postprocess_tmp",
                parent / "genomiser",
                parent / "genomiser_cadd",
                parent / "benchmark",
                parent / "ranking_stability",
                parent / "logs",
                parent / "snp_clair3" / f"{s}.longphase.vep.txt",
            ])
            paths.extend(parent / "gene_discovery" / pattern.format(s=s)
                         for pattern in GENERATED_GENE_FILES)
    if scope in ("plots", "all"):
        paths.append(root / "cohort_plots")
    if scope in ("generated", "all"):
        paths += [
            root / "cohort_analysis" / "sv_recurrence.tsv",
            root / "cohort_analysis" / "sv_recurrence_members.tsv",
            root / "cohort_analysis" / "reference" / "hon_hpo_reference.tsv",
        ]
    found = sorted({p for p in paths if (p.exists() or p.is_symlink())
                    and not p.is_symlink()}, key=str)
    for p in found:
        if any(other in p.parents for other in found):
            raise ValueError(f"Overlapping candidate targets: {p}")
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True,
                        help="Server or destination disk RESULTS folder")
    parser.add_argument("--samples",nargs="+",default=list(DEFAULT_SAMPLES))
    parser.add_argument("--scope",choices=["plots","generated","all"],default="plots")
    parser.add_argument("--delete",action="store_true",
                        help="Permanently delete selected generated outputs, not original inputs")
    parser.add_argument("--confirm",default=None,
                        help="Required with --delete: DELETE_OBSOLETE_OUTPUTS")
    parser.add_argument("--report",type=Path,default=None)
    args=parser.parse_args()
    root=args.root.expanduser().resolve()
    if not root.is_dir() or root in (Path("/"),Path.home().resolve()):
        parser.error("Provide an existing non-root RESULTS directory")
    if args.delete and args.confirm != "DELETE_OBSOLETE_OUTPUTS":
        parser.error("--delete requires --confirm DELETE_OBSOLETE_OUTPUTS")
    if not args.samples or any(not x.isdecimal() for x in args.samples):
        parser.error("Sample IDs must be numeric")
    targets=candidates(root,args.samples,args.scope)
    report=args.report or root.parent / (
        "old_thesis_cleanup_" + datetime.now().strftime("%Y%m%d_%H%M%S") + ".tsv"
    )
    if report.resolve()==root:
        parser.error("Report must be a file path")
    print(f"RESULTS root: {root}")
    print(f"Scope: {args.scope}; targets: {len(targets)}")
    for p in targets:
        print(f"  {'DIR' if p.is_dir() else 'FILE'}  {p.relative_to(root)}")
    if not args.delete:
        print("PREVIEW ONLY. No changes were made.")
        return 0
    report.parent.mkdir(parents=True,exist_ok=True)
    with report.open("x",newline="",encoding="utf-8") as handle:
        writer=csv.writer(handle,delimiter="\t")
        writer.writerow(["path_relative_to_results","status"])
        for p in targets:
            rel=p.relative_to(root)
            try:
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
                writer.writerow([str(rel),"DELETED"])
                handle.flush()
            except Exception as exc:
                writer.writerow([str(rel),f"FAILED: {exc}"])
                handle.flush()
                raise
    print(f"Deleted {len(targets)} named old outputs. Audit: {report}")
    print("Original sequencing and archived evidence folders were not selected.")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
