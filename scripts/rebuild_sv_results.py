#!/usr/bin/env python3
"""Rebuild SV evidence tables from saved outputs, without rerunning variant callers.

Use a new output directory. The original VCFs and annotations stay unchanged.
The four optional clinical tables are not required by this command.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--annotsv", required=True)
    p.add_argument("--unannotated")
    p.add_argument("--panel", required=True)
    p.add_argument("--phenotypes", help="Previously saved human gene/HPO table; no network lookup")
    p.add_argument("--hpo-seeds", help="HON/phenotype HPO seed TSV used by the current ranker")
    p.add_argument("--monarch-edges", help="Monarch edge TSV used for semantic HPO ranking")
    p.add_argument("--needlr")
    p.add_argument("--vep", help="Saved VEP text output, for a gene coverage comparison")
    p.add_argument("--annotations-dir", help="Optional installed AnnotSV resources")
    p.add_argument("--caller-vcf", action="append", default=[], metavar="CALLER=PATH")
    p.add_argument("--caller-order", default="Sniffles2,cuteSV,Delly", help="Exact original merge input order")
    p.add_argument("--sequencing-type", choices=["LRS", "SRS"], default="LRS")
    p.add_argument("--min-support", type=float, default=2)
    p.add_argument("--min-svlen", type=float, default=50)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", a.sample) or a.sample in {".", ".."}:
        p.error("Use a sample ID without slashes or spaces")
    callers = []
    for item in a.caller_vcf:
        name, sep, path = item.partition("=")
        if not sep or not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name not in a.caller_order.split(","):
            p.error("Each caller VCF needs CALLER=PATH, with CALLER in --caller-order")
        if name in [c for c, _ in callers]:
            p.error("Supply only one VCF per caller")
        callers.append((name, path))
    inputs = {k: v for k, v in vars(a).items() if k in {"vcf", "annotsv", "unannotated", "panel", "phenotypes", "hpo_seeds", "monarch_edges", "needlr", "vep"} and v}
    inputs.update({f"caller_{name}": path for name, path in callers})
    for path in inputs.values():
        if not Path(path).is_file():
            p.error(f"Input file not found: {path}")
    out = Path(a.outdir)
    if out.exists():
        p.error("Choose a new --outdir so an earlier result cannot be overwritten")
    out.mkdir(parents=True)
    scripts = Path(__file__).resolve().parent
    commands = []
    prefix = out / a.sample

    def product(suffix):
        return str(prefix) + suffix

    def run(name, *options):
        command = [sys.executable, str(scripts / name), *map(str, options)]
        commands.append(command)
        subprocess.run(command, check=True)

    provenance = {key: {"path": str(Path(path).resolve()), "sha256": sha256(path)} for key, path in inputs.items()}
    run("audit_sv_outputs.py", "--vcf", a.vcf, "--annotsv", a.annotsv,
        *(["--unannotated", a.unannotated] if a.unannotated else []),
        "--output-tsv", product("_annotation_records.tsv"), "--output-json", product("_annotation_audit.json"))
    if a.vep:
        run("audit_vep_gene_coverage.py", "--vep", a.vep, "--annotsv", a.annotsv,
            "--records", product("_annotation_records.tsv"), "--output-tsv", product("_vep_gene_differences.tsv"),
            "--output-json", product("_vep_coverage.json"))
    run("audit_annotsv_evidence.py", "--annotsv", a.annotsv,
        *(["--annotations-dir", a.annotations_dir] if a.annotations_dir else []),
        "--output-json", product("_annotation_availability.json"), "--output-tsv", product("_annotation_availability.tsv"))
    run("extract_annotsv_genes.py", "--annotsv", a.annotsv, "--output", product("_sv_genes.txt"))
    run("extract_nonpanel_genes.py", "--genes", product("_sv_genes.txt"), "--panel", a.panel,
        "--output", product("_nonpanel_genes.txt"))
    phenotypes = a.phenotypes
    if not phenotypes:
        phenotypes = product("_human_gene_phenotypes_not_supplied.tsv")
        Path(phenotypes).write_text("gene_symbol\thpo_id\toptic_neuropathy_anchor\n")
    hpo_seeds = a.hpo_seeds
    if not hpo_seeds:
        hpo_seeds = product("_hpo_seeds_not_supplied.tsv")
        Path(hpo_seeds).write_text("hpo_id\thpo_label\thon_seed_role\n")
    monarch_edges = a.monarch_edges
    if not monarch_edges:
        monarch_edges = product("_monarch_edges_not_supplied.tsv")
        Path(monarch_edges).write_text("subject\tobject\tpredicate\tcategory\n")
    run("rank_sv_gene_candidates.py", "--annotsv", a.annotsv, "--genes", product("_sv_genes.txt"),
        "--phenotypes", phenotypes, "--panel", a.panel, "--hpo-seeds", hpo_seeds,
        "--edges", monarch_edges, "--output", product("_ranked_candidates.tsv"))
    run("summarize_sv_caller_support.py", "--vcf", a.vcf, "--caller-order", a.caller_order,
        "--output", product("_caller_support_summary.tsv"))
    caller_args = []
    callers.sort(key=lambda entry: a.caller_order.split(",").index(entry[0]))
    for name, path in callers:
        parsed, filtered = product(f"_{name}_caller.tsv"), product(f"_{name}_evidence.tsv")
        run("parse_sv_caller_vcf.py", "--vcf", path, "--caller", name, "--output", parsed)
        run("filter_sv_evidence.py", "--caller-tsv", parsed, "--caller", name,
            "--min-support", a.min_support, "--min-svlen", a.min_svlen,
            "--output-tsv", filtered, "--output-ids", product(f"_{name}_passing_ids.txt"))
        # Do not filter the master with these IDs: it is the original callset.
        caller_args.extend(["--caller-tsv", filtered])
    if a.sequencing_type == "LRS" and len(callers) == len(a.caller_order.split(",")):
        caller_args.append("--jasmine-id-prefixes")
    run("build_integrated_sv_gene_tsv.py", "--vcf", a.vcf, "--annotsv", a.annotsv,
        "--sequencing-type", a.sequencing_type,
        "--annotsv-audit", product("_annotation_availability.json"),
        *(["--annotsv-unannotated", a.unannotated] if a.unannotated else []),
        *(["--needlr", a.needlr] if a.needlr else []),
        "--ranking", product("_ranked_candidates.tsv"), "--panel", a.panel,
        "--caller-summary", product("_caller_support_summary.tsv"), *caller_args,
        "--output", product("_integrated_SV_gene_analysis.tsv"))
    run("assess_sv_alleles.py", "--integrated", product("_integrated_SV_gene_analysis.tsv"),
        "--sample", a.sample, "--phenotypes", phenotypes,
        "--output", product("_allele_assessment.tsv"),
        "--hypotheses-output", product("_allele_hypotheses.tsv"), "--manifest", product("_allele_manifest.json"))
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    review_fields = ["SV_ID", "CHROM", "START", "END", "CHR2", "POS2", "SVTYPE", "SVLEN",
                     "GENES", "PANEL_STATUS", "SV_SIZE_SCOPE", "ANNOTSV_RECORD_STATUS", "ANNOTSV_UNANNOTATED_REASON",
                     "ANNOTSV_GENE_MAPPING_STATUS", "ANNOTSV_INVALID_GENE_LABELS", "CALLERS", "CALLER_COUNT",
                     "CALLER_READ_SUPPORT", "NEEDLR_AF", "NEEDLR_STATUS", "NEEDLR_AF_SOURCE_FIELD",
                     "PHENOTYPE_SCORE", "INTEGRATED_DISCOVERY_SCORE", "GENE_DISEASE_EVIDENCE_SCORE",
                     "CANDIDATE_CLASS", "OMIM", "OMIM_INHERITANCE", "GENCC_DISEASE", "GENCC_MOI",
                     "CLINGEN_HI", "CLINGEN_TS", "ACMG_CNV_CLASS", "ACMG_CNV_SCORE",
                     "SV_PATHOGENIC_DB_SOURCE", "SV_PATHOGENIC_DB_STATUS", "SV_BENIGN_DB_SOURCE", "SV_BENIGN_DB_STATUS",
                     "ALLELE_FUNCTIONAL_CONTEXT", "ALLELE_GENOTYPE", "ALLELE_DISRUPTION_STATUS",
                     "ALLELE_TECHNICAL_STATUS", "ALLELE_POPULATION_STATUS", "ALLELE_PHENOTYPE_STATUS",
                     "ALLELE_INHERITANCE_STATUS", "ALLELE_DISEASE_MECHANISM_STATUS", "ALLELE_RESEARCH_SCORE",
                     "ALLELE_UNKNOWN_DOMAINS", "ALLELE_REVIEW_FLAGS"]
    with open(product("_allele_assessment.tsv")) as src, gzip.open(product("_SV_gene_review.tsv.gz"), "wt", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, review_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            writer.writerow(row)
    manifest = {"version": "saved-output-rebuild-v1", "sample": a.sample, "inputs": provenance,
                "commands": commands, "parameters": vars(a),
                "clinical_tables": "NOT_SUPPLIED", "master_records_removed": 0,
                "limitations": ["Saved annotations can contain missing or incorrect gene names.",
                                "Gene/HPO associations are not patient-specific phenotype matches.",
                                "No mitochondrial pathway database was supplied or queried.",
                                "Filtered caller VCFs cannot recover calls discarded in the original run."],
                "scripts_sha256": {f.name: sha256(f) for f in scripts.glob("*.py")}}
    Path(product("_rebuild_manifest.json")).write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[OK] Results written to {out}. Original input files were not changed.")


if __name__ == "__main__":
    main()
