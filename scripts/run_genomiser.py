#!/usr/bin/env python3
"""Joint genome (SNV/indel + SV) prioritization with Genomiser / Exomiser CLI 15.1.

Run in a per-sample directory; do not modify the callers' original VCF files.
"""
import argparse
import hashlib
import json
import re
import shlex
import subprocess
from pathlib import Path

import yaml


def run(args, **kwargs):
    subprocess.run([str(item) for item in args], check=True, **kwargs)


def samples(vcf):
    return subprocess.check_output(["bcftools", "query", "-l", str(vcf)], text=True).splitlines()


def checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# Do not silently revert to the exome preset or exclude the noncoding genome.
NONCODING_EFFECTS = {
    "CODING_TRANSCRIPT_INTRON_VARIANT", "NON_CODING_TRANSCRIPT_INTRON_VARIANT",
    "NON_CODING_TRANSCRIPT_EXON_VARIANT", "INTERGENIC_VARIANT",
    "REGULATORY_REGION_VARIANT", "UPSTREAM_GENE_VARIANT",
    "DOWNSTREAM_GENE_VARIANT", "FIVE_PRIME_UTR_INTRON_VARIANT",
    "THREE_PRIME_UTR_INTRON_VARIANT",
}


def validate_analysis(path):
    """Check the editable YAML before launching a memory-intensive genome run."""
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        analysis = yaml.safe_load(handle)
    if not isinstance(analysis, dict):
        raise ValueError(f"{path}: expected a YAML analysis mapping")
    if analysis.get("analysisMode") not in {"PASS_ONLY", "FULL"}:
        raise ValueError("Genomiser analysisMode must be PASS_ONLY or FULL")
    for field in ("inheritanceModes", "frequencySources", "pathogenicitySources", "steps"):
        if not analysis.get(field):
            raise ValueError(f"{path}: analysis must define nonempty {field}")
    if "REMM" not in analysis["pathogenicitySources"]:
        raise ValueError(
            "Genomiser noncoding/regulatory prioritization requires REMM "
            "in pathogenicitySources"
        )
    steps = analysis["steps"]
    if not isinstance(steps, list) or any(
        not isinstance(step, dict) or len(step) != 1 for step in steps
    ):
        raise ValueError("Genomiser steps must be a list of named filter mappings")
    for step in steps:
        name, options = next(iter(step.items()))
        if name in {"intervalFilter", "genePanelFilter"}:
            raise ValueError(
                f"{name} would restrict genome-wide prioritization; "
                "remove it from the Genomiser analysis YAML"
            )
        if name == "variantEffectFilter":
            if not isinstance(options, dict):
                raise ValueError("variantEffectFilter must be a YAML mapping")
            removed = set(options.get("remove", []) or [])
            excluded = sorted(removed & NONCODING_EFFECTS)
            if excluded:
                raise ValueError(
                    "The analysis removes intronic/regulatory classes: "
                    + ", ".join(excluded)
                )
    if not any("frequencyFilter" in step for step in steps):
        raise ValueError("Include frequencyFilter to retain population-frequency evidence")
    return analysis


def prepare_vcf(vcf, sample, target, namefile):
    """Prepare an indexed, single-sample VCF without altering original variant calls."""
    observed = samples(vcf)
    if len(observed) != 1:
        raise ValueError(f"{vcf} must have exactly one genotype sample, found: {observed}")

    target = Path(target)
    if observed[0] == sample:
        # bcftools concat needs an index even if no sample rename is required.
        # Create a local indexed symlink rather than writing into caller outputs.
        if Path(str(vcf) + ".tbi").is_file() or Path(str(vcf) + ".csi").is_file():
            return str(vcf), observed[0]
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(Path(vcf).resolve())
    else:
        Path(namefile).write_text(sample + "\n", encoding="utf-8")
        run(["bcftools", "reheader", "-s", namefile, "-o", target, vcf])

    # reheader does not create a .tbi index. The index must match the *new*
    # renamed VCF, not the original caller VCF.
    run(["bcftools", "index", "-f", "-t", str(target)])
    if samples(target) != [sample]:
        raise ValueError(f"Prepared VCF does not contain sample {sample}: {target}")
    return str(target), observed[0]


def build_java_command(args, combined, options_file):
    """Build the command with an editable genome analysis, never the preset."""
    # The CLI may be invoked from an arbitrary Snakemake working directory.
    # Load the installation properties explicitly, as recommended by Exomiser.
    install_properties = Path(args.jar).resolve().parent / "application.properties"
    java_command = [args.java, f"-Xmx{args.memory_gb}g"]
    if install_properties.is_file():
        java_command.append(f"-Dspring.config.location={install_properties}")
    # Only the optional CADD comparison supplies these absolute GRCh38 paths.
    # Command-line system properties override commented application.properties
    # examples without editing the laboratory's baseline installation.
    if getattr(args, "cadd_snv", None):
        java_command.append(
            f"-Dexomiser.hg38.cadd-snv-path={Path(args.cadd_snv).resolve()}"
        )
    if getattr(args, "cadd_indel", None):
        java_command.append(
            f"-Dexomiser.hg38.cadd-in-del-path={Path(args.cadd_indel).resolve()}"
        )
    java_command += [
        f"-Dexomiser.data-directory={Path(args.data_dir).resolve()}",
        f"-Dexomiser.hg38.data-version={args.genome_version}",
        f"-Dexomiser.phenotype.data-version={args.phenotype_version}",
        f"-Dexomiser.hg38.remm-path={Path(args.remm).resolve()}",
        "-jar", str(Path(args.jar).resolve()), "analyse",
        "--sample", str(Path(args.phenopacket).resolve()),
        "--vcf", str(combined.resolve()), "--assembly", "hg38",
        "--analysis", str(Path(args.analysis).resolve()),
        "--output", str(options_file.resolve())
    ]
    return java_command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--small-vcf", required=True)
    parser.add_argument("--sv-vcf", required=True)
    parser.add_argument("--phenopacket", required=True)
    parser.add_argument("--analysis", required=True, help="Editable Genomiser genome analysis YAML")
    parser.add_argument("--jar", required=True)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--remm", required=True, help="ReMM hg38 BGZF file, with .tbi index")
    parser.add_argument("--cadd-snv", help="Optional GRCh38 CADD SNV scores .tsv.gz")
    parser.add_argument("--cadd-indel", help="Optional GRCh38 CADD indel scores .tsv.gz")
    parser.add_argument("--genome-version", default="2512")
    parser.add_argument("--phenotype-version", default="2512")
    parser.add_argument("--java", default="java")
    parser.add_argument("--memory-gb", type=int, default=32)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    if args.memory_gb < 8:
        raise ValueError("A genome-wide run needs sufficient memory (at least 8 GB configured)")
    analysis = validate_analysis(args.analysis)
    uses_cadd = "CADD" in analysis["pathogenicitySources"]
    if uses_cadd and (not args.cadd_snv or not args.cadd_indel):
        raise ValueError(
            "CADD is enabled in the analysis YAML: provide both --cadd-snv "
            "and --cadd-indel GRCh38 paths, each with a .tbi index."
        )
    if not uses_cadd and (args.cadd_snv or args.cadd_indel):
        raise ValueError("CADD score paths supplied but CADD is not in pathogenicitySources")
    inputs = [Path(args.small_vcf), Path(args.sv_vcf), Path(args.phenopacket),
              Path(args.analysis), Path(args.jar), Path(args.remm), Path(args.remm + ".tbi")]
    if uses_cadd:
        for value in (args.cadd_snv, args.cadd_indel):
            inputs.extend([Path(value), Path(value + ".tbi")])
    genome = Path(args.data_dir) / f"{args.genome_version}_hg38"
    phenotype = Path(args.data_dir) / f"{args.phenotype_version}_phenotype"
    missing = [str(path) for path in inputs if not path.is_file() or path.stat().st_size == 0]
    for folder in (genome, phenotype):
        if not folder.is_dir():
            missing.append(str(folder))
        elif not any(path.is_file() and path.stat().st_size > 0
                     for path in folder.rglob("*.mv.db")):
            missing.append(str(folder) + " (no H2 .mv.db data files)")
    if missing:
        raise FileNotFoundError("Genomiser resources are missing/empty:\n" + "\n".join(missing))
    packet = json.loads(Path(args.phenopacket).read_text(encoding="utf-8"))
    if packet.get("subject", {}).get("id") != args.sample:
        raise ValueError("Phenopacket subject must match the sample ID")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    small, old_small = prepare_vcf(args.small_vcf, args.sample, out / "small.renamed.vcf.gz",
                                   out / "small.sample_name.txt")
    sv, old_sv = prepare_vcf(args.sv_vcf, args.sample, out / "sv.renamed.vcf.gz",
                             out / "sv.sample_name.txt")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.sample):
        raise ValueError("Unsafe sample identifier")
    combined = out / f"{args.sample}.genomiser.input.vcf.gz"
    sort_tmp = out / ".sort_tmp"
    sort_tmp.mkdir(exist_ok=True)
    command = (
        "set -euo pipefail; "
        + "bcftools concat -a -Ou " + shlex.quote(small) + " " + shlex.quote(sv)
        + " | bcftools sort -T " + shlex.quote(str(sort_tmp / "sort"))
        + " -Oz -o " + shlex.quote(str(combined))
    )
    run(["bash", "-c", command])
    run(["bcftools", "index", "-f", "-t", combined])
    if samples(combined) != [args.sample]:
        raise ValueError("Merged VCF does not contain the expected sample ID")

    options_file = out / "output-options.yml"
    prefix = f"{args.sample}.genomiser"
    options_file.write_text(
        "outputContributingVariantsOnly: false\n"
        "numGenes: 0\n"
        "outputDirectory: " + json.dumps(str(out.resolve())) + "\n"
        "outputFileName: " + json.dumps(prefix) + "\n"
        "outputFormats: [TSV_GENE, TSV_VARIANT, HTML]\n",
        encoding="utf-8",
    )
    java_command = build_java_command(args, combined, options_file)
    log = out / f"{prefix}.log"
    with log.open("w", encoding="utf-8") as handle:
        run(java_command, stdout=handle, stderr=subprocess.STDOUT)
    genes = out / f"{prefix}.genes.tsv"
    variants = out / f"{prefix}.variants.tsv"
    if not genes.is_file() or not variants.is_file():
        raise RuntimeError("Genomiser did not produce both TSV files; inspect " + str(log))
    manifest = {
        "sample_id": args.sample,
        "analysis": "GENOME_PRESET_GENOMISER",
        "genome_assembly": "GRCh38",
        "small_vcf": str(Path(args.small_vcf).resolve()),
        "sv_vcf": str(Path(args.sv_vcf).resolve()),
        "original_small_vcf_sample": old_small,
        "original_sv_vcf_sample": old_sv,
        "combined_vcf": str(combined.resolve()),
        "phenopacket_sha256": checksum(args.phenopacket),
        "analysis_yaml": str(Path(args.analysis).resolve()),
        "analysis_yaml_sha256": checksum(args.analysis),
        "analysis_mode": analysis["analysisMode"],
        "pathogenicity_sources": analysis["pathogenicitySources"],
        "analysis_steps": [next(iter(step)) for step in analysis["steps"]],
        "combined_vcf_sha256": checksum(combined),
        "genome_data_version": args.genome_version,
        "phenotype_data_version": args.phenotype_version,
        "remm": str(Path(args.remm).resolve()),
        "cadd_snv": str(Path(args.cadd_snv).resolve()) if uses_cadd else None,
        "cadd_indel": str(Path(args.cadd_indel).resolve()) if uses_cadd else None,
        "genomiser_genes": str(genes.resolve()),
        "genomiser_variants": str(variants.resolve()),
        "notes": "SV interpretation is caller/representation-dependent; BND and complex SV require review."
    }
    (out / f"{prefix}.manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
