# Opt-in Genomiser + CADD sensitivity analysis using archived LRS evidence.
# Does not change the existing clinical ranking or the default 'all' target.
CADD_CONFIG = config.get("genomiser_cadd_analysis_yaml")
CADD_SNV = config.get("genomiser_cadd_snv_hg38")
CADD_INDEL = config.get("genomiser_cadd_indel_hg38")
CADD_CLINICAL = PATH + "{sample}/gene_discovery/final/{sample}_clinical_candidates.genomiser_cadd_comparison.tsv"
CADD_COHORT = PATH + "cohort_analysis/clinical_baseline_vs_cadd.tsv"

if config.get("genomiser_cadd_enabled", False):
    if not CADD_CONFIG or not CADD_SNV or not CADD_INDEL:
        raise ValueError(
            "Genomiser CADD comparison enabled: set genomiser_cadd_analysis_yaml, "
            "genomiser_cadd_snv_hg38 and genomiser_cadd_indel_hg38 to real files. "
            "Example CADD paths in application.properties are commented and not sufficient."
        )

    rule genomiser_genome_cadd:
        input:
            small=GENOMISER_SMALL_VCF,
            sv=GENOMISER_SV_VCF,
            packet=rules.prepare_patient_profile.output.packet,
            analysis=CADD_CONFIG,
            jar=config["genomiser_cli_jar"],
            remm=config["genomiser_remm_hg38"],
            remm_index=config["genomiser_remm_hg38"] + ".tbi",
            cadd_snv=CADD_SNV,
            cadd_snv_index=CADD_SNV + ".tbi",
            cadd_indel=CADD_INDEL,
            cadd_indel_index=CADD_INDEL + ".tbi",
            runner=SCRIPTS + "/run_genomiser.py",
        sv_normalizer=SCRIPTS + "/prepare_genomiser_sv_vcf.py"
        output:
            genes=PATH + "{sample}/genomiser_cadd/{sample}.genomiser.genes.tsv",
            variants=PATH + "{sample}/genomiser_cadd/{sample}.genomiser.variants.tsv",
            manifest=PATH + "{sample}/genomiser_cadd/{sample}.genomiser.manifest.json"
        params:
            outdir=PATH + "{sample}/genomiser_cadd",
            data=config["genomiser_data_dir"],
            genome=config.get("genomiser_genome_data_version", "2512"),
            phenotype=config.get("genomiser_phenotype_data_version", "2512"),
            java=config.get("genomiser_java", "java"),
            memory=int(config.get("genomiser_memory_gb", 32))
        threads: 2
        resources:
            mem_mb=lambda wc: (int(config.get("genomiser_memory_gb", 32)) + 8) * 1024
        conda:
            CONDAENV + "genomiser.yaml"
        log:
            PATH + "{sample}/logs/genomiser_cadd/{sample}.snakemake.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {log:q})"
            python {input.runner:q} --sample {wildcards.sample:q} \
              --small-vcf {input.small:q} --sv-vcf {input.sv:q} \
              --phenopacket {input.packet:q} --analysis {input.analysis:q} \
              --jar {input.jar:q} --data-dir {params.data:q} \
              --remm {input.remm:q} \
              --cadd-snv {input.cadd_snv:q} --cadd-indel {input.cadd_indel:q} \
              --genome-version {params.genome:q} \
              --phenotype-version {params.phenotype:q} \
              --java {params.java:q} --memory-gb {params.memory} \
              --outdir {params.outdir:q} > {log:q} 2>&1
            test -s {output.genes:q}
            test -s {output.variants:q}
            test -s {output.manifest:q}
            """

    rule clinical_genomiser_cadd_comparison:
        input:
            clinical=rules.finalize_candidate_ranking.output.clinical,
            baseline_manifest=rules.genomiser_genome.output.manifest,
            cadd_genes=rules.genomiser_genome_cadd.output.genes,
            cadd_variants=rules.genomiser_genome_cadd.output.variants,
            cadd_manifest=rules.genomiser_genome_cadd.output.manifest,
            script=SCRIPTS + "/merge_genomiser_cadd_comparison.py",
            helper=SCRIPTS + "/annotate_genomiser_comparison.py"
        output:
            CADD_CLINICAL
        conda:
            CONDAENV + "plots.yaml"
        shell:
            """
            python {input.script:q} --clinical {input.clinical:q} \
              --cadd-genes {input.cadd_genes:q} \
              --cadd-variants {input.cadd_variants:q} \
              --baseline-manifest {input.baseline_manifest:q} \
              --cadd-manifest {input.cadd_manifest:q} \
              --output {output:q}
            """

    rule cohort_genomiser_cadd_comparison:
        input:
            tables=expand(CADD_CLINICAL, sample=POSTPROCESS_SAMPLES),
            script=SCRIPTS + "/merge_genomiser_cadd_comparison.py"
        output:
            CADD_COHORT
        params:
            tables=lambda wc: " ".join(
                "--table " + shlex.quote(sample + "=" + CADD_CLINICAL.replace("{sample}", sample))
                for sample in POSTPROCESS_SAMPLES
            )
        conda:
            CONDAENV + "plots.yaml"
        shell:
            "python {input.script:q} {params.tables} --output {output:q}"

    # Explicit target, unlike the original default LRS downstream 'all'.
    rule cadd_comparison_all:
        input:
            CADD_COHORT,
            expand(CADD_CLINICAL, sample=POSTPROCESS_SAMPLES)
