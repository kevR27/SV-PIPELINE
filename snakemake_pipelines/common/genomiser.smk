# Required per-sample genome-wide analysis. Include from each LRS/SRS
# Snakefile after defining GENOMISER_SMALL_VCF and GENOMISER_SV_VCF.
# Missing patient profiles or ReMM resources are errors, never a skip.
rule genomiser_genome:
    input:
        small=GENOMISER_SMALL_VCF,
        sv=GENOMISER_SV_VCF,
        packet=rules.prepare_patient_profile.output.packet,
        analysis=config["genomiser_analysis_yaml"],
        jar=config["genomiser_cli_jar"],
        remm=config["genomiser_remm_hg38"],
        remm_index=config["genomiser_remm_hg38"] + ".tbi",
        runner=SCRIPTS + "/run_genomiser.py",
        sv_normalizer=SCRIPTS + "/prepare_genomiser_sv_vcf.py"
    output:
        genes=PATH + "{sample}/genomiser/{sample}.genomiser.genes.tsv",
        variants=PATH + "{sample}/genomiser/{sample}.genomiser.variants.tsv",
        manifest=PATH + "{sample}/genomiser/{sample}.genomiser.manifest.json"
    params:
        jar=config["genomiser_cli_jar"],
        data=config["genomiser_data_dir"],
        remm=config["genomiser_remm_hg38"],
        genome_version=config.get("genomiser_genome_data_version", "2512"),
        phenotype_version=config.get("genomiser_phenotype_data_version", "2512"),
        java=config.get("genomiser_java", "java"),
        memory=int(config.get("genomiser_memory_gb", 32)),
        outdir=PATH + "{sample}/genomiser"
    threads: 2
    resources:
        mem_mb=lambda wc: (int(config.get("genomiser_memory_gb", 32)) + 8) * 1024
    conda:
        CONDAENV + "genomiser.yaml"
    log:
        PATH + "{sample}/logs/genomiser/{sample}.snakemake.log"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log:q})"
        python {input.runner:q}                 --sample {wildcards.sample:q}                 --small-vcf {input.small:q} --sv-vcf {input.sv:q}                 --phenopacket {input.packet:q} --analysis {input.analysis:q}                 --jar {input.jar:q} --data-dir {params.data:q}                 --remm {input.remm:q}                 --genome-version {params.genome_version:q}                 --phenotype-version {params.phenotype_version:q}                 --java {params.java:q} --memory-gb {params.memory}                 --outdir {params.outdir:q} > {log:q} 2>&1
        test -s {output.genes:q}
        test -s {output.variants:q}
        test -s {output.manifest:q}
        """
