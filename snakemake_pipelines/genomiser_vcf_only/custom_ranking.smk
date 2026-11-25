# VCF-only interpretation uses the same scripts as the standard workflows.
# No BAM or caller producer is registered here.
import shlex

VCF_PANEL = config["candidate_genes_list"]
VCF_HON_SEEDS = hon_hpo_seed_path(config, REPO_ROOT)
VCF_NODES = config["monarch_nodes"]
VCF_EDGES = config["monarch_edges"]
VCF_ANNOTSV = config["annotsv_annotations_dir"]
VCF_VEP_CACHE = config["vep_cache_dir"]
VCF_GENE_DIR = PATH + "{sample}/gene_discovery/"
VCF_FINAL = VCF_GENE_DIR + "final/"
VCF_RARE = float(config.get("ranking_rare_af", 0.001))
VCF_MAX_AF = float(config.get("ranking_max_af", 0.01))
VCF_TOL = int(config.get("sv_gene_breakpoint_tolerance", 10000))

rule vcf_only_annotsv:
    input:
        vcf=GENOMISER_SV_VCF,
        panel=VCF_PANEL,
        annotations=VCF_ANNOTSV
    output:
        tsv=PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.tsv",
        unannotated=PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.unannotated.tsv"
    params:
        directory=PATH + "{sample}/sv/annotsv",
        benign=config.get("annotsv_benign_af", 0.01)
    conda:
        CONDAENV + "annotsv.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.directory:q}
        : > {output.unannotated:q}
        AnnotSV -SVinputFile {input.vcf:q} -annotationsDir {input.annotations:q} \
          -outputDir {params.directory:q} -outputFile {wildcards.sample}_merged_SV.annotsv.tsv \
          -genomeBuild GRCh38 -candidateGenesFile {input.panel:q} -candidateGenesFiltering 0 \
          -tx ENSEMBL -REreport 1 -benignAF {params.benign} -hpo ""
        test -s {output.tsv:q}
        """

rule vcf_only_gene_list:
    input:
        annotsv=rules.vcf_only_annotsv.output.tsv,
        script=SCRIPTS + "/extract_annotsv_genes.py"
    output:
        VCF_GENE_DIR + "{sample}_sv_genes.txt"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        "python {input.script:q} --annotsv {input.annotsv:q} --output {output:q}"

rule vcf_only_gene_phenotypes:
    input:
        genes=rules.vcf_only_gene_list.output,
        seeds=VCF_HON_SEEDS,
        nodes=VCF_NODES,
        edges=VCF_EDGES,
        script=SCRIPTS + "/monarch_human_gene_phenotypes.py"
    output:
        VCF_GENE_DIR + "{sample}_human_gene_phenotypes.tsv"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        "python {input.script:q} --genes {input.genes:q} --hpo-anchors {input.seeds:q} --nodes {input.nodes:q} --edges {input.edges:q} --output {output:q}"

rule vcf_only_gene_ranking:
    input:
        annotsv=rules.vcf_only_annotsv.output.tsv,
        genes=rules.vcf_only_gene_list.output,
        pheno=rules.vcf_only_gene_phenotypes.output,
        seeds=VCF_HON_SEEDS,
        edges=VCF_EDGES,
        panel=VCF_PANEL,
        script=SCRIPTS + "/rank_sv_gene_candidates.py",
        common=SCRIPTS + "/ranking_common.py"
    output:
        VCF_GENE_DIR + "{sample}_ranked_candidates.tsv"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        "python {input.script:q} --annotsv {input.annotsv:q} --genes {input.genes:q} --phenotypes {input.pheno:q} --hpo-seeds {input.seeds:q} --edges {input.edges:q} --panel {input.panel:q} --output {output:q}"

rule vcf_only_vep_sv:
    input:
        vcf=GENOMISER_SV_VCF,
        cache=VCF_VEP_CACHE
    output:
        PATH + "{sample}/sv/vep/{sample}.vep.tsv"
    params:
        version=config.get("vep_cache_version", 116),
        assembly=config.get("vep_assembly", "GRCh38")
    threads: 4
    conda:
        CONDAENV + "vep.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})"
        vep --input_file {input.vcf:q} --format vcf --output_file {output:q} \
            --cache --offline --dir_cache {input.cache:q} --cache_version {params.version} \
            --assembly {params.assembly:q} --tab \
            --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,BIOTYPE,EXON,INTRON,CANONICAL,PICK,VARIANT_CLASS,OverlapBP,OverlapPC,DISTANCE,STRAND" \
            --numbers --variant_class --canonical --overlaps --flag_pick --symbol --biotype \
            --max_sv_size 250000000 --buffer_size 100 --force_overwrite --fork {threads}
        """

rule vcf_only_vep_sv_summary:
    input:
        vcf=GENOMISER_SV_VCF,
        vep=rules.vcf_only_vep_sv.output,
        script=SCRIPTS + "/summarize_vep_sv_consequences.py"
    output:
        PATH + "{sample}/sv/vep/{sample}.vep.summary.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        "python {input.script:q} --vcf {input.vcf:q} --vep {input.vep:q} --output {output:q}"

rule vcf_only_annotsv_audit:
    input:
        annotsv=rules.vcf_only_annotsv.output.tsv,
        script=SCRIPTS + "/audit_annotsv_evidence.py"
    output:
        json=PATH + "{sample}/sv/annotsv/{sample}_annotation_availability.json",
        tsv=PATH + "{sample}/sv/annotsv/{sample}_annotation_availability.tsv"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        "python {input.script:q} --annotsv {input.annotsv:q} --annotations-dir {VCF_ANNOTSV:q} --output-json {output.json:q} --output-tsv {output.tsv:q}"

rule vcf_only_integrated_evidence:
    input:
        vcf=GENOMISER_SV_VCF,
        annotsv=rules.vcf_only_annotsv.output.tsv,
        unannotated=rules.vcf_only_annotsv.output.unannotated,
        audit=rules.vcf_only_annotsv_audit.output.json,
        ranking=rules.vcf_only_gene_ranking.output,
        panel=VCF_PANEL,
        vep=rules.vcf_only_vep_sv_summary.output,
        script=SCRIPTS + "/build_integrated_sv_gene_tsv.py",
        common=SCRIPTS + "/sv_evidence_common.py",
        reconciliation=SCRIPTS + "/annotsv_reconciliation.py"
    output:
        VCF_GENE_DIR + "{sample}_integrated_SV_gene_analysis.tsv"
    params:
        platform=lambda wc: VCF_INPUTS[wc.sample].get("sequencing_type", "LRS")
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} --vcf {input.vcf:q} --annotsv {input.annotsv:q} \
          --annotsv-unannotated {input.unannotated:q} --annotsv-audit {input.audit:q} \
          --ranking {input.ranking:q} --panel {input.panel:q} --vep-summary {input.vep:q} \
          --sequencing-type {params.platform:q} --output {output:q}
        """

include: "../common/allele_assessment.smk"

rule vcf_only_population:
    input:
        events=rules.assess_sv_alleles.output.tsv,
        vcf=lambda wc: [config["gnomad_sv_vcf"]] if config.get("gnomad_sv_vcf") else [],
        index=lambda wc: [config.get("gnomad_sv_index") or config["gnomad_sv_vcf"] + ".tbi"] if config.get("gnomad_sv_vcf") else [],
        script=SCRIPTS + "/annotate_gnomad_sv_exact.py"
    output:
        VCF_GENE_DIR + "{sample}.gnomad.delta.tsv.gz"
    params:
        args=lambda wc, input: ("--gnomad-vcf " + shlex.quote(str(input.vcf[0])) + " --gnomad-index " + shlex.quote(str(input.index[0]))) if input.vcf else ""
    conda:
        CONDAENV + "gnomad_sv.yaml"
    shell:
        "python {input.script:q} --input {input.events:q} {params.args} --annotation-only --output {output:q}"

rule vcf_only_events:
    input:
        evidence=rules.assess_sv_alleles.output.tsv,
        population=rules.vcf_only_population.output,
        pheno=rules.vcf_only_gene_phenotypes.output,
        seeds=VCF_HON_SEEDS,
        edges=VCF_EDGES,
        script=SCRIPTS + "/rank_sv_gene_events.py",
        common=SCRIPTS + "/ranking_common.py"
    output:
        VCF_GENE_DIR + "{sample}.events.tsv.gz"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        python {input.script:q} --input {input.evidence:q} --gnomad-delta {input.population:q} \
          --phenotypes {input.pheno:q} --hpo-seeds {input.seeds:q} --edges {input.edges:q} \
          --rare-af {VCF_RARE} --max-af {VCF_MAX_AF} --breakpoint-gene-tolerance {VCF_TOL} \
          --compact-output --output {output:q}
        """

rule vcf_only_mitocarta:
    input:
        events=rules.vcf_only_events.output,
        ranking=rules.vcf_only_gene_ranking.output,
        mito=config["mitocarta_file"],
        pathways=config["mitocarta_pathways_file"],
        script=SCRIPTS + "/annotate_mitocarta_context.py"
    output:
        VCF_GENE_DIR + "{sample}_integrated_SV_gene_analysis.final.tsv.gz"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        "python {input.script:q} --integrated {input.events:q} --ranking {input.ranking:q} --mitocarta {input.mito:q} --pathways-gmx {input.pathways:q} --output {output:q}"

rule vcf_only_candidates:
    input:
        genes=rules.vcf_only_gene_ranking.output,
        events=rules.vcf_only_mitocarta.output,
        script=SCRIPTS + "/build_candidate_tables.py",
        common=SCRIPTS + "/ranking_common.py",
        effects=SCRIPTS + "/sv_gene_effects.py"
    output:
        genes=VCF_FINAL + "{sample}_gene_candidates.pre_depth.tsv",
        sv=VCF_FINAL + "{sample}_sv_gene_candidates.pre_depth.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        "python {input.script:q} --gene-ranking {input.genes:q} --sv-gene-events {input.events:q} --breakpoint-near-bp {VCF_TOL} --gene-output {output.genes:q} --sv-output {output.sv:q}"

rule vcf_only_hon_reference:
    input:
        panel=VCF_PANEL,
        seeds=VCF_HON_SEEDS,
        nodes=VCF_NODES,
        edges=VCF_EDGES,
        pathways=config["mitocarta_pathways_file"],
        script=SCRIPTS + "/build_hon_hpo_reference.py"
    output:
        PATH + "cohort_analysis/reference/hon_hpo_reference.tsv"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        "python {input.script:q} --panel {input.panel:q} --hon-seeds {input.seeds:q} --nodes {input.nodes:q} --edges {input.edges:q} --mitopathways-gmx {input.pathways:q} --output {output:q}"

rule vcf_only_patient_hpo:
    input:
        genes=rules.vcf_only_candidates.output.genes,
        reference=rules.vcf_only_hon_reference.output,
        pheno=rules.vcf_only_gene_phenotypes.output,
        edges=VCF_EDGES,
        patient=rules.prepare_patient_profile.output.hpo,
        script=SCRIPTS + "/hpo_semantic_similarity.py",
        common=SCRIPTS + "/ranking_common.py"
    output:
        VCF_FINAL + "{sample}_hpo_semantic_similarity.tsv"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        "python {input.script:q} --genes {input.genes:q} --reference {input.reference:q} --gene-phenotypes {input.pheno:q} --edges {input.edges:q} --patient-hpo {input.patient:q} --output {output:q}"

rule vcf_only_missing_depth:
    output:
        VCF_FINAL + "{sample}_depth_not_available.tsv"
    shell:
        "mkdir -p $(dirname {output:q}); printf 'SV_ID\\tDEPTH_PATTERN\\n' > {output:q}"

rule vcf_only_patient_candidates:
    input:
        genes=rules.vcf_only_candidates.output.genes,
        sv=rules.vcf_only_candidates.output.sv,
        hpo=rules.vcf_only_patient_hpo.output,
        depth=rules.vcf_only_missing_depth.output,
        script=SCRIPTS + "/integrate_depth_into_candidates.py",
        common=SCRIPTS + "/ranking_common.py"
    output:
        genes=VCF_FINAL + "{sample}_gene_candidates.depth.tsv",
        sv=VCF_FINAL + "{sample}_sv_gene_candidates.depth.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        "python {input.script:q} --gene-candidates {input.genes:q} --sv-candidates {input.sv:q} --depth-summary {input.depth:q} --hpo-similarity {input.hpo:q} --gene-output {output.genes:q} --sv-output {output.sv:q}"

rule vcf_only_vep_small:
    input:
        vcf=GENOMISER_SMALL_VCF,
        cache=VCF_VEP_CACHE
    output:
        PATH + "{sample}/small_variants/{sample}.vep.tsv"
    params:
        version=config.get("vep_cache_version", 116),
        assembly=config.get("vep_assembly", "GRCh38")
    threads: 4
    conda:
        CONDAENV + "vep.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})"
        vep --input_file {input.vcf:q} --format vcf --output_file {output:q} \
            --cache --offline --dir_cache {input.cache:q} --cache_version {params.version} \
            --assembly {params.assembly:q} --canonical --flag_pick --symbol --biotype \
            --variant_class --regulatory --allele_number --uploaded_allele --force_overwrite --fork {threads}
        """

rule vcf_only_snv_sv:
    input:
        vep=rules.vcf_only_vep_small.output,
        small=GENOMISER_SMALL_VCF,
        sv=rules.vcf_only_patient_candidates.output.sv,
        script=SCRIPTS + "/build_snv_sv_candidates.py"
    output:
        VCF_FINAL + "{sample}_snv_sv_candidates.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        "python {input.script:q} --vep {input.vep:q} --phased-vcf {input.small:q} --sv-candidates {input.sv:q} --output {output:q}"

rule vcf_only_final_ranking:
    input:
        genes=rules.vcf_only_patient_candidates.output.genes,
        sv=rules.vcf_only_patient_candidates.output.sv,
        pairs=rules.vcf_only_snv_sv.output,
        evidence=rules.vcf_only_mitocarta.output,
        genomiser_genes=rules.genomiser_genome.output.genes,
        genomiser_variants=rules.genomiser_genome.output.variants,
        manifest=rules.genomiser_genome.output.manifest,
        script=SCRIPTS + "/finalize_candidate_ranking.py",
        clinical=SCRIPTS + "/build_diagnostic_candidates.py",
        common=SCRIPTS + "/ranking_common.py",
        comparison=SCRIPTS + "/annotate_genomiser_comparison.py",
        benchmark=SCRIPTS + "/benchmark_genomiser.py"
    output:
        genes=VCF_FINAL + "{sample}_gene_candidates.ranked.tsv",
        sv=VCF_FINAL + "{sample}_sv_gene_candidates.ranked.tsv",
        clinical=VCF_FINAL + "{sample}_clinical_candidates.ranked.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        python {input.script:q} --genes {input.genes:q} --sv {input.sv:q} \
          --snv-sv {input.pairs:q} --sample {wildcards.sample:q} \
          --clinical-evidence {input.evidence:q} --clinical-output {output.clinical:q} \
          --genomiser-genes {input.genomiser_genes:q} --genomiser-variants {input.genomiser_variants:q} \
          --gene-output {output.genes:q} --sv-output {output.sv:q}
        """
