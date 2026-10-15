# SRS-specific interpretation after genome-wide calling and annotation.
# Reuses the same research-level interpretation scripts as LRS postprocess.
# No LRS Snakefiles or shared scripts are modified.
#
# The original SRS integrated discovery VCF is never replaced by these tables.
# SV + small-variant pairs from this module are candidate reviews, not
# clinically established compound heterozygosity or SV phasing.

SRS_FINAL_DIR = PATH + "{sample}/gene_discovery/final/"
SRS_EVENT_TOL = int(config.get("srs_event_breakpoint_gene_tolerance", 10000))
SRS_LARGE_DEPTH_MIN = int(config.get("srs_large_sv_depth_min_size", 100000))
SRS_PATIENT_HPO = config.get("patient_hpo_files", {}) or {}


def srs_patient_hpo_path(sample):
    entry = SRS_PATIENT_HPO.get(str(sample))
    return str(entry) if entry else None


def srs_optional_patient_hpo(wc):
    path = srs_patient_hpo_path(wc.sample)
    return [path] if path else []


rule srs_interpret_sv_gene_events:
    input:
        integrated=diagnostic_source,
        phenotypes=rules.human_gene_phenotypes.output.tsv,
        hpo_seeds=DISCOVERY_HPO_FILE,
        edges=MONARCH_EDGES,
        script=SCRIPTS + "/rank_sv_gene_events.py",
        ranking=SCRIPTS + "/ranking_common.py"
    output:
        events=SRS_FINAL_DIR + "{sample}_sv_gene_events.tsv.gz"
    params:
        tolerance=SRS_EVENT_TOL
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.events:q})"
        python {input.script:q} --input {input.integrated:q} \
            --phenotypes {input.phenotypes:q} --hpo-seeds {input.hpo_seeds:q} \
            --edges {input.edges:q} --breakpoint-gene-tolerance {params.tolerance} \
            --compact-output --output {output.events:q}
        test -s {output.events:q}
        """


rule srs_build_final_candidate_tables:
    input:
        genes=rules.rank_genomewide_candidates.output.tsv,
        events=rules.srs_interpret_sv_gene_events.output.events,
        script=SCRIPTS + "/build_candidate_tables.py",
        effects=SCRIPTS + "/sv_gene_effects.py",
        common=SCRIPTS + "/ranking_common.py"
    output:
        genes=SRS_FINAL_DIR + "{sample}_gene_candidates.pre_depth.tsv",
        sv=SRS_FINAL_DIR + "{sample}_sv_gene_candidates.pre_depth.tsv"
    params:
        tolerance=SRS_EVENT_TOL
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --gene-ranking {input.genes:q} \
            --sv-gene-events {input.events:q} --breakpoint-near-bp {params.tolerance} \
            --gene-output {output.genes:q} --sv-output {output.sv:q}
        test -s {output.genes:q}
        test -s {output.sv:q}
        """


# Event-centered depth is additional same-BAM evidence for large DEL/DUP.
# It complements CNVpytor and never changes the integrated discovery VCF.
rule srs_large_sv_depth:
    input:
        sv=rules.srs_build_final_candidate_tables.output.sv,
        bins=rules.mosdepth_coverage.output.regions,
        script=SCRIPTS + "/summarize_large_sv_depth.py"
    output:
        summary=SRS_FINAL_DIR + "{sample}_large_sv_depth.tsv",
        bins=SRS_FINAL_DIR + "{sample}_large_sv_depth_bins.tsv"
    params:
        min_size=SRS_LARGE_DEPTH_MIN
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --candidates {input.sv:q} \
            --coverage {input.bins:q} --min-size {params.min_size} \
            --summary-output {output.summary:q} --bins-output {output.bins:q}
        test -s {output.summary:q}
        test -s {output.bins:q}
        """


rule srs_large_sv_gene_context:
    input:
        sv=rules.srs_build_final_candidate_tables.output.sv,
        depth=rules.srs_large_sv_depth.output.summary,
        gene_bed=GENE_BED,
        script=SCRIPTS + "/build_large_sv_gene_context.py"
    output:
        tsv=SRS_FINAL_DIR + "{sample}_large_sv_genes.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --candidates {input.sv:q} \
            --depth-summary {input.depth:q} --gene-bed {input.gene_bed:q} \
            --output {output.tsv:q}
        test -s {output.tsv:q}
        """


# A generic HON HPO reference is generated once per SRS results directory.
# Patient-specific HPO remains optional; absence is recorded, not inferred.
rule srs_build_hon_hpo_reference:
    input:
        panel=GENE_LIST,
        nodes=MONARCH_NODES,
        edges=MONARCH_EDGES,
        seeds=DISCOVERY_HPO_FILE,
        pathways=MITOCARTA_PATHWAYS,
        script=SCRIPTS + "/build_hon_hpo_reference.py"
    output:
        tsv=PATH + "cohort_analysis/reference/hon_hpo_reference.tsv"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv:q})"
        python {input.script:q} --panel {input.panel:q} \
            --nodes {input.nodes:q} --edges {input.edges:q} \
            --hon-seeds {input.seeds:q} --mitopathways-gmx {input.pathways:q} \
            --output {output.tsv:q}
        test -s {output.tsv:q}
        """


rule srs_hpo_semantic_similarity:
    input:
        genes=rules.srs_build_final_candidate_tables.output.genes,
        reference=rules.srs_build_hon_hpo_reference.output.tsv,
        phenotypes=rules.human_gene_phenotypes.output.tsv,
        edges=MONARCH_EDGES,
        patient=srs_optional_patient_hpo,
        script=SCRIPTS + "/hpo_semantic_similarity.py",
        ranking=SCRIPTS + "/ranking_common.py"
    output:
        tsv=SRS_FINAL_DIR + "{sample}_hpo_semantic_similarity.tsv"
    params:
        patient_arg=lambda wc: (
            "--patient-hpo " + shlex.quote(srs_patient_hpo_path(wc.sample))
            if srs_patient_hpo_path(wc.sample) else ""
        )
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --genes {input.genes:q} \
            --reference {input.reference:q} \
            --gene-phenotypes {input.phenotypes:q} --edges {input.edges:q} \
            {params.patient_arg} --output {output.tsv:q}
        test -s {output.tsv:q}
        """


rule srs_integrate_large_sv_depth:
    input:
        genes=rules.srs_build_final_candidate_tables.output.genes,
        sv=rules.srs_build_final_candidate_tables.output.sv,
        depth=rules.srs_large_sv_depth.output.summary,
        hpo=rules.srs_hpo_semantic_similarity.output.tsv,
        script=SCRIPTS + "/integrate_depth_into_candidates.py",
        ranking=SCRIPTS + "/ranking_common.py"
    output:
        genes=SRS_FINAL_DIR + "{sample}_gene_candidates.tsv",
        sv=SRS_FINAL_DIR + "{sample}_sv_gene_candidates.tsv"
    params:
        min_size=SRS_LARGE_DEPTH_MIN
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --gene-candidates {input.genes:q} \
            --sv-candidates {input.sv:q} --depth-summary {input.depth:q} \
            --hpo-similarity {input.hpo:q} --min-size {params.min_size} \
            --gene-output {output.genes:q} --sv-output {output.sv:q}
        test -s {output.genes:q}
        test -s {output.sv:q}
        """


# WhatsHap phases short variants only. There is no LongPhase SV haplotype in
# SRS, so SV+SNV cis/trans phase remains UNRESOLVED in the pair table.
rule srs_snv_sv_candidates:
    input:
        vep=rules.vep_small_variants.output.tsv,
        phased=rules.whatshap_phasing.output.phased,
        sv=rules.srs_integrate_large_sv_depth.output.sv,
        script=SCRIPTS + "/build_snv_sv_candidates.py"
    output:
        tsv=SRS_FINAL_DIR + "{sample}_snv_sv_candidates.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --vep {input.vep:q} \
            --phased-vcf {input.phased:q} --sv-candidates {input.sv:q} \
            --output {output.tsv:q}
        test -s {output.tsv:q}
        """


rule srs_finalize_candidate_ranking:
    input:
        genes=rules.srs_integrate_large_sv_depth.output.genes,
        sv=rules.srs_integrate_large_sv_depth.output.sv,
        pairs=rules.srs_snv_sv_candidates.output.tsv,
        mitochondrial=lambda wc: (
            [PATH + f"{wc.sample}/diagnostic_review/{wc.sample}.mitochondrial_gene_ranking.tsv"]
            if MITOCARTA_ENABLED else []
        ),
        evidence=diagnostic_source,
        script=SCRIPTS + "/finalize_candidate_ranking.py",
        clinical_script=SCRIPTS + "/build_diagnostic_candidates.py",
        ranking=SCRIPTS + "/ranking_common.py"
    output:
        genes=SRS_FINAL_DIR + "{sample}_gene_candidates.ranked.tsv",
        sv=SRS_FINAL_DIR + "{sample}_sv_gene_candidates.ranked.tsv",
        mito=SRS_FINAL_DIR + "{sample}_mitochondrial_gene_ranking.ranked.tsv",
        clinical=SRS_FINAL_DIR + "{sample}_clinical_candidates.ranked.tsv"
    params:
        mito_flag=lambda wc: (
            "--mitochondrial-ranking " + shlex.quote(
                PATH + f"{wc.sample}/diagnostic_review/{wc.sample}.mitochondrial_gene_ranking.tsv"
            ) if MITOCARTA_ENABLED else ""
        )
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --genes {input.genes:q} --sv {input.sv:q} \
            --snv-sv {input.pairs:q} --sample {wildcards.sample:q} \
            {params.mito_flag} --clinical-evidence {input.evidence:q} \
            --clinical-output {output.clinical:q} \
            --gene-output {output.genes:q} --sv-output {output.sv:q} \
            --mito-output {output.mito:q}
        test -s {output.genes:q}
        test -s {output.sv:q}
        test -s {output.mito:q}
        test -s {output.clinical:q}
        """
