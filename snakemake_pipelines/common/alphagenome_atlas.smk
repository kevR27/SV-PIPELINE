# Shared optional RESEARCH annotation. Source alleles come from Clair3 or
# DeepVariant, never from the restricted coding/splice SNV+SV review table.
import json
from pathlib import Path

ATLAS_ENABLED = config.get("alphagenome_atlas_enabled", False)
ATLAS_NATIVE_ENABLED = config.get("alphagenome_native_plots_enabled", False)
for atlas_setting in (ATLAS_ENABLED, ATLAS_NATIVE_ENABLED, config.get("alphagenome_native_cache_only", False)):
    if not isinstance(atlas_setting, bool):
        raise ValueError("AlphaGenome enablement settings must be true/false without quotes.")
ATLAS_AVI = config.get("alphagenome_avi_file")
ATLAS_RELEASE = config.get("alphagenome_avi_release")
if ATLAS_ENABLED and config.get("vep_assembly", "GRCh38") != "GRCh38":
    raise ValueError("AlphaGenome Atlas requires GRCh38 source/reference annotations.")
if ATLAS_ENABLED and (not ATLAS_AVI or not ATLAS_RELEASE):
    raise ValueError("Set alphagenome_avi_file and alphagenome_avi_release after the AVI download finishes. See README.md.")
if ATLAS_NATIVE_ENABLED and not ATLAS_ENABLED:
    raise ValueError("Enable the offline Atlas annotation before requesting its selected native plots.")
ATLAS_SELECTIONS = {str(s): str(p) for s, p in (config.get("alphagenome_plot_selection_files", {}) or {}).items()}
if set(ATLAS_SELECTIONS) - set(SAMPLES):
    raise ValueError("AlphaGenome plot selection contains samples absent from samples in the config.")
ATLAS_NATIVE_SAMPLES = [s for s in ATLAS_SAMPLES if s in ATLAS_SELECTIONS] if ATLAS_NATIVE_ENABLED else []
if ATLAS_NATIVE_ENABLED and not ATLAS_NATIVE_SAMPLES:
    raise ValueError("Native plots need an alphagenome_plot_selection_files entry for a sample being analysed.")
ATLAS_TISSUES = config.get("alphagenome_prediction_ontology_terms", [])
if ATLAS_NATIVE_ENABLED and (not isinstance(ATLAS_TISSUES, list) or not ATLAS_TISSUES):
    raise ValueError("Choose actual AlphaGenome tissue ontology IDs before enabling native plots.")
ATLAS_MIN_GQ = float(config.get("alphagenome_prediction_min_gq", 20))
ATLAS_MIN_DP = int(config.get("alphagenome_prediction_min_dp", 10))
ATLAS_ANNOTATION_TARGETS = expand(PATH + "{sample}/alphagenome_atlas/{sample}.summary.json", sample=ATLAS_SAMPLES) if ATLAS_ENABLED else []
ATLAS_OVERVIEW_TARGETS = expand(PATH + "{sample}/plots/17_alphagenome_atlas/atlas_overview.{ext}", sample=ATLAS_SAMPLES, ext=["pdf", "svg", "png", "source.tsv"]) if ATLAS_ENABLED else []
ATLAS_NATIVE_TARGETS = expand(PATH + "{sample}/plots/17_alphagenome_atlas/native_manifest.tsv", sample=ATLAS_NATIVE_SAMPLES)


rule alphagenome_atlas:
    input:
        vcf=ATLAS_SMALL_VCF,
        sv=ATLAS_SV_VCF,
        vep=ATLAS_VEP_TSV,
        ref=REF,
        fai=REF + ".fai",
        avi=lambda wc: ATLAS_AVI if ATLAS_ENABLED else [],
        index=lambda wc: ATLAS_AVI + ".tbi" if ATLAS_ENABLED else [],
        genes=GENE_BED,
        panel=GENE_LIST,
        script=SCRIPTS + "/annotate_alphagenome_atlas.py"
    output:
        alleles=PATH + "{sample}/alphagenome_atlas/{sample}.alleles.tsv.gz",
        review=PATH + "{sample}/alphagenome_atlas/{sample}.impact_review.tsv",
        sv=PATH + "{sample}/alphagenome_atlas/{sample}.sv_scope.tsv",
        summary=PATH + "{sample}/alphagenome_atlas/{sample}.summary.json"
    params:
        release=str(ATLAS_RELEASE or ""),
        tissue=lambda wc: str((config.get("dna_sample_tissues", {}) or {}).get(str(wc.sample), config.get("dna_sample_tissue", "blood"))),
        phred=float(config.get("alphagenome_review_phred", 20))
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        python {input.script:q} annotate --vcf {input.vcf:q} --sv-vcf {input.sv:q} \
            --vep {input.vep:q} --reference {input.ref:q} --avi {input.avi:q} \
            --gene-bed {input.genes:q} --panel {input.panel:q} --sample {wildcards.sample:q} \
            --release {params.release:q} --dna-tissue {params.tissue:q} \
            --review-phred {params.phred} --min-gq {ATLAS_MIN_GQ} --min-dp {ATLAS_MIN_DP} \
            --output {output.alleles:q} --review-output {output.review:q} \
            --sv-output {output.sv:q} --summary {output.summary:q}
        """


rule plot_alphagenome_atlas:
    input:
        summary=rules.alphagenome_atlas.output.summary,
        script=SCRIPTS + "/annotate_alphagenome_atlas.py",
        utilities=PLOTS + "/plot_utils.py"
    output:
        pdf=PATH + "{sample}/plots/17_alphagenome_atlas/atlas_overview.pdf",
        svg=PATH + "{sample}/plots/17_alphagenome_atlas/atlas_overview.svg",
        png=PATH + "{sample}/plots/17_alphagenome_atlas/atlas_overview.png",
        source=PATH + "{sample}/plots/17_alphagenome_atlas/atlas_overview.source.tsv"
    params:
        outdir=PATH + "{sample}/plots/17_alphagenome_atlas"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        "python {input.script:q} plot --summary {input.summary:q} --outdir {params.outdir:q}"


rule alphagenome_native_plots:
    input:
        context=rules.alphagenome_atlas.output.summary,
        vcf=ATLAS_SMALL_VCF,
        sv=ATLAS_SV_VCF,
        ref=REF,
        fai=REF + ".fai",
        selection=lambda wc: ATLAS_SELECTIONS[str(wc.sample)],
        script=SCRIPTS + "/annotate_alphagenome_atlas.py",
        utilities=PLOTS + "/plot_utils.py"
    output:
        manifest=PATH + "{sample}/plots/17_alphagenome_atlas/native_manifest.tsv",
        plots=directory(PATH + "{sample}/plots/17_alphagenome_atlas/native")
    params:
        cache=PATH + "{sample}/alphagenome_atlas/prediction_cache",
        terms=json.dumps(ATLAS_TISSUES),
        max_variants=int(config.get("alphagenome_native_max_variants", 10)),
        length=int(config.get("alphagenome_sequence_length", 1048576)),
        flank=int(config.get("alphagenome_min_flank", 10000)),
        model=str(config.get("alphagenome_model_version", "ALL_FOLDS")),
        cache_only="--cache-only" if config.get("alphagenome_native_cache_only", False) else ""
    conda:
        CONDAENV + "alphagenome.yaml"
    shell:
        """
        python {input.script:q} native-plots --vcf {input.vcf:q} --sv-vcf {input.sv:q} \
            --reference {input.ref:q} --sample {wildcards.sample:q} --selection {input.selection:q} \
            --outdir {output.plots:q} --cache-dir {params.cache:q} --manifest {output.manifest:q} \
            --ontology-terms-json {params.terms:q} --max-variants {params.max_variants} \
            --sequence-length {params.length} --min-flank {params.flank} --model-version {params.model:q} \
            --min-gq {ATLAS_MIN_GQ} --min-dp {ATLAS_MIN_DP} {params.cache_only}
        """


rule atlas_annotation:
    input:
        ATLAS_ANNOTATION_TARGETS


rule all_atlas_plots:
    input:
        ATLAS_OVERVIEW_TARGETS,
        ATLAS_NATIVE_TARGETS
