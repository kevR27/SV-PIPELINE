# Shared by LRS, SRS and the independent LRS postprocessing workflow.
# Included after rule all so the default target is unchanged.
import shlex

ALLELE_KEYS = ("patient_context", "disease_models", "allele_evidence", "family_genotypes")


def allele_input_path(key, sample):
    value = config.get("allele_evidence" if key == "allele_evidence" else "allele_" + key)
    if isinstance(value, dict):
        value = value.get(str(sample))
    return [str(value)] if value else []


def allele_optional_args(wc, input):
    return " ".join("--" + key.replace("_", "-") + " " + shlex.quote(str(getattr(input, key)[0]))
                    for key in ALLELE_KEYS if getattr(input, key))


rule assess_sv_alleles:
    input:
        integrated=PATH + "{sample}/gene_discovery/{sample}_integrated_SV_gene_analysis.tsv",
        phenotypes=PATH + "{sample}/gene_discovery/{sample}_human_gene_phenotypes.tsv",
        patient_context=lambda wc: allele_input_path("patient_context", wc.sample),
        disease_models=lambda wc: allele_input_path("disease_models", wc.sample),
        allele_evidence=lambda wc: allele_input_path("allele_evidence", wc.sample),
        family_genotypes=lambda wc: allele_input_path("family_genotypes", wc.sample),
        script=SCRIPTS + "/assess_sv_alleles.py",
        common=SCRIPTS + "/sv_evidence_common.py"
    output:
        tsv=PATH + "{sample}/gene_discovery/{sample}_allele_assessment.tsv",
        hypotheses=PATH + "{sample}/gene_discovery/{sample}_allele_disease_hypotheses.tsv",
        manifest=PATH + "{sample}/gene_discovery/{sample}_allele_assessment.manifest.json"
    params:
        optional=allele_optional_args,
        rare_af=config.get("allele_rare_af", 0.01),
        min_support=config.get("allele_min_support", 5),
        min_gq=config.get("allele_min_gq", 20),
        min_dp=config.get("allele_min_dp", 5)
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} \
            --integrated {input.integrated:q} --sample {wildcards.sample:q} \
            --phenotypes {input.phenotypes:q} {params.optional} \
            --rare-af {params.rare_af} --min-support {params.min_support} \
            --min-gq {params.min_gq} --min-dp {params.min_dp} \
            --output {output.tsv:q} --hypotheses-output {output.hypotheses:q} \
            --manifest {output.manifest:q}
        """
