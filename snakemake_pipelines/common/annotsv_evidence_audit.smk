# Shared annotation audit and offline Exomiser comparison for LRS and SRS.
import sys
import json
import hashlib
import shlex
from pathlib import Path
sys.path.insert(0, SCRIPTS)
from audit_annotsv_evidence import resource_paths as audit_resource_paths
from audit_annotsv_evidence import fingerprint as audit_fingerprint
from sv_evidence_common import hon_hpo_seed_path
from compare_exomiser_hon import select_hpo_query, backend_resource_paths

# The old switch is read only to migrate existing server configs. It now enables
# the offline CLI; there is no REST fallback. Prefer exomiser_enabled.
EXOMISER_ENABLED = config.get("exomiser_enabled", config.get("annotsv_exomiser_enabled", bool(config.get("hpo_terms"))))
if not isinstance(EXOMISER_ENABLED, bool):
    raise ValueError("exomiser_enabled must be true or false, without quotes.")
if "exomiser_enabled" in config and "annotsv_exomiser_enabled" in config and config["annotsv_exomiser_enabled"] != EXOMISER_ENABLED:
    raise ValueError("Exomiser enablement settings conflict. Remove the old annotsv_exomiser_enabled key.")
EXOMISER_PATIENT_HPO_FILES = {str(sample): path for sample, path in (config.get("patient_hpo_files", {}) or {}).items()}
EXOMISER_HPO_SEEDS = hon_hpo_seed_path(config, REPO_ROOT)
EXOMISER_DEFAULT_DIR = str(Path(ANNOTSV_ANNOTATIONS_DIR) / "jar" / "exomiser-cli-15.1.0")
EXOMISER_JAR = config.get("exomiser_cli_jar", str(Path(EXOMISER_DEFAULT_DIR) / "exomiser-cli-15.1.0.jar"))
EXOMISER_DATA_DIR = config.get("exomiser_data_dir", str(Path(EXOMISER_DEFAULT_DIR) / "data"))
EXOMISER_GENOME_VERSION = str(config.get("exomiser_genome_data_version", "2512"))
EXOMISER_PHENOTYPE_VERSION = str(config.get("exomiser_phenotype_data_version", "2512"))
EXOMISER_JAVA = config.get("exomiser_java", "java")
EXOMISER_MEMORY_GB = int(config.get("exomiser_memory_gb", 8))

# Scores depend on HPO terms, not the sample's variants. Score an identical
# query once, even if several samples use it, then compare per sample.
EXOMISER_SAMPLE_PROFILES = {}
EXOMISER_PROFILE_TERMS = {}
EXOMISER_PROFILE_INPUTS = {}
for sample in SAMPLES:
    patient = EXOMISER_PATIENT_HPO_FILES.get(str(sample))
    terms, scope, source, _ = select_hpo_query(EXOMISER_ENABLED, EXOMISER_HPO_SEEDS, config.get("hpo_terms"), patient)
    if EXOMISER_ENABLED and not terms:
        raise ValueError(f"Exomiser HPO query is empty for {sample}.")
    profile = "hpo_" + hashlib.sha256(json.dumps(terms).encode()).hexdigest()[:16] if EXOMISER_ENABLED else "disabled"
    if profile in EXOMISER_PROFILE_TERMS and EXOMISER_PROFILE_TERMS[profile] != terms:
        raise ValueError("Exomiser query cache identifiers collide; use a longer profile hash.")
    EXOMISER_SAMPLE_PROFILES[str(sample)] = profile
    EXOMISER_PROFILE_TERMS[profile] = terms
    EXOMISER_PROFILE_INPUTS.setdefault(profile, set())
    if source not in {".", "config:hpo_terms"}:
        EXOMISER_PROFILE_INPUTS[profile].add(source)


def exomiser_patient_hpo_inputs(wc):
    path = EXOMISER_PATIENT_HPO_FILES.get(str(wc.sample)) if EXOMISER_ENABLED else None
    return [path] if path else []


def exomiser_resource_inputs(wc):
    if not EXOMISER_ENABLED:
        return []
    # Only real reference files are tracked, never lock/trace/log files.
    # Missing files are explained by the script's extraction preflight.
    return [str(path) for path in backend_resource_paths(EXOMISER_JAR, EXOMISER_DATA_DIR,
            EXOMISER_GENOME_VERSION, EXOMISER_PHENOTYPE_VERSION) if path.is_file()]


def exomiser_profile_file(sample, filename):
    profile = EXOMISER_SAMPLE_PROFILES[str(sample)]
    return PATH + f"reference_checks/exomiser_profiles/{profile}/{filename}"


rule prepare_exomiser_hpo:
    input:
        seeds=EXOMISER_HPO_SEEDS,
        patient=exomiser_patient_hpo_inputs,
        script=SCRIPTS + "/compare_exomiser_hon.py",
        common=SCRIPTS + "/ranking_common.py"
    output:
        tsv=PATH + "{sample}/gene_discovery/{sample}_exomiser_hpo.tsv",
        terms=PATH + "{sample}/gene_discovery/{sample}_exomiser_hpo.txt"
    params:
        enabled="--enabled" if EXOMISER_ENABLED else "",
        shared_hpo=json.dumps(config.get("hpo_terms")),
        patient=lambda wc: "--patient-hpo " + shlex.quote(EXOMISER_PATIENT_HPO_FILES[str(wc.sample)]) if exomiser_patient_hpo_inputs(wc) else ""
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} prepare-hpo {params.enabled} --hpo-seeds {input.seeds:q} \
            --shared-hpo-json {params.shared_hpo:q} {params.patient} \
            --output {output.tsv:q} --terms-output {output.terms:q}
        """


rule exomiser_offline_phenotype:
    input:
        query_files=lambda wc: sorted(EXOMISER_PROFILE_INPUTS[wc.profile]),
        resources=exomiser_resource_inputs,
        script=SCRIPTS + "/compare_exomiser_hon.py",
        common=SCRIPTS + "/ranking_common.py",
        evidence_common=SCRIPTS + "/sv_evidence_common.py"
    output:
        genes=PATH + "reference_checks/exomiser_profiles/{profile}/phenotype.genes.tsv",
        job=PATH + "reference_checks/exomiser_profiles/{profile}/phenotype.job.json",
        properties=PATH + "reference_checks/exomiser_profiles/{profile}/application.properties",
        provenance=PATH + "reference_checks/exomiser_profiles/{profile}/run.json"
    params:
        enabled="--enabled" if EXOMISER_ENABLED else "",
        hpo=lambda wc: json.dumps(EXOMISER_PROFILE_TERMS[wc.profile]),
        jar=EXOMISER_JAR,
        data_dir=EXOMISER_DATA_DIR,
        genome_version=EXOMISER_GENOME_VERSION,
        phenotype_version=EXOMISER_PHENOTYPE_VERSION,
        java=EXOMISER_JAVA,
        memory_gb=EXOMISER_MEMORY_GB,
        resource_snapshot=lambda wc: json.dumps(audit_fingerprint({path: Path(path) for path in exomiser_resource_inputs(wc)}), sort_keys=True)
    # Reserve the configured cores so ordinary --cores 32 runs do not launch
    # many memory-heavy JVMs together. Java's CPU count is capped to this value.
    threads: int(config.get("threads", 32))
    resources:
        mem_mb=EXOMISER_MEMORY_GB * 1024 + 1024
    log:
        PATH + "reference_checks/exomiser_profiles/{profile}/exomiser.log"
    conda:
        CONDAENV + ("exomiser.yaml" if EXOMISER_ENABLED else "monarch.yaml")
    shell:
        """
        python {input.script:q} run-offline {params.enabled} --hpo-json {params.hpo:q} \
            --jar {params.jar:q} --data-dir {params.data_dir:q} \
            --genome-version {params.genome_version:q} --phenotype-version {params.phenotype_version:q} \
            --java {params.java:q} --memory-gb {params.memory_gb} --threads {threads} \
            --genes-output {output.genes:q} --job-output {output.job:q} \
            --properties-output {output.properties:q} --provenance-output {output.provenance:q} --log {log:q}
        """


rule compare_exomiser_hon:
    input:
        query=rules.prepare_exomiser_hpo.output.tsv,
        genes=lambda wc: exomiser_profile_file(wc.sample, "phenotype.genes.tsv"),
        provenance=lambda wc: exomiser_profile_file(wc.sample, "run.json"),
        ranking=PATH + "{sample}/gene_discovery/{sample}_ranked_candidates.tsv",
        phenotypes=PATH + "{sample}/gene_discovery/{sample}_human_gene_phenotypes.tsv",
        edges=MONARCH_EDGES,
        script=SCRIPTS + "/compare_exomiser_hon.py",
        common=SCRIPTS + "/ranking_common.py",
        evidence_common=SCRIPTS + "/sv_evidence_common.py"
    output:
        tsv=PATH + "{sample}/gene_discovery/{sample}_exomiser_comparison.tsv",
        summary=PATH + "{sample}/gene_discovery/{sample}_exomiser_comparison.summary.json"
    params:
        require="--require-exomiser" if EXOMISER_ENABLED else ""
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} compare --query {input.query:q} --exomiser-genes {input.genes:q} \
            --backend-provenance {input.provenance:q} \
            --ranking {input.ranking:q} --phenotypes {input.phenotypes:q} --edges {input.edges:q} \
            --output {output.tsv:q} --summary-output {output.summary:q} {params.require}
        """


rule exomiser_comparison:
    input:
        expand(PATH + "{sample}/gene_discovery/{sample}_exomiser_comparison.tsv", sample=SAMPLES),
        expand(PATH + "{sample}/gene_discovery/{sample}_exomiser_comparison.summary.json", sample=SAMPLES)


rule annotsv_output_reconciliation:
    input:
        vcf=globals().get("ANNOTSV_AUDIT_VCF", PATH + "{sample}/sv/merged/{sample}_merged_SV.vcf.gz"),
        annotsv=PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.tsv",
        unannotated=PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.unannotated.tsv",
        script=SCRIPTS + "/audit_sv_outputs.py",
        reconciliation=SCRIPTS + "/annotsv_reconciliation.py",
        integration=SCRIPTS + "/build_integrated_sv_gene_tsv.py",
        availability=SCRIPTS + "/audit_annotsv_evidence.py",
        common=SCRIPTS + "/sv_evidence_common.py"
    output:
        tsv=PATH + "{sample}/sv/annotsv/{sample}_annotation_records.tsv",
        json=PATH + "{sample}/sv/annotsv/{sample}_annotation_records.json"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} --vcf {input.vcf:q} --annotsv {input.annotsv:q} \
            --unannotated {input.unannotated:q} --output-tsv {output.tsv:q} --output-json {output.json:q}
        """


def annotsv_audit_resources(wc):
    return [str(path) for path in audit_resource_paths(ANNOTSV_ANNOTATIONS_DIR).values() if path.is_file()]


rule annotsv_evidence_audit:
    input:
        tsv=PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.tsv",
        resources=annotsv_audit_resources,
        script=SCRIPTS + "/audit_annotsv_evidence.py",
        common=SCRIPTS + "/sv_evidence_common.py"
    output:
        json=PATH + "{sample}/sv/annotsv/{sample}_annotation_availability.json",
        tsv=PATH + "{sample}/sv/annotsv/{sample}_annotation_availability.tsv"
    params:
        annotations=ANNOTSV_ANNOTATIONS_DIR,
        cache=PATH + "reference_checks/annotsv_bundle_inventory.json",
        resource_snapshot=lambda wc: json.dumps(audit_fingerprint(audit_resource_paths(ANNOTSV_ANNOTATIONS_DIR)), sort_keys=True)
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} --annotations-dir {params.annotations:q} \
            --genome-build GRCh38 --annotsv {input.tsv:q} --cache {params.cache:q} \
            --output-json {output.json:q} --output-tsv {output.tsv:q}
        """
