# Run after AnnotSV, which can prepare processed BED resources lazily.
# Cache bundle scans across samples; recheck file identity on every audit.
import sys
import json
import shlex
from pathlib import Path
sys.path.insert(0, SCRIPTS)
from audit_annotsv_evidence import resource_paths as audit_resource_paths
from audit_annotsv_evidence import fingerprint as audit_fingerprint
from sv_evidence_common import hon_hpo_seed_path

# One HPO query is saved before AnnotSV and reused in the separate comparison.
# Existing configs without the new switch retain their previous behaviour.
ANNOTSV_EXOMISER_ENABLED = config.get("annotsv_exomiser_enabled", bool(config.get("hpo_terms")))
if not isinstance(ANNOTSV_EXOMISER_ENABLED, bool):
    raise ValueError("annotsv_exomiser_enabled must be true or false, without quotes.")
ANNOTSV_PATIENT_HPO_FILES = {str(sample): path for sample, path in (config.get("patient_hpo_files", {}) or {}).items()}
ANNOTSV_HPO_SEEDS = hon_hpo_seed_path(config, REPO_ROOT)


def annotsv_patient_hpo_inputs(wc):
    path = ANNOTSV_PATIENT_HPO_FILES.get(str(wc.sample)) if ANNOTSV_EXOMISER_ENABLED else None
    return [path] if path else []


def exomiser_resource_inputs(wc):
    if not ANNOTSV_EXOMISER_ENABLED:
        return []
    root = Path(ANNOTSV_ANNOTATIONS_DIR)
    # Existing resources are declared so a local update reruns the comparison.
    # Missing resources produce the clear setup error in prepare-hpo.
    paths = list((root / "jar").glob("exomiser-rest-prioritiser-*.jar"))
    # Server lock/trace files are temporary activity, not new reference data.
    paths += [path for path in (root / "Annotations_Exomiser").rglob("*")
              if path.is_file() and not path.name.endswith((".lock.db", ".trace.db", ".tmp", ".log"))]
    return sorted(str(path) for path in paths)


rule prepare_annotsv_hpo:
    input:
        seeds=ANNOTSV_HPO_SEEDS,
        patient=annotsv_patient_hpo_inputs,
        resources=exomiser_resource_inputs,
        script=SCRIPTS + "/compare_exomiser_hon.py",
        common=SCRIPTS + "/ranking_common.py"
    output:
        tsv=PATH + "{sample}/gene_discovery/{sample}_exomiser_hpo.tsv",
        terms=PATH + "{sample}/gene_discovery/{sample}_exomiser_hpo.txt"
    params:
        enabled="--enabled" if ANNOTSV_EXOMISER_ENABLED else "",
        annotations=ANNOTSV_ANNOTATIONS_DIR,
        resource_snapshot=lambda wc: json.dumps(audit_fingerprint({path: Path(path) for path in exomiser_resource_inputs(wc)}), sort_keys=True),
        shared_hpo=json.dumps(config.get("hpo_terms")),
        patient=lambda wc: "--patient-hpo " + shlex.quote(ANNOTSV_PATIENT_HPO_FILES[str(wc.sample)]) if annotsv_patient_hpo_inputs(wc) else ""
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} prepare-hpo {params.enabled} --hpo-seeds {input.seeds:q} \
            --annotations-dir {params.annotations:q} \
            --shared-hpo-json {params.shared_hpo:q} {params.patient} \
            --output {output.tsv:q} --terms-output {output.terms:q}
        """


rule annotsv_exomiser_comparison_annotation:
    input:
        vcf=globals().get("ANNOTSV_AUDIT_VCF", PATH + "{sample}/sv/merged/{sample}_merged_SV.vcf.gz"),
        # Wait for the main AnnotSV run to prepare its annotation resources.
        primary=PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.tsv",
        query=rules.prepare_annotsv_hpo.output.terms
    output:
        # The gene comparison retains the scores, evidence and source checksum.
        # This second large annotation table is temporary.
        tsv=temp(PATH + "{sample}/sv/annotsv/exomiser/{sample}_exomiser.annotsv.tsv")
    params:
        outdir=PATH + "{sample}/sv/annotsv/exomiser"
    log:
        PATH + "{sample}/logs/annotsv/{sample}.exomiser.log"
    conda:
        CONDAENV + "annotsv.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir:q} $(dirname {log:q})
        hpo_terms=$(cat {input.query:q})
        if [ -n "$hpo_terms" ]; then
            java -version 2> {log:q}
            AnnotSV -SVinputFile {input.vcf:q} -annotationsDir {ANNOTSV_ANNOTATIONS_DIR:q} \
                -outputDir {params.outdir:q} -outputFile {wildcards.sample}_exomiser.annotsv.tsv \
                -genomeBuild GRCh38 -candidateGenesFile {GENE_LIST:q} -candidateGenesFiltering 0 \
                -annotationMode split -tx ENSEMBL -REreport 1 -benignAF {ANNOTSV_BENIGN_AF} \
                -rankFiltering "1 2 3 4 5 NA" -hpo "$hpo_terms" 2>&1 | tee -a {log:q}
        else
            printf 'Annotation_mode\tGene_name\tExomiser_gene_pheno_score\n' > {output.tsv:q}
            printf 'Exomiser comparison disabled in config.\n' > {log:q}
        fi
        test -s {output.tsv:q}
        """


rule compare_exomiser_hon:
    input:
        query=rules.prepare_annotsv_hpo.output.tsv,
        annotsv=rules.annotsv_exomiser_comparison_annotation.output.tsv,
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
        require="--require-exomiser" if ANNOTSV_EXOMISER_ENABLED else ""
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        python {input.script:q} compare --query {input.query:q} --annotsv {input.annotsv:q} \
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
