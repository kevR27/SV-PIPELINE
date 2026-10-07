# Run after AnnotSV, which can prepare processed BED resources lazily.
# Cache bundle scans across samples; recheck file identity on every audit.
import sys
import json
sys.path.insert(0, SCRIPTS)
from audit_annotsv_evidence import resource_paths as audit_resource_paths
from audit_annotsv_evidence import fingerprint as audit_fingerprint


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
