# One clinical table shared by LRS and SRS. The table lives outside GitHub.
PATIENT_PROFILE_SOURCE = config["patient_profiles_tsv"]
PATIENT_PROFILE_DIR = PATH + "{sample}/clinical/"

rule prepare_patient_profile:
    input:
        profiles=PATIENT_PROFILE_SOURCE,
        hpo_seeds=hon_hpo_seed_path(config, REPO_ROOT),
        script=SCRIPTS + "/prepare_patient_profiles.py"
    output:
        packet=PATIENT_PROFILE_DIR + "{sample}.phenopacket.json",
        hpo=PATIENT_PROFILE_DIR + "{sample}.patient_hpo.tsv",
        context=PATIENT_PROFILE_DIR + "{sample}.clinical_context.json",
        allele_context=PATIENT_PROFILE_DIR + "{sample}.allele_patient_context.tsv"
    params:
        profiles=PATIENT_PROFILE_SOURCE,
        seed_flag="--include-shared-hon-seeds" if config.get("hpo_include_shared_hon_seeds", False) is True else ""
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        set -euo pipefail
        test -n {params.profiles:q}
        python {input.script:q} --profiles {input.profiles:q} \
            --sample {wildcards.sample:q} \
            --hpo-seeds {input.hpo_seeds:q} \
            --phenopacket {output.packet:q} \
            --hpo-output {output.hpo:q} \
            --clinical-output {output.context:q} --allele-context-output {output.allele_context:q} {params.seed_flag}
        """
