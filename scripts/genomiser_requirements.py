"""Mandatory Genomiser configuration for both LRS and SRS analysis.

Check configuration at Snakefile load time without requiring the server-only
clinical TSV or ReMM files to be present on the GitHub Actions runner.
Actual resource and per-sample phenotype checks happen when the job runs.
"""


def require_genomiser(config, samples):
    if config.get("genomiser_enabled", True) is not True:
        raise ValueError(
            "Genomiser is required in both pipelines; do not set "
            "genomiser_enabled: false."
        )
    if not samples:
        raise ValueError("Genomiser cannot run: configure at least one sample.")
    required = (
        "patient_profiles_tsv",
        "genomiser_cli_jar",
        "genomiser_data_dir",
        "genomiser_remm_hg38",
    )
    missing = [
        name for name in required
        if not isinstance(config.get(name), str) or not config[name].strip()
    ]
    if missing:
        raise ValueError(
            "Mandatory Genomiser settings are missing: " + ", ".join(missing)
            + ". Install the resources and configure the private patient table."
        )
    if config.get("patient_hpo_files"):
        raise ValueError(
            "Use patient_profiles_tsv for both Genomiser and HPO ranking; "
            "do not also set patient_hpo_files."
        )
    memory = int(config.get("genomiser_memory_gb", 32))
    if memory < 8:
        raise ValueError("genomiser_memory_gb must be at least 8.")
    return True
