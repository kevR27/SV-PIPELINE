#!/usr/bin/env bash
# Rebuild all available BAM-free downstream thesis figures on EVERY invocation.
# Invoke with "bash run_all_plots.sh": the /DATA filesystem may be noexec.
# Reuses clinical/variant results unless missing, outdated, or incomplete.
set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

CONDA_DIR="${SNAKEMAKE_CONDA_PREFIX:-/home/casadei7/snakemake_envs/envs}"
CORES="${SNAKEMAKE_CORES:-32}"
MEM_MB="${SNAKEMAKE_MEM_MB:-45000}"

# --forcerun applies ONLY to plot rules (Python, R, filtering).
# Do not force Genomiser, SV calling, QC metrics, or clinical rankings.
snakemake -s Snakefile \
    --configfile config.yaml \
    --use-conda \
    --conda-prefix "$CONDA_DIR" \
    --cores "$CORES" \
    --resources "mem_mb=$MEM_MB" \
    --rerun-incomplete \
    -p \
    all all_saved_thesis_plots \
    --forcerun \
      archived_thesis_plots_sample \
      downstream_filtering_plots \
      r_population_gene_effect \
      r_candidate_evidence_heatmap \
      r_hpo_heatmap \
      r_large_sv_depth \
      r_sequencing_qc \
      cohort_thesis_comparison
