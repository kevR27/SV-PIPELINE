# SRS WGS pipeline — Illumina

This workflow finds structural variants throughout the genome and organises
results for hereditary optic-neuropathy research. **Nuclear genes involved in
mitochondrial function are the main biological focus; mtDNA is a secondary
analysis.** Non-panel genes remain available for review.

The layout now follows LRS: one configuration, numbered analysis stages,
shared `scripts/` and `envs/` folders, and a separate `all_thesis_plots` target.
The SRS folder still contains only its existing Snakefile, config and README.

## Start with these files

- **`config_srs_wgs.yaml`**: your samples, paths and analysis settings.
- **`Snakefile_SRS_WGS`**: the steps Snakemake runs. Usually you only edit the config.
- **This README**: preparation, commands and interpretation.

In the config, start with **sections 1–3**: replace the example BAM/result/reference
paths, check the annotation files, and choose the optional analyses. Leave
sections 4–10 at their defaults unless you intend to change a threshold.
Numeric sample IDs must be quoted. `null` means no file or optional value was
provided; `true`/`false` enable/disable a branch.

Run from this folder so paths such as `../../PANEL_OA/` resolve correctly:

```bash
cd /path/to/SV-PIPELINE/snakemake_pipelines/srs_wgs_pipeline
```

The FASTA and all annotation intervals must match the BAM alignment reference.
This workflow uses GRCh38. `nuclear_mito_candidate_bed` uses BED coordinates
(0-based start, excluded end), including the supplied `.tsv` interval file.
`gene_bed` describes gene intervals. An optional matching GTF supplies real
exon tracks in locus diagrams; without it, the diagrams show gene intervals.

## How to read the Snakefile

Each numbered section describes a biological step. Inside a rule:

| Word | Meaning |
| --- | --- |
| `input` | Files needed before the step can run |
| `output` | Files the step creates |
| `params` | Analysis settings passed to the tool |
| `threads` | CPU threads available to that step |
| `conda` | Existing environment definition in `envs/` |
| `shell` | Command that runs the tool |
| `{sample}` | Repeat the step for each configured sample |

`rule all` lists the complete analysis results. Snakemake works backwards from
those results and runs any missing steps; the order of rules in the file is
for readability. The numerical section labels are a reading guide.

| Section | Question and result |
| --- | --- |
| 1. Reference checks and QC | Are the files compatible, and is coverage/alignment adequate? |
| 2. Manta + DELLY; GRIDSS | Which breakpoint SVs are detected, and what read/assembly evidence supports them? |
| 3. DeepVariant + WhatsHap | Which small variants occur, and which nearby small variants are phased? |
| 4. CNVpytor + combined SVs | Which deletions/duplications have copy-number evidence? Keep passing depth-only CNVs too. |
| 5. AnnotSV, VEP, Monarch, population | Which genes/transcripts overlap the SV, and what annotation/frequency information is available? |
| 6. Repeats, optional MELT, secondary mtDNA | Which other variant classes need separate review? |
| 7. Evidence, MitoCarta, candidate review | Which findings fit the biological focus, including non-panel genes? |
| 8. Plots and downstream views | What do the data and review filters show? |
| 9. Optional control checks | How do results compare with known positive controls? |

The panel is used for interpretation. Manta/DELLY and DeepVariant discovery
remain genome-wide. The readable small-variant VEP table is restricted to
`nuclear_mito_candidate_bed`; the full DeepVariant VCF is preserved.
MitoCarta annotation of SVs is genome-wide. Consequently, this workflow does
not yet give equally broad annotation of non-panel small variants and SVs.

## Run the analysis and plots

After the preparation steps below, use these commands from the SRS folder.
They reuse your existing environment location.

```bash
# Check which analysis steps would run; no analysis is executed.
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 32 --dry-run

# Run the complete analysis, including the existing Samplot review branch.
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 32 all

# Generate the overview figures and the separate downstream filtering figures.
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 32 all_thesis_plots
```

To check only coverage/alignment, replace `all` with `quality_control`. To
create only the downstream filtering tables, use `further_filtering`.
The plotting target uses completed results when available; if required results
are missing, Snakemake can run the necessary analysis steps first. It does not
limit itself to drawing existing files.

## Where to start reading results

All paths below are inside the configured `path` folder.

| File/folder inside `<sample>/` | What to read |
| --- | --- |
| `qc/<sample>.wgs_qc_summary.tsv` | Coverage/alignment thresholds; `REVIEW` means inspect the sample |
| `diagnostic_review/<sample>.diagnostic_candidates.tsv` | Main ordered SV–gene review table; all rows by default |
| `diagnostic_review/<sample>.mitochondrial_gene_ranking.tsv` | Nuclear mitochondrial gene ranking, followed by secondary mtDNA context |
| `plots/00_start_here/START_HERE.txt` | Guide to the figures and an editable review TSV |
| `plots/01_qc/` | Short-read coverage figures |
| `plots/02_caller_concordance/` | Manta/DELLY overlap; GRIDSS is shown separately as supporting evidence |
| `plots/03_sv_landscape/` | SV types, sizes, chromosomes and large-event context |
| `plots/05_candidate_prioritization/` | Separate panel/non-panel figures and SRS evidence matrices |
| `plots/06_phenotype/` | Gene–HPO associations from the configured phenotype resources |
| `plots/08_phasing/` | Local small-variant phasing; this does not establish SV phase |
| `plots/12_candidate_loci/` | Diagrams for the initial review shortlist |
| `plots/14_mitochondrial_context/` | MitoCarta/ranking figures when mitochondrial annotation is enabled |
| `diagnostic_review/samplot/` | Read-evidence images; the manifest is beside the main review table |
| `downstream_filtering/` | Separate filtering views, their audit and exact thresholds |
| `downstream_filtering/plots/` | Five separate figures explaining filter outcomes and retained SVs |

Plots are exported as PNG, PDF and SVG with accompanying data tables. Counts
state whether they describe unique SVs, genes or SV–gene associations: a large
SV can overlap several genes. A colour describes the stated category or
evidence state; it is not a pathogenicity classification. Empty SV input still
produces the QC/review summaries; SV overview plots are skipped.

The SRS evidence matrix shows CNVpytor depth, GRIDSS, repeat-locus overlap and
MELT. It distinguishes both-breakpoint GRIDSS support from a partial match and
an analysis that was not run from an evaluated no-match result. Repeat-locus
overlap means the SV overlaps a catalogued locus, not that an expansion was
confirmed. CNVpytor depth requires passing quality checks and the expected
loss/gain direction to appear as computational support.

`candidate_review_top_n` selects initial rows per panel group.
`candidate_locus_top_n` and `samplot_top_n` limit image numbers.
`diagnostic_candidate_top_n: 0` keeps every candidate row; changing it to a
positive value also restricts the input available to downstream filtering.
Existing editable review notes are preserved on reruns; generated observations
are refreshed separately. Check `<sample>/logs/plots/` if plotting fails.

## Your downstream filtering views

The existing shared filter script is reused. It keeps the complete candidate
table and writes separate review views with the config section 8 thresholds:

- Main size view **strictly >1 kb**; small (<1 kb), exactly 1 kb and unresolved
  lengths are kept in separate views. Small genic SVs can still matter.
- One combined group for **measured AF ≤1% OR missing AF**. For example,
  `gt_1kb/by_af/af_le_1_percent_or_missing.tsv`. Its accompanying
  plots use different colours to retain the distinction between measured and
  unknown AF.
- Separate exon/splice and intron annotation context, caller support and
  inversion support views. Annotation context does not prove coding disruption.
- chrY/chrM events are excluded from these nuclear review views. The complete
  upstream analysis and secondary mtDNA outputs are preserved.
- Inversions require at least 3 reported supporting reads, with 3–4 and ≥5
  separated. The largest reported count from one caller is used, not a sum
  of counts across callers. Short-read callers may count fragments/read pairs
  differently; these are review thresholds, not independent molecule counts.

Missing AF may be compatible with a candidate absent from the searched
population data, but it can also result from missing resources or unsuccessful
matching. **It does not prove that the SV is new.** A 1% threshold is a broad
review filter, not a disease-specific maximum credible frequency.
The five plots explain filter outcomes, retained large and small SV context,
combined low/missing AF by caller support, and inversion support separately.

## Interpretation limits that remain

The `TECHNICAL_SUPPORT` column describes how an SV was detected:

| Label | Meaning |
| --- | --- |
| `MULTI_CALLER_PLUS_READ_DEPTH` | Manta/DELLY agreement plus passing copy-number evidence |
| `MULTI_CALLER_PLUS_GRIDSS` | Manta/DELLY agreement plus compatible GRIDSS breakpoints |
| `MULTI_CALLER` | Agreement between Manta and DELLY |
| `SINGLE_CALLER_PLUS_SUPPORTING_EVIDENCE` | One main caller plus compatible additional evidence |
| `SINGLE_CALLER` | One main caller |
| `READ_DEPTH_ONLY_CNV` | A passing CNVpytor-only deletion/duplication |
| `REVIEW_REQUIRED` | Available evidence does not establish one of these groups |

- Manta's original VCF is preserved; its official converter exposes INV3/INV5
  inversion junctions for merging. Junctions and caller agreement do not prove
  a complete two-junction inversion.
- SURVIVOR combines Manta and DELLY. Its `multi_caller` result means agreement
  between callers, not validated high confidence. GRIDSS uses the same BAM;
  it supplies a different algorithm, not an independent experiment.
- CNVpytor uses corrected read-depth calls at 10 kb and 100 kb. Its BAF layer
  supplies context; the prototype combined RD+BAF caller is not used. The
  `*.reference_info.txt` must confirm hg38 GC/mask resources.
- VEP keeps all transcript consequences with `--flag_pick`; the readable
  summary does not delete other transcripts.
- ExpansionHunter and MELT are separate variant classes and are not additional
  generic SV callers. Mutserve covers mtDNA SNVs/heteroplasmy, not a complete
  mtDNA indel analysis.
- Standard Illumina WGS does not provide native-DNA methylation for XCI.
  SRS has no native methylation/XCI branch. Local small-variant WhatsHap
  phasing does not establish an SV+SNV pair in trans.
- Nuclear mitochondrial membership, gene ranking and technical-support labels
  are research context. They do not establish the causal variant or an ACMG class.
  Missing patient/family information is not inferred.

## Optional positive-control checks

For prioritisation, fill `validation/candidate_truth.template.tsv` and set
`candidate_truth_tsv` (GENE required, SV_ID optional). For SV detection, set
`sv_truth_vcfs` and the matching `sv_truth_beds` for each control sample.
Confident regions must use the same assembly. The existing result folder is
`<sample>/benchmark/truari/`.

Truvari compares PASS events ≥50 bp by type, coordinates and size
(`--pctseq 0`, `--pctsize 0.7`). This is not sequence-resolved allele validation.
Interpret only variant classes and regions represented in the truth data;
absence outside confident regions is not a reliable negative.

## Install tools and prepare references once (manual steps)

Installation and resource downloads belong here, not inside patient-analysis
rules. Use the same GRCh38 FASTA used to align your BAMs; matching chromosome
names alone is insufficient. Record source URLs/releases and checksums of your
reference files. The configured BAMs must be coordinate-sorted, indexed,
paired-end WGS alignments.

1. **Runner and tool environments.** Activate the existing Snakemake environment.
   If there is no runner, install it once with
   `conda create -n sv-srs-runner -c conda-forge -c bioconda snakemake=8.30.0`,
   then `conda activate sv-srs-runner`. Use the environment preparation command
   in the next section. `envs/*.yaml` are version definitions, not extra installer
   scripts. After creation, run the same Snakemake command with
   `--list-conda-envs` to find the environment paths. Activate the listed QC,
   GRIDSS or CNVpytor environment for its respective command below.

2. **FASTA/BAM indexes.** In the QC environment, run `samtools faidx` on your
   configured FASTA and `samtools index` on each sorted BAM if its index is
   missing. Do not reindex a different FASTA as a substitute for the alignment
   reference. Custom BAM index paths can be set with `bam_indexes`.

3. **GRIDSS reference preparation.** In the GRIDSS environment, run this once
   for the exact configured reference, replacing the example path:

   ```bash
   SRS_REFERENCE=/DATA/Reference/hg38.fa
   SRS_GRIDSS_JAR=$(find "$CONDA_PREFIX/share" -type f -name '*gridss*jar-with-dependencies.jar' -print -quit)
   test -n "$SRS_GRIDSS_JAR"
   gridss -r "$SRS_REFERENCE" -j "$SRS_GRIDSS_JAR" -s setupreference
   ```

   Keep the generated reference/index files beside that FASTA. The analysis
   runs only `preprocess,assemble,call`; it does not install/index a reference.
   See the [official GRIDSS setup steps](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md#setupreference).

4. **CNVpytor GC/mask data.** In the CNVpytor environment, run
   `cnvpytor -download` once if its resource files are missing. This follows the
   [CNVpytor resource instructions](https://github.com/abyzovlab/CNVpytor/blob/master/GettingStarted.md).
   The analysis requires its `-ls` report to say
   `Using reference genome: hg38 [ GC: yes, mask: yes ]`.

5. **VEP release 113 human GRCh38 cache.** The pinned VEP version and cache
   release must agree. Download/unpack the matching indexed cache manually:

   ```bash
   mkdir -p /DATA/Reference/vep
   cd /DATA/Reference/vep
   curl --fail --location --remote-name https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
   tar -xzf homo_sapiens_vep_113_GRCh38.tar.gz
   ```

   Set `vep_cache_dir` to that cache root. It must contain
   `homo_sapiens/113_GRCh38/`. If you change the release, update both
   `envs/srs_vep.yaml` and `vep_cache_version` together. See the
   [official VEP cache documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html).
   Return to this workflow folder before running Snakemake.

6. **Biological annotation resources.** Install the official AnnotSV annotation
   bundle for the pinned tool release using its
   [installation instructions](https://github.com/lgmgeo/AnnotSV/blob/master/README.md),
   and set `annotsv_annotations_dir` to the real bundle. Supply the Monarch
   nodes/edges from the same KG release at the configured paths. Download
   MitoCarta inventory/pathways by following the existing
   [MitoCarta README](../../reference/mitocarta/README.md), then set both paths
   in the configuration. These files are needed to test the nuclear
   mitochondrial hypothesis rather than merely label known panel genes.

7. **Optional resources.** MELT remains off until you install its licensed
   package and the correct human GRCh38 transposon/genes files from
   [MELT](https://melt.igs.umaryland.edu/) and fill in the paths. To use population
   evidence, obtain an indexed GRCh38 gnomAD-SV VCF from
   [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) and set
   `gnomad_sv_vcf`/`gnomad_sv_index`. The analysis does not download it. With no
   resource configured, frequency is explicitly unavailable, not proven absent.

## Validate and reuse the environments

There is no separate `environment.yaml`. It duplicated the Snakemake runner
already installed on the Bologna server and was not used by the analysis rules.
Activate the existing Snakemake environment, then let `--conda-prefix` reuse
the environments already stored in `/home/casadei7/snakemake_envs/envs/`.

After editing `config_srs_wgs.yaml`, check whether every environment required
by the selected DAG is already available. Snakemake creates only a missing
environment; it does not reinstall environments that it can reuse:

```bash
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda \
  --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --conda-create-envs-only --cores 1
```

Then validate paths, indexes, optional-resource switches, and the DAG:

```bash
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda \
  --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 1 /DATA/srs_results/reference_checks/srs_preflight.tsv

snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda \
  --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 32 --dry-run
```

Replace `/DATA/srs_results` above if your configured `path` is different.

Preflight checks actual BAM header contig lengths/names against the FASTA index,
coordinate-sort metadata, index readability and a sample of paired reads. This
is a compatibility check, not verification of every aligned base. It also checks
a nonempty human GRCh38 cache of the configured VEP release. Preflight stops
before expensive callers if a BAM/index, FASTA index, GRIDSS
BWA index, annotation bundle, HPO/Monarch resource, or enabled optional resource
is missing. Each sample also receives `qc/*.wgs_qc_summary.tsv`; threshold
failures are marked `REVIEW` instead of silently discarding calls.

## Checks made for this update — 10 October 2026

The 176 regression tests passed. Workflow dry runs covered the default analyses,
disabled optional branches, enabled optional resources/control checks, and the
QC-only target. Synthetic data produced 24 overview figures and the five
downstream filtering figures; representative evidence and AF figures were
visually checked. The seven plotting/filtering steps also completed for two
synthetic samples through Snakemake using the installed Python runtime.

No real patient BAMs or full variant callers were run here. Conda is absent in
this validation environment, so creating/reusing your server environments was
not tested. All existing caller commands and configuration defaults were
preserved; this update adds plotting/review controls and changes organisation.
