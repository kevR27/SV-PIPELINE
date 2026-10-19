# SRS WGS pipeline — Illumina

This workflow finds structural variants throughout the genome and organises
results for hereditary optic-neuropathy research. **Nuclear genes involved in
mitochondrial function are the main biological focus; mtDNA is a secondary
analysis.** Non-panel genes remain available for review.

The layout now follows LRS: one configuration, numbered analysis stages,
shared `scripts/` and `envs/` folders, and a separate `all_thesis_plots` target.
The SRS folder contains the Snakefile, configuration, README and a separate `srs_interpretation.smk` module for downstream interpretation.

The existing review rule now also writes
`diagnostic_review/<sample>.clinical_candidates.tsv`. It combines the reported
AnnotSV CNV class with HON relevance and actual variant evidence; the original
biological research order and separate Exomiser comparison remain available.
Both technologies use the same clinical assessment code. Passing CNVpytor-only
CNVs remain visible as depth-supported candidates needing breakpoint/genotype
review. See the [joint interpretation rules](../../docs/SCORE_INTERPRETATION.md#joint-clinical-candidate-review)
and [exact existing-result rebuild commands](../../README.md#apply-the-clinical-update-to-existing-results).

The plotting target now also creates a genome position map, SV-spacing view,
counts/evidence UpSet, size/support comparison, population AF candidate overlay
and one-page candidate cards. Its existing mosdepth QC step saves measured depth
bins for the map and local DEL/DUP views; older QC results need that step rerun
once. Optional centromere, gap and repeat BEDs are configured in section 8.
Local-assembly views require an existing assembly alignment. SRS small-variant
phasing cannot assign SV cis/trans phase. See the
[shared plotting README](../../plots/README.md#added-genome-and-candidate-views)
for figure locations, interpretation, manual resource downloads and options.

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

## Discovery SV union and downstream interpretation

The primary Manta/DELLY SURVIVOR merge uses `survivor_min_support: 1`.
Single-caller events are retained in the complete discovery VCF, while
`multi_caller_min_support: 2` produces a separate stricter review VCF.
Caller agreement is evidence, not biological or clinical validation. GRIDSS
is additional breakpoint evidence and is not a third SURVIVOR input.
CNVpytor-only passing copy-number changes remain separately identified when
combined with the breakpoint VCF.

The SRS-only `srs_interpretation.smk` module now runs from the standard
`rule all`. It reuses LRS research scripts on SRS inputs without modifying
the LRS pipeline, the primary integrated discovery VCF, or existing
`diagnostic_review` tables.

New outputs under `<sample>/gene_discovery/final/` include:

| Output suffix | Purpose |
| --- | --- |
| `_sv_gene_events.tsv.gz` | Mechanism-aware DEL/DUP/INV/BND/INS–gene relationships, event and caller evidence |
| `_large_sv_depth.tsv`, `_large_sv_depth_bins.tsv` | Event/flank mosdepth comparisons for DEL/DUP >=100 kb |
| `_large_sv_genes.tsv` | All intersecting genes from the configured GRCh38 gene BED |
| `_hpo_semantic_similarity.tsv` | Optional patient HPO comparison; unavailable values remain explicit |
| `_snv_sv_candidates.tsv` | Coding/splice SNV/INDEL candidates in genes that also contain an SV |
| `_gene_candidates.ranked.tsv`, `_sv_gene_candidates.ranked.tsv` | Further research-ranking views |
| `_mitochondrial_gene_ranking.ranked.tsv`, `_clinical_candidates.ranked.tsv` | Separate mitochondrial and clinical-review views |

The complete genome-wide DeepVariant VCF and VEP table are kept. The SNV–SV
table is a *focused review*, not a filter of these discovery outputs.
WhatsHap can phase nearby small variants but does not, by itself, phase an
SV against an SNV. For unphased SVs, the pair table therefore records
`UNRESOLVED` rather than inferring an autosomal-recessive trans pair.
Similarly, flanking-depth changes from the same BAM are supportive context,
not independent experimental validation.

### What does downstream analysis mean here?

It means the work **after finding variants**: adding information, checking
the evidence, choosing what to review first, and making figures. These steps
are in the workflow; that does not mean they have already run on your samples.

All paths below are inside `<sample>/` in your results folder.

| Analysis | Question it helps answer | Where to look |
| --- | --- | --- |
| Variant effect and disease mechanism | Does the change remove a gene, break part of it, or change its copy number? Does that fit the known disease mechanism? | `gene_discovery/<sample>_allele_assessment.tsv` and `gene_discovery/final/<sample>_sv_gene_events.tsv.gz` |
| Population frequency | Has a similar variant been reported in a population database? | AF and population-status columns in the candidate tables; missing AF remains different from measured AF <=1% |
| HON gene and phenotype review | Which genes fit optic neuropathy and the chosen HPO terms? | `gene_discovery/final/<sample>_gene_candidates.ranked.tsv` and the separate `gene_discovery/<sample>_exomiser_comparison.tsv` |
| Read and copy-number support | Do the reads and depth support the SV? | `gene_discovery/final/<sample>_large_sv_depth.tsv`, supporting-evidence columns and `diagnostic_review/samplot/` |
| Small variant plus SV review | Is there a possible second change in the same gene? | `gene_discovery/final/<sample>_snv_sv_candidates.tsv`; check call quality and whether phase is known |
| Mitochondrial gene review | Which nuclear genes have a known mitochondrial role? | `gene_discovery/final/<sample>_mitochondrial_gene_ranking.ranked.tsv` |
| Clinical candidate review | Which findings should be reviewed first when classification and HON relevance are considered together? | Start with `gene_discovery/final/<sample>_clinical_candidates.ranked.tsv`; read flags and missing evidence, not only rank |
| Figures and separate filter views | Where are the variants, what patterns appear, and what do the review filters retain? | `plots/00_start_here/START_HERE.txt`, `plots/` and `downstream_filtering/` |

The `diagnostic_review` tables remain useful earlier views. The `final` tables
also use the later event/depth/pair review. They are two stages of the same
analysis, not independent evidence to count twice. The original HON research
ranking, clinical order and Exomiser comparison answer different questions;
agreement between them is not required to retain a candidate.

For the small-variant/SV table, `SMALL_CALL_STATUS=CALLED_ALT_PASS` means the
small variant passed the caller's filter and that exact alternative allele is
present in the genotype. It does **not** mean the variant is pathogenic.
Failed, uncalled or unmatched records remain visible as
`SMALL_VARIANT_REQUIRES_REVIEW` and cannot supply a second-hit label. A missing
filter status is also kept for review. Even a passing pair still needs review
of frequency, effect, inheritance and, where possible, family data.

### What this workflow does not yet establish

- Genome-wide small variants are found and annotated, but this is not a
  complete disease ranking of every small variant. The small-variant/SV table
  focuses on genes with an SV and selected coding/splice effects. The full VCF
  and VEP table retain the other variants, including noncoding ones.
- GRIDSS and MELT add support to existing SV rows. Findings present only in
  their own VCFs are not automatically added to the main clinical table;
  inspect those separate files too. Nearby breakpoint support alone does not
  establish the exact rearrangement structure.
- SV genotypes transferred from callers require an exact caller-ID link and
  suitable quality/depth evidence. A coordinate-only match stays unresolved;
  inspect `CALLER_EVIDENCE_MATCH` and genotype-source notes. Do not infer
  inheritance from a merged VCF genotype or a missing family record.
- AnnotSV classification and predicted effects are evidence for review, not
  proof that a finding causes HON. The joint rank is not a disease probability.

The next useful check is the existing HG002/Truvari benchmark with matching
truth files and confident regions, followed by review of missed and extra
calls. That measures performance on that control; it is not full diagnostic
validation. A separate small-variant disease review would be needed if your
goal expands from SV investigation to all possible genetic causes of HON.

To dry-run the updated complete dependency graph, including the new
interpretation stages and their upstream requirements, use:

```bash
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 -n -p all
```

To run the complete workflow, request `all` as usual. Avoid launching
SRS WGS simultaneously with a resource-intensive LRS run on the shared
server without an agreed CPU, RAM and disk-I/O allocation.

### Enabling MELT when its installation is ready

MELT is optional and currently defaults to `melt_enabled: false`.
This prevents a missing or incompatible JAR, reference family ZIP, or Java
runtime from failing the entire SRS analysis. It is not a judgment that
MEIs are biologically unimportant.

The SRS `run_melt_single.py` wrapper calls MELTv2 `Single` independently
for the four configured families, using a generated one-line reference-list
file for each family's `-t` input, then combines the outputs for review.
Check that `melt_jar`, `melt_genes_bed`, and the `melt_transposons`
archives in the config are the locally installed, matching GRCh38 files.
The MELT program and reference bundles require their own installation;
the Snakemake environment does not supply these data automatically.

After verifying these files and the required runtime, change
`melt_enabled: true` and dry-run the target:

```bash
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 -n -p all
```

Then execute a concrete MELT output file target for a single sample
when the server has available resources. Snakemake rule names containing
`{sample}` or `{family}` cannot be used directly as CLI targets: request
a concrete output path, or run `all` once ready.
The first completed MELT VCF should be inspected before enabling MELT for
the full cohort. No MELT calls should be reported as completed merely
because the option is set to `true`.

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
| Additional SRS interpretation | Which event mechanisms, local-depth patterns and candidate SNV–SV pairs merit follow-up? |
| 8. Plots and downstream views | What do the data and review filters show? |
| 9. Optional control checks | How do results compare with known positive controls? |

The panel is used for interpretation. Manta/DELLY and DeepVariant discovery
remain genome-wide. Small-variant VEP annotation now also uses the full
DeepVariant VCF. `nuclear_mito_candidate_bed` selects a focused review VCF;
it does not restrict the genome-wide annotation. MitoCarta annotation of SVs
is also genome-wide.

### Genome-wide small variants and noncoding evidence

The master annotation is
`<sample>/snp_deepvariant/<sample>.genomewide.vep.tsv.gz`. It replaces the older
`<sample>.nuclear_mito_candidates.vep.tsv`; that restricted annotation is no
longer generated. Old copies on your server are not removed automatically.
The compressed master preserves all VEP transcript consequences, including
intronic, intergenic and available regulatory/motif overlaps. VEP's
`--flag_pick` marks a preferred consequence while retaining the others.
No panel, consequence-impact, coding-only or population-frequency filter is
applied to this master. The full DeepVariant VCF retains genotype, FILTER and
quality information: annotation does not make a flagged call reliable.

VEP runs once per sample, using the matching human GRCh38 cache and
`--regulatory`. Check the log at
`<sample>/logs/<sample>.vep_small_variants.log` for warnings or skipped records.
A missing regulatory annotation does not demonstrate absence of a regulatory
effect. Intronic location alone does not establish a splice effect, and an
intergenic or regulatory overlap does not establish the affected gene.

New samples use the normal `all` command below. To annotate existing genome-wide
DeepVariant VCFs without requesting the full analysis, run from this folder:

```bash
snakemake --snakefile Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 -p small_variant_annotation
```

If a required DeepVariant VCF is missing, Snakemake will schedule its upstream
jobs. Genome-wide annotation needs more time and disk space than the former
restricted table. Small variants are retained for complementary review; this
change does not create a genome-wide small-variant clinical ranking.

### AlphaGenome Atlas

Config section **11** enables the same optional research module as LRS.
Offline genome-wide AVI lookup reads the original DeepVariant VCF, preserving
intronic/intergenic variants and unsupported indels with visible statuses.
The full VEP master supplies exact-allele consequence context; the restricted
nuclear mitochondrial VCF is not the Atlas input. Existing clinical/HON ranks
and the offline Exomiser comparison remain unchanged.

Set the official GRCh38 AVI file, index and release after downloading; the
module is disabled until then. `atlas_annotation` runs only its annotation
branch; `all_atlas_plots` makes its separate
`plots/17_alphagenome_atlas/` folder. Both join the usual analysis/plot targets
when enabled. Selected native REF/ALT plots need the official SDK, an API key,
actual model tissue IDs and a reviewed allele-selection file. Blood/muscle
DNA source is separate from the prediction tissue. Large deletions can be
explored only with reviewed endpoints and sufficient model context; other SV
types remain explicitly unsupported, without invented SV AVI scores.

Follow the [shared download, configuration, selection-file and interpretation
instructions](../../README.md#alphagenome-atlas-genome-wide-research-annotation-and-separate-plots).
These predictions cannot alter ACMG classification or support clinical decisions.

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
| `diagnostic_review/<sample>.clinical_candidates.tsv` | Separate joint AnnotSV/HON clinical review order; every integrated SV–gene row retained |
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
The clinical table always retains the complete integrated input, independent
of this top-N setting. It uses AnnotSV's reported DEL/DUP class once and keeps
its points/criteria for audit, alongside HON relevance, variant effect,
disease mechanism, genotype, support and AF. The old biological-focus order
and Exomiser comparison remain separate. Start with `CLINICAL_REVIEW_GROUP`
and its flags/unknowns; read [the joint ranking and manual CNV review guide](../../docs/SCORE_INTERPRETATION.md#joint-clinical-candidate-review)
before interpreting a candidate as diagnostic. A VUS is not upgraded by HON
matching, and a whole-CNV class is not a gene-specific causal conclusion.
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
`<sample>/benchmark/truvari/`.

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
   # The Bioconda gridss launcher locates its bundled gridss.jar automatically.
   gridss -r "$SRS_REFERENCE" -s setupreference
   ```

   Keep the generated reference/index files beside that FASTA. The analysis
   runs only `preprocess,assemble,call`; it does not install/index a reference.
   See the [official GRIDSS setup steps](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md#setupreference).

4. **CNVpytor GC/mask data.** In the CNVpytor environment, run
   `cnvpytor -download` once if its resource files are missing. This follows the
   [CNVpytor resource instructions](https://github.com/abyzovlab/CNVpytor/blob/master/GettingStarted.md).
   The analysis requires its `-ls` report to say
   `Using reference genome: hg38 [ GC: yes, mask: yes ]`.

5. **VEP release 116 human GRCh38 cache (shared with LRS).** The pinned VEP version and cache
   release must agree. Download/unpack the matching indexed cache manually:

   ```bash
   mkdir -p /DATA/Reference/vep
   cd /DATA/Reference/vep
   curl --fail --location --remote-name https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
   tar -xzf homo_sapiens_vep_116_GRCh38.tar.gz
   ```

   Set `vep_cache_dir` to that cache root. The human cache should include
   regulatory annotations for `--regulatory`. It must contain
   `homo_sapiens/116_GRCh38/`. If you change the release, update both
   `envs/srs_vep.yaml` and `vep_cache_version` together. See the
   [official VEP cache documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html).
   Return to this workflow folder before running Snakemake.

6. **Biological annotation resources.** Install the official AnnotSV annotation
   bundle for the pinned tool release using its
   [installation instructions](https://github.com/lgmgeo/AnnotSV/blob/master/README.md),
   and set `annotsv_annotations_dir` to the real bundle. Supply the Monarch
   nodes/edges from the same KG release at the configured paths; follow the
   [manual Monarch setup](../../README.md#manual-monarch-kg-setup). Download
   MitoCarta inventory/pathways by following the existing
   [MitoCarta README](../../reference/mitocarta/README.md), then set both paths
   in the configuration. These files are needed to test the nuclear
   mitochondrial hypothesis rather than merely label known panel genes.

   Exomiser is enabled as a separate OFFLINE CLI gene–phenotype comparison.
   Follow the [manual setup and interpretation guide](../../README.md#offline-exomiser-cli)
   for CLI 15.1.0, extracted hg38/phenotype data 2512 and Java 21. The shared ocular
   HON profile is used by default; its file is explicit under
   `hon_hpo_seed_terms` in the config. Patient HPO files are optional. Run the
   `exomiser_comparison` target to create the report before plotting. Its
   phenotype-only scores are not fed back into the main ranking or AnnotSV
   CNV classification. Set `exomiser_enabled: false` to write a disabled report.
   Identical HPO queries share one CLI result. No REST service is started and
   no scientific database is downloaded automatically.

   `scripts/add_supporting_evidence.py` is the active GRIDSS/MELT/ExpansionHunter
   integration. The unused older `intersect_srs_orthogonal_evidence.py` has been
   removed to avoid maintaining two versions of the same matching task.

7. **Optional resources.** MELT remains off until you install its licensed
   package and the correct human GRCh38 transposon/genes files from
   [MELT](https://melt.igs.umaryland.edu/) and fill in the paths. To use population
   evidence, obtain an indexed GRCh38 gnomAD-SV VCF from
   [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) and set
   `gnomad_sv_vcf`/`gnomad_sv_index`. The analysis does not download it. With no
   resource configured, frequency is explicitly unavailable, not proven absent.

## Bologna setup: GRIDSS, CNVpytor and VEP 116

Updated 15 October 2026. Both workflows now use VEP/cache release 116. Keep
your sample paths and results directory when updating the checkout. In your
server copy of `config_srs_wgs.yaml`, set `vep_cache_version: 116`; the matching
`envs/srs_vep.yaml` must specify `ensembl-vep=116`. Do not change only the cache
number while leaving a release-113 executable installed.

Start in the existing Snakemake environment and SRS workflow folder. Set
`exomiser_enabled: false` while its data are incomplete. Leave MELT/Atlas
disabled until their resources are ready. These switches do not install the
other tools or their data.

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/srs_wgs_pipeline

snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --conda-create-envs-only --cores 1
```

This creates missing environments for the requested analysis without running
variant callers. It includes the pinned GRIDSS 2.13.2 and CNVpytor 1.3.2
environments, and VEP 116. Do not create separate named copies if you want
Snakemake to reuse the same installations. Record the actual hashed paths:

```bash
SRS_ENV_LIST=$(mktemp /tmp/srs_envs.XXXXXX.tsv)
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --list-conda-envs > "$SRS_ENV_LIST"
cat "$SRS_ENV_LIST"
```

Keep this shell open: later commands use that temporary table. If environment
creation fails, resolve that failure before proceeding; a listed path alone
does not prove its installation completed. The table's last column is the
environment path. The commands below select it by the YAML filename rather
than guessing its hash. Stop if a selection is empty or does not name a directory.

### Manta on the server's noexec /DATA filesystem

Manta 1.6.0 uses Python 2.7. The Manta Conda environment explicitly
pins `python=2.7`; without it, `configManta.py` may stop with
`/usr/bin/env: python2: No such file or directory`. The pipeline
checks for `python2` before launching the caller, then runs
`configManta.py` through that interpreter.

Manta creates `runWorkflow.py` in its results directory on `/DATA`.
Since `/DATA` is `noexec`, the generated script is also invoked as
`python2 path/to/runWorkflow.py`, not executed directly. This avoids
a filesystem-related `Permission denied` error.

The Manta environment uses Python 2.7 and samtools 1.9. For compatibility
with this older Python, it does not force modern samtools/htslib 1.20.
The package solver still needs to be tested on the Bologna server.

After merging this Manta environment change into the server checkout,
Snakemake should create a new hashed Conda environment when the Manta rule
is requested with `--use-conda`. Check that
`python2 --version` works inside that new environment. Do not add
Python 2 to the Snakemake runner or other caller environments.

If Manta subsequently fails inside its generated workflow, inspect
`rundir/workspace/pyflow.data/logs/`; the interpreter fix alone does
not establish a successful HG002 analysis.

### GRIDSS native libraries on the Bologna server

The `/DATA` filesystem is mounted with `noexec`. GRIDSS 2.13.2 can
process reads and assemblies there, but a Snappy `.so` library extracted
there fails with `UnsatisfiedLinkError: failed to map segment`.
The SRS `gridss_sv_calling` rule now sets Java's
`org.xerial.snappy.tempdir` and `java.io.tmpdir` to
`$HOME/.cache/gridss-native` for this command only. The much larger
GRIDSS intermediate files remain under `sv/gridss/work/` on `/DATA`.
The user's home filesystem must allow executing native shared libraries.

When restarting a failed GRIDSS run, keep the existing work directory:
GRIDSS can reuse completed preprocessing and assembly outputs.
Confirm the new full log does not still extract
`snappy-*-libsnappyjava.so` under `/DATA` before starting another
whole-genome run.

### GRIDSS: activate, check, prepare the same FASTA

```bash
SRS_GRIDSS_ENV=$(awk '$1 ~ /(^|\/)gridss.yaml$/ {print $NF}' "$SRS_ENV_LIST")
test -n "$SRS_GRIDSS_ENV" && test -d "$SRS_GRIDSS_ENV" && conda activate "$SRS_GRIDSS_ENV"
command -v gridss
java -version
samtools --version

SRS_REFERENCE=/DATA/Reference/hg38.fa
# The Bioconda wrapper selects its own bundled gridss.jar.
gridss -r "$SRS_REFERENCE" -s setupreference -t 32
```

Use the exact FASTA configured for your aligned BAMs, with write access for
its indexes; do not substitute another file merely named hg38. GRIDSS's
[official quick start](https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md#setupreference)
describes this one-time step. This prepares a reference, not a patient call.
The pipeline discovers the JAR from its own active Conda environment; no
additional JAR path is needed in the configuration.

### CNVpytor: activate and prepare its data

```bash
SRS_CNVPYTOR_ENV=$(awk '$1 ~ /(^|\/)cnvpytor.yaml$/ {print $NF}' "$SRS_ENV_LIST")
test -n "$SRS_CNVPYTOR_ENV" && test -d "$SRS_CNVPYTOR_ENV" && conda activate "$SRS_CNVPYTOR_ENV"
command -v cnvpytor
python -c 'from importlib.metadata import version; print("CNVpytor", version("CNVpytor"))'
cnvpytor -download
```

The official [CNVpytor instructions](https://github.com/abyzovlab/CNVpytor/blob/master/GettingStarted.md)
separate installation from reference-data preparation. Check the hg38 files
in this same environment; testing HDF5 readability catches empty downloads or
HTML error pages saved under a `.pytor` filename:

```bash
python - <<'PY'
from pathlib import Path
from cnvpytor.genome import Genome
import h5py
for key in ("gc_file", "mask_file"):
    path = Path(Genome.reference_genomes["hg38"][key])
    assert path.is_file() and path.stat().st_size > 0, f"Missing hg38 {key}: {path}"
    with h5py.File(path, "r") as handle:
        assert len(handle), f"Empty hg38 {key}: {path}"
    print("Readable hg38 resource:", path)
PY
```

If `-download` errors or the hg38 files fail this check, retrieve just these
two files from the official pinned release, then repeat the check. The
following manual fallback replaces only hg38 GC/mask files in the activated
CNVpytor environment, using temporary downloads before replacing final files:

```bash
SRS_CNVPYTOR_DATA=$(python -c 'from cnvpytor.genome import Genome; print(Genome.pytor_datapath)')
test -d "$SRS_CNVPYTOR_DATA"
curl --fail --location --retry 3 \
  https://raw.githubusercontent.com/abyzovlab/CNVpytor/v1.3.2/cnvpytor/data/gc_hg38.pytor \
  --output "$SRS_CNVPYTOR_DATA/gc_hg38.pytor.download"
curl --fail --location --retry 3 \
  https://raw.githubusercontent.com/abyzovlab/CNVpytor/v1.3.2/cnvpytor/data/mask_hg38.pytor \
  --output "$SRS_CNVPYTOR_DATA/mask_hg38.pytor.download"
python - "$SRS_CNVPYTOR_DATA" <<'PY'
from pathlib import Path
import h5py
import sys
root = Path(sys.argv[1])
for name in ("gc_hg38.pytor", "mask_hg38.pytor"):
    with h5py.File(root / (name + ".download"), "r") as handle:
        assert len(handle), f"Empty resource: {name}"
for name in ("gc_hg38.pytor", "mask_hg38.pytor"):
    (root / (name + ".download")).replace(root / name)
print("Both hg38 resource downloads are readable HDF5 files and installed.")
PY
```

After the analysis imports read depth, the existing rule additionally checks
`cnvpytor -root <sample.pytor> -ls` for
`Using reference genome: hg38 [ GC: yes, mask: yes ]`. Resource files existing
alone do not prove the sample's reference was recognized correctly.

### VEP 116: reuse the LRS cache

```bash
SRS_VEP_ENV=$(awk '$1 ~ /(^|\/)srs_vep.yaml$/ {print $NF}' "$SRS_ENV_LIST")
test -n "$SRS_VEP_ENV" && test -d "$SRS_VEP_ENV" && conda activate "$SRS_VEP_ENV"
vep --help
ls -ld /DATA/casadei7/tools/SV-PIPELINE-main_v3/reference/vep/homo_sapiens/116_GRCh38
```

The help output should identify Ensembl VEP 116. If the LRS cache is already
complete at that location, reuse it. The standard `homo_sapiens` cache is
required by the current commands; RefSeq/merged caches require different VEP
options and are not interchangeable. Do not download a second copy merely
because SRS has a separate workflow folder.

Only if the cache is absent, download/extract the standard indexed cache into
your configured cache root, following the release-116 command in the manual
setup section above. A directory listing alone does not establish completeness;
check extraction finished and inspect the first VEP log for unavailable data.

Return to the Snakemake runner environment and run `srs_preflight`, then a dry
run. DeepVariant's Apptainer image still needs its separate preparation. These
instructions and repository checks do not claim any installation or reference
preparation has already run on Bologna.

## DeepVariant runner setup

DeepVariant now runs the official 1.10.0 image, containing the trained WGS
model and `/opt/deepvariant/bin/run_deepvariant`. The former Bioconda package
wrappers did not reliably provide this command. `envs/deepvariant.yaml`
provides Apptainer and tabix; it does not install DeepVariant models itself.
See the [official runner](https://github.com/google/deepvariant/tree/r1.10)
and [Apptainer documentation](https://apptainer.org/docs/user/latest/).

On Bologna, first create the selected Conda environments using the command
below, activate the environment containing Apptainer, then check the image:

```bash
apptainer exec docker://google/deepvariant:1.10.0 /opt/deepvariant/bin/run_deepvariant --help
```

For repeated offline runs, prepare a local image while network access is
available:

```bash
mkdir -p /DATA/Reference/deepvariant
apptainer pull /DATA/Reference/deepvariant/deepvariant_1.10.0.sif docker://google/deepvariant:1.10.0
```

Set `deepvariant_image` to the absolute `.sif` path. Local images are tracked
as rule inputs. A remote image URI requires image access or a working local
Apptainer cache. The rule binds the BAM/index, FASTA and output directories,
passes the sample name and writes a DeepVariant log under
`<sample>/logs/deepvariant/`. Server container permissions, CPU compatibility
and a real WGS run still need verification; they were not tested here.

### Optional NVIDIA GPU for DeepVariant

DeepVariant 1.10.0 has two official images: the CPU image
`docker://google/deepvariant:1.10.0` and the GPU image
`docker://google/deepvariant:1.10.0-gpu`. Only the `call_variants`
step uses the GPU (one GPU at a time); `make_examples` and
`postprocess_variants` continue on CPU.

Keep `deepvariant_use_gpu: false` until a server check passes. First run
`nvidia-smi -L` on the host. A GPU is not available to this analysis
unless the NVIDIA driver is working and you have permission to use it.
Then test the GPU image with Apptainer, preferably in a terminal separate
from any current WGS analysis:

```bash
nvidia-smi -L
apptainer exec --nv docker://google/deepvariant:1.10.0-gpu \
  python3 -c 'import tensorflow as tf; print(tf.config.list_physical_devices("GPU"))'
```

This may download the large GPU image into Apptainer's home cache.
An empty GPU list or a CUDA driver/library error means leave GPU mode
disabled and consult the server administrator. `--nv` binds the host's
NVIDIA libraries into the container; the GPU image provides the correct
DeepVariant GPU software.

If the check succeeds, set `deepvariant_use_gpu: true` in your *server*
`config_srs_wgs.yaml`. Leave `deepvariant_image` as the CPU value;
the rule automatically selects `deepvariant_gpu_image` and adds
`apptainer exec --nv`. The rule reserves one `gpus` resource. Limit
simultaneous GPU jobs with `--resources gpus=1` when starting Snakemake:

```bash
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 10 --resources gpus=1 --rerun-incomplete -p
```

A completed CPU DeepVariant VCF remains valid and is not replaced
automatically just by changing a YAML setting. To compare CPU and GPU,
use a separate output directory or deliberately rerun only after
saving the existing VCF. Do not change the configuration during an
active run, or force a rerun just to use the GPU. GPU availability
and the specific NVIDIA driver/CUDA compatibility must be verified
on Bologna; no server GPU test was performed with this code change.

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

Every BAM needs a discoverable index beside it (`sample.bam.bai`, or
`sample.bai` where supported). If using the optional `bam_indexes` mapping,
the index beside the BAM must resolve to that same file; a symlink is enough.
Preflight now rejects an isolated custom BAI that callers would not discover.

Preflight checks actual BAM header contig lengths/names against the FASTA index,
coordinate-sort metadata, index readability and a sample of paired reads. This
is a compatibility check, not verification of every aligned base. It also checks
a nonempty human GRCh38 cache of the configured VEP release. Preflight stops
before expensive callers if a BAM/index, FASTA index, GRIDSS
BWA index, annotation bundle, HPO/Monarch resource, or enabled optional resource
is missing. Each sample also receives `qc/*.wgs_qc_summary.tsv`; threshold
failures are marked `REVIEW` instead of silently discarding calls.

## Latest SRS check — 16 October 2026

The input/output review found and fixed real problems:

- Spaces in `Manta, Delly` no longer stop the workflow. Reversing the caller
  order is still rejected because it would give the wrong caller labels.
- Leaving the optional MitoCarta pathway file empty no longer breaks rule loading.
- Failed or uncalled small variants cannot become possible second hits.
  At sites with more than one alternative allele, the VEP allele number is
  used to check the right genotype allele. Unmatched alleles stay unresolved.
- Failed, reference-genotype or missing-genotype MELT calls cannot supply
  support to a candidate. The original MELT output is still kept for review.
- Plots wait for their final gene and depth tables. Several file-path quoting
  errors were fixed, and the small Python panel-filter step now uses an
  existing environment that explicitly includes Python.
- Preflight keeps an error log and shows examples of reference mismatches.
  The benchmark folder spelling is now `truvari`; old `truari` folders are
  not moved or deleted automatically.

Verification: **142 selected tests passed**. We checked **216 rule/configuration
cases** (including explicitly disabled cases), and complete `all` plus
`all_thesis_plots` dry runs passed in three configurations. These used local
test files. No patient analysis, full caller run or Bologna installation was
tested here. See [the detailed audit](../../docs/PIPELINE_INPUT_OUTPUT_AUDIT.md).

No new production script or environment file was added. Your sample paths,
VEP 116 setting, original ranking and 32 GB Exomiser setting were preserved.
The small-variant/SV quality fix is in a shared script and therefore also
applies to LRS when that downstream step is rerun; LRS calling and XCI were
not changed. Keep your server config and sample paths when updating code.
Then run the usual SRS dry run shown above and inspect the jobs before running.

## Earlier audit — 13 October 2026

See [the shared input/output audit](../../docs/PIPELINE_INPUT_OUTPUT_AUDIT.md).
245 tests passed, with follow-up SRS checks after the index-discovery fix.
459 rule/configuration cases across LRS calling, LRS postprocess and SRS
rendered successfully. Twenty environment definitions resolved in dry solves;
19 timed out and remain unverified. These checks do not certify a real server
installation or patient run. DeepVariant's runner changed as described above;
SRS used release 113 for that audit; it moved to 116 on 15 October 2026.

## Earlier checks — 10 October 2026

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
