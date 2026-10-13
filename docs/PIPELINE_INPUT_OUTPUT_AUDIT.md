# Input/output audit — 13 October 2026

The workflows had execution defects; they were not all correct before this
review. The fixes below protect the existing biological analysis rather than
introducing another ranking or another set of installation scripts.

This is a repository and synthetic-data audit. No patient BAMs, Bologna
installation, actual needLR backend, downloaded Exomiser database, Atlas API
credentials or full external-tool analyses were available. A passing test is
not evidence that these server inputs exist or that a patient result is valid.

## What changed, and where

| Area | Existing files changed | Reason and resulting behaviour |
|---|---|---|
| needLR | `snakemake_pipelines/lrs/Snakefile_LRS_update` | Old native output files could be reused without establishing their input/backend identity. Each new run now uses its own directory; reuse requires matching provenance and output checksums. Failed runs remove the public result pair and cannot select old results. The native `prefix_RESULTS.tsv` / `prefix.vcf.gz` pair is normalized to the existing public names. |
| XCI temporary pileups | `snakemake_pipelines/lrs/xci_analysis.smk` | Matching filenames alone no longer authorize temporary modkit reuse. The BAM/index/reference, haplotag summary, executable, environment and settings must match a completed manifest. Both known compressed bedMethyl suffixes remain accepted. |
| Shared reuse check | **New:** `scripts/check_run_cache.py` | One small helper serves the two costly reuse branches. Query contents and generated outputs are hashed; large backend/BAM/reference files use path, size and nanosecond modification time, with additional content hashes for small files. This is provenance, not scientific validation. Replacing a large file while deliberately preserving its metadata can evade this lightweight check. |
| SRS DeepVariant | `envs/deepvariant.yaml`, SRS Snakefile/config/README | The former package did not reliably supply the invoked runner. The official 1.10.0 image now runs through Apptainer with the WGS model, explicit sample name, bindings, log and tabix indexing. A configured local image is a tracked input. |
| VEP | `envs/vep.yaml`, LRS Snakefiles/config/extensions | LRS executable and cache release are explicitly matched at 116. SRS retains its separate 113 pin/cache. Genome-wide small-variant masters still retain noncoding records. |
| Tool environments | `envs/annotsv.yaml`, `jasmine.yaml`, `needlr.yaml`, `modkit.yaml`, `read_qc.yaml`; caller parser rules | Python use is explicit in the relevant environments; Python-only rules no longer depend on the host Python. `read_qc.yaml` lacked Bioconda for pysam and now resolves. |
| TLDR | LRS calling Snakefile | Removed a copy from a file onto itself, which would fail after TLDR had produced its result. The existing expected nonempty table remains checked. |
| Shell paths | Active LRS/SRS Snakefiles and included rules | File/executable placeholders now use Snakemake shell quoting; dirname/basename substitutions are also quoted. Optional multi-argument fragments quote their individual paths. This does not cure unquoted paths inside third-party tools: needLR input/output/backend paths containing whitespace or glob characters are rejected early. |
| SRS preflight | `scripts/preflight_srs_wgs.py`, SRS Snakefile | Reference/index inputs are tracked. A separately configured BAI must be discoverable beside its BAM, including through a symlink to the same file. Validating an isolated BAI through pysam previously did not guarantee that callers could find it. |
| Samplot | `scripts/run_samplot_candidates.py`, `plots/run_samplot_candidates.py` | Old requested images are removed before rerunning. SRS reports failed requests in its manifest and exits with an error. Unsupported balanced/complex plotting remains explicit manual review. |
| Plot completion | `plots/run_thesis_plots.py`, LRS/SRS plotting rules | A tracked `plot_run_manifest.tsv` distinguishes completed commands, missing required inputs, dry runs and failures. A completion marker alone does not prove every optional figure exists. |
| Straglr | `scripts/annotate_straglr.py` | A valid header-only TSV with an empty BED is a legitimate no-loci result and writes a header-only annotation table. Malformed or inconsistent data still fail rather than masquerading as a negative result. |
| SNV–SV review | `scripts/build_snv_sv_candidates.py` | The existing focused review streams compressed VEP data and retains only its already-defined relevant gene/consequence subset in memory. Phase lookup stores only the requested alleles. The source genome-wide master is untouched. |
| AlphaGenome Atlas | `scripts/annotate_alphagenome_atlas.py` | Exact SNV matching accepts a one-base VEP range; nonfinite genotype-quality/depth values and invalid numerical thresholds cannot silently pass quality checks. Atlas evidence remains separate from clinical classification. |

The original HON research ranking, separate offline Exomiser comparison,
clinical ranking and AF review groups remain. Missing AF and measured AF
<=1% still share a downstream review group while retaining distinct status;
missing AF is not proof of novelty. No XCI estimator, escape/PAR/XIST mask,
biological threshold or original SV filtering criterion was changed.

## What was verified

- **245 tests passed** in the available Python runtime, including new focused
  failure/reuse tests. After the final BAI-discovery change, the SRS tests were
  rerun. Fixtures cover existing classification/relevance safeguards and the
  recent Exomiser/Atlas contracts; they do not establish clinical sensitivity.
- **459 rule/configuration cases** were evaluated across LRS calling, LRS
  postprocess and SRS with default, optional-enabled and optional-disabled
  configurations. Rule inputs, outputs, parameter formatting, environment-file
  presence and rendered Bash syntax were checked. Eight unrequested optional
  cases were explicitly marked disabled. This is not a full external-tool run
  or proof that every server dependency is available.
- The **rendered needLR shell** was executed using a synthetic tool that writes
  the upstream native filename pair. Matching inputs reused results, backend
  changes reran it, and a failed rerun could not publish stale output. This
  verifies the wrapper, not needLR's population algorithm.
- All **39 environment definitions** were attempted with a time-limited
  micromamba dry solve. **20 resolved** after correcting the read-QC channel;
  **19 timed out**, rather than being demonstrated incompatible. A dry solve
  does not install/run a tool or check every pip dependency. AlphaGenome 0.10.0
  was importable in the verification runtime; native plotting uses mocked SDK
  fixtures here rather than patient API calls.

Resolved environments: cutesv, deepvariant, exomiser, expansionhunter,
jasmine, longphase, manta, melt, modkit, mosdepth, mtdna_review, mutserve,
plots, qc, read_qc, samplot, srs_samplot, survivor, tldr and ucsc_gtftogenepred.

Unverified timed-out solves: alphagenome, annotsv, benchmark, clair3, cnvpytor,
delly, gnomad_sv, gridss, monarch, needlr, r_plot, repeatmasker, sniffles2,
sniffles2_needlr, srs_vep, straglr, vep, whatshap and xci.

The no-REST/offline Exomiser configuration was retained: CLI 15.1.0, Java 21,
extracted 2512 data, phenotype-only comparison. A cohort HPO query is a search
profile, not an assertion that each patient has every finding. The AnnotSV
classification, clinical joint priority and custom HON rank remain separate
visible evidence. No diagnostic probability was added.

## Server checks required before a new analysis

These are necessary checks of actual files and runtime conditions, not extra
scientific filtering steps. Installation/download commands remain in the
[main README](../README.md) and [SRS README](../snakemake_pipelines/srs_wgs_pipeline/README.md).

| Input/runtime | What to check on Bologna |
|---|---|
| BAM and reference | Coordinate-sorted WGS BAMs; readable matching indexes beside BAMs; FASTA and `.fai` from the same assembly used in alignment; compatible contig names/lengths. Filename `hg38` alone is insufficient. SRS also needs paired-read data. |
| DeepVariant | Create its updated runtime environment and test Apptainer executing the official image. Prepare a local SIF for offline use. Check server container permissions and CPU/image compatibility. |
| VEP | Human GRCh38 indexed cache matching each workflow's executable: LRS 116, SRS 113. Include the regulatory annotations required by the current commands. A downloaded cache of another release is not interchangeable. |
| Offline Exomiser | Download complete and extract, rather than pointing at archives. Check the configured CLI JAR, Java 21, hg38/phenotype 2512 directories and HPO seed/patient files. Run a comparison before relying on scores; a download in progress is not a completed analysis. |
| AnnotSV/gene resources | Matching annotation bundle, valid gene BED/GTF/transcript identifiers, correctly formatted GenCC/OMIM/ClinGen/ClinVar resources where used. Review unavailable/conflicting evidence flags. Installed resource correctness cannot be proven by repository tests. |
| Population AF | Actual indexed needLR backend and/or gnomAD VCF, compatible build/contigs and valid allele matching. needLR paths must have no whitespace/glob characters. Missing matches remain unknown. |
| Optional SV tools | GRIDSS reference indexes/JAR; CNVpytor GC/mask resources; correct GRCh38 repeat catalog/TRF BED; licensed MELT files if enabled. Keep disabled branches disabled until their inputs exist. |
| XCI and methylation | Eligible sample, original modified-base BAM tags (MM/ML), usable chrX phase/haplotags and phase-block evidence. Promoter/CpG/escape BEDs must be real matching-build resources; the core/XIST exclusion BEDs are not substitute CpG-island/promoter BEDs. XCI R execution was not tested locally. |
| Atlas | Complete official AVI SNV resource with the expected columns/build, bgzip/tabix layout and matching alleles. Native selected plots additionally need credentials, explicit supported tissue terms and selection files. Blood sampling does not by itself identify the disease-relevant target tissue. Large-SV exploratory outputs remain distinct from validated AVI SNV lookup. |
| Plot/report outputs | Inspect the plot-run log/manifest, Samplot/native manifests and figure source tables. A valid empty table is different from missing input, failed analysis or unavailable evidence. |

Use your usual Snakemake command with `--use-conda`, first with
`--conda-create-envs-only`, then `--dry-run`. SRS also provides a preflight
TSV/JSON target described in its README. Changed environment definitions can
require a new environment despite an older environment already existing.
Do not manually rename old environments to pretend they match a new definition.

For an existing analysis, inspect the dry-run before starting: provenance and
environment changes can require affected steps to rerun. An old unsigned
needLR result is deliberately not adopted. To regenerate deleted optional
figures whose manifest/marker remains, add `-R thesis_plots_sample` to your
usual command. There is no blanket requirement to rerun every caller.

## Redundancy and remaining scientific limitations

Only one reusable cache helper was added, because both costly branches needed
the same completion/provenance check. Existing biological and plotting scripts
were repaired rather than duplicated. The manifest is a run-status record,
not another candidate rank. No new installation script or duplicate ACMG score
was introduced.

These repairs do not establish an SV's molecular effect or HON causality.
Transcript/exon predictions, mechanism compatibility and database classifications
remain evidence for review. The two OPA1-positive patients are limited controls;
they were not analysed here and cannot validate all variant classes, inheritance
models, regulatory predictions or diagnostic performance. A real small-scale
server run, inspection of actual result headers/flags and review of read/depth
plots are still needed before expanding to a cohort.
