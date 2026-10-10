# SRS WGS pipeline — Illumina

This folder contains the **short-read whole-genome sequencing pipeline** used for
Illumina data.

The pipeline is kept separate from the LRS workflow. Nothing in the LRS folder
needs to be changed to run this SRS analysis.

The files kept in this pipeline folder are:

- `Snakefile_SRS_WGS`
- `config_srs_wgs.yaml`
- this `README.md`

As in the LRS workflow, helper programs are stored in the repository-level
`scripts/` folder and Conda definitions are stored in the repository-level
`envs/` folder. Only the Snakefile, configuration and workflow documentation
remain separated here. This keeps installation simple without mixing the LRS
and SRS workflow logic.

---

## What the pipeline is trying to answer

The analysis is designed for hereditary optic neuropathy WGS and follows the
same general reasoning as the LRS pipeline:

1. Find structural variants across the genome.
2. Keep both known optic-neuropathy genes and non-panel genes.
3. Do not remove a candidate simply because only one type of evidence is
   available.
4. Add independent evidence when possible.
5. Keep mitochondrial biology separate from direct mtDNA variant calling.
6. Produce a final table that can be inspected biologically, not just a list of
   variants.

The optic-neuropathy panel is therefore used for **interpretation**, not to
restrict genome-wide discovery.

---

## Pipeline flow

| Stage | Result and biological question |
| --- | --- |
| BAM QC | Is coverage/alignment suitable for interpretation? |
| Manta + DELLY; GRIDSS support | Which genome-wide breakpoint events have supporting reads/assembly? |
| DeepVariant + WhatsHap | Which small variants occur, and which nearby small variants can be phased? |
| CNVpytor RD + BAF context | Which deletions/duplications show a copy-number change? |
| Integrated SVs → AnnotSV/VEP → population/context evidence | Which genes/transcripts may be affected, and what evidence is available? |
| Nuclear mitochondrial and optic-neuropathy interpretation | Which candidates match the study hypothesis, including non-panel genes? |
| Repeat, MEI and mtDNA branches | Which additional variant classes need separate review? |
| Candidate table + locus plots + optional benchmarks | What should be investigated next, and how does a control perform? |

---

## What each tool contributes

### Manta

Manta detects structural variants using paired-end, split-read and local
assembly evidence. It is one of the two main breakpoint-based SV callers.

Manta 1.6 reports inversion junctions as BND. The workflow preserves its original
VCF and uses Manta's supplied `convertInversion.py` to expose INV3/INV5 junctions
for type-aware merging with DELLY. These remain inversion junctions: neither
conversion nor caller agreement proves a complete two-junction inversion.

### DELLY

DELLY is the second breakpoint-based caller. It provides paired-end and
split-read evidence that can be compared with Manta.

### SURVIVOR

SURVIVOR merges the Manta and DELLY calls into one breakpoint-based SV set.

The file containing variants detected by at least two callers is named
`multi_caller`, not `high_confidence`.

This is deliberate: two callers can use related read evidence, so agreement
between Manta and DELLY is useful support but does not by itself prove that an
SV is real.

### GRIDSS

GRIDSS is kept separate from SURVIVOR.

GRIDSS performs breakpoint assembly and represents many rearrangements as
breakends. The pipeline therefore uses GRIDSS as **additional supporting
evidence from a different algorithm**. It still uses the same BAM; it is not an
independent experiment or orthogonal validation.

For events with two defined breakpoints, the final table distinguishes:

- `YES_BOTH_BREAKPOINTS`
- `PARTIAL_ONE_BREAKPOINT`
- `NO`

A one-breakpoint GRIDSS match is not treated as equivalent to support for the
complete rearrangement.

### DeepVariant

DeepVariant is run genome-wide.

This is important because the project is WGS-based. Small-variant discovery is
therefore not limited to the known optic-neuropathy panel.

The DeepVariant VCF is also used for:

- local WhatsHap phasing;
- SNP/allelic information used by CNVpytor.

The complete genome-wide VCF is preserved. A second PASS-focused VCF and VEP
table are generated for loci in `nuclear_mito_candidate_bed`, so SNV/indel
candidates in the main nuclear mitochondrial/bioenergetic hypothesis are easy
to review without pretending that the rest of the WGS callset was never made.
The repository default is the curated optic-neuropathy mitochondrial-gene BED;
replace it with a broader validated nuclear-mitochondrial BED when appropriate.
MitoCarta annotation of SVs is genome-wide; the small-variant VEP table covers
only this configured BED. This asymmetry matters when interpreting non-panel
small variants. The supplied `.tsv` intervals use 0-based, half-open BED
coordinates; the rule writes a `.bed` view so bcftools uses the correct convention.

### CNVpytor

CNVpytor adds information that the breakpoint callers do not provide directly:
**read-depth copy-number evidence**.

The workflow follows the CNVpytor documented sequence:

Read depth is imported, corrected for GC content, partitioned and called as CNVs.
DeepVariant SNPs provide a separate B-allele-frequency (BAF) context layer.

Two bin sizes are retained by default:

- 10 kb: main CNV analysis;
- 100 kb: large-CNV context.

The file `*.reference_info.txt` records the CNVpytor `-ls` output so you can
check that hg38 and the required GC/mask resources were recognized. The rule
stops if either resource is unavailable, rather than making uncorrected calls.

The pipeline does **not** use CNVpytor's prototype combined RD+BAF caller as the
main CNV callset. Standard read-depth calls are used for discovery; BAF remains
supporting information.

### Combining SV and CNV evidence

The script:

`../../scripts/combine_sv_cnv_evidence.py`

does three simple things:

1. keeps the Manta+DELLY breakpoint events;
2. adds CNVpytor read-depth evidence to overlapping DEL/DUP events;
3. keeps CNVpytor-only DEL/DUP calls if they pass the configured CNVpytor
   quality checks.

The resulting file is:

```text
<sample>/sv/integrated/<sample>_integrated_SV.vcf.gz
```

This is the structural-variant file that continues to AnnotSV and VEP.

### VEP

VEP is run with `--flag_pick`, not `--pick`.

Therefore all transcript consequences are retained. The script:

`../../scripts/summarize_vep_transcripts.py`

creates a readable per-SV/per-gene summary without deleting the other
transcripts.

### ExpansionHunter

ExpansionHunter is used for the repeat-expansion catalog.

It is a separate variant class and is **not** counted as another generic SV
caller.

### MELT

MELT is the optional mobile-element insertion branch.

It remains disabled until the local MELTv2 installation and reference files are
configured.

### Mutserve2

Mutserve2 is used for mtDNA SNVs and heteroplasmy.

This is included because hereditary optic neuropathies can be caused directly by
mitochondrial-genome variants.

The output is summarized by:

`../../scripts/summarize_mtdna_variants.py`

The current branch is deliberately described as **mtDNA SNV/heteroplasmy
analysis**, not complete mtDNA variant detection, because Mutserve2 VCF output
does not provide a complete mtDNA indel solution.

### MitoCarta

MitoCarta answers a different question from Mutserve2.

Mutserve2 evaluates variants **inside the mitochondrial genome**.

MitoCarta is used to ask whether a **nuclear gene affected by an SV encodes a
mitochondrial protein or belongs to a mitochondrial pathway**.

This is particularly useful for the optic-neuropathy hypothesis.

---

## Final technical-support labels

The final SRS table uses readable labels:

```text
MULTI_CALLER_PLUS_READ_DEPTH
MULTI_CALLER_PLUS_GRIDSS
MULTI_CALLER
SINGLE_CALLER_PLUS_SUPPORTING_EVIDENCE
SINGLE_CALLER
READ_DEPTH_ONLY_CNV
REVIEW_REQUIRED
```

These labels describe only the technical evidence supporting the call.

They are **not** pathogenicity classifications and they are **not** ACMG
classes.

---

## Files you should edit first

Open:

```text
config_srs_wgs.yaml
```

and change the real server paths for:

```yaml
samples:
path:
ref:
candidate_genes_list:
nuclear_mito_candidate_bed:
gene_bed:
exclude_bed:
expansionhunter_catalog:
vep_cache_dir:
annotsv_annotations_dir:
monarch_nodes:
monarch_edges:
```

MELT can remain disabled until its local files are ready. MitoCarta is enabled
because nuclear mitochondrial genes and bioenergetic pathways are the main
biological focus. The preflight check will stop with a clear missing-resource
message until the MitoCarta inventory and pathway files are installed.

Install the MitoCarta 3.0 inventory and pathway GMX at the configured paths.
This promotes nuclear-encoded mitochondrial genes and preserves their pathway
and subcompartment context. Mutserve remains a separate secondary mtDNA branch.

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

## Candidate review and benchmarking

The principal review output is:

```text
<sample>/diagnostic_review/<sample>.diagnostic_candidates.tsv
<sample>/diagnostic_review/<sample>.mitochondrial_gene_ranking.tsv
```

It is an ordered research/diagnostic-support table, not a pathogenicity
classification. The mitochondrial table uses the same ranking script as LRS
and places nuclear-encoded mitochondrial genes before the secondary mtDNA
section. `diagnostic_candidate_top_n: 0` keeps every ranked row, including
non-panel candidates; a positive value deliberately restricts this table.
Samplot alone limits its images to `samplot_top_n` and writes a manifest beside it. Interchromosomal BNDs are
explicitly retained for manual two-breakpoint review instead of being drawn as
false same-chromosome intervals.

To benchmark prioritization, copy `validation/candidate_truth.template.tsv`,
enter known positive-control genes/IDs, and set `candidate_truth_tsv`. To
benchmark SV detection against a sample truth VCF, add that sample under
`sv_truth_vcfs` **and the matching confident-region BED under `sv_truth_beds`**.
Truvari output is written under `<sample>/benchmark/truari/` (existing folder
name retained). The benchmark compares PASS events ≥50 bp by type, coordinates
and size (`--pctseq 0`, `--pctsize 0.7`). This is not sequence-resolved allele
validation. Interpret only classes/regions covered by the truth set; a DEL/INS
truth set cannot establish inversion/translocation performance. Absence outside
confident truth regions is not a reliable negative.

The integrated table now includes explicit gnomAD match/resource states using
the existing conservative coordinate/type matching script. No exact match can
result from representation differences. Its AF is provisional site evidence;
the allele assessment still requires verified allele-level evidence before
scoring population rarity. Missing or corrupt population indexes now fail
explicitly; they cannot silently appear as no population match. Custom `.tbi`
or `.csi` paths are passed to the reader. There is no needLR step in SRS.

For the additional large-SV review, reuse the existing
[downstream filtering and plotting guide](../../docs/DOWNSTREAM_SV_FILTERING.md),
supplying the complete `.diagnostic_candidates.tsv` as input. Missing AF and AF
≤1% enter the same retained group. Exon/intron context and inversion support
still need to pass; incomplete information goes to review. SRS has two generic
SURVIVOR callers, not three. GRIDSS, read depth and repeat/MEI tools are separate
support layers; a depth-only CNV legitimately has zero generic callers.

---

## Dry run

From this folder:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/srs_wgs_pipeline
```

first run:

```bash
snakemake \
    -s Snakefile_SRS_WGS \
    --configfile config_srs_wgs.yaml \
    --use-conda \
    --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
    --cores 32 \
    -n
```

Only after the dry run is clean should the complete workflow be launched.

---

## Main methodological references

The choices in this workflow follow the published/documented methods of the
tools rather than treating them as interchangeable callers:

- Manta — Chen et al., *Bioinformatics*, 2016.
- DELLY — Rausch et al., *Bioinformatics*, 2012.
- GRIDSS2 — Cameron et al., *Genome Biology*, 2021.
- CNVpytor — Suvakov et al., *GigaScience*, 2021.
- DeepVariant — Poplin et al., *Nature Biotechnology*, 2018.
- ExpansionHunter — Dolzhenko et al., *Genome Research*, 2017 and
  *Bioinformatics*, 2019.
- Mutserve / mtDNA-Server — Weissensteiner et al. and subsequent mtDNA-Server
  developments.

---

## Important interpretation limits

The workflow intentionally does not make conclusions that the available data
cannot support.

- Multi-caller support is not equivalent to variant truth.
- GRIDSS partial breakpoint support is kept separate from full breakpoint
  support.
- CNVpytor read depth supports copy-number change, not the exact breakpoint.
- WhatsHap here phases small variants locally. It does not phase an SV to a
  small variant. The pipeline does **not yet test an SV + SNV/indel combination
  in trans** in a recessive gene. That hypothesis requires allele-level
  integration and validated phasing/inheritance evidence.
- The genome-wide small-variant VCF is preserved, but its readable VEP candidate
  table is restricted to `nuclear_mito_candidate_bed`; genome-wide SV annotation
  does not make small-variant interpretation genome-wide.
- MitoCarta membership and phenotype similarity prioritize a hypothesis. They
  do not demonstrate altered gene function or establish a new disease gene.
- Short-read phasing is not treated as equivalent to long-read phasing.
- ExpansionHunter and MELT are separate variant classes, not extra votes for all
  SVs.
- mtDNA SNV analysis and nuclear mitochondrial-gene annotation are kept
  separate.
- Family inheritance is not inferred when trio/family data are unavailable.


## What this audit can establish

For the 10 October 2026 audit, all **169 regression tests passed**. Default and
optional-branch SRS dry runs built successfully. The XCI R plotter produced seven
separate figure families for informative and empty synthetic inputs, and the
new filtering figures were inspected.

These checks use synthetic variants, boundary cases and BAM metadata;
workflow validation checks Snakemake rule/DAG construction. These checks can
find coding/interface errors but do not validate sensitivity, specificity or
patient conclusions. Native NanoMethViz clustering against real modBAMs and the full WGS callers
were not executed in this audit. The full caller/reference stack and real WGS data must be
run on your configured server and assessed with relevant controls. Installation
instructions are manual and no new production helper scripts were added for
this audit.
