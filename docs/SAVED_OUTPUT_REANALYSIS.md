# Rechecking saved SV results

You can rebuild the evidence tables from the files you already have. The four
optional clinical TSVs can stay `null`. This does not remove the gene-level HPO
evidence or the recorded AnnotSV classes. Missing family and patient information
remains unknown in the separate allele assessment.

## Run the rebuild

From the repository root, use Python 3.11. The rebuild uses the standard library.

```bash
python scripts/rebuild_sv_results.py \
  --sample SAMPLE \
  --vcf /path/to/SAMPLE_merged_SV.vcf.gz \
  --annotsv /path/to/SAMPLE_merged_SV.annotsv.tsv \
  --unannotated /path/to/SAMPLE_merged_SV.annotsv.unannotated.tsv \
  --panel PANEL_OA/optic_neuropathy_genes.txt \
  --phenotypes /path/to/SAMPLE_human_gene_phenotypes.tsv \
  --needlr /path/to/SAMPLE_needLR_RESULTS.tsv \
  --vep /path/to/SAMPLE_SV_VEP.txt \
  --caller-vcf Sniffles2=/path/to/SAMPLE_sniffles.filtered.vcf \
  --caller-vcf cuteSV=/path/to/SAMPLE_cutesv.filtered.vcf \
  --caller-vcf Delly=/path/to/SAMPLE_delly.filtered.vcf \
  --outdir /path/to/new_reanalysis_folder
```

Replace the paths and sample name. Choose a folder that does not already exist;
the command refuses to overwrite a previous analysis. The caller order defaults
to `Sniffles2,cuteSV,Delly`, matching this repository's Jasmine workflow. If the
original merge used a different order, provide it with `--caller-order`.

Use the genome-wide AnnotSV TSV, not its panel-only view. The rejection log must
come from the unchanged VCF used in that AnnotSV run. The script checks the
recorded line number, chromosome, position and variant type, plus the allele
representation where possible. A mismatched log stops the rebuild.

The phenotype, needLR, VEP and rejection-log files, and the caller VCFs, are
optional. Omit an option when its file is unavailable. The command does not fetch
missing patient data or query Monarch again. Without the saved human phenotype
table, it cannot reproduce that part of the old gene score.

For SRS, use `--sequencing-type SRS --caller-order Manta,Delly` and the two
corresponding caller VCFs. Omit needLR.

The rebuild checks the old caller records against the current support and size
settings, but retains every master SV. Defaults are support 2 and size 50 bp;
use `--min-support` and `--min-svlen` if your study uses other values. It does not
recreate blacklist flags, COV_VAR rescue decisions or every historical filter
setting. Use the main workflow with the full study configuration for a new
callset. A filtered VCF cannot recover variants discarded by an earlier run.

## What the files contain

| File ending | What to use it for |
| --- | --- |
| `_annotation_records.tsv` | One entry per master SV, with AnnotSV's skip reason and gene-mapping status. |
| `_annotation_availability.tsv` | Which source columns and database mentions are present. |
| `_vep_gene_differences.tsv` | VEP gene annotations absent from the AnnotSV gene list, with nearby genes and unannotated records identified separately. |
| `_ranked_candidates.tsv` | Gene-level ranking, using each gene–HPO association once. |
| `_integrated_SV_gene_analysis.tsv` | The rebuilt SV–gene table, including caller and transcript evidence. |
| `_allele_assessment.tsv` | The same rows with separate evidence-domain results and predicted functional context. |
| `_SV_gene_review.tsv.gz` | A smaller review table without long sequence and JSON fields; decompress it to open the TSV in a spreadsheet. |
| `_rebuild_manifest.json` | Input hashes, script hashes, settings and commands used for the rebuild. |

The compact table is for browsing. Use the full table when inspecting
transcripts, detailed evidence or the original VCF fields.

## Changes that matter for interpretation

**Small Delly calls.** Some Delly records have no `SVLEN`. The parser now derives
their length from literal REF/ALT sequences, or from `END-POS` for interval
variants. The table records how the length was obtained. A symbolic insertion
cannot use `END-POS` as its inserted length. BND/TRA records are not filtered on
a single length, even if the caller wrote `SVLEN=0`.

**Missing AnnotSV records.** A record below AnnotSV's size threshold, a record
on an unsupported contig and an unexplained missing annotation have different
statuses. A reciprocal-breakend notice does not automatically mean the whole
record lacks annotation. The master VCF is retained.

**Incorrect gene names.** `cmpl` and `incmpl` are transcript completeness labels.
They are excluded from gene lists and rankings, while affected SVs and their
raw transcript rows remain in the integrated table. An unresolved row keeps
these annotations visible. `UNK` is retained because it is also a valid human
gene symbol. Before running AnnotSV, the main workflows check the prepared
Ensembl gene BED and stop if the incorrect labels are present.

**Inversions and breakends.** A gene entirely inside an inversion is identified
as such; overlap alone does not earn disruption points. BND transcript checks
use the appropriate chromosome and endpoint. Blacklist checks use both BND
endpoints rather than treating a remote coordinate as a local span.
`ALLELE_FUNCTIONAL_CONTEXT` describes the prediction in more detail. Expression
changes and clinical effects are not inferred from overlap alone.

**Database availability.** With no annotation directory supplied, the audit
checks the saved TSV. It distinguishes a missing source column from a source
observed elsewhere in the same file. It cannot establish everything installed
on the server, so it records `BUNDLE_NOT_INSPECTED`. Add `--annotations-dir` to
inspect the installed resources too. An empty field never becomes a benign
classification or proof of no overlap.

**VEP coverage.** Old output produced with `--pick` cannot recover the omitted
transcripts. The coverage report records that limitation. The active workflows
already use `--flag_pick`; obtaining fuller VEP output requires rerunning VEP on
the existing master VCF. Caller reruns are not needed for this change.

## If the gene resource fails its check

Run this on the server before rerunning AnnotSV:

```bash
python scripts/check_annotsv_gene_resource.py \
  --annotations-dir /path/to/AnnotSV_annotations
```

The expected Ensembl BED has ten tab-separated fields, with the gene name in
column 5. Check its preparation against the instructions for your AnnotSV
release. Keep the current files, then rebuild the resource from the matching
source annotation or install a verified replacement. Simply deleting the bad
rows would leave genes missing from the search.

The saved TSV cannot reliably tell us all the intended gene names. Once the
resource is repaired, rerun AnnotSV on the existing master VCF and rebuild the
downstream tables. Keep the resource version with the thesis methods.

Use [the gene-reference repair guide](ANNOTSV_GENE_REFERENCE_REPAIR.md) to
rebuild from the source GTF, compare transcript identities and coordinates,
and install the replacement with a backup. The repair preserves an Ensembl
gene ID when the source supplies no symbol and refreshes derived promoter
annotations. The large GTF can remain on the server.

## How this fits the optic-neuropathy study

Use the panel to examine known ON genes, then investigate nonpanel candidates
with relevant human disease/phenotype evidence and a plausible SV effect.
Review read support, rarity and exact transcript consequences for promising
alleles. A high gene score can draw attention to a gene; it does not establish
that the particular SV changes its function.

Mitochondrial pathway membership is not calculated by this rebuild. A versioned
human MitoCarta/pathway annotation would be a useful next addition for comparing
panel and nonpanel genes. That is a reference annotation task and does not
require trios or new patient questionnaires. Pathway membership would still be
biological context, not evidence that an SV caused optic neuropathy.

This study can identify plausible SV contributions and suggest new associations.
Establishing a new ON gene, a diagnosis or a functional change requires evidence
beyond these tables.

## References for the annotation checks

- [UCSC gene-model fields](https://genome.ucsc.edu/goldenPath/help/bigGenePred.html): CDS status labels and gene names are separate fields.
- [AnnotSV Ensembl gene preparation](https://github.com/lgmgeo/AnnotSV/blob/v3.5.10/share/tcl/AnnotSV/AnnotSV-genes.tcl): prepared BED layout used by the resource check.
- [Ensembl VEP options](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html): transcript selection and SV annotation options.
- [MitoCarta3.0](https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways): mitochondrial proteins and pathway annotations.
