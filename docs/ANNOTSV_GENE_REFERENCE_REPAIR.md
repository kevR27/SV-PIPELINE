# Repairing the AnnotSV Ensembl gene reference

`cmpl` and `incmpl` describe the completeness of a coding sequence. They do not
identify genes. If they appear in column 5 of `genes.ENSEMBL.sorted.bed`, the
gene reference needs rebuilding before its gene assignments can be trusted.

A likely cause is an empty `gene_name` in the source annotation. In a
genePredExt file, the name and CDS-status fields are next to each other. The
default whitespace splitting in `awk` skips empty fields; a tab-aware parser
preserves them. The UCSC converter can legitimately leave the name field empty
when `-geneNameAsName2` is used and the GTF supplies no gene name. This mechanism
reproduces the observed error, but the original GTF is needed to confirm it
for a particular installed reference.

The rebuild below uses UCSC `gtfToGenePred` for transcript coordinates and exon
structure. It reads that tool's `-infoOut` table with tabs preserved, uses the
source gene symbol when present, and otherwise keeps the Ensembl gene ID.
Unnamed coding and noncoding genes remain separate genes. `UNK` is retained as
a valid gene symbol. A missing symbol does not imply a novel disease gene.

## What you need

- The original Ensembl **GRCh38 GTF**, either plain text or gzip compressed.
- The existing `genes.ENSEMBL.sorted.bed`.
- Python 3.10 or later and UCSC `gtfToGenePred` with `-includeVersion` and
  `-infoOut` support.

Keep the GTF on the server; there is no need to upload it. To locate an existing
copy, run this from the repository root:

```bash
find reference -type f -name 'Homo_sapiens.GRCh38*.gtf*' -print
command -v gtfToGenePred
```

If the GTF was deleted after installation, identify the annotation release from
the download command, installation log or original filename, then download
that release from Ensembl. The BED alone does not reliably identify its release.
Do not select a release from the VEP version: the two resources can differ.
The [Ensembl release 113 directory](https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/)
is an example of a versioned download location, not a claim that an existing
BED came from that release. Use the chromosome GTF, not the ab initio prediction
file. If the original release cannot be recovered, document a deliberate
reference update and review the transcript changes it produces.

## Build a replacement in a new folder

Set `ANN`, `GTF` and `RELEASE` to the paths and release on your server. Replace
the example GTF path and release before running this command.

```bash
ANN="/path/to/AnnotSV_annotations"
GTF="/path/to/Homo_sapiens.GRCh38.RELEASE.chr.gtf.gz"
RELEASE="REPLACE_WITH_SOURCE_RELEASE"

python3 scripts/rebuild_annotsv_gene_resource.py \
  --gtf "$GTF" \
  --source-release "$RELEASE" \
  --existing-bed "$ANN/Annotations_Human/Genes/GRCh38/genes.ENSEMBL.sorted.bed" \
  --outdir annotsv_gene_rebuild
```

If the converter is outside your PATH, append
`--gtf-to-genepred /full/path/to/gtfToGenePred`. The script prints progress as it
runs and refuses to overwrite an existing output folder. Conversion uses a
temporary uncompressed GTF, so allow disk space for it. The installed reference
is unchanged by this build command.

The script checks the GRCh38 header, transcript and gene IDs, duplicate records,
exon and CDS bounds, transcript versions, and consistency between the converter's
two outputs. It records noncanonical contigs excluded from the AnnotSV gene BED;
chromosomes 1–22, X, Y and mitochondrial M are retained.

It also compares every existing transcript with the new reference. A missing
old transcript, changed coordinates or a changed previously valid gene name
stops the build before an installable BED is written. Added transcripts are
reported and retained: an earlier name-based filter may have discarded them.
Review these additions as well. `--allow-reference-update` permits the broader
changes only when you deliberately choose a different annotation release; it
does not bypass coordinate or identity validation.

Check these outputs:

- `reference_build.json`: build status, input hashes, source release and counts.
- `reference_changes.tsv`: transcript additions, removals, name and coordinate changes.
- `gene_identities.tsv`: gene IDs, source symbols and the label used by AnnotSV.
- `Annotations_Human/Genes/GRCh38/`: the replacement BED and transcript-version table.

The status must be `READY` before installation. A failed build leaves a report
with `FAILED`; use its error message and converter log to resolve the cause.

## Install after a successful build

Run this when no AnnotSV job is using the annotation folder:

```bash
python3 scripts/install_annotsv_gene_resource.py \
  --build-dir annotsv_gene_rebuild \
  --annotations-dir "$ANN"

python3 scripts/check_annotsv_gene_resource.py --annotations-dir "$ANN"
```

The installer verifies the generated file hashes and checks that the installed
BED still matches the original used for comparison. It backs up the old gene
files, transcript versions, temporary Ensembl BEDs, installation marker and
Ensembl promoter caches in a new `gene_reference_backup_*` folder. It installs
the BED and version table together and restores changed files if installation
raises an error.

Old Ensembl promoter caches are removed only after being backed up, because
AnnotSV derives their gene labels from the gene BED. AnnotSV regenerates them
on its next run. RefSeq and other database resources are retained. The backup
contains an installation receipt. It is much smaller than a copy of the whole
annotation bundle.

## Rerun the analysis

Use the same master VCF, genome build and study settings for the new AnnotSV
run. Refresh gene extraction, human gene–phenotype associations, panel/nonpanel
assignment, ranking, integrated tables and their downstream plots. The old
phenotype table may lack genes recovered by this repair. Refreshing those
reference associations does not require patient HPO terms or family data.

Rerunning the SV callers is unnecessary for this reference repair. Other changes
to caller filtering are a separate reason to regenerate a callset.

Do not edit old AnnotSV gene names in place: AnnotSV may have grouped multiple
transcripts under the same incorrect label. A new annotation run is needed to
recover all affected SV–gene rows. Keep the previous results for comparison.

## Sources

- [UCSC genePredExt fields](https://genome.ucsc.edu/FAQ/FAQformat.html#format9).
- [UCSC gtfToGenePred options and transcript information](https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/gtfToGenePred/gtfToGenePred.c).
- [AnnotSV 3.5.10 gene-reference format](https://github.com/lgmgeo/AnnotSV/blob/v3.5.10/share/tcl/AnnotSV/AnnotSV-genes.tcl).
- [AnnotSV 3.5.10 promoter-cache preparation](https://github.com/lgmgeo/AnnotSV/blob/v3.5.10/share/tcl/AnnotSV/AnnotSV-regulatoryelements.tcl).
