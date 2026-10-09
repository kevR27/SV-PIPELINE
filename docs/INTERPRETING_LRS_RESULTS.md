# Start interpreting completed LRS results

Start with one sample and one event. A high-ranking gene does not establish that
every SV associated with that gene is damaging. This workflow preserves the
pipeline ranks and creates an initial review workload; it does not classify
pathogenicity or change calling, analysis, or XCI outputs.

## Generate the essential review outputs

From the repository root, in the plotting environment (`envs/plots.yaml`):

```bash
for sample in 21930 25428; do
  python plots/run_thesis_plots.py \
    --root /DATA/casadei7/tools/SV-PIPELINE-main_v3/outs_test_CNTRL \
    --sample "$sample" --platform lrs --review-only \
    --gene-bed reference/gencode.genes.hg38.sorted.bed
done
```

Use the output root from your active configuration if it differs. This command
reads completed tables and creates review outputs without invoking Snakemake or
upstream analysis. It automatically prefers the final ranked SV–gene table.
Use `--candidate-table PATH` to select a different SV–gene table; do not supply
the gene-level ranking. Use `--dry-run` to inspect the planned commands.

Optionally add `--gene-gtf /path/to/matching.GRCh38.gtf.gz`. The GTF must use the
same assembly and compatible gene names. Without it, a locus figure shows gene
intervals and explicitly reports that exon detail is unavailable. With it,
one transcript is displayed: an available VEP canonical transcript, otherwise
the transcript with largest summed exon length for display only. This is not a
new consequence annotation. All selected-gene exon records and the choice are
exported for audit. Existing depth bins are added only for the exact SV ID and
chromosome; the figures do not infer a copy-number genotype.

## Open these files in order

All paths below are relative to `<output-root>/<sample>/plots/`.

| Order | Output | Purpose |
| --- | --- | --- |
| 1 | `00_start_here/START_HERE.txt` | Sample-specific navigation and source paths. |
| 2 | `01_qc/` | Coverage QC when available. Also inspect existing read-length/alignment QC before trusting calls. |
| 3 | `00_start_here/<sample>_candidate_review.tsv` | Editable review: initially up to five panel and five non-panel associations; unknown panel status remains separate. |
| 4 | `05_candidate_prioritization/<sample>_candidate_evidence_panel.*` and `..._nonpanel.*` | Separate evidence matrices with observations inside cells. Empty groups produce no figure. |
| 5 | `12_candidate_loci/top_10/` | Detailed loci from the review shortlist, with available exons/depth/methylation context and a manifest. |
| 6 | `00_start_here/evidence_availability/` | Evidence states across the full candidate table, with explicit counting units. |
| 7 | `00_start_here/mechanisms/` | Annotated geometry: loss/gain, breakpoint effects, intronic context, intact inversion-spanned genes and unresolved effects. |

Each figure has PNG, PDF and SVG versions; associated TSVs preserve exact
identifiers and values. The matrices show up to 25 associations **per panel
group** by default. `--review-top-n` and `--locus-top-n` control the smaller
initial workload. These settings do not remove candidates from your results.
Older overview figures are retained. Regenerating into an existing directory
does not delete old figures; use the current TSVs/manifests to identify the
current selection, or use `--out-dir` for a fresh plotting directory.

## Complete one review record

Identify the record by **sample + SV ID + gene**, not gene name or spreadsheet
row number. Generated columns contain observations and review flags;
the following fields deliberately require your assessment:

| Editable field | What to write |
| --- | --- |
| `WHY_GENE_MATTERS_REVIEW` | Disease association, relevant inheritance/dosage mechanism, and actual phenotype fit; record the source and date. Generic optic-neuropathy relevance is not patient-specific phenotype matching. |
| `CONTRADICTORY_EVIDENCE_REVIEW` | Evidence that weakens this particular event: incompatible mechanism, trusted population frequency, opposite depth pattern, cis pairing, or a well-supported technical concern. Separate contradictions from missing evidence. |
| `NEXT_CHECK` | One check that can change the decision: inspect breakpoint reads, verify affected coding exons/transcript, resolve phase, assess a second allele, or obtain phenotype/family information. |
| `REVIEW_DECISION` | Keep for further review, deprioritize with a stated reason, or unresolved; record reviewer and date. This is a research triage decision. |

The editable review is created once and **never overwritten on rerun**. The
shortlist and `<sample>_candidate_review.generated.tsv` are refreshed. Reconcile
changes using the sample/SV/gene key, preserving your notes.

## Read the colours critically

Reported context, computational support, review flags, evaluated no-match,
unresolved, unavailable, and not-applicable evidence are distinct. No match in
a population resource means AF remains unknown; it is not AF zero. A measured
low AF does not establish pathogenicity. The matrix uses the largest valid
explicit needLR/gnomAD event AF, excluding explicitly rejected matches.

Caller agreement and depth come from the same sequencing data. Straglr/TLDR
matches and methylation context do not by themselves validate an SV. The locus
plot is an annotation view; use existing Samplot outputs and read-level review
for breakpoint support. A partial gene overlap is not automatically coding
loss. A gene inside an inversion is not automatically disrupted.

Evidence enrichment streams selected columns from the detailed integrated
table and fills only missing values for exact SV–gene keys. It preserves final
ranks and existing observations. Ambiguous source keys are flagged without
choosing a row. Availability summaries describe evidence in the candidate
table; unavailable there can mean that an analysis exists elsewhere but was
not propagated. Event-wide domains count unique SVs; gene-specific domains
count SV–gene associations. Mechanism bars count associations, so one large
SV can contribute several genes/categories. Per-category unique-SV counts
are supplied separately and must not be summed across categories.

## Which additional plots are justified?

- **SNV–SV pairing:** inspect the existing `<sample>_snv_sv_candidates.tsv`
  for credible inheritance-compatible pairs first. A useful figure needs the
  actual SNV, SV, transcript consequence and phase evidence; label cis, trans
  or unresolved. A gene-level second-allele flag is insufficient.
- **Two-sample comparison:** use matched events, including type, coordinates,
  genotype and comparable read support. SV IDs can differ between samples;
  a shared gene is not an event match. Establish the samples' relationship
  before interpreting shared or private events. This review update does not
  create a new cross-sample matching analysis.
- **Pathway/MitoCarta:** add descriptive context after event review. Large
  intervals, gene density and shared SVs can inflate apparent enrichment.

The first unresolved questions usually concern patient phenotype, mechanism,
transcript consequence, breakpoint support, or phase. A new genome-wide bar
chart cannot answer those. Add lower-ranked direct disruptions and independent
repeat/MEI findings when justified; the initial shortlist is not exhaustive.
