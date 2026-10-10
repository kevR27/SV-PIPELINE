# Further downstream filtering of ranked SV–gene results

`scripts/filter_sv_candidates_for_review.py` makes a smaller, auditable research
review view from a **completed final ranked SV–gene table**. It leaves the input,
calling, original ranks and other pipeline analyses unchanged. It needs Python
and pandas, already available in the plotting environment. The numbered code
sections follow the biological questions: chromosome, event size, affected
features, population frequency, caller/read support, and organised exports.

## Exact default rules

| Question | Rule | Interpretation |
| --- | --- | --- |
| Main size view | **Strictly >1,000 bp**, using absolute reported SV length or a valid interval span | This is a review scope, not a pathogenicity criterion. No upper-size exclusion is added. |
| Smaller events | **<1,000 bp** saved separately, passing the same non-size filters | Keep small exon/splice events available for review. |
| Exactly 1,000 bp | Saved in `eq_1kb/` | Neither silently dropped nor labelled >1 kb. |
| Exon/intron context | Retain explicitly annotated exon/splice, intron, or mixed exon–intron context | Missing annotation goes to review; gene overlap alone does not establish exon/intron status. |
| Population AF | Largest valid explicit needLR/gnomAD-SV AF **≤0.01**, or AF unknown | Both enter one combined retained group; original AF status remains visible. Generic `AF` or sample allele balance is not used. |
| Chromosomes | Exclude Y and M/MT, including a known remote breakpoint partner | Nuclear mitochondrial genes remain. This does not remove mtDNA from the original analyses. |
| Caller priority | Additional view orders 3, 2, 1, 0 callers, then unavailable count | `retained.tsv` keeps original input order. All original rank values remain unchanged in every export. No new biological ranking is computed. |
| INV support | Reported maximum per-caller support **≥3**; distinguish 3–4 from **≥5** | A practical configurable review rule, not a validated universal inversion threshold. Missing support goes to review. |
| BND/TRA or unknown length | `length_unresolved/`, for manual review | A translocation has no meaningful single interval length. Insertion length is never inferred from END minus START. |

For read support, the script prefers `EVENT_MAX_CALLER_READ_SUPPORT`, then
`CALLER_READ_SUPPORT`/`READ_SUPPORT`, then structured caller records. It uses
the **largest reported count from one caller**, never the sum across callers.
It does not recount unique reads, assess mapping quality, or establish support
at both inversion breakpoints. Zero generic callers is a valid count for an SRS read-depth-only CNV; it must not
be labelled a one-caller event. Caller-count conflicts, ambiguous inversion
evidence joins and inconsistent event-wide fields across genes are flagged.
Original technical flags remain available; this script is not a replacement
for read-level QC.

An inversion spanning an intact gene remains interval context. Exon/intron
labels describe available transcript annotations and do not locate its exact
breakpoints. Mixed context can arise from different transcript isoforms.
Splice annotations are grouped explicitly with exonic annotations; they are
not relabelled as confirmed exon breakpoints. Whole-transcript effects have a
separate scope note. Records with only proximal/regulatory annotation fail
this requested genic view, but remain in the audit/excluded export.

## Run on the two completed samples

From the repository root, using the existing plotting environment:

```bash
for sample in 21930 25428; do
  python scripts/filter_sv_candidates_for_review.py \
    --input "outs_test_CNTRL/$sample/gene_discovery/final/${sample}_sv_gene_candidates.ranked.tsv" \
    --sample "$sample" \
    --out-dir "outs_test_CNTRL/$sample/gene_discovery/further_filtering"

  python plots/plot_downstream_sv_filtering.py \
    --filter-dir "outs_test_CNTRL/$sample/gene_discovery/further_filtering" \
    --out-dir "outs_test_CNTRL/$sample/plots/15_downstream_filtering"
done
```

If your output root differs, replace `outs_test_CNTRL` with its actual path.
Do not pass the gene-level ranking, raw AnnotSV report, or a table with repeated
SV–gene transcript rows. The script accepts a compressed TSV as well.

For a stricter inversion sensitivity comparison, use `--inv-min-reads 5` in a
**different output directory**. For another AF threshold, use `--max-af 0.001`,
again recording it as a separate analysis. These options are your study-specific
filters; they should not be presented as published clinical thresholds.

## Open these outputs first

| File or folder | What it contains |
| --- | --- |
| `START_HERE.txt` | Sample-specific reading instructions. |
| `filter_run.json` | Exact source, threshold settings and definitions. |
| `output_files.tsv` | Current files, SV–gene association counts and unique SV counts. |
| `gt_1kb/retained.tsv` | Main >1 kb view, preserving original order and rank values. |
| `gt_1kb/caller_priority.tsv` | Caller-priority presentation, preserving original order within caller groups. |
| `lt_1kb/retained.tsv` | Small candidates passing the other rules. |
| `eq_1kb/retained.tsv` | Exact-boundary candidates passing the other rules. |
| `<size>/by_sv_type/DEL/`, `DUP/`, `INV/`, `INS/`, etc. | Separate SV-type tables and exonic/splice, intronic or mixed tables where present. |
| `<size>/by_af/af_le_1_percent_or_missing.tsv` | One combined AF table; the original measured/missing AF status stays traceable in its columns. |
| `<size>/by_caller_count/` | Separate caller-count groups. |
| `<size>/by_inv_reads/` | Lower and higher inversion read-support tiers where present. |
| `<size>/review_required.tsv` and `review_by_sv_type/` | Missing/invalid annotations or support, unusable length, and conflicts. |
| `excluded/all_excluded.tsv` | Failed criteria with all recorded exclusion reasons. |
| `all_records.annotated.tsv` | Every original row and original column, plus explicit `DOWNSTREAM_*` decisions. |

AF ≤1% and missing AF are combined in `retained.tsv` and `caller_priority.tsv`.
`DOWNSTREAM_AF_REVIEW_GROUP` records this common group. Missing AF is retained
to allow review for possible new SVs, but **does not establish novelty or ultra-rarity**.
It can reflect no compatible population match, absent annotation, or an event
that the frequency procedure could not evaluate. An unconfigured resource is also missing evidence, not an absent population variant.
Explicitly rejected matches cannot supply AF. A failed population-resource FILTER
or a malformed AF value goes to review. Invalid AF units or values go to review. If one valid source
reports AF above the cutoff, the event is excluded even if another source is
missing or reports a smaller frequency.

Generated views overlap: the same event can be present in a type table, an AF
table and a caller table. Do not add their counts. Multiple affected genes are
multiple SV–gene associations, not independent SVs. A record with any known
failed criterion goes to the excluded export; otherwise unresolved information
goes to review. Exclusion from this view is not evidence of benignity.

Reruns refresh only TSVs listed as generated by the previous run. Unrelated
notes are preserved. Do not write reviewer notes into generated exports; keep
them in the editable candidate-review table created by the review launcher.

## Five useful plots

The plotting command creates separate PNG, PDF and SVG files:

| Plot folder | Figure | Counting unit |
| --- | --- | --- |
| `filter_outcomes/` | Retained, review-required and excluded records by size group | SV–gene associations; per-outcome unique SV counts are separately exported and can overlap. |
| `retained_types_and_context/` | Large retained SV types and annotation context | SV–gene associations. |
| `retained_types_and_context/` | Small retained SV types and annotation context | SV–gene associations. |
| `callers_and_population/` | One stacked bar per caller count, combining AF ≤1% and missing AF; colour preserves AF status | Unique SVs, counted once regardless of affected-gene count. |
| `inversion_support/` | Retained large inversions: 3–4 versus ≥5 reported reads, by caller count | Unique SVs. |

These plots explain what the filters did. They do not compare diagnostic yield
between the two samples or demonstrate enrichment without an appropriate
background. Keep unresolved/excluded counts visible so a smaller candidate
list cannot hide missing annotations.

For detailed review of the filtered main table, reuse the existing launcher in
a separate output folder, for example:

```bash
python plots/run_thesis_plots.py \
  --root outs_test_CNTRL --sample 21930 --platform lrs --review-only \
  --candidate-table outs_test_CNTRL/21930/gene_discovery/further_filtering/gt_1kb/retained.tsv \
  --gene-bed reference/gencode.genes.hg38.sorted.bed \
  --out-dir outs_test_CNTRL/21930/plots/16_filtered_candidate_review
```

The evidence matrices still follow the original pipeline ranks. To review a
specific caller group, supply its `by_caller_count/...tsv` as the candidate
table. This keeps caller preference distinct from the biological ranking.

## Scientific criticism of this filter design

A 1 kb cutoff is useful for focusing on larger rearrangements but can remove
small damaging exon/splice events from that view. Retaining a small-event
review branch is therefore necessary. Size alone does not establish impact.

The 1% AF rule is broad. The gnomAD team's v4 SV description reports that 96%
of its SV sites have AF below 1%; this threshold alone will not make a list
specific to a rare disease. Unknown AF is a question to investigate, not
stronger evidence than a measured low AF.

Three callers use substantially the same sequencing data and may share
alignment artefacts. Three to five supporting reads is a configurable triage
rule; reliability also depends on depth, mapping, repeats and breakpoint
consistency. Intronic events can matter through splice or regulatory mechanisms,
but an intronic annotation alone does not demonstrate either mechanism.

Keep these sources with the methodological explanation:

- [gnomAD v4 structural variants: frequency, validation and functional context](https://gnomad.broadinstitute.org/news/2023-11-v4-structural-variants).
- [Ensembl VEP structural-variant prediction and reported overlaps](https://jun2026.archive.ensembl.org/info/docs/tools/vep/script/vep_example.html#sv).

The thresholds in this script are user-selected research filters, not rules
derived from these resources. After filtering, the strongest next evidence
usually comes from event-specific locus/read review and mechanism assessment.
