# Score interpretation with existing pipeline outputs

The four optional allele TSV settings may remain `null`. This does not remove
SVs, genes, database records, the original gene-discovery score or AnnotSV's
recorded classification. They are inputs to a later, separate assessment stage.
They do not automatically populate themselves from other outputs.

| Output | Meaning | Effect of all four optional TSV inputs being absent |
| --- | --- | --- |
| `PHENOTYPE_SCORE` (0–13) | Human gene–HPO associations against configured optic-neuropathy anchors | Retained; this is gene-level relevance, not patient-specific symptom matching. |
| `INTEGRATED_DISCOVERY_SCORE` (0–19) | Gene-level phenotype, curated gene–disease evidence and a small SV-count component | Retained. It does not use the four optional files. |
| `ANNOTSV_RANKING_SCORE`, `ACMG_CNV_SCORE`, `ACMG_CNV_CLASS` | Recorded AnnotSV evidence/classification, with explicit DEL/DUP CNV scope | Retained as supplied by the current integration. No new clinical classification is calculated by the allele script. |
| `ALLELE_RESEARCH_SCORE` (0–12) | Uncalibrated sum of six evidence domains | Limited by missing data. Unknown domains earn no points; low totals do not establish benignity. |

There is no single validated "overall pathogenic score" in this workflow.
Neither a high gene-discovery score nor presence of a pathogenic database record
establishes pathogenicity of a particular allele. Large SVs can affect several
genes; an event classification is not a pathogenicity classification of every
overlapping gene. The gene-ranking table labels its aggregate AnnotSV fields as
`OVERLAPPING_SV_RECORDS_NOT_GENE_PATHOGENICITY`.

For a study without trios or patient-specific HPO terms, use the original gene
ranking and the separate allele evidence/status columns for candidate review.
Do not compare incomplete 0–12 totals with complete cases as if they were
calibrated probabilities, and do not exclude a candidate for missing data.
Disease-specific inheritance or mechanism compatibility remains unknown where
the required relationship is not available. Existing genotype fields can still
be reported without claiming segregation or de novo status.

With no optional inputs, population matches from needLR remain provisional;
they do not earn confirmed allele-rarity points. Automatically predicted
disruption, technical evidence, and gene-level dosage support can contribute
when available. The fixed maximum of 12 is not rescaled to the available data.
Generic optic-neuropathy anchors are never substituted for patient phenotypes.

## Corrected reporting behavior

- Allele assessment accepts long VCF INFO/read-name and transcript JSON fields.
  Its model identifier is now `sv-allele-evidence-v1.2`. This version also checks
  the correct BND endpoint for a gene and reports functional context separately.
- Missing OMIM tokens (`NA`, `N/A`, `None`, `NaN`, `null`, `.` and empty cells,
  case-insensitive) cannot provide positive disease-evidence points.
- `Animal Model Only` GenCC records remain visible but receive no human
  gene–disease evidence points. The gene model is now
  `phenotype13_geneDisease4_sv2_v2.1`.
- Full-row AnnotSV classification/score fields remain available in gene
  summaries as overlapping-event information. Gene-specific OMIM/GenCC data
  still require gene-specific rows.
- An expected missing `SVLEN` on BND/TRA is informational; it alone no longer
  zeros otherwise available technical-support points. Other QC flags retain
  their review behavior.
- Absent complementary source files produce `NOT_AVAILABLE`, while supplied
  files with no reported match retain `NO`. Gene summaries count unavailable
  source records separately.
- Evidence-matrix ordering uses the integrated gene-discovery score, with the
  phenotype score as an explicit legacy fallback. Panel status and optional
  context no longer silently create another numerical ranking. The exported
  matrix records `PLOT_ORDER_BASIS` and `PLOT_ORDER_SCORE`.
- The matrix displays Manta support for the SRS workflow as well as LRS callers.
- The matrix deduplicates by SV ID plus gene, retaining different events that
  share a start coordinate/type/gene. Optional match statuses that are unknown
  or unavailable are shown separately from a reported non-match.
- Context integration and gene-summary generation accept an empty callset.
  This does not guarantee that every optional plotting script can draw an
  empty dataset.

## Applying corrections to existing results

The [saved-output rebuild command](SAVED_OUTPUT_REANALYSIS.md) runs these steps
in order and writes to a new folder. It accepts the old raw tool outputs and
does not require the four clinical TSVs.

Updating repository code does not rewrite TSVs already produced on a server.
Regenerate the gene ranking, integrated table, allele assessment, downstream
context/summary tables, and selected plots from the original tool outputs in
that dependency order. Keep the original outputs and a record of the code and
database versions used. Old integrated tables lacking caller/transcript JSON
cannot recover all current evidence by running only the allele script.

Some regenerated scores can change because erroneous/duplicated/missing
evidence is handled correctly. That is separate from setting optional inputs
to `null`. In particular, deduplicating gene–HPO associations can reduce the
phenotype component of old rankings.

The AnnotSV resource-availability audit is a separate change. It does not
populate absent inversion evidence. An unmatched SV or missing database cell
must not be interpreted as an assessed negative result.

## Study questions and scope

Genome-wide SV discovery, panel/nonpanel candidate investigation, gene-level
human phenotype relevance, and supplementary repeat/MEI/phasing/methylation
context are supported analyses. Their outputs are candidate evidence.

Matched LRS–SRS concordance requires an explicit comparison of matched samples
and compatible regions/representations. The existing comparison plot consumes
an externally generated Truvari summary; it does not perform that comparison.
Gene recurrence and equivalent-SV recurrence across patients likewise require
cohort-level aggregation; within-sample caller merging does not supply these.

Diagnoses, improved diagnostic yield, new gene–disease relationships, functional
regulatory effects, and causal methylation changes require study-specific
validation beyond these scores. Lack of trios does not prevent a discovery
study; it limits the conclusions available from inheritance evidence.
