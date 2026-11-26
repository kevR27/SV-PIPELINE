# VCF-only analysis: Genomiser + custom ranking

This standalone workflow accepts existing genome-wide SNV/indel and merged
SV VCFs and performs both independent Genomiser and shared custom HON/AnnotSV
candidate ranking. It does not register any BAM-based calling or coverage rules.

Use the [full VCF-only guide](../../docs/GENOMISER_VCF_ONLY.md) for input
requirements, server configuration, execution, outputs and benchmarking.
Keep clinical profiles and confirmed truth TSVs private. The result path must
be separate from the original LRS/SRS analysis.
\n## Benchmarking and ranking stability

The normal VCF-only `all` target also evaluates fixed-rank sensitivity to
population-AF and multi-caller eligibility filters and writes
`benchmark/<sample>.benchmark_review.tsv` plus JSON provenance. This is
distinct from the optional independently confirmed-variant benchmark.
When caller counts or population AF are unavailable, results explicitly
mark those scenarios as unevaluable/unknown; no missing evidence is treated
as negative. No BAMs or external transcript-expression cohorts are required.
No storage cleanup is registered in this workflow.
