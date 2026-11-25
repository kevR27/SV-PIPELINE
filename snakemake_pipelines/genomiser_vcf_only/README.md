# VCF-only analysis: Genomiser + custom ranking

This standalone workflow accepts existing genome-wide SNV/indel and merged
SV VCFs and performs both independent Genomiser and shared custom HON/AnnotSV
candidate ranking. It does not register any BAM-based calling or coverage rules.

Use the [full VCF-only guide](../../docs/GENOMISER_VCF_ONLY.md) for input
requirements, server configuration, execution, outputs and benchmarking.
Keep clinical profiles and confirmed truth TSVs private. The result path must
be separate from the original LRS/SRS analysis.
