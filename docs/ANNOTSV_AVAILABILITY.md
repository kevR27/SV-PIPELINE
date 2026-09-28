# Automatic AnnotSV evidence availability

When only saved results are available, omit `--annotations-dir` to inspect the
TSV itself. This checks source columns and observed database mentions without
claiming that the installed resources were inspected. See the
[saved-output guide](SAVED_OUTPUT_REANALYSIS.md) for the complete rebuild command.

This report requires no trio data, patient HPO terms or manually filled TSVs. It
helps interpret database overlap evidence for the existing genome-wide SV study.
It does not add pathogenic inversion records or classify an SV clinically.

## What is checked

There are **eight source columns**, not six: P/B × loss/gain/ins/inv. Separately,
the integration summarizes six database labels: ClinVar, dbVar, gnomAD, DGV, 1000G
and ClinGen. The audit reports all 48 field/database combinations.

After AnnotSV finishes, the sibling `annotsv_evidence_audit` rule automatically:

1. Inspects the sample TSV header and counts nonempty source cells/database mentions.
2. Inspects the eight prepared SV resources for the same assembly under the configured
   `annotsv_annotations_dir`. It verifies the AnnotSV 3.5.10 BED layout, counts source
   records and database mentions, and records file paths and SHA256 checksums.
3. Writes a per-sample JSON/TSV report, then supplies it to the integrated evidence table.

The report runs after annotation because AnnotSV may prepare these resources lazily.
The expensive resource scan is cached across samples using file paths, sizes and
mtime/ctime. Changed resources invalidate the cache; per-sample columns are always
checked separately. Cache writes are atomic. Concurrent first-time audits may each
scan the bundle once. Rerun after resource preparation completes if files change
mid-scan. No source files are downloaded, altered or deleted.

Prepared paths are based on the v3.5.10 upstream implementation:

- `Annotations_Human/FtIncludedInSV/PathogenicSV/GRCh38/pathogenic_<Loss|Gain|Ins|Inv>_SV_GRCh38.sorted.bed`
- `Annotations_Human/SVincludedInFt/BenignSV/GRCh38/benign_<Loss|Gain|Ins|Inv>_SV_GRCh38.sorted.bed`

The script also accepts GRCh37/CHM13 for direct use; the current workflows explicitly
use GRCh38. Unsupported resource layouts are reported as unknown/unrecognized rather
than treated as empty. A downloaded raw database elsewhere in the bundle is not proof
that it contributed records to these processed, SV-type-specific annotation channels.

## Outputs

Per sample in `<output-root>/<sample>/sv/annotsv/`:

- `<sample>_annotation_availability.tsv`: readable field/database inventory.
- `<sample>_annotation_availability.json`: structured inventory with sample/bundle provenance.

The shared cache is `<output-root>/reference_checks/annotsv_bundle_inventory.json`.
It is an implementation cache, not a manually prepared input file. Integration verifies
that the audit's sample TSV checksum matches its actual AnnotSV input, preventing use
of a stale sample report.

Existing `SV_DB_*_OVERLAP` flags keep their YES / NOT_REPORTED / UNKNOWN /
NOT_APPLICABLE values. Additional fields explain availability without changing those
flags or breaking older summaries:

| New status | Meaning |
|---|---|
| `OVERLAP_REPORTED` | The event's source field reports relevant overlap evidence. It is not proof of the identical allele or pathogenicity. |
| `SOURCE_COLUMN_MISSING` | An expected source column is absent from this TSV; inspect output selection/version before interpreting it. |
| `PROCESSED_RESOURCE_NOT_FOUND` | The expected prepared BED was not found at the known location. This does not prove the whole database was never installed. |
| `PROCESSED_RESOURCE_EMPTY` | The prepared file exists but contains no data records. |
| `RESOURCE_SCHEMA_UNRECOGNIZED` / `RESOURCE_UNREADABLE` | The resource could not be confidently evaluated. |
| `RESOURCE_SOURCE_VALUES_EMPTY` | BED records exist but their source cells are empty. |
| `NO_OVERLAP_REPORTED` | Relevant reporting capability/resource evidence is present, but no overlap was reported for this event. This is **not** a proven biological absence or a benign classification. |
| `DATABASE_NOT_REPRESENTED_IN_PROCESSED_RESOURCES` | Both applicable processed channels were successfully inspected, but neither contains a record naming this database. Applies to this assembly, SV type and processed bundle—not every use of that database in AnnotSV. |
| `DATABASE_OBSERVED_AVAILABILITY_INCOMPLETE` / `AVAILABILITY_UNDETERMINED` | Available observations are insufficient to resolve the full resource situation. |
| `NO_ANNOTSV_MATCH` | No AnnotSV row was linked to the master SV; do not interpret this as no overlap. |
| `AVAILABILITY_NOT_AUDITED` | Direct integration was run without the optional audit JSON. |
| `NOT_APPLICABLE` | No corresponding dedicated source-field mapping is used for this SV type, such as BND. |

`SV_PATHOGENIC_DB_STATUS` and `SV_BENIGN_DB_STATUS` describe each side separately.
The six `SV_DB_<DATABASE>_AVAILABILITY` columns assess both sides conservatively.
Positive output evidence takes precedence over a missing local resource: output may
have been generated using an earlier installation. The inventory still exposes that
resource discrepancy. Output counts are TSV rows, not unique SV counts.

## Inversions and missing evidence

A P_inv field in the schema does not guarantee many pathogenic inversion records.
The upstream v3.5.10 ClinVar parser can populate inversion records; different source
loaders do not necessarily cover every SV class. Missing P_inv values can also result
from output selection, no qualifying overlap, or absent/empty prepared data. The audit
measures those distinctions for the installed bundle and the actual sample rather than
assuming that every missing inversion annotation is the same problem.

Presence of nonempty fields in another sample is useful evidence that a channel works,
but an all-empty sample column alone cannot distinguish absent resources from no
qualifying overlaps. Therefore output-column inspection is combined with prepared
resource inspection, and neither is called a complete negative database search.

## Running it

Both LRS and SRS integrated-table rules depend on the report automatically. Regenerate
the integrated table with the updated main workflow, then rerun postprocessing.
Postprocessing alone consumes completed upstream tables; it does not reannotate old TSVs.
No new configuration keys or four optional allele-evidence files are required.

For a standalone audit of an existing result:

```bash
python scripts/audit_annotsv_evidence.py \
  --annotations-dir /path/to/AnnotSV_annotations \
  --annotsv /path/to/SAMPLE_merged_SV.annotsv.tsv \
  --genome-build GRCh38 \
  --cache /path/to/output/reference_checks/annotsv_bundle_inventory.json \
  --output-json /path/to/SAMPLE_annotation_availability.json \
  --output-tsv /path/to/SAMPLE_annotation_availability.tsv
```

The script does not modify the input TSV. Direct integration can consume the result
with `--annotsv-audit /path/to/SAMPLE_annotation_availability.json`.

## Panel/nonpanel and VEP status

- Panel filtering uses exact case-insensitive gene symbols; `OPA1` does not match `OPA10`.
- The nonpanel gene list is the set difference between all overlapping genes and the panel.
- A full AnnotSV row spanning both panel and nonpanel genes remains in the panel view
  when any panel gene overlaps. This does not label every gene on that row as a panel gene.
- Plot/summary labels use exact PANEL_GENE/YES checks, so NONPANEL_GENE is not mislabeled.
- Legacy `Snakefile_LRS` includes `Snakefile_LRS_update`; both active LRS and SRS VEP rules
  use `--flag_pick`. No old `--pick` command remains in the legacy entry point.

Upstream source layout and interpretation:
[Pathogenic SV annotations](https://github.com/lgmgeo/AnnotSV/blob/v3.5.10/share/tcl/AnnotSV/AnnotSV-pathogenicsv.tcl),
[benign SV annotations](https://github.com/lgmgeo/AnnotSV/blob/v3.5.10/share/tcl/AnnotSV/AnnotSV-benignsv.tcl).
