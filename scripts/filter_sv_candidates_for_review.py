#!/usr/bin/env python3
"""Make separate downstream SV review tables without changing the input.

Use the final ranked SV-gene table, not an AnnotSV transcript report or a
gene-level ranking. The added DOWNSTREAM_* columns explain every decision.
Thresholds organise research review; they do not classify pathogenicity.
"""

import argparse
import json
import math
import re
from pathlib import Path

import pandas as pd


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL", "UNKNOWN"}
EXONIC_TERMS = {
    "coding_sequence_variant", "exon_loss_variant", "transcript_ablation",
    "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
    "missense_variant", "synonymous_variant", "inframe_deletion",
    "inframe_insertion", "protein_altering_variant",
    "5_prime_utr_variant", "3_prime_utr_variant",
    "non_coding_transcript_exon_variant", "splice_acceptor_variant",
    "splice_donor_variant", "splice_region_variant",
}
CONTEXT_NAMES = ["EXONIC_OR_SPLICE", "INTRONIC", "EXONIC_AND_INTRONIC"]


# 1. Read only explicit values. Missing information stays missing.
def value(row, columns):
    for column in columns:
        item = str(row.get(column, ".")).strip()
        if item.upper() not in MISSING:
            return item
    return "."


def number(item):
    try:
        result = float(item)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def tokens(item):
    return {x.strip().upper() for x in re.split(r"[;,|&]", str(item))
            if x.strip().upper() not in MISSING}


def chromosome(item):
    return str(item).strip().upper().removeprefix("CHR")


def chromosome_scope(row):
    endpoints = [chromosome(value(row, ["CHROM", "chrom"]))]
    partner = value(row, ["CHR2", "INFO_CHR2"])
    if partner != ".":
        endpoints.append(chromosome(partner))
    # A translocation can connect an autosome to Y or MT. Check its ALT too.
    endpoints += [chromosome(c) for c in re.findall(r"[\[\]]([^:\[\]]+):\d+", value(row, ["ALT"]))]
    if "Y" in endpoints:
        return "EXCLUDE_Y"
    if set(endpoints) & {"M", "MT"}:
        return "EXCLUDE_MITOGENOME"
    if endpoints[0] == ".":
        return "CHROMOSOME_UNKNOWN"
    return "NUCLEAR_WITHOUT_Y"


# 2. Measure event size, not gene size. BND/TRA have no single interval length.
def event_size(row):
    kind = value(row, ["SVTYPE"]).upper()
    if kind in {"BND", "TRA"}:
        return None, "NOT_APPLICABLE_BREAKEND"
    reported = number(value(row, ["SVLEN"]))
    if reported is not None and reported.is_integer() and abs(reported) > 0:
        return abs(reported), "SVLEN"
    if kind == "INS":
        # END minus START is NOT insertion length.
        ref, alt = value(row, ["REF"]), value(row, ["ALT"])
        if re.fullmatch(r"[ACGTNacgtn]+", ref) and re.fullmatch(r"[ACGTNacgtn]+", alt) and len(alt) > len(ref):
            return len(alt) - len(ref), "LITERAL_INSERTED_SEQUENCE"
        return None, "INSERTION_LENGTH_UNKNOWN"
    if kind not in {"DEL", "DUP", "INV", "CNV"}:
        return None, "UNSUPPORTED_INTERVAL_TYPE"
    for column in ["SV_SPAN_BP", "SV_EVENT_SPAN_BP"]:
        span = number(row.get(column))
        if span is not None and span.is_integer() and span > 0:
            return span, column
    start = number(value(row, ["START", "POS"]))
    end = number(value(row, ["END"]))
    remote = value(row, ["CHR2"])
    if remote != "." and chromosome(remote) != chromosome(value(row, ["CHROM"])):
        return None, "INTERCHROMOSOMAL_LENGTH_UNKNOWN"
    if start is not None and end is not None and end > start:
        return end - start, "END_MINUS_START"
    return None, "LENGTH_UNKNOWN"


# 3. Keep exon/splice and intron annotations separate from unresolved geometry.
def functional_context(row):
    status = value(row, ["VEP_MATCH_STATUS"]).upper()
    if status != "." and status not in {"GENE_TRANSCRIPT_MATCH", "MATCH", "EXACT"}:
        return "GENIC_UNRESOLVED", "VEP_GENE_MATCH_NOT_ESTABLISHED"
    regions = tokens(value(row, ["VEP_TRANSCRIPT_REGION_CLASS"]))
    consequences = tokens(value(row, ["VEP_CONSEQUENCES", "Consequence"]))
    effects = tokens(value(row, ["VEP_STRUCTURAL_EFFECT"]))
    whole = bool(regions & {"TRANSCRIPT_ABLATION"} or effects & {"WHOLE_TRANSCRIPT_LOSS", "WHOLE_TRANSCRIPT_GAIN"})
    exonic = bool(regions & {"EXONIC", "EXONIC_ONLY", "EXONIC_OR_SPLICE", "EXON_LOSS", "CODING", "UTR", "EXONIC_UTR"})
    exonic = exonic or whole or bool(consequences & {x.upper() for x in EXONIC_TERMS}) or value(row, ["VEP_EXON"]) != "."
    intronic = bool(regions & {"INTRONIC", "INTRONIC_ONLY"}) or "INTRON_VARIANT" in consequences or value(row, ["VEP_INTRON"]) != "."
    if regions & {"MIXED_EXON_INTRON", "EXONIC_AND_INTRONIC"}:
        exonic, intronic = True, True
    scope = "ANNOTATED_FEATURES_NOT_EXACT_BREAKPOINTS"
    if whole:
        scope = "WHOLE_TRANSCRIPT_EFFECT_NOT_BREAKPOINT_LOCATION"
    if exonic and intronic:
        return "EXONIC_AND_INTRONIC", scope
    if exonic:
        return "EXONIC_OR_SPLICE", scope
    if intronic:
        return "INTRONIC", scope
    if regions & {"REGULATORY_PROXIMAL", "INTERGENIC"} or consequences & {"UPSTREAM_GENE_VARIANT", "DOWNSTREAM_GENE_VARIANT", "INTERGENIC_VARIANT"}:
        return "NON_GENIC_CONTEXT", scope
    return "GENIC_UNRESOLVED", "EXON_INTRON_INFORMATION_UNAVAILABLE"


# 4. Use population AF, never sample allele balance. Unknown AF is not zero.
def population_af(row, cutoff):
    observed, sources, notes = [], [], []
    for column in ["NEEDLR_AF", "GNOMAD_SV_AF"]:
        raw = value(row, [column])
        if raw == ".":
            continue
        status = value(row, ["NEEDLR_STATUS", "POPULATION_STATUS"]).upper()
        rejected = column == "NEEDLR_AF" and any(x in status for x in ["NO_MATCH", "NO_POPULATION_MATCH", "NOT_EVALUABLE", "NOT_AVAILABLE"])
        rejected = rejected or (column == "GNOMAD_SV_AF" and value(row, ["GNOMAD_SV_EXACT_MATCH"]).upper() in {"NO", "NO_MATCH", "NOT_AVAILABLE"})
        if rejected:
            notes.append(column + "_IGNORED_REJECTED_MATCH")
            continue
        for item in re.split(r"[;,|]", raw):
            if item.strip().upper() in MISSING:
                continue
            af = number(item)
            if af is None or not 0 <= af <= 1:
                notes.append(column + "_INVALID_AF")
            else:
                observed.append(af)
                sources.append(column)
    maximum = max(observed) if observed else None
    if maximum is not None and maximum > cutoff:
        status = "MEASURED_ABOVE_CUTOFF"
    elif any("INVALID_AF" in note for note in notes):
        status = "INVALID_AF_REVIEW"
    elif maximum is not None:
        status = "MEASURED_AT_OR_BELOW_CUTOFF"
    else:
        status = "AF_UNKNOWN"
    return maximum, status, ";".join(sorted(set(sources))) or ".", ";".join(notes) or "."


# 5. Count callers separately from reads. Never add the same reads across tools.
def caller_count(row):
    count = number(value(row, ["CALLER_COUNT", "SUPP"]))
    callers = tokens(value(row, ["CALLERS"]))
    callers = {"SNIFFLES2" if c == "SNIFFLES" else c for c in callers}
    if count is not None and count.is_integer() and count >= 1:
        return int(count), "CALLER_COUNT_CONFLICT" if callers and len(callers) != count else "."
    if callers:
        return len(callers), "COUNT_FROM_CALLER_NAMES"
    return None, "CALLER_COUNT_UNKNOWN"


def max_read_support(row):
    for column in ["EVENT_MAX_CALLER_READ_SUPPORT", "CALLER_READ_SUPPORT", "READ_SUPPORT"]:
        raw = value(row, [column])
        if raw == ".":
            continue
        counts = []
        for item in re.split(r"[;,|]", raw):
            count = number(item.rsplit(":", 1)[-1])
            if count is not None and count.is_integer() and count >= 0:
                counts.append(count)
        if counts:
            return max(counts), column
    try:
        records = json.loads(value(row, ["CALLER_EVIDENCE_JSON"]))
    except (ValueError, TypeError):
        records = []
    if not isinstance(records, list):
        records = []
    counts = [number(record.get("CALLER_SUPPORT")) for record in records if isinstance(record, dict)]
    counts = [n for n in counts if n is not None and n.is_integer() and n >= 0]
    return (max(counts), "CALLER_EVIDENCE_JSON") if counts else (None, "READ_SUPPORT_UNKNOWN")


def annotate_row(row, sample, source_row, cutoff, af_cutoff, inv_min, inv_strong):
    size, size_source = event_size(row)
    if size is None:
        size_group = "LENGTH_UNRESOLVED"
    elif size > cutoff:
        size_group = "GT_CUTOFF"
    elif size < cutoff:
        size_group = "LT_CUTOFF"
    else:
        size_group = "AT_CUTOFF"
    context, context_scope = functional_context(row)
    af, af_status, af_sources, af_notes = population_af(row, af_cutoff)
    callers, caller_note = caller_count(row)
    reads, read_source = max_read_support(row)
    kind = value(row, ["SVTYPE"]).upper()
    chrom_status = chromosome_scope(row)
    excluded, review = [], []
    if value(row, ["GENE", "GENES", "Gene", "ANNotsv_Gene"]) == ".":
        review.append("GENE_IDENTITY_UNRESOLVED")
    if chrom_status.startswith("EXCLUDE_"):
        excluded.append(chrom_status)
    elif chrom_status == "CHROMOSOME_UNKNOWN":
        review.append(chrom_status)
    if af_status == "MEASURED_ABOVE_CUTOFF":
        excluded.append("AF_ABOVE_CUTOFF")
    elif af_status == "INVALID_AF_REVIEW":
        review.append("INVALID_POPULATION_AF")
    if context == "NON_GENIC_CONTEXT":
        excluded.append("NO_EXON_OR_INTRON_ANNOTATION")
    elif context == "GENIC_UNRESOLVED":
        review.append("EXON_INTRON_UNRESOLVED")
    if caller_note == "CALLER_COUNT_CONFLICT":
        review.append(caller_note)
    inv_group = "NOT_APPLICABLE"
    if kind == "INV":
        if reads is None:
            inv_group = "UNKNOWN"
        elif reads < inv_min:
            inv_group = "BELOW_MINIMUM"
        elif reads < inv_strong:
            inv_group = "MINIMUM_TO_STRONG"
        else:
            inv_group = "AT_OR_ABOVE_STRONG"
        if reads is None:
            review.append("INV_READ_SUPPORT_UNKNOWN")
        elif reads < inv_min:
            excluded.append("INV_READ_SUPPORT_BELOW_MINIMUM")
        if callers is None:
            review.append("INV_CALLER_COUNT_UNKNOWN")
        if "AMBIGUOUS" in value(row, ["CALLER_EVIDENCE_MATCH"]).upper():
            review.append("INV_AMBIGUOUS_READ_EVIDENCE_JOIN")
    inv_scope = "NOT_APPLICABLE"
    if kind == "INV":
        relation = value(row, ["SV_GENE_RELATIONSHIP", "SV_GENE_EFFECT"]).upper()
        if "SPANS_INTACT_GENE" in relation or "SPANNED_BY_INVERSION" in relation:
            inv_scope = "INTACT_GENE_SPAN_CONTEXT"
        elif "BREAKPOINT" in relation and any(term in relation for term in ["WITHIN_TRANSCRIPT", "IN_TRANSCRIPT", "IN_GENE"]):
            inv_scope = "BREAKPOINT_TRANSCRIPT_CONTEXT"
        else:
            inv_scope = "INV_MECHANISM_UNRESOLVED"
    if size is None:
        review.append("NO_SINGLE_USABLE_EVENT_LENGTH")
    return {
        "DOWNSTREAM_SAMPLE": sample, "DOWNSTREAM_SOURCE_ROW": source_row,
        "DOWNSTREAM_SV_ID": value(row, ["SV_ID", "ID"]), "DOWNSTREAM_GENE": value(row, ["GENE", "GENES", "Gene", "ANNotsv_Gene"]),
        "DOWNSTREAM_SVTYPE": kind, "DOWNSTREAM_SIZE_BP": size if size is not None else ".",
        "DOWNSTREAM_SIZE_SOURCE": size_source, "DOWNSTREAM_SIZE_GROUP": size_group,
        "DOWNSTREAM_CHROMOSOME_SCOPE": chrom_status, "DOWNSTREAM_CONTEXT": context,
        "DOWNSTREAM_CONTEXT_SCOPE": context_scope,
        "DOWNSTREAM_MAX_AF": af if af is not None else ".", "DOWNSTREAM_AF_STATUS": af_status,
        "DOWNSTREAM_AF_SOURCES": af_sources, "DOWNSTREAM_AF_NOTES": af_notes,
        "DOWNSTREAM_CALLER_COUNT": callers if callers is not None else ".", "DOWNSTREAM_CALLER_NOTE": caller_note,
        "DOWNSTREAM_MAX_READ_SUPPORT": reads if reads is not None else ".", "DOWNSTREAM_READ_SUPPORT_SOURCE": read_source,
        "DOWNSTREAM_INV_READ_GROUP": inv_group,
        "DOWNSTREAM_INV_EFFECT_SCOPE": inv_scope,
        "DOWNSTREAM_EXCLUSION_REASONS": ";".join(excluded) or ".",
        "DOWNSTREAM_REVIEW_REASONS": ";".join(review) or ".",
    }


# 6. Save organised views. Original columns, ranks and source order are kept.
def save_views(work, out, size_cutoff):
    manifest = []
    def save(frame, relative_path):
        path = out / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, sep="\t", index=False)
        manifest.append({"FILE": relative_path, "SV_GENE_ASSOCIATIONS": len(frame), "UNIQUE_SVS": frame["DOWNSTREAM_SV_ID"].nunique()})

    label = f"{size_cutoff // 1000}kb" if size_cutoff % 1000 == 0 else f"{size_cutoff}bp"
    save(work, "all_records.annotated.tsv")
    save(work[work["DOWNSTREAM_DECISION"].eq("EXCLUDED")], "excluded/all_excluded.tsv")
    for group, folder in [("GT_CUTOFF", f"gt_{label}"), ("LT_CUTOFF", f"lt_{label}"), ("AT_CUTOFF", f"eq_{label}"), ("LENGTH_UNRESOLVED", "length_unresolved")]:
        scoped = work[work["DOWNSTREAM_SIZE_GROUP"].eq(group) & ~work["DOWNSTREAM_DECISION"].eq("EXCLUDED")]
        kept = scoped[scoped["DOWNSTREAM_DECISION"].eq("RETAINED")]
        save(kept, folder + "/retained.tsv")
        save(scoped[scoped["DOWNSTREAM_DECISION"].eq("REVIEW_REQUIRED")], folder + "/review_required.tsv")
        unresolved = scoped[scoped["DOWNSTREAM_DECISION"].eq("REVIEW_REQUIRED")]
        for kind, types in unresolved.groupby("DOWNSTREAM_SVTYPE", sort=True):
            safe_kind = re.sub(r"[^A-Z0-9_-]", "_", kind) or "UNKNOWN"
            save(types, folder + f"/review_by_sv_type/{safe_kind}.tsv")
        ordered = kept.assign(_caller=pd.to_numeric(kept["DOWNSTREAM_CALLER_COUNT"], errors="coerce"))
        ordered = ordered.sort_values("_caller", ascending=False, na_position="last", kind="stable").drop(columns="_caller")
        save(ordered, folder + "/caller_priority.tsv")
        for kind, types in kept.groupby("DOWNSTREAM_SVTYPE", sort=True):
            safe_kind = re.sub(r"[^A-Z0-9_-]", "_", kind) or "UNKNOWN"
            save(types, folder + f"/by_sv_type/{safe_kind}/all.tsv")
            for context, subset in types.groupby("DOWNSTREAM_CONTEXT", sort=True):
                save(subset, folder + f"/by_sv_type/{safe_kind}/{context.lower()}.tsv")
        for af_group, subset in kept.groupby("DOWNSTREAM_AF_STATUS", sort=True):
            save(subset, folder + f"/by_af/{af_group.lower()}.tsv")
        inversions = kept[kept["DOWNSTREAM_SVTYPE"].eq("INV")]
        for read_group, subset in inversions.groupby("DOWNSTREAM_INV_READ_GROUP", sort=True):
            save(subset, folder + f"/by_inv_reads/{read_group.lower()}.tsv")
        for count, subset in ordered.groupby("DOWNSTREAM_CALLER_COUNT", sort=False):
            label_count = "unknown" if count == "." else str(count)
            save(subset, folder + f"/by_caller_count/{label_count}_callers.tsv")
    pd.DataFrame(manifest).to_csv(out / "output_files.tsv", sep="\t", index=False)
    return [x["FILE"] for x in manifest] + ["output_files.tsv"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Final ranked SV-gene TSV/TSV.GZ")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--size-cutoff", type=int, default=1000, help="Main size view is strictly greater than this many bp")
    parser.add_argument("--max-af", type=float, default=0.01, help="Population AF in 0..1 units; unknown AF is retained separately")
    parser.add_argument("--inv-min-reads", type=int, default=3)
    parser.add_argument("--inv-strong-reads", type=int, default=5, help="Higher reported read-support tier, not validation")
    args = parser.parse_args()
    if args.size_cutoff < 1 or not 0 <= args.max_af <= 1 or not 1 <= args.inv_min_reads <= args.inv_strong_reads:
        parser.error("Require size-cutoff >0, AF in 0..1, and 1 <= inv-min-reads <= inv-strong-reads")
    source, out = Path(args.input).resolve(), Path(args.out_dir).resolve()
    if source.is_relative_to(out):
        parser.error("Keep the source table outside the dedicated generated output directory")
    data = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    if any(column.startswith("DOWNSTREAM_") for column in data):
        parser.error("Use the original ranked table rather than re-filtering an already annotated export")
    if "SVTYPE" not in data or not any(c in data for c in ["SV_ID", "ID"]) or not any(c in data for c in ["GENE", "GENES", "Gene", "ANNotsv_Gene"]):
        parser.error("Expected an SV-gene table containing SVTYPE, SV_ID/ID and GENE/GENES")
    annotations = []
    for i, row in data.iterrows():
        annotations.append(annotate_row(
            row, args.sample, i+1, args.size_cutoff,
            args.max_af, args.inv_min_reads, args.inv_strong_reads,
        ))
    columns = list(annotate_row({}, args.sample, 0, args.size_cutoff, args.max_af, args.inv_min_reads, args.inv_strong_reads))
    work = pd.concat([data, pd.DataFrame(annotations, columns=columns)], axis=1)
    if work.duplicated(["DOWNSTREAM_SV_ID", "DOWNSTREAM_GENE"]).any():
        parser.error("Repeated SV-gene keys: use the final compact ranked table, not repeated transcript annotations")
    if work["DOWNSTREAM_SV_ID"].eq(".").any():
        parser.error("Every input event must have an SV ID")
    if work["DOWNSTREAM_GENE"].str.contains(r"[;,|]", regex=True).any():
        parser.error("Expected one gene per row: use the final ranked SV-gene table rather than a multi-gene interval report")
    # Population/read evidence belongs to the SV and should agree across genes.
    shared = ["DOWNSTREAM_SIZE_BP", "DOWNSTREAM_CHROMOSOME_SCOPE", "DOWNSTREAM_MAX_AF", "DOWNSTREAM_AF_STATUS", "DOWNSTREAM_CALLER_COUNT", "DOWNSTREAM_MAX_READ_SUPPORT"]
    conflicts = work.groupby("DOWNSTREAM_SV_ID")[shared].nunique().gt(1).any(axis=1)
    for idx in work.index[work["DOWNSTREAM_SV_ID"].isin(conflicts.index[conflicts])]:
        previous = work.at[idx, "DOWNSTREAM_REVIEW_REASONS"]
        work.at[idx, "DOWNSTREAM_REVIEW_REASONS"] = (previous + ";" if previous != "." else "") + "EVENT_WIDE_EVIDENCE_CONFLICT"
    work["DOWNSTREAM_DECISION"] = "RETAINED"
    work.loc[work["DOWNSTREAM_REVIEW_REASONS"].ne("."), "DOWNSTREAM_DECISION"] = "REVIEW_REQUIRED"
    work.loc[work["DOWNSTREAM_EXCLUSION_REASONS"].ne("."), "DOWNSTREAM_DECISION"] = "EXCLUDED"
    out.mkdir(parents=True, exist_ok=True)
    # Remove only files owned by a previous run of this script, never notes.
    previous = out / "filter_run.json"
    if previous.exists():
        old = json.loads(previous.read_text())
        if old.get("script") != Path(__file__).name:
            parser.error("This output directory belongs to a different generator")
        for filename in old.get("generated_files", []):
            path = (out / filename).resolve()
            if not path.is_relative_to(out) or path.suffix != ".tsv":
                parser.error("Invalid generated-file path in the previous run manifest")
        for filename in old.get("generated_files", []):
            (out / filename).unlink(missing_ok=True)
    files = save_views(work, out, args.size_cutoff)
    settings = {
        "script": Path(__file__).name,
        "sample": args.sample,
        "source": str(source),
        "size_cutoff_bp": args.size_cutoff,
        "main_size_comparison": ">",
        "max_population_af": args.max_af,
        "inv_min_reads": args.inv_min_reads,
        "inv_strong_reads": args.inv_strong_reads,
        "read_support_scope": "maximum reported count from one caller; not a sum or unique-read recount",
        "original_ranks": "preserved; retained.tsv keeps source order; caller_priority.tsv groups callers first",
        "generated_files": files,
    }
    (out / "filter_run.json").write_text(json.dumps(settings, indent=2) + "\n")
    (out / "START_HERE.txt").write_text(f"""DOWNSTREAM REVIEW: {args.sample}

Read filter_run.json for the exact thresholds and output_files.tsv for files.
Start with the gt_* / retained.tsv table: size strictly >{args.size_cutoff} bp,
annotated exon/splice or intron context, AF <={args.max_af:g} or AF unknown,
and INV maximum reported caller support >={args.inv_min_reads} reads.

AF_UNKNOWN is not demonstrated ultra-rarity. Inspect by_af/ separately.
lt_* contains small candidates passing the other filters; do not discard
small exonic/splice events. eq_* holds events exactly at the size cutoff.
review_required.tsv holds unresolved annotation/support/length or conflicts.
length_unresolved/ includes BND/TRA, which have no single interval length.
excluded/all_excluded.tsv records all failed criteria, including Y and MT.
all_records.annotated.tsv preserves every input row and original rank column.

retained.tsv follows original input order. caller_priority.tsv shows higher
caller counts first, preserving input order within each caller group.
Files by type, context, AF and callers are overlapping views, not counts to add.
Exon/splice/intron labels describe annotation, not exact breakpoint positions.
Genes spanned intact by an inversion are not automatically disrupted.
Multiple callers and 3/5-read tiers do not constitute independent validation.
Nuclear mitochondrial genes remain; only mitochondrial-chromosome events go.

Input: {source}
Full usage and scientific limits: docs/DOWNSTREAM_SV_FILTERING.md
""", encoding="utf-8")
    print(f"[OK] input associations={len(work)}; start at {out / 'START_HERE.txt'}")
    print(work["DOWNSTREAM_DECISION"].value_counts().to_string())


if __name__ == "__main__":
    main()
