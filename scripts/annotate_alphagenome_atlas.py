#!/usr/bin/env python3
"""AlphaGenome research context: offline AVI lookup and selected native plots.

The caller VCF is the source of alleles. VEP describes their consequences;
phasing describes which chromosome copy carries them. Neither supplies AVI.
This module never changes the HON, Exomiser or clinical rankings.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "plots"))

FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "ALT_INDEX", "GT", "GQ", "DP",
          "FILTER", "ALT_PRESENT", "QC_FLAGS", "OVERLAPPING_GENES", "PANEL_OVERLAP",
          "AVI_STATUS", "AVI_RAW", "AVI_PHRED", "ATLAS_LINK"]
SV_FIELDS = ["CHROM", "POS", "END", "SV_ID", "REF", "ALT", "SVTYPE", "GT", "FILTER",
             "GENOTYPE_SCOPE", "BREAKPOINT_FLAGS", "AVI_STATUS", "MODEL_REVIEW_NOTE"]
REVIEW_FIELDS = FIELDS + ["VEP_GENES", "VEP_CONSEQUENCES", "VEP_CONTEXT_STATUS"]


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def canonical(chrom):
    name = str(chrom).removeprefix("chr")
    return "chrM" if name in {"M", "MT"} else "chr" + name


def contig_name(chrom, available):
    """Allow chr1/1 aliases, but never silently choose between ambiguous contigs."""
    matches = [c for c in available if canonical(c) == canonical(chrom)]
    if len(matches) > 1:
        raise ValueError(f"Ambiguous chromosome names for {chrom}: {matches}")
    return matches[0] if matches else None


def file_info(path):
    stat = Path(path).stat()
    return {"path": str(Path(path).resolve()), "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns}


def write_table(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore", restval=".")
        writer.writeheader()
        writer.writerows(rows)


def parse_avi(line):
    # This is the official Ensembl AVI.pm layout, with 1-based positions.
    parts = line.rstrip().split("\t")
    if len(parts) != 6:
        raise ValueError("Expected six AVI columns: chrom, pos, ref, alt, raw_score, phred. "
                         "A splicing/features download is not an AVI file.")
    chrom, pos, ref, alt, raw, phred = parts
    pos, raw, phred = int(pos), float(raw), float(phred)
    if pos < 1 or ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1 or ref == alt:
        raise ValueError("AVI records must describe a 1-based, biallelic A/C/G/T SNV.")
    if not math.isfinite(raw) or not math.isfinite(phred) or phred < 0:
        raise ValueError("AVI scores must be finite; PHRED must be >=0. Negative RAW scores are valid.")
    return (canonical(chrom), pos, ref, alt), (raw, phred)


class AviLookup:
    """Read small tabix windows, rather than loading the genome-wide download."""
    def __init__(self, path, bin_size=1000):
        self.tabix = pysam.TabixFile(str(path))
        self.bin_size, self.window, self.records = bin_size, None, {}
        self.chromosomes = {}
        if not self.tabix.contigs:
            raise ValueError("AVI file/index has no chromosomes.")
        for line in self.tabix.header:
            fields = line.lstrip("#").lower().split("\t")
            if len(fields) >= 3 and fields[2] in {"ref", "reference"}:
                if len(fields) != 6 or fields[-2] not in {"raw_score", "avi_raw", "raw"} or fields[-1] not in {"phred", "avi_phred"}:
                    raise ValueError("Header does not match the official AVI six-column layout.")
        first = next(self.tabix.fetch(self.tabix.contigs[0]), None)
        if first is None:
            raise ValueError("AVI download is empty.")
        key, _ = parse_avi(first)
        # Detect an index made with zero-based coordinates for a 1-based file.
        if first not in list(self.tabix.fetch(self.tabix.contigs[0], key[1] - 1, key[1])):
            raise ValueError("AVI tabix index is not consistent with 1-based positions.")

    def get(self, chrom, pos, ref, alt):
        if chrom not in self.chromosomes:
            self.chromosomes[chrom] = contig_name(chrom, self.tabix.contigs)
        contig = self.chromosomes[chrom]
        if not contig:
            return "CONTIG_NOT_IN_ATLAS", None
        start = ((pos - 1) // self.bin_size) * self.bin_size
        window = (contig, start)
        if window != self.window:
            self.window, self.records = window, {}
            for line in self.tabix.fetch(contig, start, start + self.bin_size):
                key, score = parse_avi(line)
                if key in self.records and self.records[key] != score:
                    raise ValueError(f"Conflicting AVI records for {key}.")
                self.records[key] = score
        score = self.records.get((canonical(chrom), pos, ref, alt))
        return ("SCORED", score) if score is not None else ("ALLELE_NOT_IN_ATLAS", None)


class GeneIntervals:
    """BED overlap is genomic context, not a demonstrated regulatory target."""
    def __init__(self, path):
        self.bins = defaultdict(list)
        with open_text(path) as handle:
            for line in handle:
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                chrom, start, end, gene, *_ = line.rstrip().split("\t")
                start, end = int(start), int(end)
                if start < 0 or end <= start:
                    raise ValueError("Gene BED must use 0-based, half-open intervals.")
                for block in range(start // 100000, (end - 1) // 100000 + 1):
                    self.bins[canonical(chrom), block].append((start, end, gene))

    def at(self, chrom, pos):
        return sorted({gene for start, end, gene in self.bins[canonical(chrom), (pos - 1) // 100000]
                       if start <= pos - 1 < end})


def sample_call(record, sample, allow_caller_samples=False):
    names = list(record.samples)
    if sample in names:
        name = sample
    elif len(names) == 1:
        name = names[0]  # Single-sample caller files may use a different label.
    elif allow_caller_samples:
        # Merged SV VCFs can have one column per caller. Do not choose one as
        # the patient's reliable genotype. SV prediction needs manual review.
        return {}, (), "."
    else:
        raise ValueError(f"VCF must contain {sample} or exactly one sample; found {names}.")
    sample_record = record.samples[name]
    call = dict(sample_record)
    gt = call.get("GT") or ()
    text = ("|" if sample_record.phased else "/").join("." if a is None else str(a) for a in gt) or "."
    return call, gt, text


def allele_row(record, sample, alt_index, genes, panel, fasta, min_gq, min_dp, reference_contigs=None):
    call, gt, gt_text = sample_call(record, sample)
    ref, alt = record.ref.upper(), record.alts[alt_index - 1].upper()
    gq, dp = call.get("GQ"), call.get("DP")
    filters = list(record.filter.keys())
    present = "YES" if alt_index in gt else "UNKNOWN" if not gt or None in gt else "NO"
    flags = []
    if filters != ["PASS"]:
        flags.append("FILTER_NOT_PASS")
    if present != "YES":
        flags.append("ALT_PRESENCE_" + present)
    for value, threshold, label in [(gq, min_gq, "GQ"), (dp, min_dp, "DP")]:
        if value is None or not math.isfinite(float(value)):
            flags.append(label + "_MISSING")
        elif value < threshold:
            flags.append(label + "_LOW")
    fchrom = reference_contigs[record.contig] if reference_contigs is not None else contig_name(record.contig, fasta.references)
    if not fchrom:
        status = "CONTIG_NOT_IN_REFERENCE"
    elif not re.fullmatch("[ACGT]+", ref) or record.start + len(ref) > fasta.get_reference_length(fchrom):
        status = "INVALID_REFERENCE_ALLELE"
    elif fasta.fetch(fchrom, record.start, record.start + len(ref)).upper() != ref:
        status = "REF_MISMATCH"
    elif len(ref) != 1 or len(alt) != 1 or alt not in "ACGT" or ref == alt:
        status = "NOT_SUPPORTED_BY_OFFLINE_SNV_ATLAS"
    else:
        status = "LOOKUP"
    overlap = genes.at(record.contig, record.pos)
    return {"CHROM": record.contig, "POS": record.pos, "ID": record.id or ".", "REF": ref,
            "ALT": alt, "ALT_INDEX": alt_index, "GT": gt_text, "GQ": gq if gq is not None else ".",
            "DP": dp if dp is not None else ".", "FILTER": ";".join(filters) or ".",
            "ALT_PRESENT": present, "QC_FLAGS": ";".join(flags) or "NONE",
            "OVERLAPPING_GENES": ";".join(overlap) or ".",
            "PANEL_OVERLAP": "YES" if set(overlap) & panel else "NO", "AVI_STATUS": status,
            "AVI_RAW": ".", "AVI_PHRED": ".", "ATLAS_LINK": "."}


def add_vep_context(rows, path):
    """Join SNV consequences only when VEP identifies the exact REF/ALT allele."""
    indexed = defaultdict(list)
    for row in rows:
        key = (canonical(row["CHROM"]), int(row["POS"]), row["REF"], row["ALT"])
        indexed[key].append(row)
        row.update(VEP_GENES=set(), VEP_CONSEQUENCES=set(), VEP_CONTEXT_STATUS="EXACT_ALLELE_NOT_FOUND")
    header = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#Uploaded_variation"):
                header = line.lstrip("#").rstrip().split("\t")
            elif not line.startswith("#") and header and line.strip():
                data = dict(zip(header, line.rstrip().split("\t")))
                location = data.get("Location", "")
                if ":" not in location:
                    continue
                chrom, pos = location.rsplit(":", 1)
                if "-" in pos:
                    endpoints = pos.split("-")
                    pos = endpoints[0] if len(endpoints) == 2 and endpoints[0] == endpoints[1] else ""
                if not pos.isdigit():  # Only single-base loci can match offline SNV AVI.
                    continue
                extra = dict(item.split("=", 1) for item in data.get("Extra", "").split(";") if "=" in item)
                uploaded = extra.get("UPLOADED_ALLELE", "")
                ref = uploaded.split("/", 1)[0] if "/" in uploaded else ""
                if len(ref) != 1:
                    match = re.match(r"^.+_\d+_([ACGT])/(.+)$", data.get("Uploaded_variation", ""))
                    ref = match.group(1) if match else ""
                for row in indexed.get((canonical(chrom), int(pos), ref, data.get("Allele", "")), []):
                    symbol = extra.get("SYMBOL", data.get("Gene", "."))
                    if symbol not in {"", ".", "-"}:
                        row["VEP_GENES"].add(symbol)
                    row["VEP_CONSEQUENCES"].update(data.get("Consequence", "").split(","))
                    row["VEP_CONTEXT_STATUS"] = "EXACT_SNV_ALLELE_MATCH"
    if header is None:
        raise ValueError("Expected a VEP TSV with #Uploaded_variation header.")
    for row in rows:
        for field in ["VEP_GENES", "VEP_CONSEQUENCES"]:
            row[field] = ";".join(sorted(row[field] - {"", ".", "-"})) or "."


def sv_inventory(path, sample):
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            _, _, gt = sample_call(rec, sample, allow_caller_samples=True)
            kind = str(dict(rec.info).get("SVTYPE", "OTHER"))
            flags = [f for f in ("IMPRECISE", "CIPOS", "CIEND") if f in rec.info]
            endpoint = explicit_end(rec)
            if endpoint is None:
                flags.append("EXPLICIT_END_MISSING")
            elif rec.stop != endpoint:
                flags.append("END_SVLEN_CONFLICT")
            for alt in rec.alts or ():
                yield {"CHROM": rec.contig, "POS": rec.pos, "END": endpoint if endpoint is not None else ".", "SV_ID": rec.id or ".",
                       "REF": rec.ref, "ALT": alt, "SVTYPE": kind, "GT": gt,
                       "FILTER": ";".join(rec.filter.keys()) or ".",
                       "GENOTYPE_SCOPE": "MERGED_CALLER_COLUMNS_UNRESOLVED" if len(rec.samples) > 1 else "CALLER_REPORTED_REVIEW_REQUIRED",
                       "BREAKPOINT_FLAGS": ";".join(flags) or "PRECISION_NOT_ESTABLISHED",
                       "AVI_STATUS": "SV_NOT_SUPPORTED_BY_OFFLINE_SNV_ATLAS",
                       "MODEL_REVIEW_NOTE": "Reviewed DEL allele may be explored within model context; no validated SV impact score."
                       if kind == "DEL" else "Needs a resolved altered sequence; DUP orientation/insertion site and BND/complex junctions cannot be inferred from a type label."}


def explicit_end(record):
    # htslib may derive record.stop from SVLEN even when an explicit END
    # disagrees. Read the serialized INFO END and compare both descriptions.
    info_text = str(record).split("\t")[7]
    info = dict(item.split("=", 1) for item in info_text.split(";") if "=" in item)
    return int(info["END"]) if "END" in info else None


def annotate(args):
    if args.assembly != "GRCh38" or not args.release.strip():
        raise ValueError("Atlas requires GRCh38 and an explicit resource release/filename label.")
    lookup, genes = AviLookup(args.avi, args.bin_size), GeneIntervals(args.gene_bed)
    with open(args.panel) as handle:
        panel = {line.strip().split("\t")[0] for line in handle if line.strip() and not line.startswith("#")}
    counts, histogram, review = Counter(), Counter(), []
    with pysam.FastaFile(args.reference) as fasta, pysam.VariantFile(args.vcf) as vcf:
        reference_contigs = {}
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open_text(args.output, "wt") as handle:
            writer = csv.DictWriter(handle, FIELDS, delimiter="\t")
            writer.writeheader()
            for rec in vcf:
                if rec.contig not in reference_contigs:
                    reference_contigs[rec.contig] = contig_name(rec.contig, fasta.references)
                for index in range(1, len(rec.alts or ()) + 1):
                    row = allele_row(rec, args.sample, index, genes, panel, fasta, args.min_gq, args.min_dp, reference_contigs)
                    if row["AVI_STATUS"] == "LOOKUP":
                        row["AVI_STATUS"], score = lookup.get(rec.contig, rec.pos, row["REF"], row["ALT"])
                        if score is not None:
                            row["AVI_RAW"], row["AVI_PHRED"] = score
                            row["ATLAS_LINK"] = ("https://deepmind.google.com/science/alphagenome/atlas?q="
                                                 f"{canonical(rec.contig)}:{rec.pos}:{row['REF']}%3E{row['ALT']}&m=variant")
                            histogram[min(80, int(score[1]))] += 1
                            # This is a molecular-impact view, including poor-QC rows with visible flags.
                            # The complete table is never restricted by this threshold.
                            if score[1] >= args.review_phred:
                                review.append(row.copy())
                    counts[row["AVI_STATUS"]] += 1
                    writer.writerow(row)
    review.sort(key=lambda row: (-float(row["AVI_PHRED"]), canonical(row["CHROM"]), int(row["POS"]), row["ALT"]))
    add_vep_context(review, args.vep)
    write_table(args.review_output, review, REVIEW_FIELDS)
    write_table(args.sv_output, sv_inventory(args.sv_vcf, args.sample), SV_FIELDS)
    provenance = {"sample": args.sample, "assembly": args.assembly, "atlas_release": args.release,
                  "source_vcf": file_info(args.vcf), "avi_resource": file_info(args.avi),
                  "avi_index": file_info(args.avi + ".tbi"), "reference": file_info(args.reference),
                  "vep_context": file_info(args.vep), "sv_source": file_info(args.sv_vcf),
                  "dna_sample_tissue": args.dna_tissue, "review_phred_threshold": args.review_phred,
                  "min_gq_for_prediction": args.min_gq, "min_dp_for_prediction": args.min_dp,
                  "alleles_by_status": dict(counts), "phred_histogram": dict(histogram),
                  "review_rows": len(review), "scope": "RESEARCH_ONLY; ranks and classifications unchanged",
                  "notes": ["No panel, coding or consequence filter on the caller alleles.",
                            "AVI PHRED is a genome-wide substitution percentile, not disease probability.",
                            "AVI is not tissue-specific and combines correlated predictors; do not add separate bonuses.",
                            "No population rarity, regulatory target, inheritance or patient finding inferred.",
                            "File sizes/mtime identify local resources; complete multi-GB files are not hashed."]}
    Path(args.summary).write_text(json.dumps(provenance, indent=2) + "\n")
    lookup.tabix.close()


def summary_plot(args):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from plot_utils import save_figure, set_thesis_style
    data = json.loads(Path(args.summary).read_text())
    set_thesis_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    histogram = {int(key): value for key, value in data["phred_histogram"].items()}
    axes[0].bar(list(histogram), list(histogram.values()), color="#0072B2")
    axes[0].axvline(data["review_phred_threshold"], color="#D55E00", linestyle="--", label="Review-view threshold")
    axes[0].set(xlabel="AVI PHRED (80 includes higher scores)", ylabel="Number of caller ALT alleles", title="Predicted molecular impact; not diagnosis")
    axes[0].legend()
    statuses = data["alleles_by_status"]
    labels = {"SCORED": "SNV scored", "ALLELE_NOT_IN_ATLAS": "SNV allele absent from Atlas",
              "CONTIG_NOT_IN_ATLAS": "Chromosome absent from Atlas",
              "CONTIG_NOT_IN_REFERENCE": "Chromosome absent from FASTA", "REF_MISMATCH": "REF differs from FASTA",
              "INVALID_REFERENCE_ALLELE": "REF cannot be validated", "NOT_SUPPORTED_BY_OFFLINE_SNV_ATLAS": "Indel/other: AVI unavailable"}
    axes[1].barh([labels.get(k, k) for k in statuses], list(statuses.values()), color=["#0072B2" if k == "SCORED" else "#999999" for k in statuses])
    axes[1].set(xlabel="Number of caller ALT alleles", title="Scored and unavailable alleles")
    axes[1].tick_params(axis="y", labelsize=8)
    axes[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.suptitle(f"{data['sample']} — AlphaGenome Atlas research overview")
    fig.tight_layout()
    save_figure(fig, Path(args.outdir) / "atlas_overview", dpi=200)
    plt.close(fig)
    write_table(Path(args.outdir) / "atlas_overview.source.tsv",
                [{"measure": "PHRED_BIN", "category": k, "count": v} for k, v in sorted(histogram.items())]
                + [{"measure": "ANNOTATION_STATUS", "category": k, "count": v} for k, v in statuses.items()],
                ["measure", "category", "count"])


def selected_allele(rec, alt_index, selection, fasta, length, min_flank):
    """Resolve a caller allele; a symbolic DEL needs explicitly reviewed endpoints."""
    chrom = contig_name(rec.contig, fasta.references)
    if not chrom:
        raise ValueError("Chromosome not in reference")
    ref, alt, label = rec.ref.upper(), rec.alts[alt_index - 1].upper(), "SMALL_VARIANT_RESEARCH"
    if selection["SOURCE"] == "SV":
        info = dict(rec.info)
        if selection.get("CONFIRMED_BREAKPOINTS", "").upper() != "YES" or selection.get("CONFIRMED_ALT_PRESENT", "").upper() != "YES":
            raise ValueError("SV exploration requires CONFIRMED_BREAKPOINTS=YES and CONFIRMED_ALT_PRESENT=YES after manual review")
        if info.get("IMPRECISE", False) or any(any(v is None or int(v) != 0 for v in info.get(field, ())) for field in ("CIPOS", "CIEND")):
            raise ValueError("Imprecise/confidence-interval breakpoints cannot define an exact deletion allele")
        if str(info.get("SVTYPE", "")) != "DEL" or alt != "<DEL>":
            raise ValueError("Only reviewed symbolic DEL alleles are reconstructed; other SV types require resolved sequence")
        # VCF POS is the retained anchor. END is the last deleted base (1-based inclusive).
        endpoint = explicit_end(rec)
        if endpoint is None or endpoint <= rec.pos or endpoint - rec.start > length - 2 * min_flank:
            raise ValueError("DEL endpoints are missing/inconsistent or exceed model context with flanks")
        if rec.stop != endpoint:
            raise ValueError("DEL END and SVLEN disagree; resolve the allele before modelling")
        svlen = info.get("SVLEN")
        if svlen is not None:
            value = svlen[alt_index - 1] if isinstance(svlen, tuple) else svlen
            if abs(int(value)) != endpoint - rec.pos:
                raise ValueError("DEL END and SVLEN disagree; resolve the allele before modelling")
        ref = fasta.fetch(chrom, rec.start, endpoint).upper()
        alt, label = ref[0], "EXPLORATORY_DEL_UNVALIDATED_NO_AVI"
        if rec.ref.upper() not in {ref[0], "N"}:
            raise ValueError("DEL anchor does not match reference")
    if not re.fullmatch("[ACGT]+", ref) or not re.fullmatch("[ACGT]+", alt) or ref == alt:
        raise ValueError("Prediction requires resolved, different A/C/G/T alleles")
    if fasta.fetch(chrom, rec.start, rec.start + len(ref)).upper() != ref:
        raise ValueError("REF does not match FASTA")
    if label == "SMALL_VARIANT_RESEARCH" and max(len(ref), len(alt)) >= 50:
        label = "EXPLORATORY_SEQUENCE_SV_UNVALIDATED_NO_AVI"
    # The entire REF allele and both flanks must be inside one supported window.
    start = rec.start + len(ref) // 2 - length // 2
    end = start + length
    if start < 0 or end > fasta.get_reference_length(chrom):
        raise ValueError("Window crosses a contig boundary; choose a shorter supported context")
    if rec.start - start < min_flank or end - (rec.start + len(ref)) < min_flank:
        raise ValueError("Allele lacks required flanking context")
    return canonical(chrom), rec.pos, ref, alt, start, end, label


def native_plots(args):
    """Only manually selected alleles go to the API; saved tracks can be replotted."""
    from alphagenome.data import genome, track_data
    from alphagenome.models import dna_client
    from alphagenome.visualization import plot_components
    import importlib.metadata
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from plot_utils import save_figure, set_thesis_style
    import textwrap

    set_thesis_style()  # Includes tight export bounds for native track labels.

    terms = json.loads(args.ontology_terms_json)
    if not isinstance(terms, list) or not terms or not all(isinstance(term, str) and ":" in term for term in terms):
        raise ValueError("Choose actual model ontology IDs; DNA sample tissue is not the prediction tissue")
    with open(args.selection) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"SOURCE", "CHROM", "POS", "REF", "ALT", "SV_ID"}.issubset(reader.fieldnames or []):
            raise ValueError("Selection needs SOURCE, CHROM, POS, REF, ALT, SV_ID columns")
        selections = list(reader)
    if len(selections) > args.max_variants:
        raise ValueError(f"Selection exceeds max_variants={args.max_variants}; this is not a genome-wide API run")
    if any(s.get("SOURCE") not in {"SMALL", "SV"} for s in selections):
        raise ValueError("SOURCE must be SMALL or SV")
    wanted = defaultdict(list)
    for ix, s in enumerate(selections):
        wanted[s["SOURCE"], canonical(s["CHROM"]), int(s["POS"]), s["REF"].upper(), s["ALT"].upper()].append(ix)
    found = {}
    # 1. Check every selected allele against its original VCF and GRCh38 FASTA.
    with pysam.FastaFile(args.reference) as fasta:
        for source, path in [("SMALL", args.vcf), ("SV", args.sv_vcf)]:
            with pysam.VariantFile(str(path)) as vcf:
                for rec in vcf:
                    for alt_index, alt in enumerate(rec.alts or (), 1):
                        key = (source, canonical(rec.contig), rec.pos, rec.ref.upper(), alt.upper())
                        for ix in wanted.get(key, []):
                            s = selections[ix]
                            if source == "SV" and (s["SV_ID"] in {"", "."} or s["SV_ID"] != rec.id):
                                continue
                            if ix in found:
                                raise ValueError("Selection matches multiple caller records; resolve duplicate alleles")
                            call, gt, gt_text = sample_call(rec, args.sample, allow_caller_samples=source == "SV")
                            if list(rec.filter.keys()) != ["PASS"] or (gt and None not in gt and alt_index not in gt):
                                raise ValueError("Selected allele has failed/unknown FILTER or an explicit non-carrier genotype")
                            if source == "SMALL" and (alt_index not in gt or any(call.get(field) is None or not math.isfinite(float(call.get(field))) or call.get(field) < threshold for field, threshold in [("GQ", args.min_gq), ("DP", args.min_dp)])):
                                raise ValueError("Small-variant prediction needs observed ALT and adequate GQ/DP; the master table still retains this call")
                            found[ix] = (selected_allele(rec, alt_index, s, fasta, args.sequence_length, args.min_flank), gt_text)
    if len(found) != len(selections):
        raise ValueError("At least one selection does not exactly match its source caller VCF")
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    client, report = None, []
    sdk_version = importlib.metadata.version("alphagenome")
    for ix, selection in enumerate(selections):
        (chrom, pos, ref, alt, start, end, label), gt = found[ix]
        request = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "start": start, "end": end,
                   "ontology_terms": terms, "outputs": ["RNA_SEQ", "SPLICE_SITES"],
                   "sdk_version": sdk_version, "model_version": args.model_version,
                   "reference": file_info(args.reference), "label": label}
        digest = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()[:20]
        folder = Path(args.outdir) / digest
        folder.mkdir(exist_ok=True)
        cache_folder = Path(args.cache_dir) / digest
        cache_folder.mkdir(parents=True, exist_ok=True)
        cache = cache_folder / "tracks.npz"
        interval = genome.Interval(chrom, start, end)
        variant = genome.Variant(chrom, pos, ref, alt)
        # 2. Reuse a complete cache, or ask the model for this selected allele.
        complete_cache = all((cache_folder / name).is_file() for name in ("tracks.npz", "request.json", "track_metadata.json"))
        if not complete_cache:
            if args.cache_only:
                raise ValueError(f"No cached tracks for selected allele {digest}")
            if client is None:
                key = os.environ.get("ALPHAGENOME_API_KEY")
                if not key:
                    raise ValueError("Export ALPHAGENOME_API_KEY on the server for selected native plots; never put it in YAML")
                version = dna_client.ModelVersion[args.model_version]
                client = dna_client.create(key, model_version=version, timeout=30)
            prediction = client.predict_variant(interval, variant,
                requested_outputs=[dna_client.OutputType.RNA_SEQ, dna_client.OutputType.SPLICE_SITES],
                ontology_terms=terms)
            arrays, metadata = {}, {}
            for name in ("rna_seq", "splice_sites"):
                reference, alternate = getattr(prediction.reference, name), getattr(prediction.alternate, name)
                if reference is None or alternate is None:
                    continue
                # Skip unavailable requested tissues; never substitute an unrequested tissue.
                if reference.values.shape[1] == 0:
                    continue
                arrays[name + "_ref"], arrays[name + "_alt"] = reference.values, alternate.values
                metadata[name] = {"tracks": json.loads(reference.metadata.to_json(orient="table", index=False)),
                                  "resolution": reference.resolution,
                                  "interval": {"chromosome": reference.interval.chromosome,
                                               "start": reference.interval.start, "end": reference.interval.end}}
            if not arrays:
                raise ValueError("No requested tissue tracks are available; inspect AlphaGenome output metadata")
            np.savez_compressed(cache, **arrays)
            (cache_folder / "request.json").write_text(json.dumps(request, indent=2) + "\n")
            (cache_folder / "track_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        saved = json.loads((cache_folder / "request.json").read_text())
        if saved != request:
            raise ValueError("Prediction cache provenance differs from request")
        metadata = json.loads((cache_folder / "track_metadata.json").read_text())
        # 3. Rebuild official TrackData and export the tool's REF/ALT plots.
        with np.load(cache, allow_pickle=False) as arrays:
            for name, meta in metadata.items():
                tracks = pd.DataFrame(meta["tracks"]["data"])
                tracks.to_csv(folder / f"{name}.tracks.tsv", sep="\t", index=False)
                track_interval = genome.Interval(**meta["interval"])
                zoom_width = min(args.sequence_length, max(32768, len(ref) + 2 * args.min_flank))
                zoom_start = max(start, pos - 1 + len(ref) // 2 - zoom_width // 2)
                zoom = genome.Interval(chrom, zoom_start, min(end, zoom_start + zoom_width))
                # Four tracks per page keep tissues/strands readable and avoid
                # the native plotting library's maximum-track limit.
                for first in range(0, len(tracks), 4):
                    subset = tracks.iloc[first:first + 4].reset_index(drop=True)
                    # Keep full source names in the TSV; use readable tissue
                    # names and wrapped labels on the figure itself.
                    if "biosample_name" in subset:
                        subset["name"] = subset["biosample_name"].fillna(subset["name"])
                    subset["name"] = [textwrap.fill(str(value), width=32) for value in subset["name"]]
                    ref_data = track_data.TrackData(values=arrays[name + "_ref"][:, first:first + 4], metadata=subset, interval=track_interval, resolution=meta["resolution"])
                    alt_data = track_data.TrackData(values=arrays[name + "_alt"][:, first:first + 4], metadata=subset, interval=track_interval, resolution=meta["resolution"])
                    display_label = "Exploratory SV; unvalidated" if label.startswith("EXPLORATORY") else "Research prediction"
                    fig = plot_components.plot([plot_components.OverlaidTracks(
                        tdata={"REF": ref_data, "ALT": alt_data}, colors={"REF": "#666666", "ALT": "#CC79A7"})],
                        interval=zoom, annotations=[plot_components.VariantAnnotation([variant], alpha=0.15)], fig_width=12,
                        xlabel="GRCh38 position (0-based); VCF POS is 1-based",
                        title=f"{args.sample} — {chrom}:{pos} — {name.replace('_', ' ')}\n{display_label}; single-allele molecular signal")
                    if label.startswith("EXPLORATORY") and len(ref) > len(alt):
                        fig.text(0.5, -0.03, "Deleted positions have zero aligned ALT signal; this alone is not predicted whole-gene expression loss.", ha="center", fontsize=9)
                    prefix = name if len(tracks) <= 4 else f"{name}_page_{first // 4 + 1:02d}"
                    save_figure(fig, folder / prefix, dpi=200)
                    plt.close(fig)
        report.append({"SOURCE": selection["SOURCE"], "CHROM": chrom, "POS": pos,
                       "SV_ID": selection["SV_ID"], "GT": gt, "LABEL": label,
                       "REQUESTED_TISSUES": ";".join(terms), "AVAILABLE_OUTPUTS": ";".join(metadata),
                       "OUTPUT_NOTES": "RNA_SEQ_UNAVAILABLE_FOR_REQUESTED_TISSUES" if "rna_seq" not in metadata else "Inspect track metadata for actual tissues; splice-site class tracks may be tissue-independent.",
                       "PLOT_FOLDER": str(folder), "CACHE": str(cache)})
    write_table(args.manifest, report, ["SOURCE", "CHROM", "POS", "SV_ID", "GT", "LABEL", "REQUESTED_TISSUES", "AVAILABLE_OUTPUTS", "OUTPUT_NOTES", "PLOT_FOLDER", "CACHE"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    ann = sub.add_parser("annotate")
    for option in ("vcf", "sv-vcf", "reference", "avi", "gene-bed", "panel", "vep", "sample", "release", "output", "review-output", "sv-output", "summary"):
        ann.add_argument("--" + option, required=True)
    ann.add_argument("--assembly", default="GRCh38")
    ann.add_argument("--dna-tissue", default="UNKNOWN")
    ann.add_argument("--bin-size", type=int, default=1000)
    ann.add_argument("--review-phred", type=float, default=20)
    ann.add_argument("--min-gq", type=float, default=20)
    ann.add_argument("--min-dp", type=int, default=10)
    plot = sub.add_parser("plot")
    plot.add_argument("--summary", required=True)
    plot.add_argument("--outdir", required=True)
    native = sub.add_parser("native-plots")
    for option in ("vcf", "sv-vcf", "reference", "sample", "selection", "outdir", "cache-dir", "manifest", "ontology-terms-json"):
        native.add_argument("--" + option, required=True)
    native.add_argument("--sequence-length", type=int, choices=[16384, 131072, 524288, 1048576], default=1048576)
    native.add_argument("--min-flank", type=int, default=10000)
    native.add_argument("--max-variants", type=int, default=10)
    native.add_argument("--model-version", default="ALL_FOLDS")
    native.add_argument("--min-gq", type=float, default=20)
    native.add_argument("--min-dp", type=int, default=10)
    native.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()
    if args.command != "plot" and (not math.isfinite(args.min_gq) or args.min_gq < 0 or args.min_dp < 0):
        parser.error("Quality thresholds must be finite and nonnegative")
    if args.command == "annotate" and (args.bin_size < 1 or not math.isfinite(args.review_phred) or args.review_phred < 0):
        parser.error("bin-size must be positive and review-phred finite and >=0")
    if args.command == "native-plots" and (args.min_flank < 0 or args.max_variants < 1):
        parser.error("min-flank must be >=0 and max-variants positive")
    {"annotate": annotate, "plot": summary_plot, "native-plots": native_plots}[args.command](args)


if __name__ == "__main__":
    main()
