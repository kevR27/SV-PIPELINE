#!/usr/bin/env python3
"""Transparent SV–gene–disease research assessment, never clinical classification.

Preserve all input rows. Produce a separate disease-hypothesis table and append
one explicitly identified best hypothesis to each original SV–gene row. Missing
inputs remain unknown; no clinical label or master-callset filter is generated.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from sv_evidence_common import MISSING, number

# INFO/read-name and transcript JSON fields can exceed the csv default of 128 KiB.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

VERSION = "sv-allele-evidence-v1.2"
DOMAINS = ("POPULATION", "TECHNICAL", "DISRUPTION", "INHERITANCE", "PHENOTYPE", "DISEASE_MECHANISM")
EFFECTS = {"LOF", "COPY_GAIN", "GOF", "DOMINANT_NEGATIVE", "NO_DISRUPTION", "OTHER", "UNKNOWN"}
MOIS = {"AD", "AR", "XLD", "XLR", "MT", "UNKNOWN"}
MECHANISMS = {"LOF", "HAPLOINSUFFICIENCY", "TRIPLOSENSITIVITY", "GOF", "DOMINANT_NEGATIVE", "OTHER", "UNKNOWN"}
VALIDITY = {"DEFINITIVE", "STRONG", "MODERATE", "LIMITED", "DISPUTED", "REFUTED", "UNKNOWN"}


def known(value):
    return value is not None and str(value).strip().upper() not in MISSING


def rows(path, required=()):
    if not path:
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        return [{k: (v or ".").strip() for k, v in r.items()} for r in reader]


def unique_index(records, keys):
    result = {}
    for r in records:
        key = tuple(r.get(k, ".") for k in keys)
        if key in result:
            raise ValueError(f"Duplicate evidence key {key}; resolve the conflict explicitly")
        result[key] = r
    return result


def terms(value):
    found = set(re.split(r"[;,|\s]+", str(value))) - {"", "."}
    if any(not re.fullmatch(r"HP:\d{7}", term) for term in found):
        raise ValueError(f"Invalid HPO term list: {value}")
    return found


def genotype(value):
    if not known(value):
        return "UNKNOWN"
    parts = re.split(r"[/|]", str(value))
    if any(x not in {"0", "1"} for x in parts) or len(parts) not in {1, 2}:
        return "UNKNOWN"
    if all(x == "0" for x in parts):
        return "HOM_REF"
    if len(parts) == 1:
        return "HAPLOID_ALT"
    return "HOM_ALT" if all(x == "1" for x in parts) else "HET"


def domain(status, points=0, detail="."):
    return status, points, detail


def population(row, evidence, rare_af):
    af = number(evidence.get("population_af"))
    if known(evidence.get("population_af")) and af is None:
        raise ValueError("Invalid population_af value")
    if af is not None:
        if not 0 <= af <= 1:
            raise ValueError("population_af must be a frequency in [0, 1]")
        match = evidence.get("population_match", ".")
        source = evidence.get("population_source", ".")
        if match not in {"EXACT_ALLELE", "SEQUENCE_VALIDATED"} or not known(source):
            return domain("UNKNOWN_UNVERIFIED_EXTERNAL_AF", detail="AF requires allele match and population source")
        if af == 0:
            return domain("NOT_OBSERVED_IN_REPORTED_POPULATION", 1, f"AF=0;source={source};absence is not pathogenicity")
        if af <= rare_af:
            return domain("LOW_FREQUENCY_ALLELE_MATCH", 2, f"AF={af};source={source}")
        return domain("ABOVE_RESEARCH_AF_THRESHOLD", 0, f"AF={af};source={source};threshold={rare_af}")
    af = number(row.get("NEEDLR_AF"))
    if row.get("NEEDLR_STATUS") == "MATCHED" and af is not None and 0 <= af <= 1:
        status = "PROVISIONAL_LOW_FREQUENCY_MATCH" if af <= rare_af else "PROVISIONAL_HIGH_FREQUENCY_MATCH"
        return domain(status, 0, f"AF={af};field={row.get('NEEDLR_AF_SOURCE_FIELD', '.')};coordinate-only match;requires allele review")
    return domain("UNKNOWN", detail=f"needLR={row.get('NEEDLR_STATUS', 'NOT_SUPPLIED')}")


def technical(row, min_support):
    data = json.loads(row.get("CALLER_EVIDENCE_JSON", "[]"))
    support = [number(r.get("CALLER_SUPPORT")) for r in data]
    support = [x for x in support if x is not None and x >= 0]
    count = number(row.get("CALLER_COUNT"))
    all_flags = {flag.strip().upper() for r in data if known(r.get("EVIDENCE_FLAGS"))
                 for flag in re.split(r"[;,|]", r["EVIDENCE_FLAGS"])} - {"PASS", "NONE", ""}
    # A translocation/breakend has no meaningful single SV length. Its absence
    # is descriptive, not a technical failure. Other flags still require review.
    informational = {"NO_SVLEN"} if row.get("SVTYPE", "").upper() in {"BND", "TRA"} else set()
    flagged = sorted(all_flags - informational)
    if any(r.get("EVIDENCE_STATUS") == "FAIL" for r in data):
        flagged.append("CALLER_RECORD_FAILS_CURRENT_FILTER")
    filt = row.get("FILTER", ".")
    method = row.get("CALLER_EVIDENCE_MATCH", "NO_MATCH")
    failures = sorted({r.get("EVIDENCE_FAIL_REASONS", "UNSPECIFIED") for r in data if r.get("EVIDENCE_STATUS") == "FAIL"})
    detail = f"callers={count};max_support={max(support) if support else '.'};match={method};flags={';'.join(sorted(all_flags)) or '.'};FILTER={filt};current_filter_failures={';'.join(failures) or '.'}"
    if "AMBIGUOUS" in method or flagged or filt not in {"PASS", ".", ""}:
        return domain("REVIEW_REQUIRED", 0, detail)
    if not support:
        return domain("UNKNOWN_READ_SUPPORT", 0, detail)
    if max(support) < min_support:
        return domain("LOW_READ_SUPPORT", 0, detail)
    if count is not None and count >= 2:
        return domain("MULTICALLER_READ_SUPPORTED", 2, detail + ";shared reads are not independent validation")
    return domain("SINGLE_CALLER_READ_SUPPORTED", 1, detail)


def transcript_breakpoints(row, ann):
    """Use the correct local or remote chromosome for a BND transcript."""
    start, end = number(row.get("START")), number(row.get("END"))
    if row.get("SVTYPE") in {"BND", "TRA"}:
        contig = str(ann.get("SV_chrom", row.get("CHROM", "."))).removeprefix("chr")
        points = [number(row.get(key)) for c, key in (("CHROM", "START"), ("CHR2", "POS2"))
                  if str(row.get(c, ".")).removeprefix("chr") == contig]
    elif row.get("SVTYPE") == "INS":
        points = [start]
    else:
        points = [start, end]
    lo, hi = number(ann.get("Tx_start")), number(ann.get("Tx_end"))
    return [p for p in points if p is not None and lo is not None and hi is not None and lo < p < hi]


def functional_context(row):
    """Describe the predicted effect without turning overlap into disease proof."""
    annotations = json.loads(row.get("ANNOTSV_GENE_ROWS_JSON", "[]"))
    if not annotations:
        return "NO_GENE_SPECIFIC_ANNOTATION"
    kind = row.get("SVTYPE")
    effects = set()
    for ann in annotations:
        cds = number(ann.get("Overlapped_CDS_percent"))
        start, end = number(row.get("START")), number(row.get("END"))
        lo, hi = number(ann.get("Tx_start")), number(ann.get("Tx_end"))
        contains = None not in (start, end, lo, hi) and start <= lo and end >= hi
        if kind == "DEL" and cds is not None and 0 < cds <= 100:
            effects.add("DELETION_OVERLAPS_CDS")
        elif kind == "DUP" and contains:
            effects.add("WHOLE_TRANSCRIPT_COPY_GAIN_PREDICTED")
        elif kind == "INV" and contains:
            effects.add("TRANSCRIPT_INSIDE_INVERSION_NO_INTRAGENIC_BREAKPOINT_SHOWN")
        elif transcript_breakpoints(row, ann):
            location = str(ann.get("Location", "")).lower()
            effects.add("EXONIC_BREAKPOINT_POSSIBLE" if "exon" in location else
                        "INTRONIC_BREAKPOINT_POSSIBLE" if "intron" in location else "TRANSCRIPT_BREAKPOINT_POSSIBLE")
        else:
            effects.add("OVERLAP_WITHOUT_RESOLVED_FUNCTIONAL_EFFECT")
    return ";".join(sorted(effects))


def disruption(row, evidence):
    effect = evidence.get("effect", "UNKNOWN")
    effect = effect.upper() if known(effect) else "UNKNOWN"
    if effect not in EFFECTS:
        raise ValueError(f"Invalid effect {effect}")
    if effect != "UNKNOWN" and known(evidence.get("source")):
        level = evidence.get("effect_status", ".")
        if level not in {"PREDICTED", "VALIDATED"}:
            raise ValueError("Non-unknown effect requires PREDICTED/VALIDATED effect_status")
        return domain(f"{level}_{effect}", 0 if effect == "NO_DISRUPTION" else (2 if level == "VALIDATED" else 1),
                      evidence["source"]), effect
    annotations = json.loads(row.get("ANNOTSV_GENE_ROWS_JSON", "[]"))
    svtype = row.get("SVTYPE")
    candidates = set()
    for ann in annotations:
        cds = number(ann.get("Overlapped_CDS_percent"))
        if svtype == "DEL" and cds is not None and 0 < cds <= 100:
            candidates.add("CDS_LOSS_PREDICTED")
        start, end = number(row.get("START")), number(row.get("END"))
        txstart, txend = number(ann.get("Tx_start")), number(ann.get("Tx_end"))
        if svtype == "DUP" and None not in (start, end, txstart, txend) and start <= txstart and end >= txend:
            candidates.add("TRANSCRIPT_COPY_GAIN_PREDICTED")
        location = str(ann.get("Location", "")).lower()
        if svtype in {"INS", "INV", "BND"} and "exon" in location:
            # Location alone is insufficient for inversions encompassing a gene.
            if transcript_breakpoints(row, ann):
                candidates.add("EXONIC_BREAKPOINT_POSSIBLE")
    detail = ";".join(sorted(candidates)) or functional_context(row)
    if "CDS_LOSS_PREDICTED" in candidates:
        return domain("CDS_LOSS_PREDICTED", 1, detail + ";LOF not established"), "POSSIBLE_LOF"
    if "TRANSCRIPT_COPY_GAIN_PREDICTED" in candidates:
        return domain("TRANSCRIPT_COPY_GAIN_PREDICTED", 1, detail + ";insertion location/function unresolved"), "POSSIBLE_COPY_GAIN"
    if candidates:
        return domain("BREAKPOINT_DISRUPTION_POSSIBLE", 1, detail), "UNKNOWN"
    return domain("UNKNOWN", detail=detail), "UNKNOWN"


def patient_genotype(row, evidence, family, sample, min_gq, min_dp):
    explicit = family.get((sample, row["SV_ID"]))
    if explicit:
        gq, dp = number(explicit.get("gq")), number(explicit.get("dp"))
        gt = genotype(explicit.get("gt"))
        if not known(explicit.get("source")) or gq is None or dp is None or gq < min_gq or dp < min_dp:
            return "UNKNOWN", "FAMILY_GENOTYPE_QC_INSUFFICIENT"
        return gt, f"family_genotypes:{explicit['source']}"
    if known(evidence.get("gt")):
        gq, dp = number(evidence.get("gq")), number(evidence.get("dp"))
        if gq is not None and dp is not None and gq >= min_gq and dp >= min_dp and known(evidence.get("source")):
            return genotype(evidence["gt"]), evidence["source"]
        return "UNKNOWN", "EXTERNAL_GENOTYPE_QC_INSUFFICIENT"
    # Merged VCF GT is never treated as a reliable patient genotype.
    if row.get("CALLER_EVIDENCE_MATCH") != "IDLIST":
        return "UNKNOWN", "NO_EXACT_CALLER_ID_LINK"
    observations = []
    for r in json.loads(row.get("CALLER_EVIDENCE_JSON", "[]")):
        gq, dp = number(r.get("CALLER_GQ")), number(r.get("CALLER_DP"))
        gt = genotype(r.get("CALLER_GT"))
        if gt != "UNKNOWN" and gq is not None and dp is not None and gq >= min_gq and dp >= min_dp:
            observations.append((r.get("CALLER", "."), gt))
    values = {gt for _, gt in observations}
    if len(values) > 1:
        return "CONFLICT", json.dumps(observations, separators=(",", ":"))
    return (next(iter(values)), json.dumps(observations, separators=(",", ":"))) if values else ("UNKNOWN", "NO_QC_PASSING_CALLER_GENOTYPE")


def inheritance(row, model, patient, gt, family, min_gq, min_dp):
    moi = model.get("moi", "UNKNOWN")
    detail = f"GT_class={gt};MOI={moi}"
    if gt == "CONFLICT":
        return domain("GENOTYPE_CONFLICT", detail=detail)
    if gt == "UNKNOWN" or moi == "UNKNOWN":
        return domain("UNKNOWN", detail=detail)
    if gt == "HOM_REF":
        return domain("REFERENCE_GENOTYPE_CONFLICT", detail=detail)
    chrom = str(row.get("CHROM", "")).removeprefix("chr")
    if moi in {"AD", "AR"} and chrom not in {str(x) for x in range(1, 23)}:
        return domain("PLOIDY_OR_CHROMOSOME_REVIEW", detail=detail)
    if moi == "MT":
        return domain("UNKNOWN_MT_REQUIRES_SPECIALIZED_ANALYSIS", detail=detail)
    if moi in {"XLD", "XLR"}:
        # PAR, sex chromosome dosage and penetrance need locus-specific review.
        return domain("SEX_LINKED_REVIEW_REQUIRED", detail=detail + f";reported_sex={patient.get('sex', '.')}")
    if gt == "HAPLOID_ALT":
        return domain("PLOIDY_REVIEW_REQUIRED", detail=detail)
    if moi == "AR":
        if gt == "HOM_ALT":
            return domain("BIALLELIC_GENOTYPE_COMPATIBLE", 2, detail + ";segregation not established")
        return domain("SECOND_ALLELE_OR_TRANS_PHASE_REQUIRED", 0, detail)
    parents = []
    for relation in ("mother_id", "father_id"):
        parent = family.get((patient.get(relation, "."), row["SV_ID"]), {})
        gq, dp = number(parent.get("gq")), number(parent.get("dp"))
        valid = known(parent.get("source")) and gq is not None and dp is not None and gq >= min_gq and dp >= min_dp
        parents.append(genotype(parent.get("gt")) if valid else "UNKNOWN")
    if gt == "HET" and parents == ["HOM_REF", "HOM_REF"]:
        return domain("CANDIDATE_DE_NOVO", 2, detail + ";parentage/mosaicism/validation require review")
    return domain("DOMINANT_GENOTYPE_COMPATIBLE", 1, detail + ";segregation not established")


def transmitted_parent(row, patient, family, min_gq, min_dp):
    states = []
    for key in ("mother_id", "father_id"):
        r = family.get((patient.get(key, "."), row["SV_ID"]), {})
        gq, dp = number(r.get("gq")), number(r.get("dp"))
        if not known(r.get("source")) or gq is None or dp is None or gq < min_gq or dp < min_dp:
            return None
        states.append(genotype(r.get("gt")))
    if states == ["HET", "HOM_REF"]:
        return "maternal"
    if states == ["HOM_REF", "HET"]:
        return "paternal"
    return None


def trans_candidates(row, model, patient, family, args):
    """Family-based candidates, not validated segregation or phase conclusions.

    Only distinct, nonoverlapping same-gene SVs with compatible effects are
    considered. SNV partners require an external combined analysis.
    """
    if model.get("moi") != "AR" or row.get("SVTYPE") == "BND":
        return []
    origin = transmitted_parent(row, patient, family, args.min_gq, args.min_dp)
    if origin is None:
        return []
    result = []
    for other in getattr(args, "gene_rows", {}).get(row.get("GENES"), []):
        if other["SV_ID"] == row["SV_ID"] or other.get("CHROM") != row.get("CHROM") or other.get("SVTYPE") == "BND":
            continue
        bounds = [number(r.get(k)) for r in (row, other) for k in ("START", "END")]
        if None in bounds or not (bounds[1] < bounds[2] or bounds[3] < bounds[0]):
            continue
        opposite = transmitted_parent(other, patient, family, args.min_gq, args.min_dp)
        if opposite is None or opposite == origin:
            continue
        ev = getattr(args, "evidence_index", {}).get((args.sample, other["SV_ID"], other["GENES"], model.get("disease_id", ".")),
             getattr(args, "evidence_index", {}).get((args.sample, other["SV_ID"], other["GENES"], "."), {}))
        gt, _ = patient_genotype(other, ev, family, args.sample, args.min_gq, args.min_dp)
        _, effect = disruption(other, ev)
        mechanism = disease_mechanism(other, model, effect)[0]
        if gt == "HET" and mechanism in {"DISEASE_MECHANISM_COMPATIBLE", "DISEASE_MECHANISM_POSSIBLE"}:
            result.append(other["SV_ID"])
    return sorted(set(result))


def phenotype(gene, model, patient, gene_hpo):
    present = terms(patient.get("hpo_present", "."))
    absent = terms(patient.get("hpo_absent", "."))
    if not present:
        return domain("UNKNOWN_PATIENT_PHENOTYPES", detail="Generic optic anchors are not patient-specific evidence")
    disease_terms = terms(model.get("hpo_terms", "."))
    reference = disease_terms or gene_hpo.get(gene, set())
    scope = "DISEASE" if disease_terms else "GENE"
    if not reference:
        return domain("UNKNOWN_REFERENCE_PHENOTYPES")
    matched, negative = present & reference, absent & reference
    fraction = len(matched) / len(present)
    points = 2 if fraction >= 0.5 else (1 if matched else 0)
    status = "PATIENT_TERM_CONFLICT_REVIEW" if negative else ("EXACT_TERM_MATCH" if matched else "NO_EXACT_TERM_MATCH")
    detail = f"scope={scope};matched={';'.join(sorted(matched)) or '.'};present_count={len(present)};negative_matches={';'.join(sorted(negative)) or '.'};exact matching only, not semantic similarity"
    return domain(status, 0 if negative else points, detail)


def disease_mechanism(row, model, effect):
    validity = model.get("validity", "UNKNOWN")
    mechanism = model.get("mechanism", "UNKNOWN")
    detail = f"disease={model.get('disease_id', '.')};mechanism={mechanism};validity={validity};source={model.get('source', '.')}"
    if validity in {"DISPUTED", "REFUTED"}:
        return domain("DISEASE_ASSOCIATION_CONFLICT", detail=detail)
    compatible = {"LOF": {"LOF", "HAPLOINSUFFICIENCY"}, "COPY_GAIN": {"TRIPLOSENSITIVITY"},
                  "GOF": {"GOF"}, "DOMINANT_NEGATIVE": {"DOMINANT_NEGATIVE"}}
    if effect in compatible and mechanism in compatible[effect] and validity in {"DEFINITIVE", "STRONG", "MODERATE", "LIMITED"}:
        return domain("DISEASE_MECHANISM_COMPATIBLE", 1 if validity == "LIMITED" else 2, detail)
    if effect in {"POSSIBLE_LOF", "POSSIBLE_COPY_GAIN"} and mechanism in compatible.get(effect.removeprefix("POSSIBLE_"), set()) and validity in {"DEFINITIVE", "STRONG", "MODERATE", "LIMITED"}:
        return domain("DISEASE_MECHANISM_POSSIBLE", 1, detail)
    if effect in compatible and mechanism not in {"UNKNOWN", "OTHER"} and mechanism not in compatible[effect]:
        return domain("MECHANISM_MISMATCH_REVIEW", detail=detail)
    if model.get("disease_id", ".") == ".":
        dosage = row.get("DOSAGE_RELEVANCE", "")
        if (effect in {"LOF", "POSSIBLE_LOF"} and dosage == "HI_SUFFICIENT_EVIDENCE") or (effect in {"COPY_GAIN", "POSSIBLE_COPY_GAIN"} and dosage == "TS_SUFFICIENT_EVIDENCE"):
            return domain("GENE_DOSAGE_SUPPORT_ONLY", 1, "ClinGen HI/TS;not a disease-specific mechanism assignment")
    return domain("UNKNOWN", detail=detail)


def assess(row, model, patient, evidence, family, gene_hpo, args):
    gene = row.get("GENES", ".")
    gt, gt_source = patient_genotype(row, evidence, family, args.sample, args.min_gq, args.min_dp)
    disruption_result, effect = disruption(row, evidence)
    result = {"ALLELE_ASSESSMENT_VERSION": VERSION, "ALLELE_DISEASE_ID": model.get("disease_id", "."),
              "ALLELE_DISEASE_SOURCE": model.get("source", "."), "ALLELE_GENOTYPE": gt,
              "ALLELE_GENOTYPE_SOURCE": gt_source, "ALLELE_EFFECT": effect,
              "ALLELE_FUNCTIONAL_CONTEXT": functional_context(row)}
    components = [population(row, evidence, args.rare_af), technical(row, args.min_support), disruption_result,
                  inheritance(row, model, patient, gt, family, args.min_gq, args.min_dp),
                  phenotype(gene, model, patient, gene_hpo), disease_mechanism(row, model, effect)]
    partners = []
    if gt == "HET" and components[5][0] in {"DISEASE_MECHANISM_COMPATIBLE", "DISEASE_MECHANISM_POSSIBLE"}:
        partners = trans_candidates(row, model, patient, family, args)
        if partners:
            components[3] = domain("INFERRED_TRANS_CANDIDATE_REVIEW", 1,
                "Opposite parental transmission;partners=" + ";".join(partners) + ";parentage, allele identity and validation require review")
    result["ALLELE_TRANS_PARTNER_IDS"] = ";".join(partners) or "."
    flags, unknown = [], []
    for name, (status, points, detail) in zip(DOMAINS, components):
        result[f"ALLELE_{name}_STATUS"] = status
        result[f"ALLELE_{name}_POINTS"] = points
        result[f"ALLELE_{name}_DETAIL"] = detail
        if status.startswith("UNKNOWN"):
            unknown.append(name)
        if any(token in status for token in ("CONFLICT", "REVIEW", "REQUIRED", "PROVISIONAL", "ABOVE_", "LOW_READ")):
            flags.append(f"{name}:{status}")
    result["ALLELE_RESEARCH_SCORE"] = sum(x[1] for x in components)
    result["ALLELE_SCORE_MAXIMUM"] = 12
    result["ALLELE_UNKNOWN_DOMAINS"] = ";".join(unknown) or "."
    result["ALLELE_REVIEW_FLAGS"] = ";".join(flags) or "."
    result["ALLELE_INTERPRETATION"] = "RESEARCH_PRIORITIZATION_ONLY_NOT_ACMG"
    return result


def write_table(path, data, columns):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--integrated", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--phenotypes", help="Human gene phenotype associations")
    p.add_argument("--patient-context", help="Patient-specific HPO/family metadata TSV")
    p.add_argument("--disease-models", help="Curated one gene/disease/MOI/mechanism per row TSV")
    p.add_argument("--allele-evidence", help="Explicit source-backed allele evidence TSV")
    p.add_argument("--family-genotypes", help="Jointly genotyped exact master SV IDs, with QC/source")
    p.add_argument("--rare-af", type=float, default=0.01)
    p.add_argument("--min-support", type=int, default=5)
    p.add_argument("--min-gq", type=float, default=20)
    p.add_argument("--min-dp", type=float, default=5)
    p.add_argument("--output", required=True)
    p.add_argument("--hypotheses-output", required=True)
    p.add_argument("--manifest", required=True)
    args = p.parse_args()
    if not 0 <= args.rare_af <= 1 or min(args.min_support, args.min_gq, args.min_dp) < 0:
        p.error("Invalid research thresholds")
    data = rows(args.integrated, ("SV_ID", "GENES"))
    patients = unique_index(rows(args.patient_context, ("sample_id", "hpo_present", "hpo_absent", "source")), ("sample_id",))
    if args.patient_context and (args.sample,) not in patients:
        raise ValueError(f"No patient-context row for {args.sample}")
    patient = patients.get((args.sample,), {})
    if patient and not known(patient.get("source")):
        raise ValueError("Patient context requires source provenance")
    if terms(patient.get("hpo_present", ".")) & terms(patient.get("hpo_absent", ".")):
        raise ValueError("The same patient HPO term cannot be present and absent")
    models = rows(args.disease_models, ("gene", "disease_id", "moi", "mechanism", "validity", "source"))
    for model in models:
        model["gene"] = model["gene"].upper()
    unique_index(models, ("gene", "disease_id", "moi", "mechanism"))
    by_gene = defaultdict(list)
    for model in models:
        if model["moi"] not in MOIS or model["mechanism"] not in MECHANISMS or model["validity"] not in VALIDITY or not known(model["source"]) or not known(model["disease_id"]):
            raise ValueError(f"Invalid curated disease model: {model}")
        terms(model.get("hpo_terms", "."))
        by_gene[model["gene"].upper()].append(model)
    evidence = rows(args.allele_evidence, ("sample_id", "SV_ID", "gene", "disease_id", "source"))
    for r in evidence:
        r["gene"] = r["gene"].upper()
        if not known(r.get("source")):
            raise ValueError("Allele evidence requires source provenance")
    evidence_index = unique_index(evidence, ("sample_id", "SV_ID", "gene", "disease_id"))
    family = unique_index(rows(args.family_genotypes, ("sample_id", "SV_ID", "gt", "gq", "dp", "source")), ("sample_id", "SV_ID"))
    gene_hpo = defaultdict(set)
    for r in rows(args.phenotypes, ("gene_symbol", "hpo_id")):
        if known(r.get("hpo_id")):
            gene_hpo[r["gene_symbol"].upper()].update(terms(r["hpo_id"]))
    parameters = vars(args).copy()
    args.gene_rows = defaultdict(list)
    args.evidence_index = evidence_index
    for record in data:
        args.gene_rows[record["GENES"]].append(record)
    all_hypotheses, enriched = [], []
    for row in data:
        gene = row["GENES"].upper()
        assessments = []
        for model in by_gene.get(gene, [{}]):
            ev = evidence_index.get((args.sample, row["SV_ID"], gene, model.get("disease_id", ".")),
                                    evidence_index.get((args.sample, row["SV_ID"], gene, "."), {}))
            assessment = assess(row, model, patient, ev, family, gene_hpo, args)
            assessment["ALLELE_DISEASE_MOI"] = model.get("moi", "UNKNOWN")
            assessment["ALLELE_DISEASE_MECHANISM"] = model.get("mechanism", "UNKNOWN")
            assessments.append(assessment)
            all_hypotheses.append({**row, **assessment})
        # No evidence is combined across different diseases. Ties are deterministic.
        best = sorted(assessments, key=lambda r: (-r["ALLELE_RESEARCH_SCORE"], r["ALLELE_DISEASE_ID"], r["ALLELE_DISEASE_MOI"], r["ALLELE_DISEASE_MECHANISM"]))[0]
        enriched.append({**row, **best, "ALLELE_DISEASE_HYPOTHESIS_COUNT": len(assessments)})
    scores = sorted({r["ALLELE_RESEARCH_SCORE"] for r in enriched}, reverse=True)
    ranks = {score: i + 1 for i, score in enumerate(scores)}
    for r in enriched:
        r["ALLELE_RESEARCH_RANK"] = ranks[r["ALLELE_RESEARCH_SCORE"]]
    all_hypotheses.sort(key=lambda r: (-r["ALLELE_RESEARCH_SCORE"], r["SV_ID"], r["GENES"], r["ALLELE_DISEASE_ID"]))
    original_columns = list(data[0]) if data else list(csv.DictReader(open(args.integrated), delimiter="\t").fieldnames or [])
    # Stable schema also for an empty master callset.
    empty = assess({"SV_ID": ".", "GENES": "."}, {}, {}, {}, {}, {}, args)
    extra = list(empty) + ["ALLELE_DISEASE_MOI", "ALLELE_DISEASE_MECHANISM"]
    if set(original_columns) & set(extra):
        raise ValueError("Input already has allele assessment columns; use the original integrated table")
    write_table(args.output, enriched, original_columns + extra + ["ALLELE_DISEASE_HYPOTHESIS_COUNT", "ALLELE_RESEARCH_RANK"])
    write_table(args.hypotheses_output, all_hypotheses, original_columns + extra)
    manifest = {"model": VERSION, "parameters": parameters, "input_sha256": {}}
    for key in ("integrated", "phenotypes", "patient_context", "disease_models", "allele_evidence", "family_genotypes"):
        path = getattr(args, key)
        if path:
            digest = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1024*1024), b""):
                    digest.update(chunk)
            manifest["input_sha256"][key] = digest.hexdigest()
    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    Path(args.manifest).write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[OK] preserved_rows={len(enriched)} disease_hypotheses={len(all_hypotheses)}; research-only score")


if __name__ == "__main__":
    main()
