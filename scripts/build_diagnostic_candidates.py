#!/usr/bin/env python3
"""Prepare readable SV review tables for SRS and LRS.

The order reflects the biological focus of this project:

1. nuclear genes encoding mitochondrial proteins;
2. other established optic-neuropathy genes;
3. remaining genome-wide candidates;
4. mitochondrial-genome findings, kept as a secondary section.

The optional clinical table combines existing CNV classification with HON
relevance and variant evidence. It retains every row and never creates a
diagnostic probability or changes an ACMG class. LRS uses the same function.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from ranking_common import (
    known, mechanism_inheritance_summary, number, parse_moi,
    population_summary, summarize_gencc, technical_summary,
)


FOCUS_ORDER = {
    "NUCLEAR_MITOCHONDRIAL_PRIMARY": 1,
    "KNOWN_DISEASE_GENE": 2,
    "GENOMEWIDE_OTHER": 3,
    "MTDNA_SECONDARY": 4,
}

TECHNICAL_ORDER = {
    "MULTI_CALLER_PLUS_READ_DEPTH": 1,
    "MULTI_CALLER_PLUS_GRIDSS": 1,
    "MULTI_CALLER": 2,
    "SINGLE_CALLER_PLUS_SUPPORTING_EVIDENCE": 3,
    "SINGLE_CALLER": 4,
    "READ_DEPTH_ONLY_CNV": 5,
    "CNVPYTOR_DEPTH_ONLY": 5,
    "REVIEW_REQUIRED": 6,
}

# These are transparent review groups, not ACMG categories or learned scores.
CLINICAL_GROUP_ORDER = {
    "P_LP_HON_MECHANISM_SUPPORTED": 1,
    "P_LP_HON_INCOMPLETE_EVIDENCE": 2,
    "HON_RELEVANT_VUS_OR_UNCLASSIFIED": 3,
    "P_LP_OTHER_OR_UNRESOLVED_HON": 4,
    "OTHER_REVIEW": 5,
    "CONFLICT_OR_TECHNICAL_REVIEW": 6,
    "BENIGN_OR_LIKELY_BENIGN": 7,
}
CNV_CLASSES = {1: "BENIGN", 2: "LIKELY_BENIGN", 3: "VUS",
               4: "LIKELY_PATHOGENIC", 5: "PATHOGENIC"}


def value(row, *names, default="."):
    """Read a known value without replacing missing evidence with zero."""
    return next((str(row[name]).strip() for name in names
                 if known(row.get(name))), default)


def cnv_class(value_text):
    text = str(value_text).upper().replace(" ", "_")
    numeric = number(text)
    if numeric in CNV_CLASSES:
        return CNV_CLASSES[int(numeric)]
    if text in CNV_CLASSES.values():
        return text
    return "VUS" if text == "VARIANT_OF_UNCERTAIN_SIGNIFICANCE" else "UNCLASSIFIED"


def class_from_points(points):
    """Audit the reported class against CNV points; do not assign a new class."""
    if points >= 0.99:
        return "PATHOGENIC"
    if points >= 0.90:
        return "LIKELY_PATHOGENIC"
    if points <= -0.99:
        return "BENIGN"
    if points <= -0.90:
        return "LIKELY_BENIGN"
    return "VUS"


def clinical_row(row):
    """Assess one SV–gene row using existing evidence, without adding points."""
    flags, cautions = [], []
    kind = value(row, "SVTYPE").upper()
    is_cnv = kind in {"DEL", "DUP"}
    classification = cnv_class(value(row, "ACMG_CNV_CLASS", "ANNOTSV_ACMG_CLASS_RAW"))
    points_text = value(row, "ACMG_CNV_SCORE", "ANNOTSV_RANKING_SCORE")
    points = number(points_text)
    criteria = value(row, "ANNOTSV_RANKING_CRITERIA")
    codes = sorted(set(re.findall(r"(?<![A-Za-z0-9])([1-5][A-O](?:-\d+)?)\b", criteria)))
    source = value(row, "ANNOTSV_CLASSIFICATION_SOURCE")
    class_complete = False
    if not is_cnv:
        classification = "NOT_APPLICABLE_NON_CNV"
        audit = "NON_CNV_REQUIRES_SV_SPECIFIC_REVIEW"
    else:
        if source != "FULL_EVENT_ROW":
            if classification != "UNCLASSIFIED" or known(points_text) or codes:
                flags.append("EVENT_CLASSIFICATION_SOURCE_UNVERIFIED")
            else:
                cautions.append("EVENT_CLASSIFICATION_SOURCE_UNAVAILABLE")
        if value(row, "ANNOTSV_CLASSIFICATION_CONFLICT").upper() == "YES":
            flags.append("CONFLICTING_ANNOTSV_EVENT_RECORDS")
        if known(points_text) and points is None:
            flags.append("INVALID_ANNOTSV_POINTS")
        if points is not None and classification != "UNCLASSIFIED" and class_from_points(points) != classification:
            flags.append("ANNOTSV_CLASS_SCORE_DISAGREEMENT")
        if classification == "UNCLASSIFIED":
            cautions.append("CNV_CLASS_UNAVAILABLE")
        if points is None:
            cautions.append("CNV_POINTS_UNAVAILABLE")
        if not codes:
            cautions.append("CNV_CRITERIA_UNAVAILABLE_OR_UNPARSED")
        class_complete = not flags and classification != "UNCLASSIFIED" and points is not None and bool(codes)
        audit = "REPORTED_CLASS_POINTS_CONSISTENT_MANUAL_REVIEW_STILL_REQUIRED" if class_complete else "INCOMPLETE_OR_CONFLICTING_CLASSIFICATION_REVIEW"

    # Case/phenotype evidence can already contribute to CNV classification.
    # Its presence is disclosed; HON matching is never another ACMG point bonus.
    phenotype_codes = {
        "4A", "4B", "4C", "4D", "4E", "4F", "4G", "4H", "4I", "4J", "4L", "4M",
        "5A", "5B", "5C", "5D", "5E", "5G", "5H",
    }
    code_families = {code.split("-")[0] for code in codes}
    phenotype_in_class = (
        "CASE_OR_PHENOTYPE_CRITERIA_REPORTED" if code_families & phenotype_codes
        else ("NOT_SHOWN_IN_REPORTED_CODES" if codes else "UNKNOWN_CRITERIA")
    ) if is_cnv else "NOT_APPLICABLE_NON_CNV"

    patient_score = number(row.get("HPO_QUERY_RESNIK_NORMALIZED"))
    if value(row, "HPO_SEMANTIC_STATUS") in {"EVALUATED", "NO_REFERENCE_HPO_FOR_GENE"} and patient_score is not None:
        hon, hon_scope = patient_score, "PATIENT_OBSERVED_HPO_QUERY"
    else:
        hon = number(value(row, "HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED", "HON_SEMANTIC_SIMILARITY_NORMALIZED"))
        hon_scope = "SHARED_HON_GENE_QUERY_NOT_INDIVIDUAL_CLINICAL_FINDINGS"
    if hon is not None and not 0 <= hon <= 1:
        flags.append("INVALID_PHENOTYPE_SIMILARITY")
        hon = None
    hon_relevance = "UNKNOWN" if hon is None else (
        "STRONG_QUERY_MATCH" if hon >= 0.50 else
        "SUPPORTING_QUERY_MATCH" if hon >= 0.25 else "LIMITED_QUERY_MATCH"
    )
    if value(row, "ALLELE_PHENOTYPE_STATUS") == "PATIENT_TERM_CONFLICT_REVIEW":
        flags.append("REPORTED_PATIENT_PHENOTYPE_CONFLICT")

    disease = summarize_gencc(
        [value(row, "GENCC"), value(row, "GENE_DISEASE_EVIDENCE_LEVEL", "GENE_DISEASE_EVIDENCE")],
        has_omim=known(row.get("OMIM")),
    )
    if disease["conflict"] == "YES" or value(row, "GENE_DISEASE_EVIDENCE_CONFLICT").upper() == "YES":
        flags.append("GENE_DISEASE_EVIDENCE_CONFLICT")
    reviewed_model = known(row.get("ALLELE_DISEASE_ID")) and known(row.get("ALLELE_DISEASE_SOURCE"))
    model_validity_match = re.search(r"(?:^|;)validity=([^;]+)", value(row, "ALLELE_DISEASE_MECHANISM_DETAIL"))
    model_validity = model_validity_match.group(1) if reviewed_model and model_validity_match else "UNKNOWN"
    disease_scope = "REVIEWED_DISEASE_MODEL_HON_FIT_REQUIRES_REVIEW" if reviewed_model else "GENE_LEVEL_ASSOCIATION_NOT_DISEASE_SPECIFIC"
    if not reviewed_model:
        cautions.append("DISEASE_SPECIFIC_HON_MODEL_NOT_AVAILABLE")
    if disease["best"] in {"DISPUTED", "REFUTED"} or model_validity in {"DISPUTED", "REFUTED"}:
        flags.append("DISPUTED_OR_REFUTED_DISEASE_ASSOCIATION")

    # For SRS, allele assessment already describes predicted CDS loss/copy gain.
    # Only use that functional prediction if no geometric relationship exists.
    # Transcript overlap alone does not establish loss of function.
    effect = value(row, "ALLELE_EFFECT")
    mechanism_row = dict(row)
    if not known(row.get("SV_GENE_RELATIONSHIP")):
        if kind in {"DEL", "DUP"} and effect in {"LOF", "POSSIBLE_LOF"}:
            mechanism_row["SV_GENE_RELATIONSHIP"] = "PARTIAL_GENE_OVERLAP"
        elif kind == "DUP" and effect in {"COPY_GAIN", "POSSIBLE_COPY_GAIN"}:
            mechanism_row["SV_GENE_RELATIONSHIP"] = "WHOLE_GENE_DOSAGE_CONTEXT"
    mechanism = mechanism_inheritance_summary(mechanism_row)
    mechanism_class = value(row, "FINAL_INHERITANCE_MECHANISM_CLASS", "INHERITANCE_MECHANISM_CLASS", default=mechanism["category"])
    disease_mechanism = value(row, "ALLELE_DISEASE_MECHANISM_STATUS")
    supported_mechanism = mechanism_class in {
        "STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY", "SUPPORTED_DISEASE_MECHANISM", "AR_TRANS_SECOND_ALLELE_SUPPORTED",
    } or (reviewed_model and disease_mechanism in {"DISEASE_MECHANISM_COMPATIBLE", "DISEASE_MECHANISM_POSSIBLE"})
    if mechanism_class == "DOSAGE_MECHANISM_CONFLICT_REVIEW" or disease_mechanism in {
        "MECHANISM_MISMATCH_REVIEW", "DISEASE_ASSOCIATION_CONFLICT",
    }:
        flags.append("DISEASE_MECHANISM_CONFLICT")
    if not known(effect) or effect in {"UNKNOWN", "NO_DISRUPTION"}:
        supported_mechanism = False
        cautions.append("GENE_SPECIFIC_FUNCTIONAL_EFFECT_UNRESOLVED")
    if not supported_mechanism:
        cautions.append("DISEASE_MECHANISM_NOT_ESTABLISHED")

    # Only ALLELE_GENOTYPE has the caller/external QC used by allele assessment.
    # Merged VCF GT and unrecorded parental genotypes are never substituted.
    gt = value(row, "ALLELE_GENOTYPE", default="UNKNOWN")
    genotype_source_known = known(row.get("ALLELE_GENOTYPE_SOURCE"))
    if not genotype_source_known:
        cautions.append("QC_GENOTYPE_SOURCE_UNAVAILABLE")
    inheritance_status = value(row, "ALLELE_INHERITANCE_STATUS")
    mois = parse_moi(value(row, "ALLELE_DISEASE_MOI") if reviewed_model else ".",
                     value(row, "GENCC_MOI"), value(row, "OMIM_INHERITANCE"))
    if reviewed_model and value(row, "ALLELE_DISEASE_MOI") in {"AD", "AR", "XLD", "XLR", "MT"}:
        mois = {value(row, "ALLELE_DISEASE_MOI")}
    inheritance = "UNKNOWN_GENOTYPE_OR_DISEASE_MODEL"
    if gt in {"CONFLICT", "HOM_REF"} or "CONFLICT" in inheritance_status:
        flags.append("GENOTYPE_OR_INHERITANCE_CONFLICT")
        inheritance = "CONFLICT_REVIEW"
    elif gt == "HET" and mois == {"AR"}:
        trans_candidate = (value(row, "RECESSIVE_PAIR_STATUS") == "AR_TRANS_SNV_SV_CANDIDATE"
                           or known(row.get("ALLELE_TRANS_PARTNER_IDS")))
        inheritance = "TRANS_PARTNER_CANDIDATE_REVIEW" if trans_candidate else "SECOND_ALLELE_REQUIRED"
    elif (gt == "HET" and mois == {"AD"}) or (gt == "HOM_ALT" and mois == {"AR"}):
        inheritance = "GENOTYPE_COMPATIBLE_SEGREGATION_NOT_ESTABLISHED"
    elif gt == "HOM_ALT" and mois == {"AD"}:
        flags.append("DOMINANT_HOMOZYGOUS_GENOTYPE_REVIEW")
    elif mois & {"XLD", "XLR", "MT"} or gt == "HAPLOID_ALT":
        inheritance = "SEX_PLOIDY_OR_MTDNA_SPECIALIST_REVIEW"
    elif {"AD", "AR"}.issubset(mois):
        inheritance = "MIXED_DISEASE_MODELS_REVIEW"
    if "PLOIDY" in inheritance_status or "SEX_LINKED" in inheritance_status:
        inheritance = "SEX_PLOIDY_OR_MTDNA_SPECIALIST_REVIEW"
    if inheritance == "GENOTYPE_COMPATIBLE_SEGREGATION_NOT_ESTABLISHED" and not genotype_source_known:
        inheritance = "GENOTYPE_SOURCE_UNVERIFIED_REVIEW"
    inheritance_supported = inheritance == "GENOTYPE_COMPATIBLE_SEGREGATION_NOT_ESTABLISHED"
    if not inheritance_supported:
        cautions.append(inheritance)

    technical = technical_summary(row)
    allele_technical = value(row, "ALLELE_TECHNICAL_STATUS")
    technical_supported = allele_technical in {"MULTICALLER_READ_SUPPORTED", "SINGLE_CALLER_READ_SUPPORTED"}
    if not known(allele_technical):
        technical_supported = technical["tier"] in {"STRONG", "MODERATE"} and number(technical["max_read_support"]) is not None
    if allele_technical in {"REVIEW_REQUIRED", "LOW_READ_SUPPORT"} or technical["review"] == "YES" or value(row, "FILTER", default="PASS") != "PASS":
        technical_supported = False
        flags.append("CALLER_QC_OR_READ_SUPPORT_REVIEW")
    if value(row, "DEPTH_SUPPORT_CLASS") == "DOES_NOT_SHOW_EXPECTED_COPY_CHANGE" or value(row, "DEPTH_PATTERN") == "NOT_CLEAR":
        technical_supported = False
        flags.append("COPY_CHANGE_DEPTH_CONFIRMATION_REQUIRED")
    if not technical_supported:
        flags.append("TECHNICAL_CONFIRMATION_REQUIRED")

    population_row = dict(row)
    for column in ("GNOMAD_SV_AF", "NEEDLR_AF"):
        af = number(row.get(column))
        if known(row.get(column)) and (af is None or not 0 <= af <= 1):
            population_row[column] = "."
            flags.append("INVALID_POPULATION_AF")
    population = population_summary(population_row)
    af = number(population["max_af"])
    population_status = "UNKNOWN_AF" if af is None else (
        "REPORTED_AF_LE_1_PERCENT" if af <= 0.01 else "REPORTED_AF_GT_1_PERCENT_REVIEW"
    )
    if af is None:
        cautions.append("MISSING_AF_IS_NOT_PROOF_OF_NOVELTY")
    if af is not None and af > 0.01:
        flags.append("COMMON_OR_PROVISIONAL_POPULATION_MATCH_REVIEW")
    if population["conflict"] == "YES":
        flags.append("POPULATION_AF_SOURCE_CONFLICT")

    # Database region overlap is useful context, but may describe another allele.
    # It never adds another ClinVar/pathogenicity bonus on top of AnnotSV.
    pathogenic_region = known(row.get("SV_PATHOGENIC_DB_SOURCE"))
    benign_region = known(row.get("SV_BENIGN_DB_SOURCE"))
    database_audit = "NO_REPORTED_REGION_MATCH_NOT_PROOF_OF_ABSENCE"
    if pathogenic_region or benign_region:
        database_audit = "REGION_OVERLAP_VERIFY_ALLELE_IDENTITY_AND_DISEASE"
        cautions.append("DATABASE_REGION_OVERLAP_IS_NOT_AN_EXACT_VARIANT_CLASSIFICATION")
    if pathogenic_region and benign_region:
        database_audit = "PATHOGENIC_AND_BENIGN_REGION_RECORDS_REQUIRE_REVIEW"
        cautions.append("MIXED_DATABASE_REGION_EVIDENCE_REVIEW_ALLELE_IDENTITY")

    relevant = hon_relevance in {"STRONG_QUERY_MATCH", "SUPPORTING_QUERY_MATCH"}
    positive_class = classification in {"PATHOGENIC", "LIKELY_PATHOGENIC"}
    disease_specific_support = (reviewed_model and model_validity in {"DEFINITIVE", "STRONG", "MODERATE"}
                               and disease_mechanism in {"DISEASE_MECHANISM_COMPATIBLE", "DISEASE_MECHANISM_POSSIBLE"})
    if classification in {"BENIGN", "LIKELY_BENIGN"} and class_complete:
        group = "BENIGN_OR_LIKELY_BENIGN"
    elif flags:
        group = "CONFLICT_OR_TECHNICAL_REVIEW"
    elif positive_class and relevant:
        group = "P_LP_HON_MECHANISM_SUPPORTED" if (
            class_complete and hon_relevance == "STRONG_QUERY_MATCH"
            and supported_mechanism and inheritance_supported and disease_specific_support
        ) else "P_LP_HON_INCOMPLETE_EVIDENCE"
    elif relevant:
        group = "HON_RELEVANT_VUS_OR_UNCLASSIFIED"
    elif positive_class:
        group = "P_LP_OTHER_OR_UNRESOLVED_HON"
    else:
        group = "OTHER_REVIEW"

    return {
        "CLINICAL_REVIEW_GROUP": group,
        "CLINICAL_CNV_CLASS": classification,
        "CLINICAL_CNV_CLASS_AUDIT": audit,
        "CLINICAL_CNV_FRAMEWORK": "ACMG_CLINGEN_COPY_NUMBER_LOSS" if kind == "DEL" else "ACMG_CLINGEN_COPY_NUMBER_GAIN" if kind == "DUP" else "NOT_APPLICABLE_NON_CNV",
        "CLINICAL_CLASSIFICATION_UNIT": "WHOLE_CNV_EVENT_NOT_EACH_GENE" if is_cnv else "SV_SPECIFIC_INTERPRETATION_REQUIRED",
        "CLINICAL_CNV_CRITERIA_CODES": ";".join(codes) if is_cnv and codes else ".",
        "CLINICAL_CNV_FRAMEWORK_AUDIT": ";".join(
            f"section_{section}:" + (",".join(c for c in codes if c.startswith(str(section))) or "NO_REPORTED_CODE")
            for section in range(1, 6)
        ) if is_cnv else "NOT_APPLICABLE_NON_CNV",
        "CLINICAL_CNV_MANUAL_OR_UNAVAILABLE": (
            "1:breakpoint/transcript/regulatory content;2:loss-HI or gain-TS and exact region identity;"
            "3:event-wide coding gene count;4:independent case/control evidence and phenotype specificity;"
            "5:patient/family segregation and callable genotypes;NO_REPORTED_CODE_DOES_NOT_MEAN_CRITERION_EXCLUDED"
        ) if is_cnv else "Balanced/other SVs need mechanism-specific interpretation; do not apply the DEL/DUP point framework",
        "CLINICAL_PHENOTYPE_IN_CNV_CRITERIA": phenotype_in_class,
        "CLINICAL_HON_RELEVANCE": hon_relevance,
        "CLINICAL_HON_QUERY_SCOPE": hon_scope if hon is not None else "NO_USABLE_PHENOTYPE_SCORE",
        "CLINICAL_HON_QUERY_SIMILARITY": hon if hon is not None else ".",
        "CLINICAL_DISEASE_EVIDENCE_LEVEL": disease["best"],
        "CLINICAL_DISEASE_EVIDENCE_SCOPE": disease_scope,
        "CLINICAL_DISEASE_MODEL_VALIDITY": model_validity,
        "CLINICAL_DATABASE_OVERLAP_AUDIT": database_audit,
        "CLINICAL_MECHANISM_CLASS": mechanism_class,
        "CLINICAL_MECHANISM_SUPPORTED": "YES" if supported_mechanism else "NO_OR_UNRESOLVED",
        "CLINICAL_INHERITANCE_STATUS": inheritance,
        "CLINICAL_TECHNICAL_STATUS": "READ_SUPPORTED_NOT_INDEPENDENT_VALIDATION" if technical_supported else "REVIEW_REQUIRED",
        "CLINICAL_POPULATION_STATUS": population_status,
        "CLINICAL_MAX_REPORTED_AF": af if af is not None else ".",
        "CLINICAL_REVIEW_FLAGS": ";".join(sorted(set(flags))) or ".",
        "CLINICAL_MISSING_OR_CONTEXT_NOTES": ";".join(sorted(set(cautions))) or ".",
        "CLINICAL_MANUAL_REVIEW": "Confirm breakpoints/copy change, transcript effect, disease-specific mechanism, phenotype and segregation; audit AnnotSV automated CNV criteria",
        "CLINICAL_RANKING_METHOD": "JOINT_REVIEW_GROUPS_V1_NO_ADDITIVE_ACMG_OR_PHENOTYPE_POINTS",
        "CLINICAL_USE_NOTE": "Candidate review only; no diagnostic probability, new ACMG class or validated diagnosis; a multi-gene CNV class is not proof for each gene",
    }


def build_clinical_table(evidence: pd.DataFrame, context: pd.DataFrame | None = None) -> pd.DataFrame:
    """Keep every evidence row; add final LRS phase/HPO/depth context if present."""
    # Large transcript/caller audit payloads stay in the original evidence file.
    # They are not used here; drop them before making per-row dictionaries.
    result = evidence.drop(columns=[c for c in evidence if c.endswith("_JSON")]).copy()
    gene_column = find_column(result, ["GENE", "GENES", "ANNotsv_Gene"])
    if gene_column is None or "SV_ID" not in result:
        raise ValueError("Clinical evidence requires SV_ID and a gene column")
    result["GENE"] = [value(row, gene_column, "GENE", "GENES", "ANNotsv_Gene") for row in result.to_dict("records")]
    if result.duplicated(["SV_ID", "GENE"]).any():
        raise ValueError("Clinical evidence has duplicate SV–gene keys; resolve them before ranking")
    if context is not None and not context.empty:
        columns = [c for c in context if c.startswith(("HPO_", "FINAL_", "DEPTH_", "RECESSIVE_"))
                   or c in {"SV_GENE_RELATIONSHIP", "SV_GENE_EFFECT", "MEDIAN_DEPTH_INSIDE", "MEDIAN_DEPTH_FLANKS"}]
        context_view = context[["SV_ID", "GENE", *columns]]
        result = result.merge(context_view, on=["SV_ID", "GENE"], how="left", validate="one_to_one", suffixes=("", "_context"))
        for column in columns:
            if column + "_context" in result:
                use_context = result[column + "_context"].map(known)
                result[column] = result[column].astype(object)
                result.loc[use_context, column] = result.loc[use_context, column + "_context"]
                result = result.drop(columns=column + "_context")
    assessments = [clinical_row(row) for row in result.to_dict("records")]
    field_names = list(clinical_row({}))
    for column in field_names:
        result[column] = [assessment[column] for assessment in assessments]
    def sort_key(row):
        # The class is used once. Points/criterion count, panel and MitoCarta
        # membership are not extra bonuses. No candidate is filtered out.
        hon_order = {"STRONG_QUERY_MATCH": 3, "SUPPORTING_QUERY_MATCH": 2, "LIMITED_QUERY_MATCH": 1}
        class_order = {"PATHOGENIC": 3, "LIKELY_PATHOGENIC": 2, "VUS": 1}
        disease_order = {"DEFINITIVE": 3, "STRONG": 2, "MODERATE": 1}
        disease_level = row["CLINICAL_DISEASE_MODEL_VALIDITY"] if known(row.get("ALLELE_DISEASE_ID")) else row["CLINICAL_DISEASE_EVIDENCE_LEVEL"]
        return (CLINICAL_GROUP_ORDER[row["CLINICAL_REVIEW_GROUP"]],
                -hon_order.get(row["CLINICAL_HON_RELEVANCE"], 0),
                -class_order.get(row["CLINICAL_CNV_CLASS"], 0),
                -int(row["CLINICAL_MECHANISM_SUPPORTED"] == "YES"),
                -int(row["CLINICAL_INHERITANCE_STATUS"] == "GENOTYPE_COMPATIBLE_SEGREGATION_NOT_ESTABLISHED"),
                -disease_order.get(disease_level, 0),
                value(row, "SV_ID"), value(row, "GENE"))
    order = sorted(range(len(result)), key=lambda i: sort_key(result.iloc[i]))
    result = result.iloc[order].reset_index(drop=True)
    result.insert(0, "CLINICAL_REVIEW_RANK", range(1, len(result) + 1))
    first = ["CLINICAL_REVIEW_RANK", "SV_ID", "GENE", *field_names]
    return result[first + [c for c in result if c not in first]]


def find_column(table: pd.DataFrame, possible_names: list[str]) -> str | None:
    """Return the first matching column, ignoring capitalization."""
    lower_names = {str(column).lower(): column for column in table.columns}

    for name in possible_names:
        if name in table.columns:
            return name
        if name.lower() in lower_names:
            return lower_names[name.lower()]

    return None


def positive_panel_status(values: pd.Series) -> pd.Series:
    """Recognize the panel labels used by the shared LRS/SRS scripts."""
    normalized = values.fillna("").astype(str).str.upper().str.strip()
    return normalized.isin({"PANEL_GENE", "YES", "TRUE", "1"})


def add_biological_focus(table: pd.DataFrame) -> pd.DataFrame:
    """Label each row according to the study's biological priority."""
    result = table.copy()
    result["DIAGNOSTIC_FOCUS"] = "GENOMEWIDE_OTHER"

    encoding = result.get(
        "MITOCARTA_ENCODING",
        pd.Series(".", index=result.index),
    ).fillna(".")

    panel_values = result.get(
        "PANEL_STATUS",
        result.get("panel_gene", pd.Series(".", index=result.index)),
    )

    is_panel_gene = positive_panel_status(panel_values)
    is_nuclear_mito = encoding.eq("NUCLEAR_MITOCHONDRIAL_GENE")
    is_mtdna = encoding.eq("MTDNA_ENCODED_GENE")

    result.loc[is_panel_gene, "DIAGNOSTIC_FOCUS"] = "KNOWN_DISEASE_GENE"
    result.loc[is_nuclear_mito, "DIAGNOSTIC_FOCUS"] = (
        "NUCLEAR_MITOCHONDRIAL_PRIMARY"
    )
    result.loc[is_mtdna, "DIAGNOSTIC_FOCUS"] = "MTDNA_SECONDARY"

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--top-n", type=int, default=0, help="0 keeps all rows; a positive value creates a limited review view")
    parser.add_argument("--output", required=True)
    parser.add_argument("--clinical-output", help="Additional ACMG/HON joint review table; always retains every row")
    args = parser.parse_args()
    if args.top_n < 0:
        parser.error("top-n must be zero or positive")

    candidates = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    gene_column = find_column(
        candidates,
        ["GENE", "GENES", "Gene", "gene", "ANNotsv_Gene"],
    )
    if gene_column is None:
        raise ValueError("Candidate input has no recognized gene column.")

    candidates = add_biological_focus(candidates)
    if args.clinical_output:
        clinical = build_clinical_table(candidates)
        clinical_path = Path(args.clinical_output)
        clinical_path.parent.mkdir(parents=True, exist_ok=True)
        clinical.to_csv(clinical_path, sep="\t", index=False)
    candidates["_FOCUS_ORDER"] = (
        candidates["DIAGNOSTIC_FOCUS"].map(FOCUS_ORDER).fillna(9)
    )

    score_column = find_column(
        candidates,
        [
            "ALLELE_RESEARCH_SCORE",
            "EVENT_GENE_RELEVANCE_SCORE",
            "INTEGRATED_DISCOVERY_SCORE",
            "PHENOTYPE_SCORE",
        ],
    )
    if score_column:
        candidates["_BIOLOGICAL_SCORE"] = pd.to_numeric(
            candidates[score_column],
            errors="coerce",
        ).fillna(0)
    else:
        candidates["_BIOLOGICAL_SCORE"] = 0

    technical_support = candidates.get(
        "TECHNICAL_SUPPORT",
        pd.Series("REVIEW_REQUIRED", index=candidates.index),
    )
    candidates["_TECHNICAL_ORDER"] = (
        technical_support.map(TECHNICAL_ORDER).fillna(9)
    )

    candidates = candidates.sort_values(
        ["_FOCUS_ORDER", "_BIOLOGICAL_SCORE", "_TECHNICAL_ORDER"],
        ascending=[True, False, True],
        kind="stable",
    )
    if args.top_n:
        candidates = candidates.head(args.top_n)

    candidates = candidates.copy()
    candidates.insert(
        0,
        "DIAGNOSTIC_REVIEW_RANK",
        range(1, len(candidates) + 1),
    )
    candidates["DIAGNOSTIC_USE_NOTE"] = (
        "Research/diagnostic-support candidate; orthogonal validation and "
        "clinical interpretation are required."
    )

    candidates = candidates.drop(
        columns=["_FOCUS_ORDER", "_BIOLOGICAL_SCORE", "_TECHNICAL_ORDER"]
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output, sep="\t", index=False)

    print(f"[OK] diagnostic_candidates={len(candidates)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
