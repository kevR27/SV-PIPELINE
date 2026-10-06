#!/usr/bin/env python3
"""Shared ranking helpers for SV/gene prioritization.

The functions in this module keep evidence axes separate:
- phenotype relevance;
- curated gene-disease evidence;
- SV-gene mechanism and inheritance compatibility;
- technical support;
- explicit population frequency.

They are research-prioritization helpers, not ACMG/AMP pathogenicity rules.
"""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from functools import lru_cache

import pandas as pd


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL", "UNKNOWN"}

DIRECT_RELATIONSHIPS = {
    "WHOLE_GENE_DOSAGE_CONTEXT",
    "PARTIAL_GENE_OVERLAP",
    "INSERTION_WITHIN_TRANSCRIPT",
    "BREAKPOINT_WITHIN_TRANSCRIPT",
}

PROXIMAL_RELATIONSHIPS = {
    "INSERTION_PROXIMAL_TO_GENE",
    "BREAKPOINT_PROXIMAL_TO_GENE",
    "GENE_PROXIMAL_INTERVAL",
}

CONTEXT_RELATIONSHIPS = {
    "INVERSION_SPANS_INTACT_GENE",
    "INTERVAL_CONTEXT_ONLY",
}

GENCC_SCORES = {
    "DEFINITIVE": 4.0,
    "STRONG": 3.5,
    "MODERATE": 2.5,
    "SUPPORTIVE": 1.5,
    "LIMITED": 0.5,
    "ANIMAL MODEL ONLY": 0.0,
    "DISPUTED": 0.0,
    "REFUTED": 0.0,
    "NO KNOWN DISEASE RELATIONSHIP": 0.0,
}

POSITIVE_GENCC = {
    "DEFINITIVE", "STRONG", "MODERATE", "SUPPORTIVE", "LIMITED",
}
CONTRADICTORY_GENCC = {
    "DISPUTED", "REFUTED", "NO KNOWN DISEASE RELATIONSHIP",
}

GENE_TIER_PRIORITY = {
    "HIGH": 3,
    "MODERATE": 2,
    "SUPPORTING": 1,
    "LIMITED": 0,
}

TECHNICAL_PRIORITY = {
    "STRONG": 3,
    "MODERATE": 2,
    "SUPPORTING": 1,
    "REVIEW": 0,
}

POPULATION_PRIORITY = {
    "VERY_RARE_MAX_AF_LE_0.001": 3,
    "LOW_FREQUENCY_MAX_AF_0.001_TO_0.01": 2,
    "UNKNOWN": 1,
    "TOO_COMMON_MAX_AF_GT_0.01": 0,
}

MECHANISM_PRIORITY = {
    "STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY": 5,
    "AR_TRANS_SECOND_ALLELE_SUPPORTED": 5,
    "SUPPORTED_DISEASE_MECHANISM": 4,
    "AD_LOF_GENOTYPE_UNRESOLVED": 3,
    "DIRECT_EFFECT_RECESSIVE_SECOND_ALLELE_REQUIRED": 3,
    "MIXED_AD_AR_DISEASE_MODEL_REVIEW": 3,
    "DIRECT_EFFECT_MOI_UNKNOWN": 3,
    "SEX_LINKED_REVIEW": 3,
    "AD_HOMOZYGOUS_ALT_REVIEW": 2,
    "PROXIMAL_CONTEXT": 2,
    "DIRECT_COPY_GAIN_WITHOUT_TS_SUPPORT": 2,
    "INTERVAL_CONTEXT_ONLY": 1,
    "DOSAGE_MECHANISM_CONFLICT_REVIEW": 0,
    "UNRESOLVED": 0,
    "MTDNA_SPECIALIZED_REVIEW": 0,
}


def known(value) -> bool:
    return value is not None and str(value).strip().upper() not in MISSING


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except Exception:
        return None


def first_existing(df: pd.DataFrame, names: list[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        hit = lower.get(name.lower())
        if hit is not None:
            return hit
    return None


def split_values(value) -> list[str]:
    if not known(value):
        return []
    return [
        token.strip()
        for token in re.split(r"[;|,]", str(value))
        if known(token)
    ]


def normalize_gencc_term(value: str) -> str | None:
    text = str(value or "").strip().upper()
    if not text:
        return None
    if "NO KNOWN" in text:
        return "NO KNOWN DISEASE RELATIONSHIP"
    if "ANIMAL" in text:
        return "ANIMAL MODEL ONLY"
    for term in (
        "DEFINITIVE", "STRONG", "MODERATE", "SUPPORTIVE", "LIMITED",
        "DISPUTED", "REFUTED",
    ):
        if term in text:
            return term
    return None


def summarize_gencc(
    values: list[str],
    has_omim: bool,
    hon_semantic_normalized: float | None = None,
) -> dict:
    """Summarize curated evidence without letting generic disease evidence dominate.

    Full disease-evidence weight is used when the gene has clear HON phenotype
    similarity. Otherwise the evidence is retained at reduced weight because the
    GenCC/OMIM disease may concern a phenotype unrelated to optic neuropathy.

    This is currently a gene-level HON-context proxy, not disease-ID-specific
    HPO matching, because the integrated AnnotSV/GenCC rows do not provide a
    stable disease identifier that can always be joined to Monarch.
    """
    terms = set()
    for raw in values:
        for piece in split_values(raw) or [raw]:
            term = normalize_gencc_term(piece)
            if term:
                terms.add(term)

    positive = sorted(
        (t for t in terms if t in POSITIVE_GENCC),
        key=lambda t: GENCC_SCORES[t],
        reverse=True,
    )
    contradictory = bool(terms & CONTRADICTORY_GENCC)
    conflict = bool(positive and contradictory)

    if positive:
        best = positive[0]
        raw_score = GENCC_SCORES[best]
        source = "GenCC"
    elif "ANIMAL MODEL ONLY" in terms:
        best = "ANIMAL MODEL ONLY"
        raw_score = 0.0
        source = "GenCC"
    elif terms:
        best = sorted(terms)[0]
        raw_score = 0.0
        source = "GenCC"
    elif has_omim:
        best = "OMIM_DISEASE_EVIDENCE"
        raw_score = 1.0
        source = "OMIM"
    else:
        best = "NO_CURATED_EVIDENCE"
        raw_score = 0.0
        source = "NONE"

    semantic = hon_semantic_normalized if hon_semantic_normalized is not None else 0.0
    context_factor = 1.0 if semantic >= 0.25 else 0.4
    conflict_factor = 0.5 if conflict else 1.0
    adjusted = raw_score * context_factor * conflict_factor

    return {
        "best": best,
        "raw_score": raw_score,
        "adjusted_score": adjusted,
        "source": source,
        "conflict": "YES" if conflict else "NO",
        "all_terms": ";".join(sorted(terms)) if terms else ".",
        "hon_context_factor": context_factor,
        "conflict_factor": conflict_factor,
        "scope": (
            "GENE_LEVEL_HON_SEMANTIC_PROXY_NOT_DISEASE_SPECIFIC"
            if source != "NONE"
            else "NO_CURATED_DISEASE_EVIDENCE"
        ),
    }


def gene_relevance(
    hon_semantic_normalized: float | None,
    disease_adjusted_score: float | None,
) -> dict:
    """Balanced phenotype/disease relevance used for tiering.

    The continuous score is for display/tie-breaking only. Primary ordering
    should use the discrete tier so tiny score differences cannot dominate
    SV mechanism or technical evidence.
    """
    hon = min(max(hon_semantic_normalized or 0.0, 0.0), 1.0)
    disease_signal = min(max((disease_adjusted_score or 0.0) / 4.0, 0.0), 1.0)

    display_score = 5.0 * hon + 5.0 * disease_signal

    if hon >= 0.50 and disease_signal >= 0.50:
        tier = "HIGH"
    elif hon >= 0.50 or disease_signal >= 0.50:
        tier = "MODERATE"
    elif hon >= 0.15 or disease_signal > 0:
        tier = "SUPPORTING"
    else:
        tier = "LIMITED"

    return {
        "tier": tier,
        "tier_rank": GENE_TIER_PRIORITY[tier],
        "display_score": round(display_score, 3),
        "hon_signal": round(hon, 6),
        "disease_signal": round(disease_signal, 6),
    }


def parse_moi(*values) -> set[str]:
    """Normalize GenCC/OMIM inheritance text into broad review categories."""
    text = " ; ".join(str(v) for v in values if known(v)).upper()
    result = set()

    patterns = [
        ("AD", r"AUTOSOMAL\s+DOMINANT|\bAD\b|MONOALLELIC"),
        ("AR", r"AUTOSOMAL\s+RECESSIVE|\bAR\b|BIALLELIC"),
        ("XLD", r"X[- ]LINKED\s+DOMINANT|\bXLD\b"),
        ("XLR", r"X[- ]LINKED\s+RECESSIVE|\bXLR\b"),
        ("MT", r"MITOCHONDRIAL|MATERNAL|\bMT\b"),
    ]
    for label, pattern in patterns:
        if re.search(pattern, text):
            result.add(label)
    return result


def inheritance_class(mois: set[str]) -> str:
    if not mois:
        return "UNKNOWN"
    if mois == {"AD"}:
        return "AUTOSOMAL_DOMINANT"
    if mois == {"AR"}:
        return "AUTOSOMAL_RECESSIVE"
    if mois <= {"XLD", "XLR"}:
        return "X_LINKED"
    if mois == {"MT"}:
        return "MITOCHONDRIAL"
    return "MIXED_OR_MULTIPLE_MOI"


def parse_clingen_dosage(value) -> int | None:
    """Extract ClinGen dosage score while preserving special 30/40 codes."""
    if not known(value):
        return None
    text = str(value).upper()
    # Prefer explicit special codes when present.
    if re.search(r"(^|\D)40(\D|$)", text):
        return 40
    if re.search(r"(^|\D)30(\D|$)", text):
        return 30
    for score in (3, 2, 1, 0):
        if re.search(rf"(^|\D){score}(\D|$)", text):
            return score
    if "SUFFICIENT" in text:
        return 3
    if "EMERGING" in text:
        return 2
    if "LITTLE" in text:
        return 1
    if "UNLIKELY" in text:
        return 40
    return None


def extract_loeuf(row) -> float | None:
    available = set(row.index) if hasattr(row, "index") else set(row)
    for name in (
        "GNOMAD_LOEUF", "LOEUF", "LOEUF_score", "ANNOTSV_LOEUF",
    ):
        if name in available:
            value = number(row.get(name))
            if value is not None:
                return value
    return None


def constraint_tiebreak(row) -> float:
    """Higher is more LoF constrained; only a modest late tie-break."""
    loeuf = extract_loeuf(row)
    if loeuf is None:
        return 0.0
    return max(0.0, min(1.0, 1.0 - loeuf))


def relationship_priority(value: str) -> int:
    value = str(value or "").upper()
    if value in DIRECT_RELATIONSHIPS:
        return 3
    if value in PROXIMAL_RELATIONSHIPS:
        return 2
    if value in CONTEXT_RELATIONSHIPS:
        return 1
    return 0


def inferred_effect_class(svtype: str, relationship: str) -> str:
    svtype = str(svtype or "").upper()
    relationship = str(relationship or "").upper()

    if relationship in CONTEXT_RELATIONSHIPS:
        return "INTERVAL_CONTEXT"
    if relationship in PROXIMAL_RELATIONSHIPS:
        return "PROXIMAL_CONTEXT"
    if relationship not in DIRECT_RELATIONSHIPS:
        return "UNRESOLVED"

    if svtype in {"DEL"}:
        return "LOSS_OR_LOF"
    if svtype in {"DUP"}:
        return "COPY_GAIN"
    if svtype in {"INS", "INV", "BND", "TRA"}:
        return "DISRUPTION_OR_LOF"
    if svtype == "CNV":
        return "DOSAGE_CHANGE"
    return "DIRECT_OTHER"


def mechanism_inheritance_summary(row) -> dict:
    relationship = str(row.get("SV_GENE_RELATIONSHIP", ".") or ".").upper()
    svtype = str(row.get("SVTYPE", ".") or ".").upper()
    effect = inferred_effect_class(svtype, relationship)

    gene = str(
        row.get("GENES", row.get("GENE", row.get("ANNotsv_Gene", ".")))
    ).upper()
    chrom = str(row.get("CHROM", "") or "").removeprefix("chr").upper()

    mois = parse_moi(row.get("GENCC_MOI"), row.get("OMIM_INHERITANCE"))
    hi = parse_clingen_dosage(row.get("CLINGEN_HI"))
    ts = parse_clingen_dosage(row.get("CLINGEN_TS"))
    if hi == 30:
        mois.add("AR")

    moi_class = inheritance_class(mois)
    gt = str(row.get("ALLELE_GENOTYPE", "UNKNOWN") or "UNKNOWN").upper()

    if chrom in {"M", "MT"} or gene.startswith("MT-") or "MT" in mois:
        category = "MTDNA_SPECIALIZED_REVIEW"
        detail = "mtDNA/mitochondrial inheritance requires heteroplasmy and mtDNA-specific interpretation"
    elif effect == "INTERVAL_CONTEXT":
        category = "INTERVAL_CONTEXT_ONLY"
        detail = "gene is interval/inversion context without a direct transcript breakpoint"
    elif effect == "PROXIMAL_CONTEXT":
        category = "PROXIMAL_CONTEXT"
        detail = "breakpoint/insertion is proximal but not demonstrated to disrupt the transcript"
    elif effect in {"LOSS_OR_LOF", "DISRUPTION_OR_LOF"}:
        if hi == 40:
            category = "DOSAGE_MECHANISM_CONFLICT_REVIEW"
            detail = "ClinGen dosage sensitivity unlikely for haploinsufficiency"
        elif gt == "HOM_ALT" and "AR" in mois:
            category = "STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY"
            detail = "biallelic SV genotype is compatible with an AR disease model; segregation still requires review"
        elif hi == 3:
            category = "STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY"
            detail = "ClinGen HI=3 supports monoallelic loss sensitivity"
        elif hi == 2:
            category = "SUPPORTED_DISEASE_MECHANISM"
            detail = "ClinGen HI=2 provides emerging haploinsufficiency support"
        elif {"AD", "AR"}.issubset(mois):
            category = "MIXED_AD_AR_DISEASE_MODEL_REVIEW"
            detail = (
                "gene has both dominant and recessive disease models; without a "
                "disease-specific mechanism match this SV is not promoted as a "
                "simple dominant or recessive event"
            )
        elif "AD" in mois:
            if gt == "HET":
                category = "SUPPORTED_DISEASE_MECHANISM"
                detail = (
                    "heterozygous direct loss/disruption is compatible with a "
                    "reported dominant disease model; disease-specific molecular "
                    "mechanism and segregation still require review"
                )
            elif gt == "HOM_ALT":
                category = "AD_HOMOZYGOUS_ALT_REVIEW"
                detail = (
                    "gene has a dominant disease model but the SV genotype is "
                    "homozygous alternate; review disease mechanism, dosage, "
                    "viability and caller genotype before prioritizing"
                )
            else:
                category = "AD_LOF_GENOTYPE_UNRESOLVED"
                detail = (
                    "direct loss/disruption is compatible with a reported dominant "
                    "disease model, but the patient SV genotype is unresolved"
                )
        elif "AR" in mois:
            category = "DIRECT_EFFECT_RECESSIVE_SECOND_ALLELE_REQUIRED"
            detail = (
                "direct SV effect in an AR gene; a second pathogenic allele or "
                "biallelic SV state is required"
            )
        elif mois & {"XLD", "XLR"}:
            category = "SEX_LINKED_REVIEW"
            detail = "X-linked disease model requires sex/ploidy and locus-specific review"
        else:
            category = "DIRECT_EFFECT_MOI_UNKNOWN"
            detail = "direct structural effect but inheritance/mechanism is not established"
    elif effect == "COPY_GAIN":
        if ts == 40:
            category = "DOSAGE_MECHANISM_CONFLICT_REVIEW"
            detail = "ClinGen dosage sensitivity unlikely for triplosensitivity"
        elif ts == 3:
            category = "STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY"
            detail = "ClinGen TS=3 supports copy-gain sensitivity"
        elif ts == 2:
            category = "SUPPORTED_DISEASE_MECHANISM"
            detail = "ClinGen TS=2 provides emerging triplosensitivity support"
        else:
            category = "DIRECT_COPY_GAIN_WITHOUT_TS_SUPPORT"
            detail = "direct copy gain; dominance alone does not establish triplosensitivity"
    else:
        category = "UNRESOLVED"
        detail = "SV-gene mechanism could not be resolved"

    return {
        "effect_class": effect,
        "moi_set": ";".join(sorted(mois)) if mois else ".",
        "inheritance_class": moi_class,
        "clingen_hi_score": hi if hi is not None else ".",
        "clingen_ts_score": ts if ts is not None else ".",
        "category": category,
        "priority": MECHANISM_PRIORITY[category],
        "detail": detail,
    }


def parse_read_support(value) -> float | None:
    if not known(value):
        return None
    numbers = []
    for token in re.split(r"[;|]", str(value)):
        if ":" in token:
            token = token.rsplit(":", 1)[-1]
        val = number(token)
        if val is not None:
            numbers.append(val)
    return max(numbers) if numbers else None


def technical_summary(row) -> dict:
    callers = number(row.get("CALLER_COUNT"))
    if callers is None:
        callers = number(row.get("SUPP")) or 0.0
    read_support = parse_read_support(row.get("CALLER_READ_SUPPORT"))
    flags = str(row.get("CALLER_EVIDENCE_FLAGS", "") or "").upper()
    match = str(row.get("CALLER_EVIDENCE_MATCH", "") or "").upper()

    review_tokens = (
        "RESCUED_COV_VAR", "LOW_GQ", "IMPRECISE", "BLACKLIST_REGION",
        "CALLER_RECORD_FAILS_CURRENT_FILTER",
    )
    flagged = "AMBIGUOUS" in match or any(t in flags for t in review_tokens)

    if flagged:
        tier = "REVIEW"
    elif callers >= 2 and read_support is not None and read_support >= 5:
        tier = "STRONG"
    elif callers >= 2 or (read_support is not None and read_support >= 5):
        tier = "MODERATE"
    elif read_support is not None and read_support >= 2:
        tier = "SUPPORTING"
    else:
        tier = "REVIEW"

    return {
        "tier": tier,
        "priority": TECHNICAL_PRIORITY[tier],
        "caller_count": callers,
        "max_read_support": read_support if read_support is not None else ".",
        "review": "YES" if flagged else "NO",
    }


def population_summary(
    row,
    rare_af: float = 0.001,
    max_af: float = 0.01,
) -> dict:
    """Use the maximum explicit AF across needLR and exact gnomAD-SV.

    The maximum is conservative for rare-disease prioritization: one common
    population source is enough to trigger review. AnnotSV benign-region AFmax
    is not mixed into this exact/provisional event AF because it can describe
    an overlapping but non-equivalent allele.
    """
    values = []
    for source, col in (
        ("gnomAD-SV_exact", "GNOMAD_SV_AF"),
        ("needLR_provisional", "NEEDLR_AF"),
    ):
        value = number(row.get(col))
        if value is not None:
            values.append((source, value))

    if not values:
        return {
            "tier": "UNKNOWN",
            "priority": POPULATION_PRIORITY["UNKNOWN"],
            "max_af": ".",
            "min_af": ".",
            "sources": ".",
            "conflict": "NO",
        }

    afs = [v for _, v in values]
    maximum = max(afs)
    minimum = min(afs)
    conflict = minimum <= rare_af and maximum > max_af

    if maximum > max_af:
        tier = "TOO_COMMON_MAX_AF_GT_0.01"
    elif maximum <= rare_af:
        tier = "VERY_RARE_MAX_AF_LE_0.001"
    else:
        tier = "LOW_FREQUENCY_MAX_AF_0.001_TO_0.01"

    return {
        "tier": tier,
        "priority": POPULATION_PRIORITY[tier],
        "max_af": maximum,
        "min_af": minimum,
        "sources": ";".join(source for source, _ in values),
        "conflict": "YES" if conflict else "NO",
    }


def recurrence_tiebreak(row) -> int:
    """Late cautionary tie-break only; recurrence is not an exclusion rule."""
    for col in ("COHORT_SAMPLE_COUNT", "RECURRENT_SAMPLE_COUNT", "SAMPLE_COUNT"):
        value = number(row.get(col))
        if value is not None:
            # event_sort_tuple is sorted ascending. A smaller cohort count is
            # therefore preferred only after all major evidence axes.
            return int(value)
    return 0


def event_sort_tuple(row) -> tuple:
    gene_tier = str(row.get("GENE_RELEVANCE_TIER", "LIMITED")).upper()
    mechanism = str(row.get("INHERITANCE_MECHANISM_CLASS", "UNRESOLVED"))
    technical = str(row.get("EVENT_TECHNICAL_TIER_V2", "REVIEW"))
    population = str(row.get("EVENT_POPULATION_TIER_V2", "UNKNOWN"))

    return (
        -GENE_TIER_PRIORITY.get(gene_tier, 0),
        -MECHANISM_PRIORITY.get(mechanism, 0),
        -TECHNICAL_PRIORITY.get(technical, 0),
        -POPULATION_PRIORITY.get(population, 1),
        -constraint_tiebreak(row),
        recurrence_tiebreak(row),
        -float(number(row.get("GENE_RELEVANCE_DISPLAY_SCORE")) or 0.0),
        str(row.get("SV_ID", row.get("ID", "."))),
        str(row.get("GENES", row.get("GENE", "."))),
    )


def load_hpo_background(edges_path: str):
    parents = defaultdict(set)
    all_gene_terms = defaultdict(set)
    with open(edges_path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for edge in reader:
            subject = edge.get("subject", "")
            obj = edge.get("object", "")
            predicate = edge.get("predicate", "")
            category = edge.get("category", "")
            if (
                subject.startswith("HP:")
                and obj.startswith("HP:")
                and "subclass" in predicate.lower()
            ):
                parents[subject].add(obj)
            if (
                subject.startswith("HGNC:")
                and obj.startswith("HP:")
                and "GeneToPhenotypicFeatureAssociation" in category
            ):
                all_gene_terms[subject].add(obj)
    return parents, all_gene_terms


def load_hpo_seed_terms(path: str) -> set[str]:
    terms = set()
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            value = str(
                row.get("hpo_id")
                or row.get("HPO_ID")
                or row.get("term")
                or ""
            ).strip()
            if value.startswith("HP:"):
                terms.add(value)
    return terms


def load_gene_hpo_terms(path: str) -> dict[str, set[str]]:
    terms = defaultdict(set)
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = str(
                row.get("gene_symbol")
                or row.get("GENE")
                or row.get("gene")
                or ""
            ).strip().upper()
            value = str(
                row.get("hpo_id")
                or row.get("HPO_ID")
                or ""
            ).strip()
            if gene and value.startswith("HP:"):
                terms[gene].add(value)
    return terms


def generic_hon_semantic_map(
    phenotype_path: str,
    seed_path: str,
    edges_path: str,
) -> dict[str, float]:
    engine = semantic_engine(edges_path)
    seeds = load_hpo_seed_terms(seed_path)
    gene_terms = load_gene_hpo_terms(phenotype_path)
    return {
        gene: engine["normalized_bma"](terms, seeds)
        for gene, terms in gene_terms.items()
    }


def semantic_engine(edges_path: str):
    parents, all_gene_terms = load_hpo_background(edges_path)

    @lru_cache(maxsize=None)
    def ancestors(term):
        result = {term}
        stack = list(parents.get(term, ()))
        while stack:
            current = stack.pop()
            if current in result:
                continue
            result.add(current)
            stack.extend(parents.get(current, ()))
        return frozenset(result)

    propagated_count = defaultdict(int)
    for terms in all_gene_terms.values():
        propagated = set()
        for term in terms:
            propagated.update(ancestors(term))
        for term in propagated:
            propagated_count[term] += 1

    n_genes = max(len(all_gene_terms), 1)
    max_ic = -math.log(1 / (n_genes + 1)) if n_genes > 0 else 1.0

    def ic(term):
        return -math.log((propagated_count.get(term, 0) + 1) / (n_genes + 1))

    @lru_cache(maxsize=None)
    def pair_similarity(a, b):
        common = ancestors(a) & ancestors(b)
        if not common:
            return 0.0
        return max(ic(term) for term in common)

    def bma(left, right):
        left = tuple(sorted(set(left)))
        right = tuple(sorted(set(right)))
        if not left or not right:
            return 0.0
        left_best = [max(pair_similarity(a, b) for b in right) for a in left]
        right_best = [max(pair_similarity(b, a) for a in left) for b in right]
        return (
            sum(left_best) / len(left_best)
            + sum(right_best) / len(right_best)
        ) / 2.0

    def normalized_bma(left, right):
        if not left or not right or max_ic <= 0:
            return 0.0
        return min(1.0, bma(left, right) / max_ic)

    return {
        "parents": parents,
        "all_gene_terms": all_gene_terms,
        "max_ic": max_ic,
        "pair_similarity": pair_similarity,
        "bma": bma,
        "normalized_bma": normalized_bma,
    }
