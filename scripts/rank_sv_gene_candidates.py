#!/usr/bin/env python3
"""Prioritize genes intersected by the genome-wide SV callset.

The primary ranking is tier-based rather than a fine-grained additive score.
Generic hereditary-optic-neuropathy (HON) phenotype relevance is calculated
with IC-weighted Resnik best-match-average similarity to the HON seed terms.
Curated gene-disease evidence is kept separate and down-weighted when the gene
has little HON phenotype similarity.

Panel membership is reported but is not used as a discovery score. Panel and
non-panel ranks are therefore directly comparable within their own groups.

This is research prioritization, not variant pathogenicity classification.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from ranking_common import (
    MISSING,
    combine_hon_semantic,
    gene_relevance,
    inheritance_class,
    known,
    parse_moi,
    semantic_engine,
    summarize_gencc,
)
from sv_evidence_common import INVALID_GENE_LABELS, gene_symbols

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


def read_list(path: str) -> set[str]:
    out = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            value = line.strip()
            if value and not value.startswith("#"):
                out.add(value.upper())
    return out


def read_hpo_seeds(path: str) -> dict[str, set[str]]:
    out = {"core": set(), "context": set()}
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="	")
        for row in reader:
            hp = str(row.get("hpo_id", "")).strip()
            if not hp.startswith("HP:"):
                continue
            role = str(row.get("hon_seed_role", "") or "").upper()
            if "MITO_SYNDROMIC_CONTEXT" in role:
                out["context"].add(hp)
            else:
                out["core"].add(hp)
    return out


def first(row: dict, names: list[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if known(value):
            return str(value)
    lower = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        value = lower.get(name.lower(), "")
        if known(value):
            return str(value)
    return ""


def is_positive_flag(value: str) -> bool:
    return str(value or "").strip().upper() in {"YES", "TRUE", "1", "Y"}


def first_numeric_text(values: list[str]) -> str:
    cleaned = sorted({str(v).strip() for v in values if known(v)})
    return ";".join(cleaned) if cleaned else "."


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotsv", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--phenotypes", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--hpo-seeds", required=True)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    genes = read_list(args.genes) - INVALID_GENE_LABELS
    panel = read_list(args.panel)
    hon_seed_groups = read_hpo_seeds(args.hpo_seeds)
    semantic = semantic_engine(args.edges)

    pheno = defaultdict(
        lambda: {
            "hpos": set(),
            "anchors": set(),
            "sources": set(),
        }
    )
    with open(args.phenotypes, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="	"):
            gene = str(row.get("gene_symbol", "")).upper().strip()
            if gene not in genes:
                continue
            hpo_id = str(row.get("hpo_id", "")).strip()
            if hpo_id.startswith("HP:"):
                pheno[gene]["hpos"].add(hpo_id)
            if (
                hpo_id.startswith("HP:")
                and str(row.get("optic_neuropathy_anchor", "0")) == "1"
            ):
                pheno[gene]["anchors"].add(hpo_id)
            source = str(row.get("source", "")).strip()
            if source:
                pheno[gene]["sources"].add(source)

    annotsv = defaultdict(
        lambda: {
            "sv_ids": set(),
            "ranking_score": [],
            "ranking_criteria": [],
            "acmg_class": [],
            "omim": [],
            "omim_morbid": [],
            "omim_inheritance": [],
            "gencc_classification": [],
            "gencc_disease": [],
            "gencc_moi": [],
            "clinvar": [],
            "constraint": [],
            "hi": [],
            "ts": [],
            "pli": [],
            "loeuf": [],
            "loeuf_bin": [],
            "sv_records": {},
        }
    )

    with open(args.annotsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="	")
        for row in reader:
            gene_field = first(
                row,
                ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL"],
            )
            if not gene_field:
                continue

            row_genes = set(gene_symbols(gene_field))
            sv_id = first(row, ["SV_ID", "ID", "AnnotSV_ID"])

            for gene in row_genes:
                if gene not in genes:
                    continue

                data = annotsv[gene]
                if sv_id:
                    data["sv_ids"].add(sv_id)
                    info_text = first(row, ["INFO"])
                    info_svtype = ""
                    if info_text:
                        match = re.search(
                            r"(?:^|;)SVTYPE=([^;]+)",
                            info_text,
                            flags=re.IGNORECASE,
                        )
                        if match:
                            info_svtype = match.group(1).upper()
                    svtype = (
                        info_svtype
                        or first(row, ["SVTYPE", "SV_type", "svtype"]).upper()
                    )
                    raw_svlen = first(row, ["SV_length", "SVLEN", "svlen"])
                    try:
                        abs_svlen = abs(float(raw_svlen))
                    except (TypeError, ValueError):
                        abs_svlen = None

                    previous = data["sv_records"].get(sv_id)
                    current = {"svtype": svtype or ".", "abs_svlen": abs_svlen}
                    if (
                        previous is None
                        or (
                            previous.get("abs_svlen") is None
                            and abs_svlen is not None
                        )
                    ):
                        data["sv_records"][sv_id] = current

                for key, aliases in {
                    "ranking_score": [
                        "AnnotSV_ranking_score", "AnnotSV ranking score"
                    ],
                    "ranking_criteria": [
                        "AnnotSV_ranking_criteria", "AnnotSV ranking criteria"
                    ],
                    "acmg_class": ["ACMG_class", "ACMG class"],
                }.items():
                    value = first(row, aliases)
                    if value:
                        data[key].append(value)

                # Gene-specific evidence must come from the gene-specific rows.
                if (
                    str(row.get("Annotation_mode", "")).lower() == "full"
                    or len(row_genes) != 1
                ):
                    continue

                aliases = {
                    "omim": [
                        "OMIM_phenotype", "OMIM",
                        "AnnotSV_OMIM", "AnnotSV_OMIM_evidence",
                    ],
                    "omim_morbid": [
                        "OMIM_morbid", "OMIM_morbid_candidate"
                    ],
                    "omim_inheritance": [
                        "OMIM_inheritance", "OMIM inheritance"
                    ],
                    "gencc_classification": [
                        "GenCC_classification", "GENCC_classification",
                        "GENCC", "GenCC",
                    ],
                    "gencc_disease": [
                        "GenCC_disease", "GENCC_disease"
                    ],
                    "gencc_moi": ["GenCC_moi", "GENCC_moi"],
                    "clinvar": ["ClinVar", "ClinVar_SV"],
                    "constraint": [
                        "HI", "TS", "pLI", "LOEUF", "LOEUF_bin",
                        "GnomAD_pLI",
                    ],
                    "hi": ["HI"],
                    "ts": ["TS"],
                    "pli": ["pLI", "GnomAD_pLI"],
                    "loeuf": ["LOEUF", "LOEUF_score"],
                    "loeuf_bin": ["LOEUF_bin"],
                }
                for key, names in aliases.items():
                    value = first(row, names)
                    if value:
                        data[key].append(value)

    rows = []
    for gene in sorted(genes):
        p = pheno[gene]
        a = annotsv[gene]

        hon_core_similarity = semantic["normalized_bma"](
            p["hpos"],
            hon_seed_groups["core"],
        )
        hon_context_similarity = (
            semantic["normalized_bma"](
                p["hpos"],
                hon_seed_groups["context"],
            )
            if hon_seed_groups["context"]
            else 0.0
        )
        hon_similarity = combine_hon_semantic(
            hon_core_similarity,
            hon_context_similarity,
        )
        phenotype_score_0_10 = 10.0 * hon_similarity

        has_omim = bool(a["omim"]) or any(
            is_positive_flag(value) for value in a["omim_morbid"]
        )
        gencc = summarize_gencc(
            a["gencc_classification"],
            has_omim=has_omim,
            hon_semantic_normalized=hon_similarity,
        )
        relevance = gene_relevance(
            hon_similarity,
            float(gencc["adjusted_score"]),
        )

        moi_set = parse_moi(
            ";".join(a["gencc_moi"]),
            ";".join(a["omim_inheritance"]),
        )
        inheritance = inheritance_class(moi_set)

        sv_records = list(a["sv_records"].values())
        sized = [
            record
            for record in sv_records
            if record.get("abs_svlen") is not None
        ]
        sv_types = sorted({
            record["svtype"]
            for record in sv_records
            if record.get("svtype") not in ("", ".")
        })

        sv_count_lt100kb = sum(r["abs_svlen"] < 100_000 for r in sized)
        sv_count_100kb_1mb = sum(
            100_000 <= r["abs_svlen"] < 1_000_000 for r in sized
        )
        sv_count_1mb_10mb = sum(
            1_000_000 <= r["abs_svlen"] < 10_000_000 for r in sized
        )
        sv_count_ge10mb = sum(r["abs_svlen"] >= 10_000_000 for r in sized)
        breakpoint_defined_count = sum(
            r.get("svtype") in {"INV", "BND", "TRA"} for r in sv_records
        )
        max_sv_size_bp = max(
            (r["abs_svlen"] for r in sized),
            default=None,
        )

        has_disease_evidence = float(gencc["adjusted_score"]) > 0
        has_hon_context = hon_similarity > 0

        if gene in panel:
            candidate_group = "PANEL_GENE"
        elif has_hon_context and has_disease_evidence:
            candidate_group = "NONPANEL_HPO_AND_DISEASE_EVIDENCE"
        elif has_hon_context:
            candidate_group = "NONPANEL_HPO_OVERLAP"
        elif has_disease_evidence:
            candidate_group = "NONPANEL_HUMAN_DISEASE_GENE"
        else:
            candidate_group = "OTHER_NONPANEL_CANDIDATE"

        rows.append({
            "gene": gene,
            "panel_gene": "YES" if gene in panel else "NO",
            "PANEL_STATUS": "PANEL_GENE" if gene in panel else "NONPANEL_GENE",
            "SV_count": len(a["sv_ids"]),
            "SV_types": ";".join(sv_types) if sv_types else ".",
            "SV_count_lt100kb": sv_count_lt100kb,
            "SV_count_100kb_to_1Mb": sv_count_100kb_1mb,
            "SV_count_1Mb_to_10Mb": sv_count_1mb_10mb,
            "SV_count_ge10Mb": sv_count_ge10mb,
            "breakpoint_defined_INV_BND_count": breakpoint_defined_count,
            "max_SV_size_bp": (
                int(max_sv_size_bp) if max_sv_size_bp is not None else "."
            ),
            "human_HPO_count": len(p["hpos"]),
            "optic_neuropathy_anchor_HPO_count": len(p["anchors"]),
            "HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED": round(
                hon_core_similarity, 6
            ),
            "HON_MITO_SYNDROMIC_CONTEXT_SIMILARITY_NORMALIZED": round(
                hon_context_similarity, 6
            ),
            "HON_SEMANTIC_SIMILARITY_NORMALIZED": round(hon_similarity, 6),
            "HON_SEMANTIC_SCORE_0_10": round(phenotype_score_0_10, 3),
            "phenotype_score": round(phenotype_score_0_10, 3),
            "phenotype_score_scope": (
                "GENERIC_HON_RESNIK_BMA_NOT_PATIENT_SPECIFIC"
            ),
            "gene_disease_evidence_raw_score": round(
                float(gencc["raw_score"]), 3
            ),
            "gene_disease_evidence_score": round(
                float(gencc["adjusted_score"]), 3
            ),
            "gene_disease_evidence_level": gencc["best"],
            "gene_disease_evidence_source": gencc["source"],
            "gene_disease_evidence_conflict": gencc["conflict"],
            "GENE_DISEASE_HON_CONTEXT_FACTOR": gencc["hon_context_factor"],
            "GENE_DISEASE_CONTEXT_SCOPE": gencc["scope"],
            "GENE_RELEVANCE_TIER": relevance["tier"],
            "GENE_RELEVANCE_TIER_RANK": relevance["tier_rank"],
            "GENE_RELEVANCE_DISPLAY_SCORE": relevance["display_score"],
            # Backward-compatible numeric column. It is a display/tie-break
            # score only; tier is the primary ranking axis.
            "integrated_discovery_score": relevance["display_score"],
            "integrated_discovery_score_scope": (
                "BALANCED_HON_SEMANTIC_AND_CURATED_DISEASE_DISPLAY_SCORE_"
                "PRIMARY_ORDER_USES_GENE_RELEVANCE_TIER"
            ),
            "GenCC_classifications": gencc["all_terms"],
            "GenCC_disease": ";".join(sorted(set(a["gencc_disease"]))) or ".",
            "GenCC_moi": ";".join(sorted(set(a["gencc_moi"]))) or ".",
            "OMIM_inheritance": (
                ";".join(sorted(set(a["omim_inheritance"]))) or "."
            ),
            "GENE_MOI_SET": ";".join(sorted(moi_set)) if moi_set else ".",
            "GENE_INHERITANCE_CLASS": inheritance,
            "ClinGen_HI": first_numeric_text(a["hi"]),
            "ClinGen_TS": first_numeric_text(a["ts"]),
            "gnomAD_pLI": first_numeric_text(a["pli"]),
            "gnomAD_LOEUF": first_numeric_text(a["loeuf"]),
            "gnomAD_LOEUF_bin": first_numeric_text(a["loeuf_bin"]),
            "SV_evidence_score": 0.0,
            "AnnotSV_ranking_scores": (
                ";".join(sorted(set(a["ranking_score"]))) or "."
            ),
            "AnnotSV_ranking_criteria": (
                ";".join(sorted(set(a["ranking_criteria"]))) or "."
            ),
            "AnnotSV_ACMG_classes": (
                ";".join(sorted(set(a["acmg_class"]))) or "."
            ),
            "AnnotSV_classification_scope": (
                "OVERLAPPING_SV_RECORDS_NOT_GENE_PATHOGENICITY"
            ),
            "AnnotSV_OMIM_evidence": (
                ";".join(sorted(set(a["omim"]))) or "."
            ),
            "AnnotSV_GENCC_evidence": (
                ";".join(sorted(set(a["gencc_classification"]))) or "."
            ),
            "AnnotSV_ClinVar_evidence": (
                ";".join(sorted(set(a["clinvar"]))) or "."
            ),
            "AnnotSV_constraint_evidence": (
                ";".join(sorted(set(a["constraint"]))) or "."
            ),
            "CANDIDATE_CLASS": candidate_group,
            "OMIM": ";".join(sorted(set(a["omim"]))) or ".",
            "GENCC": gencc["all_terms"],
            "ANNOTSV_GENERAL_CLASSIFICATION": (
                ";".join(sorted(set(a["acmg_class"]))) or "."
            ),
            "candidate_group": candidate_group,
            "classification": candidate_group,
            "ranking_model": (
                "HON_corePlusCappedMitoContext_resnikBMA__balancedDiseaseEvidence__tiered_v6"
            ),
            "interpretation": (
                "Research prioritization only. Primary gene ordering uses a "
                "discrete relevance tier derived from continuous HON semantic "
                "similarity and curated gene-disease evidence. Panel membership, "
                "SV count and event size do not add gene-relevance points. "
                "Inheritance, dosage mechanism, genotype, population frequency "
                "and technical support are evaluated at the SV-gene event level."
            ),
        })

    tier_order = {"HIGH": 3, "MODERATE": 2, "SUPPORTING": 1, "LIMITED": 0}
    rows.sort(
        key=lambda row: (
            -tier_order.get(row["GENE_RELEVANCE_TIER"], 0),
            -float(row["GENE_RELEVANCE_DISPLAY_SCORE"]),
            -float(row["HON_SEMANTIC_SIMILARITY_NORMALIZED"]),
            -float(row["gene_disease_evidence_score"]),
            row["gene"],
        )
    )

    counters = defaultdict(int)
    for index, row in enumerate(rows, start=1):
        row["GENE_RANK_GLOBAL"] = index
        group = row["PANEL_STATUS"]
        counters[group] += 1
        row["GENE_RANK_WITHIN_PANEL_STATUS"] = counters[group]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else [
        "gene", "panel_gene", "PANEL_STATUS", "GENE_RELEVANCE_TIER",
        "GENE_RELEVANCE_DISPLAY_SCORE", "GENE_INHERITANCE_CLASS",
        "candidate_group", "interpretation",
    ]
    with output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="	",
            lineterminator="
",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[OK] ranked_genes={len(rows)} "
        f"panel={sum(r['PANEL_STATUS'] == 'PANEL_GENE' for r in rows)} "
        f"nonpanel={sum(r['PANEL_STATUS'] == 'NONPANEL_GENE' for r in rows)} "
        f"output={output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
