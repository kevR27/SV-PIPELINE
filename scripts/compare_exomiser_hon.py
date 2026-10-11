#!/usr/bin/env python3
"""Compare AnnotSV's Exomiser phenotype scores with the custom HON method.

Two steps share this script:
  prepare-hpo: save the exact HPO query that AnnotSV will receive;
  compare: score the same query with the custom semantic method and compare ranks.

The original gene ranking is read, never rewritten. These are gene–phenotype
scores, not probabilities that an SV causes disease.
"""
from __future__ import annotations

import argparse
import csv
import json
import hashlib
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from ranking_common import semantic_engine
from sv_evidence_common import gene_symbols

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}
QUERY_FIELDS = ["hpo_id", "hpo_label", "query_scope", "query_source"]
COMPARISON_FIELDS = [
    "GENE", "PANEL_STATUS", "QUERY_SCOPE", "QUERY_HPO_IDS", "QUERY_SOURCE",
    "CUSTOM_GENE_RANK", "CUSTOM_GENE_TIER", "CUSTOM_GENE_DISPLAY_SCORE",
    "CUSTOM_HON_PHENOTYPE_SCORE_0_10", "CUSTOM_QUERY_STATUS",
    "CUSTOM_QUERY_PHENOTYPE_SCORE_0_1", "EXOMISER_STATUS",
    "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_RAW_SCORES",
    "EXOMISER_GENE_SPECIFIC_ROWS", "EXOMISER_HUMAN_PHENO_EVIDENCE",
    "EXOMISER_MOUSE_PHENO_EVIDENCE", "EXOMISER_FISH_PHENO_EVIDENCE",
    "COMPARISON_STATUS", "CUSTOM_QUERY_RANK", "EXOMISER_QUERY_RANK",
    "RANK_DIFFERENCE_EXOMISER_MINUS_CUSTOM",
]


def read_table(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_table(path, rows, fields):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hpo_ids(values):
    """Accept comma/space separated IDs or a YAML list; reject mistyped IDs."""
    if values is None:
        return []
    if isinstance(values, str):
        values = re.split(r"[,;\s]+", values.strip())
    if not isinstance(values, list):
        raise ValueError("hpo_terms must be a list of HPO IDs or a comma-separated string.")
    cleaned = sorted({str(value).strip() for value in values if str(value).strip()})
    bad = [value for value in cleaned if not re.fullmatch(r"HP:\d{7}", value)]
    if bad:
        raise ValueError(f"Invalid HPO IDs: {', '.join(bad)}. Expected HP: followed by seven digits.")
    return cleaned


def patient_terms(path):
    """A patient file is a complete query, overriding the shared profile."""
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    first = lines[0] if lines else ""
    if first.split("\t")[0].strip().lower() in {"hpo_id", "hpo", "term"}:
        reader = csv.DictReader(lines, delimiter="\t")
        rows, fields = list(reader), list(reader.fieldnames or [])
        column = next((field for field in fields if field.lower() in {"hpo_id", "hpo", "term"}), None)
        if column is None:
            raise ValueError(f"Patient HPO table needs hpo_id or HPO_ID: {path}")
        return hpo_ids([row[column] for row in rows])
    return hpo_ids([line.split("\t")[0].strip() for line in lines])


def check_backend_files(directory):
    """Check local files only; installation stays in the README."""
    root = Path(directory)
    jars = list((root / "jar").glob("exomiser-rest-prioritiser-*.jar"))
    expected = "exomiser-rest-prioritiser-14.1.0.jar"
    if len(jars) != 1 or jars[0].name != expected or not jars[0].stat().st_size:
        raise ValueError(f"AnnotSV 3.5.10 needs one compatible REST JAR: {root / 'jar' / expected}. See README.md: Manual Exomiser setup.")
    releases = [path for path in (root / "Annotations_Exomiser").glob("*") if path.is_dir() and path.name.isdigit()]
    if not releases:
        raise ValueError(f"No Exomiser phenotype-data release found under {root / 'Annotations_Exomiser'}. See README.md: Manual Exomiser setup.")
    selected = max(releases, key=lambda path: int(path.name))
    databases = list(selected.glob("**/*phenotype*.mv.db"))
    if selected.name != "2406" or not any(path.stat().st_size for path in databases):
        raise ValueError(f"AnnotSV selects {selected}, but this pinned comparison requires unzipped phenotype data 2406. Keep other data releases in a separate bundle. See README.md: Manual Exomiser setup.")


def prepare_query(args):
    labels = {}
    if not args.enabled:
        terms, scope, source = [], "EXOMISER_DISABLED", "."
    elif args.patient_hpo:
        terms = patient_terms(args.patient_hpo)
        scope, source = "PATIENT_HPO_PROFILE", args.patient_hpo
    else:
        terms = hpo_ids(json.loads(args.shared_hpo_json))
        scope, source = "SHARED_HON_PROFILE", "config:hpo_terms"
        if not terms:
            seeds, _ = read_table(args.hpo_seeds)
            core = [row for row in seeds if row.get("hon_seed_role") == "CORE_OCULAR_HON"]
            terms = hpo_ids([row.get("hpo_id", "") for row in core])
            labels = {row["hpo_id"]: row.get("hpo_label", ".") for row in core}
            source = args.hpo_seeds
    if args.enabled and not terms:
        raise ValueError("Exomiser is enabled but the HPO query is empty. Set hpo_terms or use a seed file with CORE_OCULAR_HON rows.")
    if args.enabled and args.annotations_dir:
        check_backend_files(args.annotations_dir)
    rows = [{"hpo_id": term, "hpo_label": labels.get(term, "."),
             "query_scope": scope, "query_source": source} for term in terms]
    if not rows:
        rows = [{"hpo_id": ".", "hpo_label": ".", "query_scope": scope, "query_source": source}]
    write_table(args.output, rows, QUERY_FIELDS)
    target = Path(args.terms_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(",".join(terms), encoding="utf-8")
    print(f"[OK] query_scope={scope} hpo_terms={len(terms)} output={args.output}")


def value(row, name):
    """Column spelling is case-insensitive, but gene-specific scope is strict."""
    return next((str(v or "").strip() for key, v in row.items()
                 if key.lower() == name.lower()), "")


def valid_score(raw):
    try:
        number = float(raw)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and 0 <= number <= 1 else None


def exomiser_by_gene(rows):
    grouped = defaultdict(list)
    for row in rows:
        genes = gene_symbols(value(row, "Gene_name"))
        # A full row may summarize the best gene in a multi-gene SV. It must
        # never donate that score to every gene covered by the event.
        if value(row, "Annotation_mode").lower() == "split" and len(genes) == 1:
            grouped[genes[0].upper()].append(row)
    return grouped


def extract_exomiser(rows, has_column, enabled):
    result = {"EXOMISER_GENE_PHENO_SCORE": ".", "EXOMISER_RAW_SCORES": ".",
              "EXOMISER_GENE_SPECIFIC_ROWS": len(rows)}
    for organism in ["HUMAN", "MOUSE", "FISH"]:
        evidence = sorted({value(row, organism.title() + "_pheno_evidence") for row in rows}
                          - MISSING)
        result[f"EXOMISER_{organism}_PHENO_EVIDENCE"] = " | ".join(evidence) or "."
    raw = [value(row, "Exomiser_gene_pheno_score") for row in rows]
    result["EXOMISER_RAW_SCORES"] = ";".join(sorted(set(raw) - {""})) or "."
    scores = set()
    for item in raw:
        score = valid_score(item)
        if score is not None:
            scores.add(score)
    if not enabled:
        status = "DISABLED"
    elif not has_column:
        status = "COLUMN_MISSING"
    elif not rows:
        status = "NO_GENE_SPECIFIC_ROW"
    elif not scores:
        status = "NO_VALID_SCORE"
    elif len(scores) > 1:
        status = "CONFLICTING_GENE_SCORES"
    else:
        result["EXOMISER_GENE_PHENO_SCORE"] = next(iter(scores))
        status = "AVAILABLE" if all(valid_score(item) is not None for item in raw) else "AVAILABLE_WITH_INVALID_ROWS"
    result["EXOMISER_STATUS"] = status
    return result


def compare(args):
    ranking, _ = read_table(args.ranking)
    annotation, annotation_fields = read_table(args.annotsv)
    query, _ = read_table(args.query)
    scopes = {row["query_scope"] for row in query}
    if len(scopes) != 1:
        raise ValueError("Saved HPO query must describe exactly one profile.")
    scope = scopes.pop()
    terms = hpo_ids([row["hpo_id"] for row in query if row["hpo_id"] != "."])
    enabled = scope != "EXOMISER_DISABLED"
    if enabled and not terms:
        raise ValueError("Saved Exomiser HPO query is empty.")
    engine = semantic_engine(args.edges) if enabled else None
    known_terms = set()
    if engine:
        known_terms.update(engine["parents"])
        for parents in engine["parents"].values():
            known_terms.update(parents)
        for background_terms in engine["all_gene_terms"].values():
            known_terms.update(background_terms)
    missing_query_terms = sorted(set(terms) - known_terms) if enabled else []
    phenotypes, _ = read_table(args.phenotypes)
    gene_terms = defaultdict(set)
    for row in phenotypes:
        if str(row.get("hpo_id", "")).startswith("HP:"):
            gene_terms[str(row.get("gene_symbol", "")).upper()].add(row["hpo_id"])
    grouped = exomiser_by_gene(annotation)
    has_score = any(field.lower() == "exomiser_gene_pheno_score" for field in annotation_fields)
    output, seen = [], set()
    for original in ranking:
        gene = str(original.get("gene", original.get("GENE", ""))).strip().upper()
        if not gene or gene in seen:
            raise ValueError("The original gene ranking must contain one row per gene.")
        seen.add(gene)
        row = dict.fromkeys(COMPARISON_FIELDS, ".")
        row.update({"GENE": gene, "PANEL_STATUS": original.get("PANEL_STATUS", "."),
                    "QUERY_SCOPE": scope, "QUERY_HPO_IDS": ";".join(terms) or ".",
                    "QUERY_SOURCE": " | ".join(sorted({r["query_source"] for r in query})),
                    "CUSTOM_GENE_RANK": original.get("GENE_RANK_GLOBAL", "."),
                    "CUSTOM_GENE_TIER": original.get("GENE_RELEVANCE_TIER", "."),
                    "CUSTOM_GENE_DISPLAY_SCORE": original.get("GENE_RELEVANCE_DISPLAY_SCORE", "."),
                    "CUSTOM_HON_PHENOTYPE_SCORE_0_10": original.get("HON_SEMANTIC_SCORE_0_10", ".")})
        if not enabled:
            row["CUSTOM_QUERY_STATUS"] = "DISABLED"
        elif missing_query_terms:
            row["CUSTOM_QUERY_STATUS"] = "QUERY_TERMS_MISSING_FROM_MONARCH"
        elif sum(engine["pair_similarity"](term, term) for term in terms) <= 0:
            row["CUSTOM_QUERY_STATUS"] = "QUERY_HAS_NO_INFORMATION_CONTENT"
        elif not gene_terms[gene]:
            row["CUSTOM_QUERY_STATUS"] = "NO_HUMAN_HPO_ANNOTATIONS"
        else:
            row["CUSTOM_QUERY_STATUS"] = "AVAILABLE"
            row["CUSTOM_QUERY_PHENOTYPE_SCORE_0_1"] = round(
                engine["normalized_query_coverage"](terms, gene_terms[gene]), 6)
        row.update(extract_exomiser(grouped[gene], has_score, enabled))
        row["COMPARISON_STATUS"] = "PAIRED" if (
            row["CUSTOM_QUERY_STATUS"] == "AVAILABLE"
            and row["EXOMISER_STATUS"].startswith("AVAILABLE")) else "NOT_COMPARABLE"
        output.append(row)
    scored = sum(row["EXOMISER_STATUS"].startswith("AVAILABLE") for row in output)
    if args.require_exomiser and enabled and output and not scored:
        raise ValueError("Exomiser was requested but no valid gene-specific scores were returned. Check the AnnotSV log, Java 17+, phenotype database, REST JAR, application.properties and gene-ID mapping. Negative scores such as -1 are unavailable, not zero. See README.md: Exomiser through AnnotSV.")
    paired = pd.DataFrame([row for row in output if row["COMPARISON_STATUS"] == "PAIRED"])
    statistics = {"paired_genes": len(paired), "spearman_rank_correlation": None, "top_gene_overlap": {}}
    if not paired.empty:
        custom = pd.to_numeric(paired["CUSTOM_QUERY_PHENOTYPE_SCORE_0_1"])
        exomiser = pd.to_numeric(paired["EXOMISER_GENE_PHENO_SCORE"])
        custom_ranks = custom.rank(ascending=False, method="min")
        exomiser_ranks = exomiser.rank(ascending=False, method="min")
        ranks = {gene: (c, e) for gene, c, e in zip(paired.GENE, custom_ranks, exomiser_ranks)}
        for row in output:
            if row["GENE"] in ranks:
                c, e = ranks[row["GENE"]]
                row.update({"CUSTOM_QUERY_RANK": c, "EXOMISER_QUERY_RANK": e,
                            "RANK_DIFFERENCE_EXOMISER_MINUS_CUSTOM": e - c})
        if len(paired) >= 3 and custom.nunique() > 1 and exomiser.nunique() > 1:
            correlation = custom.rank(method="average").corr(exomiser.rank(method="average"))
            statistics["spearman_rank_correlation"] = round(float(correlation), 6)
        for k in [10, 20]:
            cset = set(paired.loc[custom_ranks <= k, "GENE"])
            eset = set(paired.loc[exomiser_ranks <= k, "GENE"])
            statistics["top_gene_overlap"][str(k)] = {
                "custom_genes": len(cset), "exomiser_genes": len(eset),
                "shared_genes": sorted(cset & eset),
                "jaccard": round(len(cset & eset) / len(cset | eset), 6)}
    summary = {"query_scope": scope, "query_hpo_ids": terms,
               "source_files": {name: {"path": str(path), "sha256": checksum(path)} for name, path in
                                {"query": args.query, "annotsv_exomiser": args.annotsv, "original_ranking": args.ranking,
                                 "human_phenotypes": args.phenotypes, "monarch_edges": args.edges}.items()},
               "candidate_genes": len(output), "exomiser_scored_genes": scored,
               "custom_query_missing_from_monarch": missing_query_terms,
               "exomiser_status_counts": {str(k): int(v) for k, v in pd.Series([r["EXOMISER_STATUS"] for r in output], dtype=str).value_counts().items()},
               **statistics, "rank_definition": "Descending phenotype score among the same paired genes; ties share the best rank and can enlarge top-k sets.",
               "interpretation": "Comparison of gene–phenotype prioritization only. Original gene ranks are unchanged. Scores are not averaged. Agreement is not independent validation or proof of SV pathogenicity."}
    write_table(args.output, output, COMPARISON_FIELDS)
    target = Path(args.summary_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(f"[OK] genes={len(output)} exomiser_scored={scored} paired={len(paired)} output={args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare-hpo", help="Save the shared HON or patient HPO query")
    prepare.add_argument("--enabled", action="store_true")
    prepare.add_argument("--hpo-seeds", required=True)
    prepare.add_argument("--shared-hpo-json", default="null")
    prepare.add_argument("--patient-hpo")
    prepare.add_argument("--annotations-dir", help="Check compatible local Exomiser files when enabled; never download")
    prepare.add_argument("--output", required=True)
    prepare.add_argument("--terms-output", required=True)
    comparison = commands.add_parser("compare", help="Write a separate comparison; preserve the original ranking")
    for name in ["ranking", "annotsv", "query", "phenotypes", "edges", "output", "summary-output"]:
        comparison.add_argument("--" + name, required=True)
    comparison.add_argument("--require-exomiser", action="store_true")
    args = parser.parse_args()
    try:
        prepare_query(args) if args.command == "prepare-hpo" else compare(args)
    except ValueError as error:
        parser.exit(2, f"[ERROR] {error}\n")


if __name__ == "__main__":
    main()
