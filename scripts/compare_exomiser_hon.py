#!/usr/bin/env python3
"""Run offline Exomiser phenotype matching and compare it with the HON method.

Three steps share this script:
  prepare-hpo: save the exact HPO query and its shared/patient scope;
  run-offline: run the installed CLI's PHENOTYPE_ONLY preset on a query;
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
import subprocess
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

from sv_evidence_common import gene_symbols

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}
QUERY_FIELDS = ["hpo_id", "hpo_label", "query_scope", "query_source"]
COMPARISON_FIELDS = [
    "GENE", "PANEL_STATUS", "QUERY_SCOPE", "QUERY_HPO_IDS", "QUERY_SOURCE",
    "CUSTOM_GENE_RANK", "CUSTOM_GENE_TIER", "CUSTOM_GENE_DISPLAY_SCORE",
    "CUSTOM_HON_PHENOTYPE_SCORE_0_10", "CUSTOM_QUERY_STATUS",
    "CUSTOM_QUERY_PHENOTYPE_SCORE_0_1", "EXOMISER_STATUS",
    "EXOMISER_BACKEND", "EXOMISER_CLI_VERSION", "EXOMISER_PHENOTYPE_DATA_VERSION",
    "EXOMISER_GENOME_DATA_VERSION",
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


def backend_resource_paths(jar, data_dir, genome_version, phenotype_version):
    """Explicit files used by CLI 15.1.0; never select the newest folder silently."""
    root = Path(data_dir)
    genome = root / f"{genome_version}_hg38"
    phenotype = root / f"{phenotype_version}_phenotype"
    prefix = f"{genome_version}_hg38"
    paths = [Path(jar), genome / f"{prefix}_transcripts_ensembl.ser",
             genome / f"{prefix}_genome.mv.db", genome / f"{prefix}_variants.mv.db",
             genome / f"{prefix}_clinvar.mv.db", phenotype / f"{phenotype_version}_phenotype.mv.db",
             phenotype / "rw_string_10.mv", phenotype / "rw_string_9_05_id2index.gz"]
    # Keep the distribution's lib folder beside its JAR. Some distributions
    # depend on these libraries; updates must invalidate the cached result too.
    paths.extend(sorted((Path(jar).parent / "lib").glob("*.jar")))
    return paths


def check_backend_files(jar, data_dir, genome_version, phenotype_version):
    """Check extraction and file identity, not scientific database integrity."""
    if Path(jar).name != "exomiser-cli-15.1.0.jar":
        raise ValueError("This offline interface is verified for exomiser-cli-15.1.0.jar. Set exomiser_cli_jar to the CLI JAR, not a REST JAR.")
    for version in (genome_version, phenotype_version):
        if not re.fullmatch(r"\d{4}", str(version)):
            raise ValueError("Exomiser data versions must be four-digit release IDs, for example '2512'.")
    paths = backend_resource_paths(jar, data_dir, genome_version, phenotype_version)
    missing = [str(path) for path in paths if not path.is_file() or not path.stat().st_size]
    if missing:
        raise ValueError("Offline Exomiser files are missing or empty. Finish downloading AND extracting the configured data releases. Missing: "
                         + "; ".join(missing) + ". See README.md: Offline Exomiser CLI.")
    if not zipfile.is_zipfile(jar):
        raise ValueError(f"The CLI JAR is not a valid ZIP/JAR file: {jar}. Check the completed download.")
    return paths


def select_hpo_query(enabled, hpo_seeds, shared_hpo=None, patient_hpo=None):
    """One selection method for both the saved query and shared CLI cache."""
    labels = {}
    if not enabled:
        return [], "EXOMISER_DISABLED", ".", labels
    if patient_hpo:
        return patient_terms(patient_hpo), "PATIENT_HPO_PROFILE", str(patient_hpo), labels
    terms = hpo_ids(shared_hpo)
    source = "config:hpo_terms"
    if not terms:
        seeds, _ = read_table(hpo_seeds)
        core = [row for row in seeds if row.get("hon_seed_role") == "CORE_OCULAR_HON"]
        terms = hpo_ids([row.get("hpo_id", "") for row in core])
        labels = {row["hpo_id"]: row.get("hpo_label", ".") for row in core}
        source = str(hpo_seeds)
    return terms, "SHARED_HON_PROFILE", source, labels


def prepare_query(args):
    terms, scope, source, labels = select_hpo_query(
        args.enabled, args.hpo_seeds, json.loads(args.shared_hpo_json), args.patient_hpo)
    if args.enabled and not terms:
        raise ValueError("Exomiser is enabled but the HPO query is empty. Set hpo_terms or use a seed file with CORE_OCULAR_HON rows.")
    rows = [{"hpo_id": term, "hpo_label": labels.get(term, "."),
             "query_scope": scope, "query_source": source} for term in terms]
    if not rows:
        rows = [{"hpo_id": ".", "hpo_label": ".", "query_scope": scope, "query_source": source}]
    write_table(args.output, rows, QUERY_FIELDS)
    target = Path(args.terms_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(",".join(terms), encoding="utf-8")
    print(f"[OK] query_scope={scope} hpo_terms={len(terms)} output={args.output}")


def write_json(path, data):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def run_offline(args):
    """Run a phenotype query once, without a VCF, invented sex or pedigree."""
    terms = hpo_ids(json.loads(args.hpo_json))
    output = Path(args.genes_output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    log = Path(args.log)
    log.parent.mkdir(parents=True, exist_ok=True)
    if not args.enabled:
        write_table(output, [], ["GENE_SYMBOL", "EXOMISER_GENE_PHENO_SCORE"])
        write_json(args.job_output, {"status": "DISABLED"})
        Path(args.properties_output).write_text("# Offline Exomiser disabled.\n", encoding="utf-8")
        write_json(args.provenance_output, {"status": "DISABLED", "backend": "OFFLINE_CLI_PHENOTYPE_ONLY", "query_hpo_ids": []})
        log.write_text("Offline Exomiser disabled in config.\n", encoding="utf-8")
        return
    if not terms:
        raise ValueError("Offline Exomiser is enabled but its HPO query is empty.")
    if args.memory_gb < 1 or args.threads < 1:
        raise ValueError("Exomiser memory_gb and threads must be positive integers.")
    resources = check_backend_files(args.jar, args.data_dir, args.genome_version, args.phenotype_version)
    try:
        java_version = subprocess.run([args.java, "-version"], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"Cannot run the configured Java executable {args.java}: {error}") from error
    version_text = java_version.stdout + java_version.stderr
    match = re.search(r'version\s+"(\d+)(?:\.(\d+))?', version_text)
    major = int(match.group(2)) if match and match.group(1) == "1" else int(match.group(1)) if match else 0
    if major < 21:
        raise ValueError(f"Offline Exomiser 15.1.0 requires Java 21 or newer; configured Java reports {major}. Use --use-conda with envs/exomiser.yaml or set exomiser_java to a Java 21 executable.")
    properties = Path(args.properties_output).resolve()
    properties.parent.mkdir(parents=True, exist_ok=True)
    # A small run-specific file avoids accidentally enabling hg19 or inheriting
    # the REST template's port placeholders. External properties are untouched.
    data_dir = str(Path(args.data_dir).resolve()).replace("\\", "\\\\")
    if any(char in data_dir for char in "\n\r"):
        raise ValueError("Exomiser data path must not contain a line break.")
    properties.write_text(
        "# Generated from pipeline config; hg38 only. No server or downloads.\n"
        f"exomiser.data-directory={data_dir}\n"
        f"exomiser.hg38.data-version={args.genome_version}\n"
        f"exomiser.phenotype.data-version={args.phenotype_version}\n", encoding="utf-8")
    job = {"sample": {"proband": "HPO_QUERY", "genomeAssembly": "hg38", "hpoIds": terms},
           "preset": "PHENOTYPE_ONLY",
           "outputOptions": {"numGenes": 0, "outputContributingVariantsOnly": False,
                             "outputDirectory": str(output.parent),
                             "outputFileName": "phenotype", "outputFormats": ["TSV_GENE"]}}
    # Sample here identifies the query; it is not a record of observed patient
    # findings. The per-sample saved query records shared versus patient scope.
    write_json(args.job_output, job)
    if output.name != "phenotype.genes.tsv":
        raise ValueError("Offline CLI genes-output must be named phenotype.genes.tsv.")
    command = [args.java, f"-Xmx{args.memory_gb}g", f"-XX:ActiveProcessorCount={args.threads}",
               "-Dspring.config.location=" + properties.as_uri(), "-jar", str(Path(args.jar).resolve()),
               "analyse", "--analysis", str(Path(args.job_output).resolve())]
    # Remove only this generated result before rerunning, so a stale successful
    # TSV cannot hide a failed CLI job. Reference data are never modified here.
    output.unlink(missing_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write(version_text + "\n")
        handle.flush()
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise ValueError(f"Offline Exomiser exited with code {result.returncode}; see {log}.")
    if not output.is_file():
        raise ValueError(f"Offline Exomiser did not create {output}; see {log}.")
    rows, fields = read_table(output)
    if not {"GENE_SYMBOL", "EXOMISER_GENE_PHENO_SCORE"}.issubset(fields):
        raise ValueError(f"Unexpected CLI TSV_GENE header in {output}; expected GENE_SYMBOL and EXOMISER_GENE_PHENO_SCORE.")
    if not any(valid_score(row.get("EXOMISER_GENE_PHENO_SCORE")) is not None for row in rows):
        raise ValueError(f"Offline CLI returned no valid phenotype scores; see {log}.")
    write_json(args.provenance_output, {
        "status": "COMPLETED", "backend": "OFFLINE_CLI_PHENOTYPE_ONLY", "cli_version": "15.1.0",
        "genome_assembly": "GRCh38", "genome_data_version": args.genome_version,
        "phenotype_data_version": args.phenotype_version, "query_hpo_ids": terms,
        "java_version": version_text.strip(), "command": command,
        "cli_jar_sha256": checksum(args.jar), "job_sha256": checksum(args.job_output),
        "properties_sha256": checksum(properties), "genes_sha256": checksum(output),
        "resources": [{"path": str(path.resolve()), "size_bytes": path.stat().st_size,
                       "mtime_ns": path.stat().st_mtime_ns} for path in resources],
        "resource_check": "Required nonempty files checked; full databases are not hashed or clinically validated.",
        "interpretation": "Phenotype-only query across reference genes. No patient variants, genotype, inheritance or diagnostic probability were analysed."})
    print(f"[OK] offline_exomiser_genes={len(rows)} query_terms={len(terms)} output={output}")


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
        # Native CLI TSV_GENE rows are already gene-specific. MOI-specific
        # duplicates remain visible; inconsistent phenotype scores are flagged.
        symbol = value(row, "GENE_SYMBOL")
        if symbol:
            genes = gene_symbols(symbol)
            if len(genes) == 1:
                grouped[genes[0].upper()].append(row)
            continue
        # Read-only compatibility for old saved AnnotSV comparisons. No REST
        # backend is started; event-wide summaries still cannot donate scores.
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
    # Only this step needs pandas/semantic scoring. Query selection can also
    # be imported by Snakemake without installing them in its own environment.
    import pandas as pd
    from ranking_common import semantic_engine
    ranking, _ = read_table(args.ranking)
    score_file = args.exomiser_genes or args.annotsv
    annotation, annotation_fields = read_table(score_file)
    query, _ = read_table(args.query)
    scopes = {row["query_scope"] for row in query}
    if len(scopes) != 1:
        raise ValueError("Saved HPO query must describe exactly one profile.")
    scope = scopes.pop()
    terms = hpo_ids([row["hpo_id"] for row in query if row["hpo_id"] != "."])
    enabled = scope != "EXOMISER_DISABLED"
    if enabled and not terms:
        raise ValueError("Saved Exomiser HPO query is empty.")
    provenance = json.loads(Path(args.backend_provenance).read_text()) if args.backend_provenance else {}
    if args.exomiser_genes and not args.backend_provenance:
        raise ValueError("CLI comparison requires --backend-provenance to identify the query, versions and completed run.")
    if args.exomiser_genes and provenance.get("backend") != "OFFLINE_CLI_PHENOTYPE_ONLY":
        raise ValueError("CLI comparison expects an OFFLINE_CLI_PHENOTYPE_ONLY result, not a combined variant analysis.")
    if provenance:
        expected_status = "COMPLETED" if enabled else "DISABLED"
        if provenance.get("status") != expected_status or provenance.get("query_hpo_ids") != terms:
            raise ValueError("Offline CLI result does not match the saved HPO query or enablement state. Rerun exomiser_comparison.")
        if enabled and provenance.get("genes_sha256") != checksum(score_file):
            raise ValueError("CLI gene TSV checksum differs from the completed run. Rerun exomiser_comparison.")
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
                    "CUSTOM_HON_PHENOTYPE_SCORE_0_10": original.get("HON_SEMANTIC_SCORE_0_10", "."),
                    "EXOMISER_BACKEND": provenance.get("backend", "LEGACY_ANNOTSV_IMPORT"),
                    "EXOMISER_CLI_VERSION": provenance.get("cli_version", "."),
                    "EXOMISER_PHENOTYPE_DATA_VERSION": provenance.get("phenotype_data_version", "."),
                    "EXOMISER_GENOME_DATA_VERSION": provenance.get("genome_data_version", ".")})
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
        raise ValueError("Exomiser was requested but no valid gene-specific scores were returned. Check the offline CLI log, Java 21+, extracted data versions and gene symbols. Negative scores such as -1 are unavailable, not zero. See README.md: Offline Exomiser CLI.")
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
                                {"query": args.query, "exomiser_gene_scores": score_file, "original_ranking": args.ranking,
                                 "human_phenotypes": args.phenotypes, "monarch_edges": args.edges}.items()},
               "backend_provenance": provenance,
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
    prepare.add_argument("--output", required=True)
    prepare.add_argument("--terms-output", required=True)
    offline = commands.add_parser("run-offline", help="Run the installed CLI's phenotype-only preset; never download")
    offline.add_argument("--enabled", action="store_true")
    offline.add_argument("--hpo-json", required=True)
    offline.add_argument("--java", default="java", help="Java 21+ executable, supplied by the dedicated conda environment")
    offline.add_argument("--jar", required=True)
    offline.add_argument("--data-dir", required=True)
    offline.add_argument("--genome-version", required=True)
    offline.add_argument("--phenotype-version", required=True)
    offline.add_argument("--memory-gb", type=int, default=8)
    offline.add_argument("--threads", type=int, default=1)
    for name in ["genes-output", "job-output", "properties-output", "provenance-output", "log"]:
        offline.add_argument("--" + name, required=True)
    comparison = commands.add_parser("compare", help="Write a separate comparison; preserve the original ranking")
    for name in ["ranking", "query", "phenotypes", "edges", "output", "summary-output"]:
        comparison.add_argument("--" + name, required=True)
    sources = comparison.add_mutually_exclusive_group(required=True)
    sources.add_argument("--exomiser-genes", help="Native CLI TSV_GENE result")
    sources.add_argument("--annotsv", help="Read-only import of an older saved AnnotSV result")
    comparison.add_argument("--backend-provenance", help="Completed offline CLI run JSON")
    comparison.add_argument("--require-exomiser", action="store_true")
    args = parser.parse_args()
    try:
        {"prepare-hpo": prepare_query, "run-offline": run_offline, "compare": compare}[args.command](args)
    except (ValueError, OSError) as error:
        parser.exit(2, f"[ERROR] {error}\n")


if __name__ == "__main__":
    main()
