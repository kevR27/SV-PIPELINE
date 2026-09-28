#!/usr/bin/env python3
"""Compare saved VEP gene annotations with AnnotSV, keeping nearby genes separate."""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from sv_evidence_common import gene_symbols, invalid_gene_labels

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vep", required=True)
    p.add_argument("--annotsv", required=True)
    p.add_argument("--records", required=True, help="Record audit from audit_sv_outputs.py")
    p.add_argument("--output-tsv", required=True)
    p.add_argument("--output-json", required=True)
    a = p.parse_args()
    with open(a.records) as f:
        records = {r["SV_ID"]: r for r in csv.DictReader(f, delimiter="\t")}
    genes, invalid = defaultdict(set), set()
    with open(a.annotsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            identifier = row.get("SV_ID", row.get("ID", "."))
            genes[identifier].update(gene_symbols(row.get("Gene_name")))
            if invalid_gene_labels(row.get("Gene_name")):
                invalid.add(identifier)
    counts, ids, pairs, commands = Counter(), set(), set(), []
    fields = ["SV_ID", "VEP_GENE", "VEP_TRANSCRIPT", "VEP_CONSEQUENCE", "RELATION",
              "ANNOTSV_RECORD_STATUS", "ANNOTSV_GENE_LIST_HAS_INVALID_LABEL"]
    Path(a.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(a.vep) as f, open(a.output_tsv, "w", newline="") as output:
        def data_lines():
            for line in f:
                if line.startswith("##"):
                    if "command-line:" in line:
                        commands.append(line.strip())
                    continue
                yield line
        writer = csv.DictWriter(output, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in csv.DictReader(data_lines(), delimiter="\t"):
            identifier = row.get("#Uploaded_variation", row.get("Uploaded_variation", "."))
            ids.add(identifier); counts["vep_rows"] += 1
            extra = dict(x.split("=", 1) if "=" in x else (x, "") for x in row.get("Extra", "").split(";"))
            symbols = gene_symbols(row.get("SYMBOL", extra.get("SYMBOL")))
            for gene in symbols:
                pairs.add((identifier, gene))
                if gene in genes[identifier]:
                    continue
                consequence = row.get("Consequence", ".")
                nearby = bool(set(consequence.split(",")) & {"upstream_gene_variant", "downstream_gene_variant"})
                status = records.get(identifier, {}).get("ANNOTSV_RECORD_STATUS", "ID_NOT_IN_MASTER")
                relation = ("NEARBY_GENE" if nearby else "ANNOTSV_RECORD_NOT_ANNOTATED" if not status.startswith("ANNOTATED")
                            else "GENE_ANNOTATION_DIFFERENCE_REQUIRES_REVIEW")
                counts[relation] += 1
                writer.writerow(dict(zip(fields, [identifier, gene, row.get("Feature", "."), consequence,
                    relation, status, "YES" if identifier in invalid else "NO"])))
    pick = any(re.search(r"(?:^|\s)--pick(?:\s|$)", c) for c in commands)
    flag_pick = any(re.search(r"(?:^|\s)--flag_pick(?:\s|$)", c) for c in commands)
    summary = {"counts": dict(counts), "master_ids": len(records), "vep_ids": len(ids), "vep_gene_pairs": len(pairs),
               "master_ids_without_vep": len(set(records) - ids), "vep_ids_not_in_master": len(ids - set(records)),
               "transcript_selection": "PICK_RESTRICTED" if pick else "FLAG_PICK" if flag_pick else "NOT_DETERMINED",
               "interpretation": "Differences can reflect nearby genes, rejected records, transcript releases or incorrect gene resources. They are not automatic evidence of lost disease genes.",
               "input_sha256": {k: hashlib.file_digest(open(v, "rb"), "sha256").hexdigest()
                                for k,v in {"vep":a.vep,"annotsv":a.annotsv,"records":a.records}.items()}}
    Path(a.output_json).parent.mkdir(parents=True, exist_ok=True)
    Path(a.output_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k:v for k,v in summary.items() if k != "input_sha256"}))


if __name__ == "__main__":
    main()
