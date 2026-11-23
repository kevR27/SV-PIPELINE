#!/usr/bin/env python3
"""Convert one private patient table into HPO and Phenopacket inputs.

The clinical hypotheses and onset estimates are intentionally NOT treated as
confirmed diagnoses or converted into phenotype codes.
"""
import argparse
import csv
import json
import re
from pathlib import Path
from datetime import datetime, timezone

HPO_PATTERN = re.compile(r"HP:\d{7}$")


def terms(value):
    ids = [item.strip() for item in re.split(r"[,;]", value or "") if item.strip()]
    bad = [item for item in ids if not HPO_PATTERN.fullmatch(item)]
    if bad:
        raise ValueError("Invalid HPO identifiers: " + ", ".join(bad))
    return sorted(set(ids))


def years(value, field):
    if not str(value or "").strip():
        return None
    number = float(value)
    if not (0 <= number <= 120):
        raise ValueError(f"{field} must be between 0 and 120 years")
    return number


def load_profile(path, sample):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise ValueError("Patient table must be tab-delimited and include sample_id")
        matches = [row for row in reader if row["sample_id"].strip() == sample]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one clinical row for {sample}, found {len(matches)}")
    return matches[0]


def make_profile(row, sample, labels):
    present = terms(row.get("hpo_present", ""))
    absent = terms(row.get("hpo_absent", ""))
    overlap = sorted(set(present) & set(absent))
    if overlap:
        raise ValueError("HPO terms cannot be both present and absent: " + ", ".join(overlap))
    if not present:
        raise ValueError(f"{sample}: no confirmed HPO terms; phenotype ranking cannot run")
    minimum = years(row.get("onset_min_years"), "onset_min_years")
    maximum = years(row.get("onset_max_years"), "onset_max_years")
    if minimum is not None and maximum is not None and minimum > maximum:
        raise ValueError("onset_min_years exceeds onset_max_years")
    onset_status = row.get("onset_status", "").strip().lower()
    if onset_status and onset_status not in {"confirmed", "estimated", "unknown"}:
        raise ValueError("onset_status must be confirmed, estimated or unknown")
    if (minimum is not None or maximum is not None) and onset_status == "unknown":
        raise ValueError("An unknown onset cannot have a numerical estimate")
    sex = row.get("sex", "").strip().upper()
    if sex and sex not in {"MALE", "FEMALE", "UNKNOWN_SEX"}:
        raise ValueError("sex must be MALE, FEMALE or UNKNOWN_SEX")

    features = []
    for hpo_id in present:
        feature = {"type": {"id": hpo_id}}
        if hpo_id in labels:
            feature["type"]["label"] = labels[hpo_id]
        features.append(feature)
    # Absent findings remain in clinical metadata, not the score-generating
    # query. Both ranking methods receive the same positive HPO terms.

    # Exomiser 15.1 accepts a GA4GH v1 Phenopacket; keep a minimal, stable
    # representation. Do not pretend that suspected diagnoses are confirmed.
    packet = {
        "id": sample,
        "subject": {"id": sample},
        "phenotypicFeatures": features,
        "metaData": {
            "created": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
            "createdBy": "SV-PIPELINE patient profile importer",
            "resources": [{
                "id": "hp",
                "name": "human phenotype ontology",
                "url": "http://purl.obolibrary.org/obo/hp.owl",
                "version": "not-specified",
                "namespacePrefix": "HP",
                "iriPrefix": "http://purl.obolibrary.org/obo/HP_"
            }],
            "phenopacketSchemaVersion": "1.0"
        }
    }
    if sex:
        packet["subject"]["sex"] = sex
    context = {
        "sample_id": sample,
        "hpo_present": present,
        "hpo_absent": absent,
        "onset_min_years": minimum,
        "onset_max_years": maximum,
        "onset_status": onset_status or "unknown",
        "diagnostic_hypothesis": row.get("diagnostic_hypothesis", "").strip(),
        "hypothesis_status": row.get("hypothesis_status", "").strip(),
        "interpretation": "Hypothesis and onset are for clinical review only, not scored as known disease"
    }
    return packet, context


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--hpo-seeds")
    parser.add_argument("--phenopacket", required=True)
    parser.add_argument("--hpo-output", required=True)
    parser.add_argument("--clinical-output", required=True)
    args = parser.parse_args()

    labels = {}
    if args.hpo_seeds:
        with open(args.hpo_seeds, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                labels[row["hpo_id"]] = row.get("hpo_label", "")
    packet, context = make_profile(load_profile(args.profiles, args.sample), args.sample, labels)
    for filename in (args.phenopacket, args.hpo_output, args.clinical_output):
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
    Path(args.phenopacket).write_text(json.dumps(packet, indent=2) + "\n", encoding="utf-8")
    with open(args.hpo_output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["hpo_id", "hpo_label"])
        for identifier in context["hpo_present"]:
            writer.writerow([identifier, labels.get(identifier, "")])
    Path(args.clinical_output).write_text(json.dumps(context, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
