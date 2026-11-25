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


def make_profile(row, sample, labels, hon_seed_terms=()):
    # Clinical findings and shared HON references must remain distinguishable.
    observed = terms(row.get("hpo_present", ""))
    absent = terms(row.get("hpo_absent", ""))
    overlap = sorted(set(observed) & set(absent))
    if overlap:
        raise ValueError("HPO terms cannot be both present and absent: " + ", ".join(overlap))
    if not observed:
        raise ValueError(f"{sample}: no confirmed HPO terms; phenotype ranking cannot run")

    # Add the shared HON profile only for patients with documented optic atrophy.
    # It must not be imputed into unrelated controls (for example HG002).
    # Documented absences take precedence over reference seed terms.
    has_optic_atrophy = "HP:0000648" in observed
    reference = terms(";".join(hon_seed_terms)) if has_optic_atrophy else []
    skipped_absent = sorted(set(reference) & set(absent))
    reference_added = sorted(set(reference) - set(observed) - set(absent))
    scoring_terms = sorted(set(observed) | set(reference_added))
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
    for hpo_id in scoring_terms:
        feature = {"type": {"id": hpo_id}}
        if hpo_id in labels:
            feature["type"]["label"] = labels[hpo_id]
        features.append(feature)
    # The Phenopacket and the custom-ranking HPO TSV use the same combined
    # scoring terms. Provenance in clinical_context distinguishes actual
    # observed features from cohort-level reference assumptions.

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
        "hpo_present": observed,
        "hpo_absent": absent,
        "hpo_scoring_terms": scoring_terms,
        "hpo_hon_seed_profile_applied": has_optic_atrophy,
        "hpo_hon_reference_terms": reference,
        "hpo_hon_reference_added": reference_added,
        "hpo_hon_reference_excluded_as_absent": skipped_absent,
        "hpo_scoring_term_sources": {
            item: ("PATIENT_OBSERVED" if item in observed else "SHARED_HON_REFERENCE")
            for item in scoring_terms
        },
        "hpo_scoring_note": (
            ("Shared HON reference terms have been added for an explicit "
             "sensitivity analysis; only hpo_present represents patient-observed findings.")
            if reference_added else
            "Only explicitly patient-observed HPO terms are used for scoring; "
            "blank clinical fields remain unknown."
        ),
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
    parser.add_argument("--hpo-seeds", required=True)
    parser.add_argument(
        "--include-shared-hon-seeds", action="store_true",
        help="Opt-in legacy analysis: add cohort HON reference terms to scoring"
    )
    parser.add_argument("--phenopacket", required=True)
    parser.add_argument("--hpo-output", required=True)
    parser.add_argument("--clinical-output", required=True)
    parser.add_argument("--allele-context-output")
    args = parser.parse_args()

    labels = {}
    shared_hon_terms = []
    with open(args.hpo_seeds, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "hpo_id" not in reader.fieldnames:
            raise ValueError("HON seed table must contain an hpo_id column")
        for row in reader:
            hpo_id = str(row["hpo_id"] or "").strip()
            if not HPO_PATTERN.fullmatch(hpo_id):
                raise ValueError(f"Invalid HON seed term: {hpo_id}")
            labels[hpo_id] = row.get("hpo_label", "")
            if str(row.get("hon_seed_role") or "").strip().upper() in {"", "CORE_OCULAR_HON"}:
                shared_hon_terms.append(hpo_id)
    if not shared_hon_terms:
        raise ValueError("HON seed table is empty")
    packet, context = make_profile(
        load_profile(args.profiles, args.sample), args.sample, labels,
        hon_seed_terms=shared_hon_terms if args.include_shared_hon_seeds else ()
    )
    for filename in (args.phenopacket, args.hpo_output, args.clinical_output):
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
    Path(args.phenopacket).write_text(json.dumps(packet, indent=2) + "\n", encoding="utf-8")
    with open(args.hpo_output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["hpo_id", "hpo_label"])
        for identifier in context["hpo_scoring_terms"]:
            writer.writerow([identifier, labels.get(identifier, "")])
    Path(args.clinical_output).write_text(json.dumps(context, indent=2) + "\n", encoding="utf-8")
    if args.allele_context_output:
        row = load_profile(args.profiles, args.sample)
        path = Path(args.allele_context_output)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            fields = ["sample_id", "hpo_present", "hpo_absent", "sex", "mother_id", "father_id", "source"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerow({"sample_id": args.sample, "hpo_present": ";".join(context["hpo_present"]) or ".",
                             "hpo_absent": ";".join(context["hpo_absent"]) or ".",
                             "sex": row.get("sex") or ".", "mother_id": row.get("mother_id") or ".",
                             "father_id": row.get("father_id") or ".", "source": "private_patient_profiles_observed"})


if __name__ == "__main__":
    main()
