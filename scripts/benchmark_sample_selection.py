"""Restrict optional benchmark runs to independently confirmed samples."""
import csv


def confirmed_benchmark_samples(path, allowed_samples):
    if not path:
        return []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample_id", "truth_status"}.issubset(reader.fieldnames or []):
            raise ValueError("Benchmark truth requires sample_id and truth_status")
        confirmed = {
            str(row["sample_id"]).strip()
            for row in reader
            if str(row["truth_status"]).strip().upper() == "CONFIRMED"
        }
    return [str(sample) for sample in allowed_samples if str(sample) in confirmed]
