#!/usr/bin/env python3
"""Install a validated gene rebuild, backing up the old files and promoter cache."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
from pathlib import Path

from rebuild_annotsv_gene_resource import BAD_NAMES, read_bed, sha256


GENES = Path("Annotations_Human/Genes/GRCh38")
REGULATORY = Path("Annotations_Human/FtIncludedInSV/RegulatoryElements/GRCh38")
BED = GENES / "genes.ENSEMBL.sorted.bed"
VERSIONS = GENES / "transcript_version.ENSEMBL.tsv"


def install(build, root):
    build, root = Path(build).resolve(), Path(root).resolve()
    metadata = json.loads((build / "reference_build.json").read_text())
    if metadata.get("status") != "READY" or metadata.get("genome_build") != "GRCh38":
        raise ValueError("The rebuild is not marked READY for GRCh38")
    expected = {BED: metadata["bed_sha256"], VERSIONS: metadata["transcript_versions_sha256"]}
    for relative, digest in expected.items():
        if sha256(build / relative) != digest:
            raise ValueError(f"The rebuilt file changed after validation: {relative}")
    rows = read_bed(build / BED)
    if any(row[4].upper() in BAD_NAMES for row in rows.values()):
        raise ValueError("The rebuilt BED still contains incorrect gene labels")
    versions = {}
    for line in (build / VERSIONS).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 2 or fields[0] in versions:
            raise ValueError("Invalid or duplicated transcript version entry")
        versions[fields[0]] = fields[1]
    if set(versions) != set(rows):
        raise ValueError("Transcript version entries do not match the rebuilt BED")
    if not root.is_dir() or not (root / BED).is_file():
        raise ValueError("The installed annotation folder does not contain the original BED")
    if sha256(root / BED) != metadata["existing_bed_sha256"]:
        raise ValueError("The installed BED differs from the file used in this rebuild")
    lock = root / ".gene-reference-update.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise ValueError(f"Another installation may be running: {lock}") from error
    os.close(descriptor)
    try:
        if sha256(root / BED) != metadata["existing_bed_sha256"]:
            raise ValueError("The installed BED changed before the installation lock was acquired")
        return install_locked(build, root, metadata)
    finally:
        lock.unlink(missing_ok=True)


def install_locked(build, root, metadata):
    backup = Path(tempfile.mkdtemp(prefix="gene_reference_backup_", dir=root))
    originals = {}
    targets = [root / BED, root / VERSIONS]
    caches = sorted((root / REGULATORY).glob("promoter_*bp_ENSEMBL_GRCh38.sorted.bed"))
    caches += sorted((root / GENES).glob("genes.ENSEMBL.sorted.tmp*.bed"))
    caches += [root / ".installed"]
    for target in targets + caches:
        if target.is_symlink():
            raise ValueError(f"Resolve this symlink before installing: {target}")
        if target.exists():
            relative = target.relative_to(root)
            saved = backup / relative
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(target, saved)
            originals[target] = saved
    changes = []
    try:
        for relative in (BED, VERSIONS):
            target = root / relative
            descriptor, temporary = tempfile.mkstemp(prefix=".new_gene_resource_", dir=target.parent)
            os.close(descriptor)
            temporary = Path(temporary)
            try:
                shutil.copy2(build / relative, temporary)
                os.replace(temporary, target)
                changes.append(target)
            finally:
                temporary.unlink(missing_ok=True)
        for cached in caches:
            if cached in originals:
                cached.unlink()
                changes.append(cached)
        receipt = {"status": "INSTALLED", "source_release": metadata["source_release"],
                   "bed_sha256": metadata["bed_sha256"],
                   "backup_directory": str(backup),
                   "promoter_files_to_regenerate": [str(p.relative_to(root)) for p in caches
                                                    if p in originals and "promoter_" in p.name]}
        (backup / "installation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    except Exception:
        for target in reversed(changes):
            if target in originals:
                shutil.copy2(originals[target], target)
            else:
                target.unlink(missing_ok=True)
        raise
    print(f"Installed the rebuilt gene BED and transcript versions in {root / GENES}")
    print(f"Original files saved in {backup}")
    print("AnnotSV will regenerate the affected promoter files on its next run.")
    print("Rerun AnnotSV and the downstream gene analysis on the existing master VCF.")
    return backup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--annotations-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        install(args.build_dir, args.annotations_dir)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"[ERROR] {error}\n")


if __name__ == "__main__":
    main()
