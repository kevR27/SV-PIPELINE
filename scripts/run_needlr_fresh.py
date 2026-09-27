#!/usr/bin/env python3
"""Run needLR in a fresh directory; never reuse unverified native results."""
import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    for name in ['vcf', 'backend', 'outdir', 'tsv', 'output_vcf', 'manifest']:
        p.add_argument('--' + name.replace('_', '-'), required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--backend-version', default='UNSPECIFIED')
    a = p.parse_args()
    root = Path(a.outdir)
    root.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix='run-', dir=root))
    command = ['needLR', 'annotate', '-O', str(run / 'native'), '-T', str(a.threads), '-B', a.backend, '--all', a.vcf]
    subprocess.run(command, check=True)
    candidates = list(run.rglob('*_RESULTS.tsv'))
    pairs = [(path, path.with_name(path.name.removesuffix('_RESULTS.tsv') + '.vcf.gz')) for path in candidates]
    pairs = [(t, v) for t, v in pairs if t.stat().st_size and v.exists() and v.stat().st_size]
    if len(pairs) != 1:
        raise RuntimeError(f'Expected one complete needLR result pair in {run}, found {len(pairs)}')
    digest = hashlib.sha256()
    with open(a.vcf, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b''):
            digest.update(chunk)
    manifest = {'command': command, 'query_sha256': digest.hexdigest(), 'backend': a.backend,
                'backend_version': a.backend_version, 'native_run_directory': str(run)}
    for source, dest in zip(pairs[0], [a.tsv, a.output_vcf]):
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, dest)
    Path(a.manifest).write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
