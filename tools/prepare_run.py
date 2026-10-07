#!/usr/bin/env python3
"""Assemble an external workflow run tree from code and archived data."""
import argparse
import hashlib
import shutil
from pathlib import Path

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.digest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--code-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--data-root', type=Path, required=True,
                        help='External data directory laid out like this repository')
    parser.add_argument('--output', type=Path, required=True, help='New directory outside both Git repositories')
    args = parser.parse_args()
    code, data, output = (p.expanduser().resolve() for p in (args.code_root, args.data_root, args.output))
    if not code.is_dir() or not data.is_dir():
        parser.error('Both source directories must exist')
    if output.exists():
        parser.error('Output already exists; choose a new run directory')
    for parent in (output, *output.parents):
        if (parent / '.git').exists():
            parser.error('Output must be outside Git working trees')
    if output.is_relative_to(code) or output.is_relative_to(data):
        parser.error('Output must be outside the source directories')
    mapping = {}
    for root in (code, data):
        for source in sorted(root.rglob('*')):
            relative = source.relative_to(root)
            if '.git' in relative.parts or source.name == '.DS_Store':
                continue
            if source.is_symlink():
                parser.error('Refusing symlink source: ' + str(source))
            if not source.is_file():
                continue
            if relative in mapping and sha256(source) != sha256(mapping[relative]):
                parser.error('Different code/data files share a path: ' + str(relative))
            mapping[relative] = source
    for relative, source in mapping.items():
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if sha256(source) != sha256(target):
            raise RuntimeError('Copy verification failed: ' + str(target))
    print(f'Prepared and verified {len(mapping)} files in {output}')
    print('Edit config.local.sh in this external run tree to select inputs, outputs and installed software.')

if __name__ == '__main__':
    main()
