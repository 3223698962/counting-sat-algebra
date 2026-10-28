#!/usr/bin/env python3
"""Fail unless tracked text is UTF-8, without NUL bytes or a UTF-16/UTF-32 BOM."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXTENSIONS = {
    '.py',
    '.md',
    '.toml',
    '.yml',
    '.yaml',
    '.json',
    '.txt',
    '.cnf',
    '.lean',
}
SKIP_DIRS = {'.git', '__pycache__', '.pytest_cache', '.ruff_cache', '.mypy_cache'}
BOMS = {
    b'\xff\xfe': 'UTF-16 LE',
    b'\xfe\xff': 'UTF-16 BE',
    b'\xff\xfe\x00\x00': 'UTF-32 LE',
    b'\x00\x00\xfe\xff': 'UTF-32 BE',
}


def _iter_files():
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix.lower() not in EXTENSIONS:
            continue
        yield path


def main() -> int:
    errors: list[str] = []
    checked = 0
    for path in _iter_files():
        checked += 1
        data = path.read_bytes()
        relative = path.relative_to(ROOT)
        if b'\x00' in data:
            errors.append(f'{relative} contains a NUL byte')
            continue
        for bom, label in BOMS.items():
            if data.startswith(bom):
                errors.append(f'{relative} starts with a {label} BOM')
                break
        else:
            try:
                data.decode('utf-8')
            except UnicodeDecodeError as exc:
                errors.append(f'{relative} is not UTF-8 ({exc.reason} at byte {exc.start})')
    if checked == 0:
        errors.append('encoding gate found no text files')
    if errors:
        print('encoding gate failed:', file=sys.stderr)
        for error in errors:
            print(f'- {error}', file=sys.stderr)
        return 1
    print(f'encoding gate ok: {checked} UTF-8 text files, no NUL bytes')
    return 0


if __name__ == '__main__':
    sys.exit(main())
