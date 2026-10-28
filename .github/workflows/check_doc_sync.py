#!/usr/bin/env python3
"""Fail when the documented package version drifts from pyproject.toml.

docs/STATUS.md and references/lean/README.md must both mention the current
version. The version string in src/counting_sat_algebra/__init__.py must
match pyproject.toml. A longer version such as 0.22.10 does not satisfy a
check for 0.22.1.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _pyproject_version() -> str:
    text = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')
    match = re.search(r'(?m)^version = "([^"]+)"', text)
    if match is None:
        raise SystemExit('pyproject.toml is missing a version')
    return match.group(1)


def _init_version() -> str:
    text = (ROOT / 'src' / 'counting_sat_algebra' / '__init__.py').read_text(encoding='utf-8')
    match = re.search(r'(?m)^__version__ = "([^"]+)"', text)
    if match is None:
        raise SystemExit('__init__.py is missing __version__')
    return match.group(1)


def _mentions(text: str, version: str) -> bool:
    return re.search(rf'(?<![\w.]){re.escape(version)}(?![\w.])', text) is not None


def main() -> int:
    package = _pyproject_version()
    init = _init_version()
    errors: list[str] = []
    if package != init:
        errors.append(f'pyproject.toml version {package} != __version__ {init}')
    required = {
        'docs/STATUS.md': ROOT / 'docs' / 'STATUS.md',
        'references/lean/README.md': ROOT / 'references' / 'lean' / 'README.md',
    }
    for label, path in required.items():
        if not path.is_file():
            errors.append(f'missing {label}')
            continue
        if not _mentions(path.read_text(encoding='utf-8'), package):
            errors.append(f'{label} does not mention package version {package}')
    if errors:
        print('doc sync gate failed:', file=sys.stderr)
        for error in errors:
            print(f'- {error}', file=sys.stderr)
        return 1
    print(f'doc sync gate ok: STATUS and Lean README mention {package}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
