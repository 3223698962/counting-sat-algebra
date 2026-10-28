#!/usr/bin/env python3
"""Fail when a Python module under src/counting_sat_algebra grows past the limit.

The hard limit is 65536 bytes (64 KiB). After the reference_verify split, the
largest module is verify.py, which is under that cap. A softer limit of 49152
bytes (48 KiB) applies to every file except the three modules grandfathered
below. Those three may sit between the two limits until a later, tested split.
New files over the soft limit fail. Any file over the hard limit fails.

docs/ARCHITECTURE.md must name both limits and the grandfathered modules.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'src' / 'counting_sat_algebra'
ARCHITECTURE = ROOT / 'docs' / 'ARCHITECTURE.md'
HARD_LIMIT = 65536
SOFT_LIMIT = 49152
GRANDFATHERED = frozenset(
    {
        'verify.py',
        'parameter_scan.py',
        'equation_selection.py',
    }
)


def main() -> int:
    errors: list[str] = []
    if not ARCHITECTURE.is_file():
        errors.append('docs/ARCHITECTURE.md is missing')
    else:
        documented = ARCHITECTURE.read_text(encoding='utf-8')
        for token in (str(HARD_LIMIT), str(SOFT_LIMIT), *sorted(GRANDFATHERED)):
            if token not in documented:
                errors.append(f'docs/ARCHITECTURE.md does not mention {token}')
    for path in sorted(PACKAGE.rglob('*.py')):
        size = path.stat().st_size
        relative = path.relative_to(ROOT)
        if size > HARD_LIMIT:
            errors.append(f'{relative} is {size} bytes, over the {HARD_LIMIT}-byte hard limit')
        elif size > SOFT_LIMIT and path.name not in GRANDFATHERED:
            errors.append(f'{relative} is {size} bytes, over the {SOFT_LIMIT}-byte soft limit')
    if errors:
        print('module size gate failed:', file=sys.stderr)
        for error in errors:
            print(f'- {error}', file=sys.stderr)
        return 1
    print(
        f'module size gate ok: hard limit {HARD_LIMIT} bytes, '
        f'soft limit {SOFT_LIMIT} bytes, grandfathered {", ".join(sorted(GRANDFATHERED))}'
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
