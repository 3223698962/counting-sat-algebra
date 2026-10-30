#!/usr/bin/env python3
"""Fail when a standing honesty flag is true.

The flags are theorem_1_1_complete, lean_proof_complete, and
recursive_contraction. They must stay false in the committed reference
report, the step verification artifacts, and in Python assignments under
src/counting_sat_algebra. Absence of recursive_contraction is a failure on
the reference report and on the Step 20 through Step 30 artifacts,
because those checks are the ones that record the flag.
"""

from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path

FLAGS = ('theorem_1_1_complete', 'lean_proof_complete', 'recursive_contraction')
ROOT = Path(__file__).resolve().parents[2]
STEP_ARTIFACTS = tuple(
    ROOT / 'artifacts' / f'step{step}-verification.json' for step in range(15, 31)
)
REFERENCE_REPORT = ROOT / 'benches' / 'reference_alignment.json'
RECURSIVE_REQUIRED = {
    REFERENCE_REPORT,
    ROOT / 'artifacts' / 'step20-verification.json',
    ROOT / 'artifacts' / 'step21-verification.json',
    ROOT / 'artifacts' / 'step22-verification.json',
    ROOT / 'artifacts' / 'step23-verification.json',
    ROOT / 'artifacts' / 'step24-verification.json',
    ROOT / 'artifacts' / 'step25-verification.json',
    ROOT / 'artifacts' / 'step26-verification.json',
    ROOT / 'artifacts' / 'step27-verification.json',
    ROOT / 'artifacts' / 'step28-verification.json',
    ROOT / 'artifacts' / 'step29-verification.json',
    ROOT / 'artifacts' / 'step30-verification.json',
}


def _version() -> str:
    text = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')
    match = re.search(r'(?m)^version = "([^"]+)"', text)
    if match is None:
        raise SystemExit('pyproject.toml is missing a version')
    return match.group(1)


def _walk(value, path: str, found: list[tuple[str, object]]) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            child = f'{path}.{key}'
            if key in FLAGS:
                found.append((child, item))
            _walk(item, child, found)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _walk(item, f'{path}[{index}]', found)


def _check_json(path: Path, version: str, errors: list[str]) -> None:
    if not path.is_file():
        errors.append(f'missing honesty artifact: {path.relative_to(ROOT)}')
        return
    report = json.loads(path.read_text(encoding='utf-8'))
    relative = path.relative_to(ROOT)
    if report.get('version') != version:
        errors.append(f'{relative} version is {report.get("version")!r}, expected {version}')
    for flag in ('theorem_1_1_complete', 'lean_proof_complete'):
        if report.get(flag) is not False:
            errors.append(f'{relative} top-level {flag} must be false')
    if report.get('recursive_contraction') not in (None, False):
        errors.append(f'{relative} top-level recursive_contraction must be false')
    found: list[tuple[str, object]] = []
    _walk(report, '$', found)
    for pointer, value in found:
        if value is not False:
            errors.append(f'{relative} {pointer} is {value!r}, expected false')
    if path in RECURSIVE_REQUIRED and not any(
        pointer.endswith('recursive_contraction') for pointer, _ in found
    ):
        errors.append(f'{relative} does not record recursive_contraction=false')


class _FlagVisitor(ast.NodeVisitor):
    def __init__(self, path: Path):
        self.path = path
        self.errors: list[str] = []

    def _true(self, node: ast.AST | None) -> bool:
        return isinstance(node, ast.Constant) and node.value is True

    def _name_of(self, node: ast.AST) -> str | None:
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Attribute):
            return node.attr
        if isinstance(node, ast.Subscript):
            index = node.slice
            if isinstance(index, ast.Constant) and isinstance(index.value, str):
                return index.value
        return None

    def _reject(self, node: ast.AST, name: str) -> None:
        self.errors.append(f'{self.path.relative_to(ROOT)}:{node.lineno} sets {name}=True')

    def visit_Assign(self, node: ast.Assign) -> None:
        if self._true(node.value):
            for target in node.targets:
                name = self._name_of(target)
                if name in FLAGS:
                    self._reject(node, name)
        self.generic_visit(node)

    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        name = self._name_of(node.target)
        if name in FLAGS and self._true(node.value):
            self._reject(node, name)
        self.generic_visit(node)

    def visit_keyword(self, node: ast.keyword) -> None:
        if node.arg in FLAGS and self._true(node.value):
            self._reject(node, node.arg)
        self.generic_visit(node)

    def visit_Dict(self, node: ast.Dict) -> None:
        for key, value in zip(node.keys, node.values):
            if isinstance(key, ast.Constant) and key.value in FLAGS and self._true(value):
                self._reject(key, key.value)
        self.generic_visit(node)


def _check_sources(errors: list[str]) -> None:
    root = ROOT / 'src' / 'counting_sat_algebra'
    for path in sorted(root.rglob('*.py')):
        tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
        visitor = _FlagVisitor(path)
        visitor.visit(tree)
        errors.extend(visitor.errors)


def main() -> int:
    version = _version()
    errors: list[str] = []
    for path in (REFERENCE_REPORT, *STEP_ARTIFACTS):
        _check_json(path, version, errors)
    _check_sources(errors)
    if errors:
        print('honesty gate failed:', file=sys.stderr)
        for error in errors:
            print(f'- {error}', file=sys.stderr)
        return 1
    print(
        f'honesty gate ok: version {version}; '
        'theorem_1_1_complete, lean_proof_complete, and recursive_contraction stay false'
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
