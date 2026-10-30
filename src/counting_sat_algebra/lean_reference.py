"""Read-only source audit of the supplied Lean reference, without kernel claims.

The audit follows local OAI imports from Main.lean, ignores nested comments
and strings, and records source hashes and code-level proof placeholders.
Mathlib and other Lake dependencies still need the actual Lean toolchain.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

ENTRY = 'OAI.Algorithms.SubexponentialSAT.Main'
DECLARATIONS = (
    ('OAI/Algorithms/SubexponentialSAT/SubexpOrderedRawMainTheorem.lean', 'mainTheorem_explicit'),
    ('OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean', 'compiled_exact_count_and_sat'),
    ('OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean', 'exact_count_and_sat'),
    ('OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean', 'exact_count_eth'),
)


def lean_code(source: str) -> str:
    """Erase strings, line comments and nested block comments, retaining lines."""
    output = list(source)
    depth, position, quoted, line, escaped = 0, 0, False, False, False
    while position < len(source):
        char = source[position]
        pair = source[position : position + 2]
        if line:
            if char == '\n':
                line = False
            else:
                output[position] = ' '
        elif depth:
            if pair == '/-':
                depth += 1
                output[position : position + 2] = [' ', ' ']
                position += 2
                continue
            if pair == '-/':
                depth -= 1
                output[position : position + 2] = [' ', ' ']
                position += 2
                continue
            if char not in '\r\n':
                output[position] = ' '
        elif quoted:
            if char not in '\r\n':
                output[position] = ' '
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == '"':
                quoted = False
        elif pair == '--':
            line = True
            output[position : position + 2] = [' ', ' ']
            position += 2
            continue
        elif pair == '/-':
            depth = 1
            output[position : position + 2] = [' ', ' ']
            position += 2
            continue
        elif char == '"':
            quoted = True
            output[position] = ' '
        position += 1
    return ''.join(output)


def audit_lean_reference(root: str | Path) -> dict:
    """Check source availability and syntax markers; never label this a Lean build."""
    root = Path(root).resolve()
    entries = sorted({ENTRY, *(file[:-5].replace('/', '.') for file, _ in DECLARATIONS)})
    pending, visited, records, missing, external = list(entries), set(), {}, set(), set()
    placeholders, axioms = [], []
    declaration_code = {}
    source_bytes = 0
    while pending:
        module = pending.pop()
        if module in visited:
            continue
        visited.add(module)
        path = root.joinpath(*module.split('.')).with_suffix('.lean')
        if not path.is_file():
            missing.add(module)
            continue
        data = path.read_bytes()
        source_bytes += len(data)
        code = lean_code(data.decode('utf-8'))
        relative = path.relative_to(root).as_posix()
        records[relative] = hashlib.sha256(data).hexdigest()
        if any(relative == name for name, _ in DECLARATIONS):
            declaration_code[relative] = code
        for match in re.finditer(r'\b(?:sorry|admit)\b', code):
            placeholders.append(
                dict(
                    file=relative, line=code.count('\n', 0, match.start()) + 1, token=match.group()
                )
            )
        for match in re.finditer(
            r'(?m)^\s*(?:(?:private|protected|noncomputable|unsafe)\s+)*axiom\b', code
        ):
            axioms.append(dict(file=relative, line=code.count('\n', 0, match.end()) + 1))
        for match in re.finditer(r'(?m)^\s*import\s+([^\r\n]+)', code):
            for imported in match.group(1).split():
                if imported.startswith('OAI.'):
                    pending.append(imported)
                else:
                    external.add(imported.split('.')[0])
    declarations = [
        dict(
            file=file,
            declaration=name,
            body_present=bool(
                re.search(
                    rf'\btheorem\s+{re.escape(name)}\b[\s\S]*?:=\s*\S',
                    declaration_code.get(file, ''),
                )
            ),
        )
        for file, name in DECLARATIONS
    ]
    digest = hashlib.sha256()
    for file, file_hash in sorted(records.items()):
        digest.update(f'{file}\0{file_hash}\n'.encode('utf-8'))
    toolchain = root / 'lean-toolchain'
    return dict(
        schema_version=1,
        scope='Static import-closure and source-marker audit; no Lean execution or kernel verification.',
        entry_module=ENTRY,
        audited_entry_modules=entries,
        toolchain=toolchain.read_text(encoding='utf-8').strip() if toolchain.is_file() else None,
        local_modules=len(records),
        source_bytes=source_bytes,
        source_sha256=digest.hexdigest(),
        missing_local_imports=sorted(missing),
        external_dependency_roots=sorted(external),
        proof_placeholders=sorted(placeholders, key=lambda row: (row['file'], row['line'])),
        axiom_declarations=sorted(axioms, key=lambda row: (row['file'], row['line'])),
        declarations=declarations,
        static_checks_passed=bool(records)
        and not missing
        and not placeholders
        and not axioms
        and all(row['body_present'] for row in declarations),
        kernel_verified=False,
        python_algorithm_complete=False,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--root',
        type=Path,
        default=Path(__file__).resolve().parents[2] / 'references' / 'lean' / 'math' / 'lean',
    )
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    report = audit_lean_reference(args.root)
    payload = json.dumps(report, ensure_ascii=False, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding='utf-8')
    print(payload, end='')
    return 0 if report['static_checks_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
