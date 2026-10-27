"""Strict DIMACS CNF input preserving the number of declared variables."""

from pathlib import Path

from .cnf import CNF


def parse_dimacs(text: str) -> CNF:
    """Read comments, multiline clauses, empty clauses and a required header."""
    header = None
    clauses = []
    pending = []
    for line_number, raw in enumerate(text.lstrip('\ufeff').splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith('c'):
            continue
        if line == '%':
            break
        tokens = line.split()
        if tokens[0] == 'p':
            if header is not None or len(tokens) != 4 or tokens[1] != 'cnf':
                raise ValueError(f"line {line_number}: expected one 'p cnf n m' header")
            try:
                n, m = int(tokens[2]), int(tokens[3])
            except ValueError as exc:
                raise ValueError(f"line {line_number}: header counts must be integers") from exc
            if n < 0 or m < 0:
                raise ValueError(f"line {line_number}: header counts must be nonnegative")
            header = n, m
            continue
        if header is None:
            raise ValueError(f"line {line_number}: literals appear before the DIMACS header")
        n, m = header
        for token in tokens:
            try:
                value = int(token)
            except ValueError as exc:
                raise ValueError(f"line {line_number}: invalid literal {token!r}") from exc
            if not value:
                clauses.append(pending)
                pending = []
                if len(clauses) > m:
                    raise ValueError(f"line {line_number}: more clauses than declared")
            else:
                if abs(value) > n:
                    raise ValueError(f"line {line_number}: literal outside the declared variables")
                pending.append((abs(value) - 1, value < 0))
    if header is None:
        raise ValueError("missing DIMACS header")
    if pending:
        raise ValueError("the last clause is missing its terminating zero")
    if len(clauses) != header[1]:
        raise ValueError(f"expected {header[1]} clauses, found {len(clauses)}")
    return CNF(header[0], clauses)


def read_dimacs(path: str | Path) -> CNF:
    return parse_dimacs(Path(path).read_text(encoding='utf-8-sig'))


def format_dimacs(cnf: CNF) -> str:
    lines = [f"p cnf {cnf.n} {cnf.m}"]
    for clause in cnf.clauses:
        literals = [str(-(var + 1) if neg else var + 1) for var, neg in clause]
        lines.append(' '.join(literals + ['0']))
    return '\n'.join(lines) + '\n'
