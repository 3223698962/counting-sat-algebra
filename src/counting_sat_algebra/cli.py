"""Count a supplied DIMACS formula and emit a reproducible exact result."""

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .component_count import ComponentCountResult, count_cnf_components
from .dimacs import read_dimacs
from .tensor_count import ContractionBudgets


def count_document(result: ComponentCountResult) -> dict:
    data = asdict(result)
    count = result.count
    data['count'] = None if count is None else (str(count) if count.bit_length() <= 4096 else hex(count))
    data['count_encoding'] = None if count is None else ('decimal' if count.bit_length() <= 4096 else 'hex')
    data['satisfiable'] = None if count is None else count > 0
    data['algorithm'] = ('independent components with finite Taylor tensor contraction' if result.backend == 'tensor'
                         else 'independent components with sparse shared-measurement local counts')
    return data


def count_main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog='counting_sat_algebra count', description='Exactly count a DIMACS CNF within explicit component and enumeration budgets.')
    parser.add_argument('path', type=Path, help='DIMACS .cnf input')
    parser.add_argument('--backend', choices=('common', 'tensor'), default='common', help='local algebraic counting backend')
    parser.add_argument('--max-component-variables', type=int, default=20)
    parser.add_argument('--max-tensor-width', type=int, default=12, help='maximum variables in an elimination bucket')
    parser.add_argument('--max-tensor-degree', type=int, default=64, help='maximum Taylor degree of a support component')
    parser.add_argument('--vector-budget', type=int, default=500_000)
    parser.add_argument('--total-vector-budget', type=int, default=5_000_000)
    parser.add_argument('--no-memo', action='store_true', help='disable reuse of identical local formulas')
    parser.add_argument('--json', action='store_true', help='print a JSON result')
    parser.add_argument('--output', type=Path, help='save the JSON result to this path')
    args = parser.parse_args(argv)
    try:
        cnf = read_dimacs(args.path)
        result = count_cnf_components(cnf, max_component_variables=args.max_component_variables,
                                     vector_budget=args.vector_budget, total_vector_budget=args.total_vector_budget,
                                     memo_entries=0 if args.no_memo else 4096, backend=args.backend,
                                     tensor_budgets=ContractionBudgets(max_bucket_variables=args.max_tensor_width,
                                                                       max_component_degree=args.max_tensor_degree))
        document = count_document(result)
        document.update(input=str(args.path.resolve()), n=cnf.n, m=cnf.m)
        payload = json.dumps(document, ensure_ascii=False, indent=2) + '\n'
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(payload, encoding='utf-8')
    except (OSError, UnicodeError, ValueError) as exc:
        parser.error(str(exc))
    if args.json:
        print(payload, end='')
    else:
        print(result.status.upper())
        if result.counted:
            print(f"count = {document['count']}")
            print(f"components = {result.components}; declared free variables = {result.free_variables}")
            print(f"local evaluations = {result.local_evaluations}; memo hits = {result.memo_hits}; local Boolean vectors = {result.boolean_vectors}")
        else:
            print(result.reason)
        if args.output is not None:
            print(f"JSON result saved to {args.output.resolve()}")
    return 0 if result.counted else (2 if result.status == 'refused' else 1)
