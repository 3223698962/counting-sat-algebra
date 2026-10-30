"""Count a supplied DIMACS formula and emit a reproducible exact result."""

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .component_count import ComponentCountResult, count_cnf_components
from .dimacs import read_dimacs
from .lean_greedy import GreedyBudgets
from .lean_sparsify import StreamingBudgets, StreamingCountResult, count_cnf_streaming
from .lean_ring import LeanAlgebraBudgets
from .lean_taylor import TaylorBudgets
from .tensor_count import ContractionBudgets


def count_document(result: ComponentCountResult | StreamingCountResult) -> dict:
    data = asdict(result)
    count = result.count
    data['count'] = None if count is None else (str(count) if count.bit_length() <= 4096 else hex(count))
    data['count_encoding'] = None if count is None else ('decimal' if count.bit_length() <= 4096 else 'hex')
    data['satisfiable'] = None if count is None else count > 0
    if isinstance(result, StreamingCountResult):
        data['algorithm'] = (
            'Lean-guided streaming disjoint sparsification with sparse-code greedy hash leaves'
            if result.leaf_backend == 'hash' else
            'Lean-guided streaming disjoint sparsification with pointwise correction-pair leaves'
            if result.leaf_backend == 'stream' else
            'Lean-guided streaming disjoint sparsification with native ternary squarefree queries'
            if result.leaf_backend == 'native' else
            'Lean-guided streaming disjoint sparsification with finite algebraic leaves'
        )
        data['count_word_encoding'] = 'little-endian binary, zero represented by one zero bit'
    else:
        data['algorithm'] = ('independent components with finite Taylor tensor contraction' if result.backend == 'tensor'
                             else 'independent components with sparse shared-measurement local counts')
    return data


def count_main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog='counting_sat_algebra count', description='Exactly count a DIMACS CNF within explicit component and enumeration budgets.')
    parser.add_argument('path', type=Path, help='DIMACS .cnf input')
    parser.add_argument('--backend', choices=('common', 'tensor', 'lean'), default='common', help='algebraic backend or streaming disjoint sparsification')
    parser.add_argument('--leaf-backend', choices=('common', 'tensor', 'native', 'hash', 'stream'), default='common', help='algebraic leaf evaluator for --backend lean; native uses ternary squarefree sums, hash runs sparse-code greedy selection with a separate accepted query, stream uses pointwise correction pairs without support tables')
    parser.add_argument('--max-stream-pairs', type=int, default=2_000_000, help='global correction-pair budget across all stream leaves')
    parser.add_argument('--max-stream-evaluations', type=int, default=20_000_000, help='global point-evaluation budget across all stream leaves')
    parser.add_argument('--hash-depth', type=int, default=1, help='candidate depth of the hash leaf; the Lean schedule depth is distanceDenom, above the Python budgets')
    parser.add_argument('--max-hash-candidates', type=int, default=300_000, help='maximum canonical ternary candidates per hash round and leaf')
    parser.add_argument('--high-loop', action='store_true', help='route the hash leaf marked and accepted queries through the streamed high-variable rank loop')
    parser.add_argument('--max-ring-precision', type=int, default=64, help='maximum t for native residues modulo 3^t')
    parser.add_argument('--max-algebra-terms', type=int, default=16_384, help='maximum terms per native squarefree polynomial')
    parser.add_argument('--max-algebra-products', type=int, default=2_000_000, help='global native monomial product budget across all leaves')
    parser.add_argument('--max-algebra-operations', type=int, default=8_000_000, help='global native scalar operation budget across all leaves')
    parser.add_argument('--max-algebra-bits', type=int, default=16_000_000, help='native logical encoding storage budget')
    parser.add_argument('--sparsify-width', type=int, help='fixed width bound k for --backend lean; inferred if omitted')
    parser.add_argument('--sparsify-M', type=int, help='explicit branching parameter; omitted selects the first passing dyadic M')
    parser.add_argument('--sparsify-sigma', type=str, default='1', help='positive rational exponent budget for the dyadic search')
    parser.add_argument('--max-sparsify-nodes', type=int, default=100_000)
    parser.add_argument('--max-sparsify-candidates', type=int, default=2_000_000)
    parser.add_argument('--max-frontier-bits', type=int, default=16_000_000)
    parser.add_argument('--max-formula-bits', type=int, default=16_000_000)
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
        tensor_limits = ContractionBudgets(max_bucket_variables=args.max_tensor_width,
                                          max_component_degree=args.max_tensor_degree)
        if args.backend == 'lean':
            from fractions import Fraction

            sigma = Fraction(args.sparsify_sigma)
            if sigma <= 0:
                raise ValueError('sparsify sigma must be positive')
            result = count_cnf_streaming(
                cnf, k=args.sparsify_width, M=args.sparsify_M,
                sigma_numerator=sigma.numerator, sigma_denominator=sigma.denominator,
                budgets=StreamingBudgets(
                    max_nodes=args.max_sparsify_nodes,
                    max_search_candidates=args.max_sparsify_candidates,
                    max_frontier_bits=args.max_frontier_bits, max_formula_bits=args.max_formula_bits,
                ), leaf_backend=args.leaf_backend,
                max_component_variables=args.max_component_variables, vector_budget=args.vector_budget,
                total_vector_budget=args.total_vector_budget, memo_entries=0,
                tensor_budgets=tensor_limits,
                native_budgets=LeanAlgebraBudgets(
                    max_precision=args.max_ring_precision, max_polynomial_terms=args.max_algebra_terms,
                    max_monomial_products=args.max_algebra_products,
                    max_scalar_operations=args.max_algebra_operations,
                    max_algebra_bits=args.max_algebra_bits,
                ),
                greedy_budgets=GreedyBudgets(
                    depth=args.hash_depth,
                    max_candidates=args.max_hash_candidates,
                ),
                high_loop=args.high_loop,
                taylor_budgets=TaylorBudgets(
                    max_pairs=args.max_stream_pairs,
                    max_assignments=args.max_stream_evaluations,
                ),
            )
        else:
            result = count_cnf_components(cnf, max_component_variables=args.max_component_variables,
                                         vector_budget=args.vector_budget, total_vector_budget=args.total_vector_budget,
                                         memo_entries=0 if args.no_memo else 4096, backend=args.backend,
                                         tensor_budgets=tensor_limits)
        document = count_document(result)
        document.update(input=str(args.path.resolve()), n=cnf.n, m=cnf.m)
        payload = json.dumps(document, ensure_ascii=False, indent=2) + '\n'
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(payload, encoding='utf-8')
    except (OSError, UnicodeError, ValueError, ZeroDivisionError) as exc:
        parser.error(str(exc))
    if args.json:
        print(payload, end='')
    else:
        print(result.status.upper())
        if result.counted:
            print(f"count = {document['count']}")
            if isinstance(result, StreamingCountResult):
                print(f"sparsification leaves = {result.leaves_counted}; nodes = {result.nodes}; maximum pending nodes = {result.max_pending_nodes}")
            else:
                print(f"components = {result.components}; declared free variables = {result.free_variables}")
            print(f"local evaluations = {result.local_evaluations}; memo hits = {result.memo_hits}; local Boolean vectors = {result.boolean_vectors}")
            if isinstance(result, StreamingCountResult) and result.leaf_backend == 'native':
                print(f"native monomial products = {result.monomial_products}; scalar operations = {result.scalar_operations}; maximum polynomial terms = {result.peak_polynomial_terms}")
            if isinstance(result, StreamingCountResult) and result.leaf_backend == 'hash':
                print(f"hash rounds = {result.hash_rounds}; marked queries = {result.hash_marked_queries}; accepted queries = {result.hash_accepted_queries}; final violation mass = {result.hash_final_mass}")
            if isinstance(result, StreamingCountResult) and result.leaf_backend == 'stream':
                print(f"stream pairs = {result.stream_pairs}; point evaluations = {result.stream_point_evaluations}")
        else:
            print(result.reason)
        if args.output is not None:
            print(f"JSON result saved to {args.output.resolve()}")
    return 0 if result.counted else (2 if result.status == 'refused' else 1)
