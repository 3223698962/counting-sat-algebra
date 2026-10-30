"""Step 23 finite checks for the Lean-guided Python counting frontier."""

from __future__ import annotations

from itertools import product

from .protocol import _row, stamp_completion
from ..cnf import CNF, is_satisfied
from ..lean_encoding import (
    count_word,
    decode_cnf,
    decode_count_word,
    encode_cnf,
    formula_bit_length,
)
from ..lean_parameters import choose_sparse_parameter, parameter_test
from ..lean_sparsify import (
    STEP23_SCOPE,
    StreamingBudgets,
    count_cnf_streaming,
    normalize_lean_cnf,
    select_priority,
    split_priority,
)


def _models(cnf: CNF) -> set[tuple[int, ...]]:
    return {bits for bits in product((0, 1), repeat=cnf.n) if is_satisfied(bits, cnf)}


def run_lean_streaming_suite() -> list[dict]:
    print('=== (Y) Lean-guided streaming disjoint counting ===', flush=True)
    rows = []
    star = CNF(4, [[(0, False), (variable, False)] for variable in range(1, 4)])
    bits = encode_cnf(star)
    rows.append(
        _row(
            'lean-canonical-encoding',
            star.n,
            star.m,
            len(bits) == formula_bit_length(star) == 138
            and _models(decode_cnf(bits)) == _models(star),
            'CHECKED',
            encoded_bits=len(bits),
            declared_n=decode_cnf(bits).n,
            count=None,
        )
    )

    repeated = CNF(
        5, [[(0, False)], [(0, False), (1, False)], [(1, False), (1, True)], [(0, False)]]
    )
    normalized = normalize_lean_cnf(repeated)
    rows.append(
        _row(
            'lean-normalization-keeps-declared-n',
            repeated.n,
            repeated.m,
            normalized.n == repeated.n
            and normalized.clauses == [[(0, False)]]
            and _models(normalized) == _models(repeated),
            'CHECKED',
            normalized_clauses=normalized.m,
            declared_n=normalized.n,
            unit_propagation=False,
            count=None,
        )
    )

    parameter = choose_sparse_parameter(3)
    rows.append(
        _row(
            'lean-first-passing-dyadic-M',
            0,
            0,
            parameter.M == 512
            and parameter.exponent == 9
            and parameter_test(3, 1, 1, 9)
            and all(not parameter_test(3, 1, 1, e) for e in range(9)),
            'CHECKED',
            M=parameter.M,
            exponent=parameter.exponent,
            occurrence_bound=parameter.occurrence_bound,
            count=None,
        )
    )

    clauses = tuple(frozenset(clause) for clause in star.clauses)
    choice = select_priority(clauses, 2)
    if choice is None:
        raise RuntimeError('the reference star has no selected core')
    left, right = split_priority(clauses, choice)
    a, b = (CNF(star.n, [list(clause) for clause in child]) for child in (left, right))
    models_a, models_b = _models(a), _models(b)
    rows.append(
        _row(
            'lean-all-selected-disjoint-children',
            star.n,
            star.m,
            choice.degree == 3
            and choice.threshold == 2
            and not (models_a & models_b)
            and models_a | models_b == _models(star),
            'CHECKED',
            selected_clauses=choice.degree,
            threshold=choice.threshold,
            left_models=len(models_a),
            right_models=len(models_b),
            count=None,
        )
    )

    for backend in ('common', 'tensor'):
        counted = count_cnf_streaming(star, M=2, leaf_backend=backend)
        rows.append(
            _row(
                f'lean-streamed-{backend}-leaf-count',
                star.n,
                star.m,
                counted.counted
                and counted.count == len(_models(star)) == 9
                and counted.leaves_counted == 2
                and counted.nodes == 3
                and counted.max_pending_nodes <= counted.depth + 1
                and counted.same_declared_n
                and counted.completed_leaves_retained == 0,
                'COUNTED',
                count=counted.count,
                leaf_backend=backend,
                leaves=counted.leaves_counted,
                nodes=counted.nodes,
                depth=counted.depth,
                max_pending_nodes=counted.max_pending_nodes,
                peak_frontier_bits=counted.peak_frontier_bits,
                boolean_vectors=counted.boolean_vectors,
                completed_leaves_retained=0,
            )
        )

    forest = CNF(
        9, [[(3 * i, False), (3 * i + offset, False)] for i in range(3) for offset in (1, 2)]
    )
    counted = count_cnf_streaming(forest, M=2, leaf_backend='tensor')
    rows.append(
        _row(
            'lean-frontier-discards-completed-leaves',
            forest.n,
            forest.m,
            counted.counted
            and counted.count == len(_models(forest)) == 125
            and counted.leaves == 8
            and counted.max_pending_nodes == 4
            and counted.depth == 3
            and counted.completed_leaves_retained == 0,
            'COUNTED',
            count=counted.count,
            leaves=counted.leaves,
            nodes=counted.nodes,
            depth=counted.depth,
            max_pending_nodes=counted.max_pending_nodes,
            completed_leaves_retained=counted.completed_leaves_retained,
        )
    )

    for formula, expected, name in (
        (CNF(8, []), 256, 'free-bits'),
        (CNF(8, [[]]), 0, 'empty-clause'),
    ):
        counted = count_cnf_streaming(formula)
        rows.append(
            _row(
                f'lean-{name}-without-enumeration',
                formula.n,
                formula.m,
                counted.count == expected
                and counted.boolean_vectors == 0
                and counted.count_word_bits <= formula.n + 1,
                'COUNTED',
                count=counted.count,
                boolean_vectors=counted.boolean_vectors,
                count_word_bits=counted.count_word_bits,
            )
        )

    rows.append(
        _row(
            'lean-binary-count-word',
            8,
            0,
            count_word(0) == b'\x00'
            and count_word(196) == b'\x00\x00\x01\x00\x00\x00\x01\x01'
            and decode_count_word(count_word(196)) == 196,
            'CHECKED',
            count=None,
            zero_word=[0],
            count_word=list(count_word(196)),
        )
    )

    for name, budgets in (
        ('node', StreamingBudgets(max_nodes=1)),
        ('frontier', StreamingBudgets(max_frontier_bits=1)),
        ('search', StreamingBudgets(max_search_candidates=0)),
    ):
        refused = count_cnf_streaming(star, M=2, budgets=budgets)
        rows.append(
            _row(
                f'lean-{name}-preflight-refusal',
                star.n,
                star.m,
                refused.status == 'refused'
                and refused.count is None
                and refused.boolean_vectors == refused.leaves_counted == refused.nodes == 0,
                'REFUSED',
                count=None,
                boolean_vectors=refused.boolean_vectors,
                leaves_counted=refused.leaves_counted,
                reason=refused.reason,
            )
        )

    refused = count_cnf_streaming(star, M=2, total_vector_budget=12)
    rows.append(
        _row(
            'lean-global-leaf-budget-no-partial-total',
            star.n,
            star.m,
            refused.status == 'refused'
            and refused.count is None
            and refused.leaves_counted == 1
            and refused.boolean_vectors == 10,
            'REFUSED',
            count=None,
            boolean_vectors=refused.boolean_vectors,
            leaves_counted=refused.leaves_counted,
            reason=refused.reason,
        )
    )
    for row in rows:
        stamp_completion(row)
        print(f'  [Y] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP23_SCOPE}', flush=True)
    return rows
