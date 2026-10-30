"""Step 25 finite witnesses for scheduled sparse codes and greedy selection."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product

from .protocol import _row, stamp_completion
from ..cnf import CNF, clause_poly_bool
from ..lean_code import (
    CodeBudgets,
    SparseTannerCode,
    canonical_hash_candidates,
    certify_sparse_code,
    construct_sparse_code,
    search_local_ternary_matrix,
    search_walsh_shifts,
    sparse_code_parameters,
)
from ..lean_greedy import (
    STEP25_SCOPE,
    GreedyBudgets,
    count_cnf_greedy,
    occurrence_bound,
    select_history_rows,
)
from ..lean_ring import LeanAlgebraBudgets
from ..lean_sparsify import count_cnf_streaming

PARAMETER_CASES = ((4, 1), (16, 1), (16384, 1), (256, 2))
CERTIFICATE_CASES = ((4, 1), (16, 1), (64, 1), (256, 2))


def _brute_models(cnf: CNF) -> int:
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if all(clause_poly_bool(bits, clause) for clause in cnf.clauses):
            total += 1
    return total


def _accepted(cnf: CNF, rows, bits) -> bool:
    bad = [1 - clause_poly_bool(bits, clause) for clause in cnf.clauses]
    return all(sum(weight * bad[index] for index, weight in row.terms) % 3 == 0 for row in rows)


def _brute_mass(cnf: CNF, rows) -> int:
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if _accepted(cnf, rows, bits):
            total += sum(1 - clause_poly_bool(bits, clause) for clause in cnf.clauses)
    return total


def _brute_trajectory(cnf: CNF, initial, depth: int, max_rounds: int):
    rows = list(initial)
    trail = [_brute_mass(cnf, rows)]
    while trail[-1] > 0 and len(trail) <= max_rounds:
        best = best_score = None
        for _rank, candidate in canonical_hash_candidates(cnf.m, depth):
            score = _brute_mass(cnf, [candidate, *rows])
            if best is None or score < best_score:
                best, best_score = candidate, score
        rows = [best, *rows]
        trail.append(best_score)
    return rows, trail


def _local_distance_ok(matrix, columns: int, distance: int) -> bool:
    for weight in range(1, distance):
        for support in combinations(range(columns), weight):
            for coefficients in product((1, 2), repeat=weight):
                if all(
                    sum(row[i] * c for i, c in zip(support, coefficients)) % 3 == 0
                    for row in matrix
                ):
                    return False
    return True


def _radix_first_local(columns: int, rows: int, distance: int):
    digits = [0] * (rows * columns)
    tried = 0
    while True:
        tried += 1
        matrix = tuple(tuple(digits[i * columns : (i + 1) * columns]) for i in range(rows))
        if _local_distance_ok(matrix, columns, distance):
            return matrix, tried
        for position in range(len(digits)):
            digits[position] += 1
            if digits[position] < 3:
                break
            digits[position] = 0
        else:
            return None, tried


def _walsh_ok(dimension: int, shifts, denominator: int) -> bool:
    for character in range(1, 1 << dimension):
        value = sum(-1 if (character & shift).bit_count() % 2 else 1 for shift in shifts)
        if denominator * abs(value) > len(shifts):
            return False
    return True


def _first_walsh(dimension: int, degree: int, denominator: int):
    radix = 1 << dimension
    digits = [0] * degree
    tried = 0
    while True:
        tried += 1
        if _walsh_ok(dimension, tuple(digits), denominator):
            return tuple(digits), tried
        for position in range(degree):
            digits[position] += 1
            if digits[position] < radix:
                break
            digits[position] = 0
        else:
            return None, tried


def _slot_rows(m: int, dimension: int, shifts, local):
    degree = len(shifts)
    count = 2 * (1 << dimension) * len(local)
    output = []
    for index in range(count):
        chronological = count - 1 - index
        vertex, within = divmod(chronological, 2 * len(local))
        side, local_index = divmod(within, len(local))
        terms = []
        for port, coefficient in enumerate(local[local_index]):
            source = (vertex ^ shifts[port]) if side else vertex
            column = source * degree + port
            if column < m:
                terms.append((column, coefficient))
        output.append(tuple(terms))
    return output


def _scheduled_parameters_match(m: int, B: int) -> bool:
    parameters = sparse_code_parameters(m, B)
    logarithm = m.bit_length()
    gap = 200 * B * logarithm
    bound = 8 * gap * gap * logarithm
    degree = min(m, bound)
    dimension = 0
    while degree * (1 << dimension) < m:
        dimension += 1
    return (
        parameters.log_size == logarithm
        and parameters.gap_denominator == gap
        and parameters.degree_bound == bound
        and parameters.degree == degree
        and parameters.dimension == dimension
        and parameters.local_rows == degree // (8 * B) + 1
        and parameters.local_distance == degree // (100 * B * logarithm) + 1
        and parameters.distance_denominator == 20000 * B * B * logarithm * logarithm
        and parameters.total_rows == 2 * (1 << dimension) * parameters.local_rows
    )


def run_lean_greedy_suite() -> list[dict]:
    print('=== (AA) Scheduled sparse codes and violation-mass greedy selection ===', flush=True)
    rows = []
    scheduled = all(_scheduled_parameters_match(m, B) for m, B in PARAMETER_CASES)
    rows.append(
        _row(
            'lean-code-scheduled-integer-parameters',
            0,
            0,
            scheduled,
            'CHECKED',
            count=None,
            cases=[f'm={m},B={B}' for m, B in PARAMETER_CASES],
        )
    )

    identity = all(
        Fraction(1, 100 * B * m.bit_length())
        * (Fraction(1, 100 * B * m.bit_length()) - Fraction(1, 200 * B * m.bit_length()))
        == Fraction(1, sparse_code_parameters(m, B).distance_denominator)
        for m, B in PARAMETER_CASES
    )
    certificate_ok = True
    certificate_details = []
    for m, B in CERTIFICATE_CASES:
        result = construct_sparse_code(m, B)
        parameters, certificate = result.parameters, result.certificate
        valid = (
            result.checked
            and certificate is not None
            and certificate.relative_distance >= Fraction(1, parameters.distance_denominator)
            and parameters.distance_denominator * certificate.minimum_weight >= m
            and certificate.depth <= parameters.distance_denominator
        )
        certificate_ok = certificate_ok and valid
        certificate_details.append(
            dict(
                m=m, B=B, minimum_weight=certificate.minimum_weight, kernel_depth=certificate.depth
            )
        )
    rows.append(
        _row(
            'lean-code-relative-distance-and-kernel-bound',
            0,
            0,
            identity and certificate_ok,
            'CHECKED',
            count=None,
            certificates=certificate_details,
        )
    )

    witness = search_local_ternary_matrix(2, 1, 2)
    first, tried = _radix_first_local(2, 1, 2)
    exhausted, exhausted_tried = _radix_first_local(2, 1, 3)
    rows.append(
        _row(
            'lean-code-local-distance-first-witness',
            0,
            0,
            witness.status == 'checked'
            and witness.witness == first
            and witness.candidates == tried
            and _local_distance_ok(first, 2, 2)
            and exhausted is None
            and search_local_ternary_matrix(2, 1, 3).status == 'not_found',
            'CHECKED',
            count=None,
            witness=[list(row) for row in first],
            candidates=tried,
            exhausted_candidates=exhausted_tried,
        )
    )

    walsh = search_walsh_shifts(1, 2, 2)
    walsh_first, walsh_tried = _first_walsh(1, 2, 2)
    wide = search_walsh_shifts(2, 4, 2)
    wide_first, wide_tried = _first_walsh(2, 4, 2)
    rows.append(
        _row(
            'lean-code-walsh-first-witness',
            0,
            0,
            walsh.status == 'checked'
            and walsh.witness == walsh_first
            and walsh.candidates == walsh_tried
            and wide.status == 'checked'
            and wide.witness == wide_first
            and wide.candidates == wide_tried,
            'CHECKED',
            count=None,
            small_witness=list(walsh_first),
            small_candidates=walsh_tried,
            wide_witness=list(wide_first),
            wide_candidates=wide_tried,
        )
    )

    hand_shifts = (0, 1, 2, 3)
    hand_local = ((1, 2, 1, 0), (0, 1, 2, 1))
    hand = SparseTannerCode(10, 2, hand_shifts, hand_local, 600, 1)
    oracle = _slot_rows(10, 2, hand_shifts, hand_local)
    regenerated = tuple(hand.row(i).terms for i in range(hand.row_count))
    rebuild = SparseTannerCode(10, 2, hand_shifts, hand_local, 600, 1)
    rebuild_rows = tuple(rebuild.row(i).terms for i in range(rebuild.row_count))
    scheduled_regenerate = True
    for m, B in ((4, 1), (16, 1)):
        code = construct_sparse_code(m, B).certificate.code
        scheduled_regenerate = (
            scheduled_regenerate
            and code.row_count == sparse_code_parameters(m, B).total_rows
            and all(code.row(i).terms == code.rows[i].terms for i in range(code.row_count))
        )
    rows.append(
        _row(
            'lean-code-slot-rows-regenerate',
            0,
            0,
            regenerated == tuple(oracle)
            and rebuild_rows == regenerated
            and len(regenerated) == hand.row_count == 2 * (1 << 2) * len(hand_local)
            and all(
                len(terms) <= hand.degree
                and all(
                    0 <= column < 10 and coefficient in (0, 1, 2) for column, coefficient in terms
                )
                for terms in regenerated
            )
            and scheduled_regenerate,
            'CHECKED',
            count=None,
            hand_rows=hand.row_count,
            hand_terms=[len(terms) for terms in regenerated],
        )
    )

    arithmetic_ok = True
    for m, B in CERTIFICATE_CASES:
        certificate = construct_sparse_code(m, B).certificate
        relative = certificate.relative_distance
        weight = max(1, (m * relative.numerator + relative.denominator - 1) // relative.denominator)
        arithmetic_ok = (
            arithmetic_ok
            and weight == certificate.minimum_weight
            and certificate.depth == (m + weight - 1) // weight
            and certificate.row_term_bound == certificate.code.degree
            and certificate.column_incidence_bound == 2 * len(certificate.code.local_matrix)
        )
    rows.append(
        _row(
            'lean-code-certificate-arithmetic',
            0,
            0,
            arithmetic_ok,
            'CHECKED',
            count=None,
        )
    )

    good = SparseTannerCode(2, 1, (1, 0), ((1, 1),), 600, 2)
    good_certificate = certify_sparse_code(good)
    rejected = 0
    for bad in (
        SparseTannerCode(2, 1, (1, 0), ((1, 1),), 600, 3),
        SparseTannerCode(2, 1, (1, 1), ((1, 1),), 600, 2),
    ):
        try:
            certify_sparse_code(bad)
        except ValueError:
            rejected += 1
    rows.append(
        _row(
            'lean-code-certify-accept-and-reject',
            0,
            0,
            good_certificate.relative_distance == Fraction(599, 600)
            and good_certificate.minimum_weight == 2
            and good_certificate.depth == 1
            and rejected == 2,
            'CHECKED',
            count=None,
            rejected=rejected,
        )
    )

    scheduled_refused = construct_sparse_code(16384, 1)
    big = construct_sparse_code(1 << 20, 1)
    small = search_local_ternary_matrix(4, 2, 2, budgets=CodeBudgets(max_candidates=1))
    shifts = search_walsh_shifts(3, 4, 2, budgets=CodeBudgets(max_candidates=1))
    rows.append(
        _row(
            'lean-code-budget-refusals',
            0,
            0,
            scheduled_refused.status == 'refused'
            and scheduled_refused.certificate is None
            and scheduled_refused.candidates == 0
            and big.status == 'refused'
            and big.certificate is None
            and big.candidates == 0
            and small.status == 'refused'
            and small.witness is None
            and shifts.status == 'refused'
            and shifts.witness is None,
            'REFUSED',
            count=None,
            scheduled_reason=scheduled_refused.reason,
            large_reason=big.reason,
            local_reason=small.reason,
            walsh_reason=shifts.reason,
        )
    )

    stream = list(canonical_hash_candidates(3, 2))
    width = (3 * 3).bit_length()
    expected_terms = {
        tuple((digit % 3, digit // 3) for digit in digits) for digits in product(range(9), repeat=2)
    }
    ranks = [rank for rank, _ in stream]
    rows.append(
        _row(
            'lean-code-canonical-candidate-stream',
            3,
            0,
            len(stream) == 81
            and ranks == sorted(ranks)
            and len(set(ranks)) == 81
            and {row.terms for _, row in stream} == expected_terms
            and all(
                rank
                == sum(
                    (index + 3 * weight) << (width * position)
                    for position, (index, weight) in enumerate(row.terms)
                )
                for rank, row in stream
            ),
            'CHECKED',
            count=None,
            candidates=len(stream),
        )
    )

    star = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
    trajectory = count_cnf_greedy(star, depth=1)
    initial = trajectory.rows[trajectory.rounds :]
    brute_rows, brute_trail = _brute_trajectory(star, initial, 1, trajectory.max_rounds)
    rows.append(
        _row(
            'lean-greedy-star-matches-brute-trajectory',
            star.n,
            star.m,
            trajectory.counted
            and trajectory.count == _brute_models(star) == 9
            and trajectory.final_mass == 0
            and trajectory.accepted_queries == 1
            and trajectory.rows == tuple(brute_rows)
            and brute_trail[0] == trajectory.initial_mass
            and len(brute_trail) - 1 == trajectory.rounds
            and _brute_mass(star, trajectory.rows) == 0,
            'COUNTED',
            count=trajectory.count,
            rounds=trajectory.rounds,
            initial_mass=trajectory.initial_mass,
            final_mass=trajectory.final_mass,
            marked_queries=trajectory.marked_queries,
            occurrence_bound=occurrence_bound(star),
            brute_trail=brute_trail,
        )
    )

    selected = select_history_rows(star, depth=1)
    accepted_set = {
        bits for bits in product((0, 1), repeat=star.n) if _accepted(star, selected.rows, bits)
    }
    model_set = {
        bits
        for bits in product((0, 1), repeat=star.n)
        if all(clause_poly_bool(bits, clause) for clause in star.clauses)
    }
    rows.append(
        _row(
            'lean-greedy-zero-mass-accepts-exactly-models',
            star.n,
            star.m,
            selected.selected
            and selected.final_mass == 0
            and accepted_set == model_set
            and len(model_set) == 9,
            'CHECKED',
            count=None,
            accepted=len(accepted_set),
            models=len(model_set),
        )
    )

    demo = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
    forest = CNF(
        6,
        [
            [(0, False), (1, False)],
            [(0, False), (2, False)],
            [(3, False), (4, False)],
            [(3, False), (5, False)],
        ],
    )
    streamed_ok = True
    streamed_details = []
    for name, cnf, expected in (('star', star, 9), ('demo', demo, 196), ('forest', forest, 25)):
        streamed = count_cnf_streaming(cnf, M=2, leaf_backend='hash')
        valid = (
            streamed.counted
            and streamed.count == expected == _brute_models(cnf)
            and streamed.hash_final_mass == 0
            and streamed.hash_accepted_queries == streamed.leaves_counted
            and streamed.boolean_vectors == 0
            and streamed.same_declared_n
        )
        streamed_ok = streamed_ok and valid
        streamed_details.append(
            dict(
                name=name,
                count=streamed.count,
                leaves=streamed.leaves_counted,
                marked_queries=streamed.hash_marked_queries,
                rounds=streamed.hash_rounds,
            )
        )
    rows.append(
        _row(
            'lean-greedy-streamed-hash-leaves',
            0,
            0,
            streamed_ok,
            'COUNTED',
            count=None,
            cases=streamed_details,
        )
    )

    budget_capped = count_cnf_greedy(star, depth=1, budgets=GreedyBudgets(max_candidates=1))
    round_capped = count_cnf_greedy(forest, depth=1, budgets=GreedyBudgets(max_rounds=1))
    algebra_capped = count_cnf_greedy(
        star, depth=1, native_budgets=LeanAlgebraBudgets(max_scalar_operations=2000)
    )
    rows.append(
        _row(
            'lean-greedy-refusals-without-partial-count',
            0,
            0,
            budget_capped.status == 'refused'
            and budget_capped.count is None
            and budget_capped.marked_queries == 0
            and budget_capped.candidates_scanned == 0
            and budget_capped.code_status == ''
            and round_capped.status == 'refused'
            and round_capped.count is None
            and round_capped.rounds == 1
            and round_capped.final_mass > 0
            and round_capped.marked_queries > 0
            and algebra_capped.status == 'refused'
            and algebra_capped.count is None
            and algebra_capped.count_word_bits == 0
            and algebra_capped.marked_queries > 0,
            'REFUSED',
            count=None,
            budget_reason=budget_capped.reason,
            round_reason=round_capped.reason,
            algebra_reason=algebra_capped.reason,
            round_marked_queries=round_capped.marked_queries,
        )
    )

    free = count_cnf_greedy(CNF(8, []))
    empty = count_cnf_greedy(CNF(8, [[]]))
    single_cnf = CNF(3, [[(0, False)]])
    single = count_cnf_greedy(single_cnf)
    rows.append(
        _row(
            'lean-greedy-terminal-cases',
            0,
            0,
            free.counted
            and free.count == 256
            and free.rounds == 0
            and free.initial_mass == 0
            and free.marked_queries == 0
            and empty.counted
            and empty.count == 0
            and empty.accepted_queries == 1
            and single.counted
            and single.count == _brute_models(single_cnf) == 4,
            'COUNTED',
            count=None,
            free_bits=free.count,
            empty_clause=empty.count,
            single_clause=single.count,
        )
    )

    shared = count_cnf_streaming(
        star, M=2, leaf_backend='hash', greedy_budgets=GreedyBudgets(max_marked_queries=50)
    )
    rows.append(
        _row(
            'lean-greedy-shared-budget-no-partial-total',
            star.n,
            star.m,
            shared.status == 'refused'
            and shared.count is None
            and shared.count_word_bits == 0
            and shared.hash_marked_queries == 50
            and shared.leaves_counted == 1,
            'REFUSED',
            count=None,
            marked_queries=shared.hash_marked_queries,
            leaves_counted=shared.leaves_counted,
            reason=shared.reason,
        )
    )

    for row in rows:
        stamp_completion(row)
        print(f'  [AA] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP25_SCOPE}', flush=True)
    return rows
