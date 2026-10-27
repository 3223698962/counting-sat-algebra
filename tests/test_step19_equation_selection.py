"""Finite Section 6 checks. Not Appendix D, not a rho bound, and not Theorem 1.1."""

import json
from fractions import Fraction

from counting_sat_algebra import __version__
from counting_sat_algebra.__main__ import main
from counting_sat_algebra.algebra_count import count_replacement_sum
from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.equation_selection import (
    STEP19_SCOPE, Section6Budgets, build_disjoint_sparsification,
    choose_branching_M, clause_violation_polynomial, construct_initial_matrix,
    count_via_disjoint_sparsification, count_violations_common_zeros,
    expand_binary_row, greedy_select_equations, kernel_support_check, largest_s,
    lemma_6_1_parameters, lemma_6_3_occurrence_bound, normalize_cnf,
    proposition_6_4_schedule, search_local_code, variable_occurrences,
)
from counting_sat_algebra.reference_verify import run_section6_suite
from counting_sat_algebra.sparse_systems import (
    SparseSystem, brute_common_zeros, count_common_zeros,
)
from counting_sat_algebra.sparsify import sparsify


def test_violation_polynomials_and_duplicate_cancellation():
    positive = clause_violation_polynomial(1, [(0, False)])
    assert [set(monomial) for monomial in positive.monomials] == [set(), {0}]
    mixed = clause_violation_polynomial(2, [(0, False), (1, True)])
    assert [set(monomial) for monomial in mixed.monomials] == [{1}, {0, 1}]
    assert clause_violation_polynomial(1, [(0, False), (0, True)]).is_zero
    assert clause_violation_polynomial(2, []).is_one
    repeated = clause_violation_polynomial(1, [(0, False), (0, False)])
    assert [set(monomial) for monomial in repeated.monomials] == [set(), {0}]
    violations = [positive]
    assert expand_binary_row(1, violations, [(0, 1), (0, 1)]).is_zero
    restored = expand_binary_row(1, violations, [(0, 1), (0, 1), (0, 1)])
    assert [set(monomial) for monomial in restored.monomials] == [set(), {0}]
    assert 'Theorem 1.1' in STEP19_SCOPE and 'Appendix D' in STEP19_SCOPE


def test_lemma_6_1_parameters_and_tiny_kernel():
    parameters = lemma_6_1_parameters(3, 3, 1)
    assert (parameters.mu, parameters.P, parameters.D0) == (2, 80_000, 2_560_000)
    assert parameters.q_G == 1 and parameters.h == 0 and parameters.local_search_bits == 0
    built = construct_initial_matrix(3, 3, 1)
    assert built.constructed and built.matrix is not None and built.rows_allocated
    assert built.search_started is False
    matrix = built.matrix
    assert matrix.nrows == 2 and matrix.nrows <= 3 * 3
    assert matrix.max_row_weight == 3 <= parameters.D0
    kernel = kernel_support_check(
        matrix, numerator=parameters.m, denominator=parameters.P, max_boolean_vectors=64,
    )
    assert kernel.ok and kernel.vectors_scanned == 8
    assert kernel.min_nonzero_weight == 2 and kernel.nonzero_kernel_vectors == 3
    assert kernel.min_nonzero_weight * parameters.P >= parameters.m
    refused = kernel_support_check(
        matrix, numerator=parameters.m, denominator=parameters.P, max_boolean_vectors=1,
    )
    assert refused.refused and refused.vectors_scanned == 0 and refused.count is None


def test_lemma_6_1_search_refuses_before_allocation():
    refused = construct_initial_matrix(2048, 2048, 1)
    assert refused.refused and refused.matrix is None
    assert refused.rows_allocated is False and refused.search_started is False
    assert refused.parameters is not None and refused.parameters.local_search_bits == 524_288
    assert 'no local matrix allocated' in refused.reason


def test_local_code_search_is_budgeted_and_deterministic():
    found = search_local_code(3, 2, 2, max_bits=16, max_words=20)
    assert found == ((0, 1, 2), ())
    refused = search_local_code(3, 2, 2, max_bits=2, max_words=20)
    assert refused is None


def test_greedy_count_matches_brute_and_common_zeros():
    cnf = CNF(2, [[(0, False)], [(1, False)]])
    greedy = greedy_select_equations(cnf, list_length=1, B=1)
    assert greedy.counted and greedy.count == brute_count(cnf) == 1
    assert greedy.rounds_completed == 6 and greedy.counts_performed == 24
    assert greedy.boolean_vectors_scanned == 96 and greedy.equations is not None
    system = SparseSystem(cnf.n, 2, 1, greedy.equations)
    common = count_common_zeros(system)
    assert common.counted and common.count == brute_common_zeros(system) == 1
    assert [set(monomial) for monomial in greedy.equations[2].monomials] == [set(), {0}]
    assert greedy.equations[3].is_zero

    def common_counter(equations):
        result = count_common_zeros(SparseSystem(cnf.n, 2, 1, equations))
        return None if result.refused else result.count

    through_step17 = greedy_select_equations(cnf, list_length=1, B=1, counter=common_counter)
    assert through_step17.counted and through_step17.count == 1
    assert through_step17.counts_performed == 24
    assert through_step17.boolean_vectors_scanned == 0
    replacement = count_replacement_sum(SparseSystem(
        cnf.n, 1, 1, tuple(equation for equation in greedy.equations if not equation.is_zero),
    ))
    assert replacement.counted and replacement.count == 1


def test_greedy_budgets_return_no_count():
    cnf = CNF(2, [[(0, False)], [(1, False)]])
    candidate = greedy_select_equations(
        cnf, list_length=1, B=1, budgets=Section6Budgets(max_candidate_lists=1),
    )
    assert candidate.refused and candidate.count is None
    assert candidate.candidates_enumerated == 0 and candidate.boolean_vectors_scanned == 0
    assert candidate.matrix_allocated is False and candidate.equations is None
    cube = greedy_select_equations(
        cnf, list_length=1, B=1, budgets=Section6Budgets(max_boolean_vectors=1),
    )
    assert cube.refused and cube.count is None and cube.boolean_vectors_scanned == 0
    assert cube.candidates_enumerated == 0 and 'no assignments scanned' in cube.reason
    paper_P = lemma_6_1_parameters(cnf.n, cnf.m, 1).P
    assert paper_P == 80_000
    paper = greedy_select_equations(cnf, list_length=paper_P, B=1)
    assert paper.refused and paper.count is None and paper.candidates_enumerated == 0
    assert paper.matrix_allocated is False


def test_normalization_keeps_n_and_does_not_propagate_units():
    normalized = normalize_cnf(CNF(3, [
        [(0, False), (0, False), (1, True)],
        [(0, False), (0, True)],
        [(0, False), (1, True)],
        [(0, False), (1, True), (2, False)],
        [(2, False)],
    ]))
    assert normalized.n == 3
    assert normalized.clauses == [[(2, False)], [(0, False), (1, True)]]
    kept = normalize_cnf(CNF(2, [[(0, False)], [(0, True), (1, False)]]))
    assert kept.clauses == [[(0, False)], [(0, True), (1, False)]]
    assert brute_count(kept) == 1


def test_disjoint_sparsification_partitions_a_small_cnf():
    cnf = CNF(3, [
        [(0, False), (1, False)],
        [(0, False), (2, False)],
        [(1, False), (2, False)],
    ])
    counted = count_via_disjoint_sparsification(
        cnf, k=2, M=2, verify_disjoint=True,
    )
    assert counted.counted and counted.count == brute_count(cnf) == 4
    assert counted.leaf_counts == (3, 1)
    assert counted.leaves == 2 and counted.nodes == 3 and counted.depth == 1
    assert counted.same_declared_n and counted.occurrence_bound == 6
    assert counted.max_occurrence == 1
    assert counted.tree is not None
    assert all(leaf.n == cnf.n for leaf in counted.tree.leaves)
    assert all(max(variable_occurrences(leaf)) <= counted.occurrence_bound for leaf in counted.tree.leaves)
    common = sum(count_violations_common_zeros(leaf) for leaf in counted.tree.leaves)
    assert common == 4
    seen = set()
    for leaf in counted.tree.leaves:
        models = set()
        for mask in range(1 << cnf.n):
            bits = [(mask >> var) & 1 for var in range(cnf.n)]
            if all(
                any((bits[var] == 0) if neg else (bits[var] == 1) for var, neg in clause)
                for clause in leaf.clauses
            ):
                models.add(mask)
        assert seen.isdisjoint(models)
        seen |= models
    assert len(seen) == 4
    direct_empty = count_via_disjoint_sparsification(CNF(3, [[]]), k=2, M=2)
    assert direct_empty.counted and direct_empty.count == 0 and direct_empty.boolean_vectors_scanned == 0
    direct_full = count_via_disjoint_sparsification(CNF(3, []), k=2, M=2)
    assert direct_full.counted and direct_full.count == 8 and direct_full.boolean_vectors_scanned == 0


def test_step5_copying_increases_n_and_section6_does_not():
    source = CNF(2, [
        [(0, False), (1, False)],
        [(0, False), (1, True)],
        [(0, True), (1, False)],
        [(0, True), (1, True)],
    ])
    copied = sparsify(source)
    assert copied.n_after > source.n
    tree = build_disjoint_sparsification(source, k=2, M=2)
    assert tree.built and tree.leaves
    assert all(leaf.n == source.n for leaf in tree.leaves)
    counted = count_via_disjoint_sparsification(source, k=2, M=2, verify_disjoint=True)
    assert counted.counted and counted.count == brute_count(source) == 0
    assert counted.same_declared_n


def test_tree_budget_returns_no_leaves_and_no_count():
    cnf = CNF(3, [
        [(0, False), (1, False)],
        [(0, False), (2, False)],
        [(1, False), (2, False)],
    ])
    refused = count_via_disjoint_sparsification(
        cnf, k=2, M=2, budgets=Section6Budgets(max_tree_nodes=1),
    )
    assert refused.refused and refused.count is None
    assert refused.leaves == 0 and refused.tree is None
    assert refused.boolean_vectors_scanned == 0 and refused.leaf_counts == ()


def test_proposition_6_4_integers_do_not_claim_the_rho_bound():
    schedule = proposition_6_4_schedule(20)
    assert (schedule.lambda_n, schedule.nu_n, schedule.s) == (5, 3, 0)
    assert schedule.exhaustive_search and schedule.sigma is None
    assert schedule.rho_bound_claimed is False
    assert schedule.theorem_1_1_complete is False and schedule.lean_proof_complete is False
    assert largest_s(5, 3) == 0
    assert largest_s(521_284, 19) == 2
    choice = choose_branching_M(2, Fraction(1))
    assert choice.chosen and choice.M == 64 and choice.C == 4 and choice.h == 10
    assert lemma_6_3_occurrence_bound(2, 2) == 6
    earlier = choose_branching_M(2, Fraction(1), Section6Budgets(max_branching_bits=5))
    assert earlier.refused and earlier.M is None


def test_step19_reference_rows_and_cli(tmp_path, capsys):
    rows = run_section6_suite()
    assert [row['label'] for row in rows] == [
        'CHECKED', 'REFUSED', 'COUNTED', 'REFUSED', 'REFUSED', 'REFUSED', 'COUNTED', 'REFUSED', 'CHECKED',
    ]
    assert all(row['ok'] for row in rows)
    assert rows[2]['count'] == 1 and rows[6]['count'] == 4
    assert rows[6]['leaf_counts'] == [3, 1]
    assert all(row['theorem_1_1_complete'] is False and row['lean_proof_complete'] is False for row in rows)
    assert all(row['appendix_d_implemented'] is False for row in rows)
    output = tmp_path / 'step19.json'
    assert main(['verify', '--step19', '--output', str(output)]) == 0
    report = json.loads(output.read_text(encoding='utf-8'))
    assert report['version'] == __version__
    assert report['steps'] == [19] and report['all_checks_passed'] is True
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert 'Appendix D' in report['scope'] and 'Section 6' in report['scope']
    assert 'rho' in report['scope']
    text = capsys.readouterr().out
    assert 'ALL TESTS PASSED' in text
    assert 'Theorem 1.1 is not implemented' in text
    assert '9/9 checks passed' in text
