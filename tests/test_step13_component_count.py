"""Independent component counting, normalization and exact declared-bit counts."""

import random
from dataclasses import replace

import pytest

from counting_sat_algebra.cnf import CNF, brute_count, random_3cnf
from counting_sat_algebra.component_count import count_cnf_components
from counting_sat_algebra.sparse_placement import SparsePlacementBudgets
from counting_sat_algebra.structured_auxiliary import build_structured_auxiliary


@pytest.mark.parametrize('cnf', [
    CNF(0, []), CNF(0, [[]]), CNF(5, []),
    CNF(3, [[(0, False)], [(0, True)]]),
    CNF(8, [[(0, False), (1, False), (2, True)], [(3, False), (4, True), (5, False)]]),
    CNF(4, [[(0, False), (0, True), (1, True)], [(2, False)], [(2, False), (2, False)]]),
])
def test_normalization_components_and_free_bits_match_brute_force(cnf):
    result = count_cnf_components(cnf)
    assert result.counted and result.count == brute_count(cnf)


def test_random_connected_and_disconnected_formulas_match_brute_force():
    rng = random.Random(13)
    for n in range(3, 7):
        cnf = random_3cnf(n, n + 2, rng=rng)
        cnf.n += 2
        result = count_cnf_components(cnf)
        assert result.counted and result.count == brute_count(cnf)


def test_large_disconnected_3cnf_is_counted_without_a_cube_oracle(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('a full-cube oracle was used')

    monkeypatch.setattr('counting_sat_algebra.cnf.brute_count', boom)
    monkeypatch.setattr('counting_sat_algebra.hermite.eval_H_hasse_jet', boom)
    monkeypatch.setattr('counting_sat_algebra.replacements.eval_H_at', boom)
    clauses = [[(3 * i, False), (3 * i + 1, True), (3 * i + 2, False)] for i in range(10923)]
    cnf = CNF(32771, clauses)  # 32769 active variables plus two declared free bits.
    result = count_cnf_components(cnf)
    assert result.counted and result.count == 4 * 7 ** 10923
    assert result.components == result.components_counted == 10923
    assert result.free_variables == 2 and result.largest_component == 3
    assert result.local_evaluations == 1 and result.memo_hits == 10922
    assert result.boolean_vectors <= 500_000


def test_memo_disabled_and_taylor_local_backend_preserve_the_count():
    cnf = CNF(6, [[(0, False), (1, True), (2, False)], [(3, False), (4, True), (5, False)]])
    auxiliary = build_structured_auxiliary(2, 3, 2)
    cached = count_cnf_components(cnf, auxiliary=auxiliary)
    uncached = count_cnf_components(cnf, auxiliary=auxiliary, memo_entries=0)
    assert cached.count == uncached.count == brute_count(cnf) == 49
    assert cached.local_evaluations == 1 and cached.memo_hits == 1
    assert uncached.local_evaluations == 2 and uncached.memo_hits == 0
    assert uncached.boolean_vectors == 2 * cached.boolean_vectors


def test_budget_refusals_do_not_report_a_partial_count(monkeypatch):
    cnf = CNF(5, [[(v, False), (v + 1, False)] for v in range(4)])
    result = count_cnf_components(cnf, max_component_variables=4)
    assert result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    small = CNF(2, [[(0, False), (1, False)]])
    for kwargs in ({'vector_budget': 0}, {'total_vector_budget': 0}, {'max_variables': 1},
                   {'placement_budgets': replace(SparsePlacementBudgets(), max_trials=0)}):
        result = count_cnf_components(small, **kwargs)
        assert result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    with pytest.raises(ValueError, match='nonnegative'):
        count_cnf_components(small, memo_entries=-1)


def test_global_vector_budget_counts_only_executed_local_work():
    cnf = CNF(2, [[(0, False)], [(1, True)]])
    first = count_cnf_components(CNF(1, [[(0, False)]]))
    result = count_cnf_components(cnf, total_vector_budget=2 * first.boolean_vectors)
    assert result.status == 'refused' and result.count is None
    assert result.components_counted == 1 and result.local_evaluations == 2
    assert result.boolean_vectors == first.boolean_vectors
