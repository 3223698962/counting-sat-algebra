"""Structured full-parameter coordinate queries against dense small slices."""

from fractions import Fraction
import itertools
from types import SimpleNamespace

import pytest

from counting_sat_algebra.replacement_matrix import build_auxiliary_slice
from counting_sat_algebra.replacements import eval_poly
from counting_sat_algebra.structured_auxiliary import build_structured_auxiliary


@pytest.mark.parametrize('r,Q,R', [(2, 1, 1), (2, 2, 2), (3, 1, 1), (2, 3, 2), (3, 2, 1)])
def test_all_small_matrix_entries_and_polynomial_rows_match(r, Q, R):
    dense = build_auxiliary_slice(r, Q, R)
    structured = build_structured_auxiliary(r, Q, R)
    assert structured.dimension == dense.dimension
    assert structured.measurement_mode_count == len(dense.C_rows)
    for u, row in enumerate(dense.labels):
        assert structured.label_at(u) == row
        for y, column in enumerate(dense.labels):
            assert structured.matrix_entry(row, column) == dense.M0[u][y]
            assert structured.replacement_entry(row, column) == dense.rows[u][y]
    # Descriptors can have a different order from the dense plus-then-minus basis.
    vectors = [tuple(structured.measurement_entry(mode, label) for label in dense.labels) for mode in structured.measurement_modes()]
    assert set(vectors) == set(dense.C_rows)


def test_full_parameters_can_be_queried_without_enumerating_coordinates(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('coordinates or modes were enumerated')

    monkeypatch.setattr('counting_sat_algebra.structured_auxiliary.itertools', SimpleNamespace(product=boom))
    operator = build_structured_auxiliary(2, 192, 32)
    assert operator.dimension == 2 ** 384
    assert operator.A0 == 32 and len(operator.f_coefficients) == 32
    row = ((0,) * 192, (1,) * 192)
    swapped = (row[1], row[0])
    assert operator.replacement_entry(row, row) == (Fraction(1),)
    assert operator.replacement_entry(row, swapped) == (Fraction(0),)
    one = ((0,) + (1,) * 191, (1,) + (0,) * 191)
    polynomial = operator.replacement_entry(row, one)
    assert len(polynomial) <= 32 and eval_poly(list(polynomial), Fraction(0)) == 0
    assert eval_poly(list(polynomial), Fraction(1)) == -operator.matrix_entry(row, one)
    for j in operator.interval:
        assert operator.f_values[j] == (-1) ** j
    with pytest.raises(ValueError, match='no modes enumerated'):
        next(operator.measurement_modes())


def test_invalid_labels_modes_and_intervals_are_rejected():
    operator = build_structured_auxiliary(2, 2, 1)
    with pytest.raises(ValueError, match='length-Q'):
        operator.matrix_entry(((0,), (1,)), ((0,), (1,)))
    with pytest.raises(ValueError, match='index'):
        operator.label_at(-1)
    with pytest.raises(ValueError, match='inside I'):
        operator.measurement_entry(((0, 0, 0),) * 2, ((0, 0), (0, 0)))
    with pytest.raises(ValueError, match='interval'):
        build_structured_auxiliary(2, 2, 1, -1)
    with pytest.raises(ValueError, match='1 <= Q'):
        build_structured_auxiliary(2, 0, 1, 0)


def test_few_exceptional_modes_do_not_walk_the_full_eigenbasis(monkeypatch):
    generated = []

    def counted_product(*args, **kwargs):
        for value in itertools.product(*args, **kwargs):
            generated.append(value)
            yield value

    monkeypatch.setattr('counting_sat_algebra.structured_auxiliary.itertools',
                        SimpleNamespace(product=counted_product, combinations=itertools.combinations))
    operator = build_structured_auxiliary(2, 8, 8, A0=0)
    modes = list(operator.measurement_modes(max_modes=1))
    assert operator.dimension == 65536 and operator.measurement_mode_count == 1
    assert len(modes) == len(generated) == 1
    assert modes[0] == ((-1, 0, 1),) * 8
    empty = build_structured_auxiliary(2, 8, 9, A0=0)
    assert list(empty.measurement_modes(max_modes=0)) == []
    assert len(generated) == 1
