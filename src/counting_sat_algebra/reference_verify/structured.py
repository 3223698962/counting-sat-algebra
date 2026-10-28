"""Step 11 finite checks: structured auxiliary entries, including full-parameter queries.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction

from .protocol import _row
from ..replacement_matrix import build_auxiliary_slice
from ..structured_auxiliary import build_structured_auxiliary


def run_structured_auxiliary_suite() -> list[dict]:
    print('=== (M) Structured auxiliary operator, including R=32, Q=192 ===', flush=True)
    rows = []
    for r, Q, R in ((2, 2, 2), (2, 3, 2)):
        dense = build_auxiliary_slice(r, Q, R)
        structured = build_structured_auxiliary(r, Q, R)
        ok = all(
            structured.matrix_entry(row, col) == dense.M0[u][y]
            and structured.replacement_entry(row, col) == dense.rows[u][y]
            for u, row in enumerate(dense.labels)
            for y, col in enumerate(dense.labels)
        )
        rows.append(
            _row(
                f'structured-d{dense.dimension}',
                0,
                0,
                ok,
                'CHECKED',
                dimension=dense.dimension,
                entry_comparisons=dense.dimension**2,
            )
        )
        print(
            f'  [M] D={dense.dimension}: all matrix and replacement entries {"CHECKED" if ok else "FAIL"}',
            flush=True,
        )
    full = build_structured_auxiliary(2, 192, 32)
    row = ((0,) * 192, (1,) * 192)
    ok = (
        full.dimension == 2**384
        and full.replacement_entry(row, row) == (Fraction(1),)
        and full.replacement_entry(row, (row[1], row[0])) == (Fraction(0),)
        and all(full.f_values[j] == (-1) ** j for j in full.interval)
    )
    rows.append(
        _row(
            'full-parameter-coordinate-queries',
            0,
            0,
            ok,
            'CHECKED',
            dimension=str(full.dimension),
            measurement_modes=str(full.measurement_mode_count),
            dense_matrix_allocated=False,
            modes_enumerated=False,
            count_checked=False,
        )
    )
    print(
        f'  [M] R=32, Q=192, D=2**384: coordinate queries {"CHECKED" if ok else "FAIL"}', flush=True
    )
    return rows
