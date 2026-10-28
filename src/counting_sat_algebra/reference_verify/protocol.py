"""Shared labels and standing honesty flags for finite reference checks.

Success labels are COUNTED, CHECKED, and REFUSED. A check that does not hold
is recorded as FAIL. ``clean_refusal`` is the budget-refusal shape used by
the suites: status ``refused``, no count, and nothing scanned.

``theorem_1_1_complete``, ``lean_proof_complete``, and ``recursive_contraction``
stay false. This module does not implement Theorem 1.1 and does not refute ETH.
There is no general 2^{o(n)} algorithm here.
"""

from __future__ import annotations

COUNTED = 'COUNTED'
CHECKED = 'CHECKED'
REFUSED = 'REFUSED'
FAIL = 'FAIL'

THEOREM_1_1_COMPLETE = False
LEAN_PROOF_COMPLETE = False
RECURSIVE_CONTRACTION = False

# Insertion order is part of the reference JSON. Do not reorder these keys.
COMPLETION_FLAGS = {
    'theorem_1_1_complete': THEOREM_1_1_COMPLETE,
    'lean_proof_complete': LEAN_PROOF_COMPLETE,
}


def _row(name, n, m, ok, label, **details):
    """One finite reference row. The success label is kept only when the check holds."""
    return dict(
        name=name,
        n=n,
        m=m,
        brute=details.pop('brute', -1),
        ok=ok,
        label=label if ok else FAIL,
        **details,
    )


def clean_refusal(status: str, count, scanned: int) -> bool:
    """True when a budget refusal recorded no count and scanned nothing."""
    return status == 'refused' and count is None and scanned == 0


def completion_fields() -> dict:
    """Section 6 row fields for the two completion flags. Both stay false."""
    return dict(COMPLETION_FLAGS)


def stamp_completion(record: dict) -> dict:
    """Force the Appendix D honesty flags false on one row.

    The update order matches the historical call: ``recursive_contraction``,
    then the two completion flags. Keys that are already present keep their
    position; missing keys are appended in that order.
    """
    record.update(
        recursive_contraction=RECURSIVE_CONTRACTION,
        theorem_1_1_complete=THEOREM_1_1_COMPLETE,
        lean_proof_complete=LEAN_PROOF_COMPLETE,
    )
    return record
