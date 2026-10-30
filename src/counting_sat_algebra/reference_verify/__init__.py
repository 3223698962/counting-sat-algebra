"""Executable finite alignment checks for reference Sections 3–6 and Appendix D.

Section 5 checks are the finite monomial, Lemma 5.1, and truncated-product
objects, plus the finite Section 5.3–5.5 replacement sum in A_ell. Section 6
checks are budgeted equation selection and disjoint sparsification on a fixed
n. Step 20 adds Appendix D coefficient rings, networks and local recovery.
Step 21 adds finite collision packing and prospective child topologies.
Step 22 adds finite D.4.5 weighted modes. Step 23 checks Lean-compatible
syntax and a streamed disjoint frontier. Step 24 checks constructed ternary
indicator digits, shared marker truncation and native accepted/marked
queries. Step 25 checks the scheduled sparse ternary code, its regenerated
rows and the violation-mass greedy selector. Step 26 checks the early
high-incidence threshold, the high-variable frontend and the streamed
conditioned rank loop. Step 27 checks the counting extension ring, the fixed
nodes and the Q*J > degree Hermite condition. Step 28 checks the streamed
correction pairs and the Taylor readback. The active roadmap follows the
supplied SqrtLogSpace program; packing and the whole-machine resource refinement remain open. These checks do not establish Python asymptotic bounds or replace
Lean kernel verification. Appendix D packet tensors and recursive modular
contraction remain open.

Each ``run_*_suite`` family lives in its own module. This package re-exports
the entry points imported by ``verify`` and the tests.
"""

from .appendix_d import run_appendix_d_suite
from .asymptotic import run_asymptotic_family_suite
from .common import run_common_count_suite
from .component import run_component_count_suite
from .design import run_reference_design_suite
from .lean_ext import run_lean_ext_suite
from .lean_greedy import run_lean_greedy_suite
from .lean_high import run_lean_high_suite
from .lean_streaming import run_lean_streaming_suite
from .lean_native import run_lean_native_suite
from .lean_taylor import run_lean_taylor_suite
from .packet_packing import run_packet_packing_suite
from .report import write_reference_report
from .section5 import run_section5_suite
from .section5_replacement import run_section5_replacement_suite
from .section6 import run_section6_suite
from .sparse_placement import run_sparse_placement_suite
from .structured import run_structured_auxiliary_suite
from .support import run_support_count_suite
from .tensor import run_tensor_count_suite
from .weighted_modes import run_weighted_modes_suite

REFERENCE_SUITES = (
    (9, run_reference_design_suite),
    (10, run_common_count_suite),
    (11, run_structured_auxiliary_suite),
    (12, run_sparse_placement_suite),
    (13, run_component_count_suite),
    (14, run_support_count_suite),
    (15, run_tensor_count_suite),
    (16, run_asymptotic_family_suite),
    (17, run_section5_suite),
    (18, run_section5_replacement_suite),
    (19, run_section6_suite),
    (20, run_appendix_d_suite),
    (21, run_packet_packing_suite),
    (22, run_weighted_modes_suite),
    (23, run_lean_streaming_suite),
    (24, run_lean_native_suite),
    (25, run_lean_greedy_suite),
    (26, run_lean_high_suite),
    (27, run_lean_ext_suite),
    (28, run_lean_taylor_suite),
)

__all__ = [
    'REFERENCE_SUITES',
    'run_appendix_d_suite',
    'run_asymptotic_family_suite',
    'run_common_count_suite',
    'run_component_count_suite',
    'run_lean_streaming_suite',
    'run_lean_ext_suite',
    'run_lean_greedy_suite',
    'run_lean_high_suite',
    'run_lean_native_suite',
    'run_lean_taylor_suite',
    'run_packet_packing_suite',
    'run_reference_design_suite',
    'run_section5_replacement_suite',
    'run_section5_suite',
    'run_section6_suite',
    'run_sparse_placement_suite',
    'run_structured_auxiliary_suite',
    'run_support_count_suite',
    'run_tensor_count_suite',
    'run_weighted_modes_suite',
    'write_reference_report',
]
