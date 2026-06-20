"""Validation harness: conservation laws and benchmark curves as first-class data.

Two complementary disciplines from the handoff:

* **Conservation as an invariant** (section 1.6): carbon, nitrogen, and mass
  balances must hold to tolerance. A model that quietly creates carbon is broken
  regardless of how good its curves look. :mod:`conservation` provides
  kinetics-agnostic checkers; a model supplies the conserved-quantity function.

* **Benchmark curves** (section 2.2): the acceptance criteria for Milestone 1,
  encoded as data *before* the model is tuned. :mod:`benchmarks` holds the
  declarative specs and a comparator that also ingests real measured series when
  such datasets become available.
"""

from fermentation.validation.benchmarks import (
    BENCHMARKS,
    BenchmarkSpec,
    ReferenceSeries,
    compare_series,
)
from fermentation.validation.conservation import (
    assert_conserved,
    assert_nonnegative,
    max_drift,
)

__all__ = [
    "BENCHMARKS",
    "BenchmarkSpec",
    "ReferenceSeries",
    "assert_conserved",
    "assert_nonnegative",
    "compare_series",
    "max_drift",
]
