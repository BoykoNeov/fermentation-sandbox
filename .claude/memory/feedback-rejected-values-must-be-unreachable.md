---
name: feedback-rejected-values-must-be-unreachable
description: "When rejecting a value as unphysical, check every machine-readable field it could still enter through — a green suite doesn't prove it's unreachable"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 4c5542b6-994d-42ff-9b5f-a6dbc7d14d50
  modified: 2026-07-20T09:18:09.319Z
---

If a value is rejected because it breaks a physical law or invariant, the
rejection is **not complete until the value is unreachable by any code path**.
Writing it into a neighbouring machine-readable field with a prose "do not use
this" warning does not count.

**Why:** at D-117 the measured `f_non_ehrlich_phenylalanine = 0.975` was proven
to breach carbon conservation, so the parameter shipped a sourced lower bound
(0.531) instead — but the *honest* interval was written as
`uncertainty: {low: 0.53, high: 0.975}` with a `DO NOT ENSEMBLE-SAMPLE` note.
The full 1184-test suite passed and it was pushed. `uncertainty` is not a
comment: `runtime/ensemble.py::sample_parameters` draws **every** parameter over
`[low, high]`, so ~1 draw in 900 would land in the breach. **Tests exercise the
default value; only the sampler ever visits band edges** — so no test could have
caught it, and it was found by re-reading `ensemble.py`.

**How to apply:** ask "what else consumes this field?" before treating a
rejection as done — samplers, ensembles, schema validators, serialisers, config
readers. Prefer *not expressing* an unusable value in a consumed field over
expressing it with a warning. In this repo the fix was a zero-width band
(`hi <= lo` ⇒ the sampler pins it) plus the real value held as a test constant
that asserts the breach, which keeps it recorded and executable but unreachable.
Do **not** instead cap the field at "the largest value that happens to run" —
that is a magic number with no provenance and it falsely claims the truth might
be that value. Related: [[feedback-never-pipe-checks-to-tail]] — both are cases
where a green signal was not evidence of correctness.
