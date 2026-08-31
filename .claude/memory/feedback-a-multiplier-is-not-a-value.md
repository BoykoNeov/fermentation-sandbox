---
name: feedback-a-multiplier-is-not-a-value
description: "A sweep keyed on multipliers of a shipped parameter silently becomes a different sweep the day that parameter moves, and nothing goes red"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ab06c4ba-e5e4-4ff7-a3de-6478c4835471
  modified: 2026-08-31T13:01:32.416Z
---

A test grid whose arms are **multipliers** of a compiled parameter reads identically to absolute
values for exactly as long as that parameter is 1.0. The moment it moves, every absolute value in
the file's prose is wrong by that factor, and no assertion notices — the arms still exist, the
ratios between them are unchanged, and the tests stay green.

**Why:** at D-253 `test_assimilable_nitrogen_uptake`'s capacity sweep was written as
`scale={"amino_acid_uptake_capacity_ratio": ratio}` where `_run` multiplies. Moving the shipped
value 1.0 → 2.6 would have turned every "`r = 10`" in that file's docstrings into 26 and every
"200× sweep" into a differently-centred one. The prose would have been false and the suite green.
The hazard is invisible while the parameter sits at 1.0, which is exactly when such a sweep gets
written — a value of 1.0 is what a "this is a bound, not a level" parameter starts at, so the
files most likely to carry this are the ones whose parameter is most likely to be re-landed later.

**How to apply:** when a docstring names an **absolute** value, the code must set an absolute
value. Offer both affordances explicitly and name them for what they do — `override=` replaces,
`scale=` multiplies — rather than having one and meaning the other. Key the arm that means "what
ships" off a constant read from the parameter file, not off a literal and not off `1.0`, and
assert the grid still contains it. A commensurability probe genuinely means "times that ratio" and
should keep multiplying; a calibration sweep never does.

Related: [[feedback-a-parameter-can-be-pinned-and-drawn]], [[feedback-check-the-blocker-is-still-blocking]].
