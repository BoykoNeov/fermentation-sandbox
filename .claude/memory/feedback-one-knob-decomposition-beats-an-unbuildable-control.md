---
name: feedback-one-knob-decomposition-beats-an-unbuildable-control
description: D-245 - when the clean control needs a source you do not have, hold the fixture fixed and sweep the ONE parameter the change moved; the confound becomes a constant instead of an excuse
metadata:
  node_type: memory
  type: feedback
---

**A confound you cannot remove can still be held constant.** D-244 left six red guards with an
asserted cause and a known confound: the fixture runs at 405.4 mg N/L, 2.25x its source paper's
must, and the clean control — the same probe on a commensurate medium — cannot be built, because
that paper's synthetic-medium composition is not in this repo and composing one from a generic
grape-must partition is the category error D-244 §6 already killed. The temptation is a sweep
across musts, which re-opens the same invention. What worked instead: hold the fixture **exactly**
as it is and sweep **only `biomass_N_fraction`**, from the value the old evaluation point produced
to the shipped one. All six thresholds crossed monotonically inside that one-parameter sweep, so
the nitrogen richness is a *constant* in the comparison and cannot be the cause. Six asserted
attributions became one measured one, and the sweep also exposed which leg of each ratio moved.

**Why:** attribution needs a comparison where exactly one thing differs, and people reach for that
by building the ideal fixture. When the ideal fixture needs a number nobody published, the beat
stalls or invents. But the *change under test* is usually a single resolved parameter, and that
parameter is directly settable in the running param map — so the counterfactual can be constructed
from the shipped fixture instead of from a medium. It is cheaper than the control you wanted and
strictly less inventive. Related: [[feedback-measure-which-side-before-building]] and
[[feedback-a-control-needs-mechanical-reach]] — the knob must be the channel the change actually
travelled down, or the sweep proves nothing.

**How to apply:** identify the ONE value the change resolved differently, set it directly in the
resolved parameter dict handed to the runtime (never on a property that re-resolves —
[[feedback-a-parameter-can-be-pinned-and-drawn]]), and sweep it end to end. Report the crossings,
and report every quantity that did **not** move: constants across the sweep are the strongest part
of the result — here, three precursor pools identical to five decimal places said the numerator was
supply-limited and pinned, which is what turned "the draw's share rose" into "the denominator fell".
