---
name: feedback-pair-the-arms-before-comparing-spreads
description: "Two ensembles drawn over different name sets draw different members, so comparing their spreads measures the draw order as much as the change"
metadata:
  node_type: memory
  type: feedback
---

Comparing the spread of one ensemble against another that was drawn over a **different set of
names** is not a comparison. `_resolve_sample_names` returns a sorted tuple and the sampler walks
it, so adding a name reorders every subsequent draw: the two arms integrate different members, and
`max − min` over a dozen of them is far too noisy to carry the difference.

**Why:** at D-240 I priced the never-drawn banded parameters by running the ensemble over the
default set and again over `default ∪ undrawn`. It reported the wine ferment's ethanol spread
widening **1.672×** — a big, quotable, wrong number. The same run reported `h2s` *narrowing* to
**0.42×**, which is the identical noise wearing the other sign and is what gave it away. Paired on
**identical draws** — same `only=` set both arms, the only difference being whether `y0` follows
the member — the ferment's real figure is **0.323**. The unpaired number was wrong by ~5× and in
the direction that flattered the beat.

**How to apply:** hold the draw fixed and vary the **one** thing under test. Force the same
`only=` set into both arms and change only the mechanism (here `y0_for_member`), so the members are
the same trajectories with one channel switched. The shipped arm then doubles as the negative
control: at D-240 it reported **exactly 0.000000** spread on every slot, which is what makes the
other column mean something. If the two arms genuinely cannot share a draw, say the number is
noise-dominated rather than quoting it — and never quote a widening measured at n≈12 without
showing what the same statistic does to a channel that should not have moved.
See [[feedback-a-one-parameter-sweep-is-not-the-band]] for the other half of this trap: an
edge-to-edge sweep is an upper bound, not what a triangular sampler actually reaches.
