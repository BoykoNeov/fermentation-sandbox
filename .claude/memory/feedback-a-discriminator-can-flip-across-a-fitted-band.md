---
name: feedback-a-discriminator-can-flip-across-a-fitted-band
description: "A test that decides between two hypotheses can fail by reversing its verdict across a fitted parameter's own band — a different failure from proving nothing"
metadata:
  node_type: memory
  type: feedback
---

At D-232 I designed a test to separate two explanations for beer's growth-extent gap: ask what
share of the yeast crop must already have settled out by day 1 versus by day 3. If it needs MORE
gone early than late, flocculation-timed settling is impossible.

I predicted it would fail by being **inert** — unmoved across `mu_max`'s uncertainty band, the way
the extent fold is (0.010 across the whole band). It is the opposite. Day 1 moves **0.348** across
that band and day 3 moves **0.0028**, and the movement **reverses the verdict**: at the low edge
settling looks possible, at the high edge impossible. Day 3 is inert only because it is the peak
the growth fit normalises on.

**Why:** the two failures are not the same and the difference decides what to do next. *Inert*
means the observable carries no information about anything. *Sign-flipping* means it carries
information about the FIT rather than about the hypotheses — the parameter it moves with was
fitted on the very curve being compared against, so the "discriminator" is a readout of where you
chose to stand inside a band. An inert test is dead; a flipping test becomes usable the moment the
band narrows, and is worth pinning rather than forgetting.

**How to apply:** before believing any test that decides between hypotheses, **run it at both
edges of every fitted parameter it touches** and check the verdict is the same at both. Report the
span, not just the central answer. If the verdict flips, say "sign-indeterminate across <param>'s
band", never "inconclusive" — and pin it in a guard that goes red when the edges start agreeing,
because that is exactly when the test becomes evidence. See [[feedback-pin-the-band-not-the-nominal]]
and [[feedback-a-red-set-firing-wider-is-a-finding]].
