---
name: feedback-a-parameter-can-be-pinned-and-drawn
description: "A parameter used at compile time AND at runtime is only half-sampled, and the surviving half can have the wrong sign"
metadata:
  node_type: memory
  type: feedback
---

When one parameter is read at **compile time** (to build the initial state) *and* at **runtime**
(inside a rate law), sampling it moves only the runtime half. The compile-time half stays pinned at
the nominal, because an ensemble re-draws parameters but **not** `y0`. If the two halves were
designed to cancel, sampling breaks the cancellation in one direction — and the surviving term can
carry the **opposite sign** to what the parameter's name promises.

**Why:** at D-206, `must_aa_fraction_methionine` splits the amino-acid dose at compile time and
scales the depletion gate at runtime. Those cancel exactly (`gate = (dose/Σf)/(K + dose/Σf)`), which
is the documented design. But a *drawn* fraction moved only `K·f`, so "this must has less
methionine" reached the model as a **looser gate** and produced **more** methional: **+7.14 %** at
the band floor, **13.4 %** across the band, while the honest change — recompiling so the pool really
moves, 1.605 → 6.381 mg/L, a **4×** span — moved the output **0.06 %**. The channel that should
matter is silent; the one that moves it runs backwards.

**How to apply:** before attributing anything to a sampled parameter, ask **where else it is read**,
and specifically whether anything consumes it *before* `y0` is built. Then measure the two arms
separately and pair them: the **sampler** arm (fixed `y0`, patched parameter map — what the ensemble
actually does) against a **control** arm that recompiles so every role moves together. Same
parameter, same Process, same scenario, one step earlier — the same channel
[[feedback-a-control-needs-mechanical-reach]]. If they disagree in size or sign, the parameter is
not doing what its name says, and a census run on an *un-dosed* scenario cannot see it (D-166 listed
this whole set as structurally inert because the input was never dosed —
[[feedback-check-the-blocker-is-still-blocking]]). Generalise from the **registry**, never from the
one member you noticed [[feedback-a-magnitude-defence-can-be-scoped-to-the-smallest-member]].
