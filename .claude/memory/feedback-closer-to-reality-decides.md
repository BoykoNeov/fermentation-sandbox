---
name: feedback-closer-to-reality-decides
description: "When options trade fidelity against smaller blast radius, pick the one closer to reality — do not offer \"hold it unchanged\" as a peer option"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 99326074-907e-4a86-84c2-8cc80588bab0
  modified: 2026-08-11T11:44:45.956Z
---

When I present a choice between an option that is **closer to reality** and one that is
merely **safer / smaller blast radius / leaves current behaviour untouched**, the owner picks
closer to reality. Stated by them at D-178: *"1 is closer to reality, right? we are always
aiming at what is closer to reality."*

**Why:** the project's bar is correspondence with reality, not convenience — it is prime
directive material, already in [[user-boykoneov]]. A "keep it as it is for now" option is
conservatism about the *simulation's* comfort, not about the *physics*. Deferring a change I
have already measured and believe is right just leaves known-wrong behaviour shipped.

**How to apply:** do not offer "hold it unchanged" as a peer option when I have evidence the
change is more faithful. Ship the faithful version and record the size of the consequence.
Reserve a genuine question for cases where **which option is more faithful is itself unclear**
— e.g. a fitted term whose calibration frame does not bind (D-178's peptide buffer), or a
conflict between two printed sources. Those are real forks; "it would move a number a lot" is
not.

This does **not** license skipping measurement or provenance. Closer-to-reality still has to
be *demonstrated* — measured, sourced, tiered — and the size of the change recorded. It also
does not license silently widening scope: state the consequence plainly, then proceed.

A surprising magnitude, or my own framing turning out backwards, is a reason to **re-check the
reasoning**, not a reason to retreat to the status quo. At D-178 I predicted beer's ethyl
acetate would age *faster* under acid catalysis and measured a 5–20× *slowdown*; the chemistry
was right (the rate is anchored to a wine at pH 3.3, and beer is less acidic than that) and
only my framing was inverted.

Related: [[feedback-discuss-disagreements]], [[feedback-mutate-the-premise-before-building-the-guard]].
