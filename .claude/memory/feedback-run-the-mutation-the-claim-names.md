---
name: feedback-run-the-mutation-the-claim-names
description: A coverage claim names a specific mutation -- running a DIFFERENT one and inferring is not evidence; and an upper-bound-only pin goes green on the change you declined
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 76396680-20bf-4e86-af29-a50e3b00ae04
  modified: 2026-08-14T07:38:46.705Z
---

**A "nothing in the suite sees this" claim names a mutation. Run *that* one.** D-203 shipped the
sentence *"if a later beat inflates fermentation acetaldehyde or the α-KB residual, this term grows
silently and no test notices"* — inferred from a mutation that **deleted** the term (a
counterfactual gate). Deletion and inflation are different directions through different asserts, so
the sentence was an extrapolation wearing a measurement's clothes. D-204 ran the named one:
`k_alpha_kb_excretion` ×3, **inside its own declared band**, moved the offset exactly **3.00×** and
**both** existing sotolon tests stayed GREEN. The claim was right — and it cost one probe to stop
being a guess. [[feedback-mutate-the-premise-before-building-the-guard]]

**Two design rules the same beat produced, both about a pin going green when it should not:**

- **Assert BOTH edges.** An upper bound alone catches inflation and goes **GREEN on deletion** —
  i.e. on the exact change the record considered and *declined*. The lower edge is what makes a
  **decision** visible to the suite, so a future reversal argues with a red test instead of slipping
  past a green one. [[feedback-pin-the-band-not-the-nominal]]
- **Never share one threshold across arms of different magnitude.** Four arms spanned 0.0253 →
  0.2847; a ceiling loose enough for the largest sits ~10× above the smallest and would be blind to
  a 3× move there. **The tight arms become decorative** and the test reads as coverage it does not
  have. One pin per arm. [[feedback-grep-finds-claims-not-guards]]

**Why:** a guard's whole value is the mutation it would catch, and both failure modes above produce
a green test on the change you most wanted to see. Inferring across mutations is how a suite
acquires assertions that are true and useless.

**How to apply:** before writing the guard, run three arms — the named inflation (expect RED), the
deletion/gate (expect RED on the *other* edge), and an unmutated baseline (expect GREEN,
[[feedback-verify-the-restore-between-mutation-arms]]). Check each RED is the guard's own `assert`
failing on the number it names: the first gate arm here red-ed on a `NameError` at import and was
**discarded, not counted**. Keep the mutation **in the parameter's declared band** so it is a legal
model position, not a broken one. Choose the tolerance for **reach** — state the smallest ratio it
catches — and label it CONSTRUCTED; these are pins on an output, not an uncertainty band. Related:
[[feedback-a-tautology-can-smuggle-an-attribution]] (a value equal across arms *by construction* is
not a finding — do not assert it), [[feedback-pair-the-red-with-an-ordering-preserving-baseline]].
