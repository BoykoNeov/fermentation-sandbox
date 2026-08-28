---
name: feedback-a-park-slot-inherits-its-physics
description: "A stand-in parked in an existing slot inherits everything that slot participates in; D-248 parked intracellular nitrogen in `N` and it kept titrating the must"
metadata:
  node_type: memory
  type: feedback
---

D-248 needed somewhere to hold assimilable nitrogen the yeast had transported in ahead of demand.
It parked it in the `N` slot — and `N` is read by the acid-base charge balance at the must's mean
charge per mole of nitrogen (D-209). So nitrogen already **inside a cell** went on titrating the
liquid around it: `N` 300 → 436.8 mg N/L and pH 3.030 → 3.216 on a 2 g/L amino-acid-dosed wine,
**+0.215** against the same run with uptake off.

**Mass was never the issue.** Nitrogen conserved to 1e-8 in every arm; the slot exceeded its own
starting *share*, never the must's declared total. The defect was **charge** — and "the pool went
up" reads like creation, so say which conserved quantity is intact before anyone chases it.

**The invariant was already written down, in a place nobody re-reads.**
`nitrogen_charge_excess`'s docstring (D-210) says the excess is *"constant except at dose events …
only an addition of differently-charged nitrogen moves it. No Process touches this slot."* Uptake
**is** such an addition. The prohibition existed; the beat that broke it never opened the file.

**How to apply:** before parking a stand-in quantity in an existing slot, enumerate **every reader
of that slot**, not just the writers you are joining — an equilibrium solve, a Monod gate, a
conservation weight, a charge balance. Grep the slot NAME across `core/`, and read the docstring
of any slot that exists to correct the one you are about to write. If a reader would give the
stand-in physics it should not have, the slot is wrong and a new one is the honest cost.
Related: [[feedback-a-declared-quantity-can-have-a-second-channel]], [[feedback-a-lump-is-right-in-total-and-wrong-in-fate]].
