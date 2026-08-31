---
name: feedback-solver-noise-size-is-not-portable
description: "The magnitude AND the sign of a solver-noise excursion are properties of one trajectory, never of the model; guards that pin either will fire on an unrelated change"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: ab06c4ba-e5e4-4ff7-a3de-6478c4835471
  modified: 2026-08-31T13:02:26.874Z
---

A tiny negative excursion in a state variable, or a tiny non-zero residue in one that "should" be
empty, is a property of **that integration** — that must, that schedule, that output grid, that
tolerance. It is not a property of the model, so neither its size nor its sign travels.

**Why:** D-253 hit both halves in one beat.

- **Size.** D-251 measured `stored_nitrogen` dipping to −2.3e-9 on Crépin's must and recorded
  `assert_nonnegative`'s 1e-9 default as "within a factor of two of firing". Landing the uptake
  capacity moved that must's dip *inward* to −1.4e-10 — and fired the guard on a **different**
  must (the 2 g/L dosed one) at −2.36e-9. The brittleness did not go away, it moved. Measuring it
  on one scenario licensed no statement about another.
- **Sign.** A guard asserted three Strecker aldehydes were *identically* `0.0` with no substrate.
  It was reading which side of zero the precursor's own residue landed on: leucine ends at
  −6.5e-11, so its gate clamps and gives an exact zero; isoleucine ends at +3.0e-15 and leaks
  3.6e-17. The capacity move flipped one, with nothing about the mechanism changed.

**How to apply:** never pin `== 0.0` on a `solve_ivp` output — bound it by a floor chosen from a
**physical** scale (the sensory threshold, the assay's detection limit) and say how many orders of
margin that is, then add a non-vacuity arm on the premise so the floor cannot pass while the
premise quietly stops holding. When a noise excursion does cross a guard, the discriminator is
never the size: tighten the solver. Noise collapses with tolerance while the quantity's real peak
stays put; a genuine defect is a property of the derivative and survives any tolerance. If a bound
must be loosened, pair it in the same test with a tightened run checked at the *un*-loosened
default — it is the pair that forbids the defect, not either number.

Related: [[feedback-a-new-state-slot-moves-every-tight-pin]], [[feedback-a-guard-that-hardcodes-an-input-cannot-price-it]].
