---
name: prohibitions-osmotic-high-sugar
description: "D-192 — the high-sugar osmotic/substrate brake: what is BUILT, what is deliberately forfeited, and which shapes are refuted by measurement"
metadata: 
  node_type: memory
  type: project
  originSessionId: 2619e25d-f5ec-4f22-8c72-48d1de2dc5a3
  modified: 2026-08-12T08:23:01.114Z
---

# Osmotic / substrate inhibition at high sugar — D-192

**SETTLED. Do not re-propose any of this as unbuilt.** Read before touching substrate
inhibition, the Coleman envelope, or the sweet-must scenarios.

## BUILT

`OsmoticSubstrateInhibition` — a **wine-only** `RateModifier` in
`core/kinetics/osmotic.py`, wired via `_OSMOTIC_INHIBITION_MODIFIERS`. Scales **uptake AND
growth AND the amino-acid swap** (the swap must ride with growth — D-32's coupling).

```
f(S) = 1                                        for S <= 300 g/L   (literal 1.0, early return)
f(S) = 1 / (1 + ((S - 300)/74.6)**2)             above
```

Effect: an 881 g/L must goes 19.76 % → **5.59 % ABV** (shipped n=2; 2.00 % at n=6), +199 g/L
residual, and is **genuinely STUCK** — unchanged from 1 y to 20 y, viable biomass ~0.
**At 32-40 °Brix the brake changes the PATH, not the DESTINATION**: mid-run sugar differs by
15/94/241 g/L while final ABV moves 0.0006/0.005/0.030 %. Never call it "negligible" there —
it is engaged; only the endpoint is unmoved (a supply-limited flux still arrives).

## NEVER re-litigate

- **The threshold is 300 g/L, NOT the Handbook's printed 200 g/L onset.** 200 is *below every
  normal must* (24 °Brix = 245 g/L). The model's `mu_max`/`K_n`/`k'_d` are Coleman 2007's,
  fitted over 265-300 g/L **with no substrate term**, so inhibition below 300 is already inside
  them. Lowering the threshold double-counts it and destroys the byte-for-byte inertness.
- **The Handbook's "300 g/L can yield less alcohol than 200 g/L" is DELIBERATELY FORFEITED.**
  Measured: reaching it needs a Haldane `K_i ≈ 17`, which is a **92 % brake at 200 g/L** (a
  global rate cut, not high-sugar inhibition) and puts the Coleman RMSE at **170.9 against a
  2.0 threshold**. Pinned as its own test so it reads as a decision.
- **The Haldane form (Ghose & Tyagi eq. 23.12, the literature's own) is REFUTED ON SHAPE.** Its
  group `S²/(K_S+S)` grows only **1.51×** from 200→300 g/L; flat-at-200-and-biting-at-300 needs
  **~19×** (steepness n ≈ 7.3). One smooth constant cannot do both. Don't rebuild it.
- **Never make the far anchor a hard zero.** A wall makes a sweet must an **absorbing state**
  (uptake 0, growth nitrogen-capped, nothing removes sugar — it could never ferment, ever).
  Real must at that concentration ferments over *years*.
- **`n` has a hard floor of 2** (C¹ smoothness; below it there's a derivative corner where the
  brake engages and the modifier *raises*). Not an uncertainty bound.
- **`K` and `n` are a DERIVED PAIR, not independent bands.** `K = 325/((1/0.05-1)**(1/n))`. The
  admissible set is a *curve*, not the box. Sweeping either at the other's edge gives an
  off-curve pair that misses the 625 g/L anchor — a different model, not a weaker one.
- **Beer is inert by NOT BEING WIRED**, which is stronger than a parameter value. It also
  *cannot* be wired by accident: `ProcessSet` refuses a modifier naming an absent target and
  beer has no `AminoAcidAssimilation`.
- **All three constants sit on band edges on purpose** and are in the D-166 census — and unlike
  the copper ones they ARE sampler-reachable, so "harmless because unread" does not apply.

## Corrections this record shipped

- **D-129's "late, not stuck" is FALSE** at high enough sugar: measured to 20 years the arm
  converges and STAYS there, +199 g/L residual, because ethanol kills the biomass before the
  brake lifts. Always quote the SHIPPED both-sites numbers (5.59/2.00), not an uptake-only arm.
- **D-14's Coleman margin was overstated 4×**: the docstring claimed "~1.3 g/L ... ~50 % margin";
  measured **1.7959 / 1.5956** against 2.0 = **~11 %**. Corrected in place. There is not 0.2 g/L
  of room inside that envelope for any new term.

## OPEN — owner's call, deliberately not folded in

**`_SWEET_BRIX = 70.0` in `test_aging_scenario.py` is not a botrytis must** (947 g/L raw;
Sauternes is ~35-40 °Brix). All eleven scenarios still pass, but they now age a **1.4 % ABV**
wine instead of 19.8 % — eleven tests whose names say sweet WINE. Re-anchoring it is a *scenario* decision, not physics.

**Tokay Aszu figures are NOT a validation anchor** — confounded by cold cellars and
botrytis antifungals. Context only.

See [[feedback-measure-which-side-before-building]], [[feedback-a-band-is-per-parameter-a-claim-is-joint]],
[[feedback-a-form-can-be-too-gradual-for-its-own-anchors]].
