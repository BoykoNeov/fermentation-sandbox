---
name: feedback-a-margin-is-a-claim-about-what-holds-it-open
description: "Recording \"X is inert, margin M%\" is incomplete until you find which parameters hold M open — D-165's 0.5% ethanol headroom turned out to rest on two byproduct-yield bands, and a documented toggle closes it"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 10a4c2c5-2d8d-4887-b416-58ee04aef9d0
  modified: 2026-08-09T13:34:51.232Z
---

D-165 recorded that `EthanolToleranceDeath` is exactly inert on a 24-Brix must, with the
peak ensemble ethanol 119.4 g/L against an `ethanol_tolerance` band low of 120.0 — "a 0.5 %
margin", framed as a property of a normal must.

D-166 asked what holds that 0.5 % open, and the answer was not the sugar load. The must
carries 245.3 g/L sugar, so the **Gay-Lussac ceiling is 125.4 g/L — already above the band
low of 120**. The gap exists only because the D-16 realised-yield carbon diversion
(`Y_glycerol_sugar` 0.035, `Y_byproduct_sugar` 0.012) carves carbon out of the ethanol flux.
Both are sampled in every wine ensemble:

- shipped nominal → peak E 117.280, **+2.27 %** margin, DEAD
- both yield bands at their **low** edges → 119.696, **+0.25 %**, DEAD (tighter than the
  figure D-165 quoted, and reachable by two ordinary band-edge draws)
- diversion off, which `uptake.py` advertises as "the theoretical Gay-Lussac core (togglable
  off)" → **122.959**, −2.47 %, **the gate straddles at the reference must**

**Why:** a margin is a difference between two quantities, so it is a claim about *everything
that feeds either side*. Attributing it to the obvious axis (here: Brix, or "a normal must")
leaves the real lever unnamed — and the real lever was a documented toggle, not an exotic
scenario. It also makes the margin look more robust than it is: 0.5 % across an ensemble
sounds like noise, 0.25 % from two named band edges sounds like a design constraint.

**How to apply:** when you record an inertness margin, do not stop at the number. Ask what
sets the quantity on each side, compute its *theoretical* bound (the ceiling the mechanism
could reach if the diverting terms were off), and check whether that bound already crosses
the threshold. If it does, the margin belongs to the diverting parameters — name them, state
their bands and tiers, and test the ordering so a widened band cannot close the gap silently.

Related: [[feedback-pin-the-band-not-the-nominal]],
[[feedback-validate-calibrations-in-the-frame-that-binds]],
[[feedback-measure-which-side-before-building]].
