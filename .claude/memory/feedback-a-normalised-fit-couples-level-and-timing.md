---
name: feedback-a-normalised-fit-couples-level-and-timing
description: "When a rate is fitted on a curve normalised by its own endpoint, the level and the timing are ONE observable — so a repair that lands the level moves the timing, and a target can turn out to contradict itself"
metadata:
  node_type: memory
  type: feedback
---

A repair aimed at a measured **level** can be refuted by the **timing** measured on the *same*
data, and the coupling that makes this happen is easy to miss because it lives in the
normalisation rather than in the model.

**The case (D-258).** Beer's growth `mu_max` is fitted by least squares on the growth fraction at
days 1-2, **normalised on each curve's own peak** — a deliberate choice (D-211), because cells/mL
and g/L differ by a conversion nothing sources and the normalisation cancels it exactly. Beer's
growth *extent* separately overshoots Tyrell's counted crop by 1.55x, and the owner re-opened an
oxygen→growth coupling to bring the ceiling down.

Lowering the ceiling to land the counted fold pushed the day-1 **normalised** fraction to
0.494-0.604 against the **0.235-0.448 measured on the same panel of the same figure**. The
target was unreachable not because the mechanism was wrong but because **landing the level
broke the shape**, and both come from one measurement. At a *different* extent target (a printed
4-5x from a textbook) the same ceiling change kept the timing inside the measured spread — so
the obstacle was the target, not the repair.

**Why:** dividing by the endpoint makes the curve's shape a function of where the endpoint is.
Move the ceiling and the same rate reaches its (now nearer) ceiling sooner, so the normalised
fraction at any fixed hour rises. Level and timing stop being two observables you can satisfy
independently; they are one, and a fit made on the normalised curve has already spent the
freedom. The trap is that the fit's own docstring reads like a *timing* calibration — nothing
in it mentions the extent — so the coupling is invisible from the rate's side.
Cf. [[feedback-a-normalisation-is-a-free-parameter]], which is the same mechanism seen from the
reference-point side, and [[feedback-two-measurements-one-figure-can-disagree]], which is the
two-*series* version; this is the one-series version and it is nastier, because there is no
second series to alert you.

**How to apply:** before building anything that moves a **ceiling, endpoint or extent**, ask what
was fitted on a curve normalised by that quantity, and **run the target through the fit's own
normalised observable first** — it is one integration, far cheaper than the build. State the
result *per target*: "the refit re-opens" is usually not a property of the repair but of which
number the repair is aimed at. And when a target fails this check, that is a finding about the
**target's internal consistency**, worth recording as a constraint on *any* future repair of the
same quantity — not a verdict on the mechanism you happened to be testing.
