---
name: feedback-fit-the-observable-not-the-consequence
description: "Calibrate a parameter on the observable it IS, and keep the downstream consequence out of the fit as an out-of-sample check — the value that scores best downstream is usually the wrong one"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: dcd1ccd2-410e-4bfa-9d67-bd5e911c9dec
  modified: 2026-08-17T14:46:19.347Z
---

When two measured quantities both constrain one parameter, **fit on the one the parameter
directly is, and hold the downstream one out of the fit** as an independent check. Do not ship
the value that maximises the downstream score.

**The case (D-211).** Beer's growth rate could be scored two ways off one figure: measured cell
counts (what a growth rate *is*) and the measured pH course (three mechanisms downstream). The
counts admitted 0.031–0.040 and least-squares picked **0.034**; the pH course was fully inside
only for 0.038–0.048. So 0.040 scored **8 of 8** pH days and 0.034 scored **7 of 8**. Shipping
0.040 was tempting and wrong: D-209 §8 had already recorded that its charge term is a **lower
bound**, because assimilation removes the nitrogen pool's charge but not its buffering, and that
unbuilt half pushes pH **down**. A day-1 residual sitting slightly too *alkaline* is exactly
what the missing term would close. Tuning the rate to erase it would have consumed that headroom
and booked a missing buffer term as a growth rate — an upstream parameter absorbing a downstream
omission, where it can never be found again.

**Why:** a downstream observable is a *sum* of the parameter you are fitting and every term you
have not built. Fitting to it silently assigns the unbuilt terms' share to your parameter, and
the fit **looks better** for it — the score improves precisely because the error is being
absorbed rather than exposed. Held out instead, that same residual is evidence: its **sign**
tells you which unbuilt term it belongs to. This is the identifiability shape of
[[feedback-a-normalisation-is-a-free-parameter]] and [[feedback-count-the-anchors-before-adding-a-parameter]],
seen from the calibration side rather than the design side.

**How to apply:** name which observable the parameter *is* before scoring anything. Fit on that
alone, with a stated objective, and report the downstream check separately with the word
"out-of-sample" attached. If the downstream check is worse at the fitted value than at some
other admissible value, that gap is a **finding to record with its sign**, not a reason to move.
Check the archive first for an already-identified unbuilt term with the matching sign — if one
exists, the residual is corroboration, not error. And when the two admissible sets merely
overlap, say so: an intersection is not a derivation, and picking its midpoint is fitting by eye
([[feedback-a-units-fork-is-not-a-band]]).
