# D-268 findings — every growth-anchored arm, re-priced on the sourced composition

The pre-registration is `PREREGISTRATION.md`, written before the first number. The probe is
`probe.py` (arms) and `lambda.py` (the over-draw dial); their raw output is `findings.json`.
Every arm below carries the per-species quadrature closure the D-259 fixture asserts, and every
one read 1.0000.

## 0. The reproduction, which licenses everything else

Run on the composition D-259 *stated*, the harness reproduces the archive exactly:

| arm | published | measured here |
|---|---|---|
| growth sink alone, Crépin leucine split | D-260: 27.58 % | **27.58 %** |
| growth sink alone, leucine tracer | D-260: 5.900 % | **5.900 %** |
| joint arm, Crépin leucine split | D-266: 54.50 % | **54.50 %** |
| joint arm, leucine tracer | D-266: 3.617 % | **3.617 %** |
| joint arm, ile / val / thr splits | D-266: 66.2 / 65.2 / 81.2 | **66.22 / 65.25 / 81.19** |

Both compositions are kept live in the suite for exactly this reason: the superseded arms stay
runnable, so the re-pricing is auditable rather than asserted.

## 1. The weights themselves

`w_i`, g precursor drawn per g of biomass built. Stated -> sourced, and the multiplier:

| | lo | mid | hi |
|---|---|---|---|
| leucine | 0.02400 -> 0.03856 (**1.607x**) | 0.03375 -> 0.04338 (1.285x) | 0.04500 -> 0.04820 (1.071x) |
| isoleucine | 0.01600 -> 0.02835 (**1.772x**) | 0.02250 -> 0.03190 (1.418x) | 0.03000 -> 0.03544 (1.181x) |
| valine | 0.01800 -> 0.03151 (1.751x) | 0.02475 -> 0.03545 (1.432x) | 0.03250 -> 0.03939 (1.212x) |
| threonine | 0.01600 -> 0.02435 (1.522x) | 0.02250 -> 0.02739 (1.217x) | 0.03000 -> 0.03043 (**1.014x**) |
| phenylalanine | 0.01400 -> 0.02279 (1.628x) | 0.02025 -> 0.02564 (1.266x) | 0.02750 -> 0.02849 (1.036x) |

The multiplier is largest at the low edge and smallest at the high one, because the stated
bracket's high edge was its closest guess to the measurement. **That is why the re-priced
bracket is narrower as well as higher.**

## 2. P1 and P2 — every split rises, every bracket narrows

D-259's uncorrected (pre-modifier) splits, % of consumed precursor reaching the lump:

| | lo | mid | hi | span |
|---|---|---|---|---|
| leucine | 13.07 -> **19.44** | 17.45 -> **21.35** | 21.97 -> **23.17** | 8.90 -> **3.73** |
| isoleucine | 19.83 -> 30.43 | 25.79 -> 32.97 | 31.64 -> 35.33 | 11.81 -> 4.90 |
| valine | 19.33 -> 29.60 | 24.78 -> 32.12 | 30.20 -> 34.46 | 10.87 -> 4.86 |
| threonine | 35.72 -> 45.74 | 43.81 -> 48.65 | 50.91 -> 51.26 | 15.19 -> 5.52 |
| phenylalanine | 97.20 -> 98.25 | 98.03 -> 98.43 | 98.54 -> 98.59 | 1.34 -> 0.34 |

D-260's corrected (post-modifier) leucine bracket: **21.3-33.7 -> 30.3-35.2 %**, span 12.4 -> 4.9.

## 3. P3 — D-104's 20.9 % moves from the top of the bracket to inside it

D-259 read D-104's unrecorded 20.9 % at the **top** of its uncorrected bracket (13.1-22.0). On
the sourced composition the bracket is 19.4-23.2 and **20.9 falls inside it**, between the low
edge and the middle. D-260's modifier correction still puts the corrected bracket entirely above
it — by 9.4 points now (30.3) rather than 0.4 (21.3). Both of D-260's conclusions survive.

## 4. P4 and P5 — the refusal and the positive finding both survive

* Corrected leucine tops out at **35.22 %** against Crépin's 77-86 %. D-259's and D-260's refusal
  of the growth-anchored sink is untouched.
* Phenylalanine lands **98.25-98.59 %** against its sourced 97.5 %, inside D-259's 96-99.5 window.
* The control that makes that a finding rather than an artefact of supply limitation: leucine
  reads 19.4-23.2 %, so the margin narrowed by 1.2 points and is still 75 points wide.

## 5. P6 — the two dials, treated differently on purpose

* `_D260_MATCHED_F`, the shipped-form `f` **defined** as the one matching the growth-anchored
  split: **0.174 -> 0.2135**. Re-derived, not retuned. Its two asserts are not both evidence:
  the split-agreement one (reads 0.001 points) compares a hardcoded literal against the number
  that literal was set from, so it is a **staleness tripwire on the literal** — this beat is what
  made it fire — and not a statement about the mechanisms. The load-bearing assert is the other:
  the two arms' **tracers** agree to **1.0001x** at a matched split, which has teeth because they
  run different mechanisms.
* `_D260_LAMBDA`, the over-draw **priced at a number**: **left at 5.0**, and what it buys is
  re-reported. **This produced a correction to D-260's own prose.** On the stated composition
  λ=5 landed the split at **73.29 %** — *below* Crépin's 77 % floor — while D-260's docstring
  said it "lands the split in Crépin's band"; the guard that should have caught that was pinned
  at 70, four points below her floor, so prose and assert could disagree. On the sourced
  composition λ=5 lands **77.83 %**, genuinely inside 77-86, and the crossing now sits at
  λ≈4.8. The guard is re-pinned against her floor.

  Sweep on the sourced composition (kappa = 0.01, modifiers attached, 20001 points, closure
  1.0000 throughout): λ = 1.5 / 2.0 / 2.5 / 3.0 / 4.0 / 4.5 / 4.75 / 5.0 gives split =
  52.08 / 58.93 / 64.06 / 68.04 / 73.82 / 75.99 / 76.95 / 77.83 % and tracer =
  3.904 / 3.346 / 2.928 / 2.604 / 2.132 / 1.956 / 1.878 / 1.806 %.

## 6. P7 — THE HEADLINE: D-266's central result is retracted by its own composition

D-266 pre-registered that the joint repair would push the leucine tracer below Rollero's 3.4 %
floor, measured 3.617 %, and scored its own prediction **wrong**. On the sourced composition:

| | leucine split | leucine tracer | verdict |
|---|---|---|---|
| joint, stated composition (D-266) | 54.50 % | 3.617 % | inside Rollero's 3.4-8.2 |
| joint, sourced, **lo** | 57.77 % | **3.356 %** | below the floor |
| joint, sourced, **mid** | **60.61 %** | **3.131 %** | below the floor |
| joint, sourced, **hi** | 63.09 % | **2.934 %** | below the floor |

**The tracer is below Rollero's floor across the whole of what is left of the bracket.** The
pre-registered mechanism was right; what was wrong was the composition it had been priced at.
The denominator is untouched (isoamyl 2176.6 µM in every joint arm, +2.5 % on the shipped
2123.5, level-preserved at the Wang anchor by construction), so this is a different point on
D-260's line and not an escape from it — the joint arm now misses **both** leucine observables
it was offered as satisfying one of.

## 7. P8 — the costs, all confirmed, and one D-266 statement changed

* The three splits Crépin's own numbers are imposed on over-shoot further: ile **66.2 -> 73.5**,
  val **65.2 -> 73.0**, thr **81.2 -> 84.0**, against her 51 / 41 / 38.
* **D-266's "the growth sink alone under-shoots two of them" no longer holds.** On the sourced
  composition the sink alone reads ile 46.9 (under 51), **val 46.0 (over 41)**, thr 62.6 (over
  38). It under-shoots one, not two.
* The D-117 sparing-credit refund rises **0.236 -> 0.276x** growth's draw — the refund scales
  with the draw — against the shipped 0.584x and the hard ceiling of 1.0.
* Carbon and nitrogen ledgers close to 0.0 relative drift; net dS/dt never positive.

## 8. The one axis that improves, and why it is not a win

Rollero's 13C-leucine enrichment, joint arm, against Table S2:

| must | Table S2 | D-266 (stated) | D-268 (sourced) |
|---|---|---|---|
| SM70 | 3.4-3.5 % | 1.78 (under) | **1.54 (under)** |
| SM250 | 4.2-4.7 % | 5.44 (over) | **4.69 (inside)** |
| SM425 | 6.8-8.2 % | 8.33 (top edge) | **7.18 (inside)** |

Two of three move from outside the band to inside it. **Scored apart from the total it is
divided by, this is two errors cancelling more closely, not agreement**: the labelled AMOUNT is
86.7 µM against his 55.2-64.2 at SM250 (1.35x) and 147.3 against 65.0-70.3 at SM425 (2.10x),
while the total is ~2x his (2051.7 against 793-1034 µM at SM425).

Valine's tracer stays wrong in the same direction and further out (0.68 / 1.62 / 2.08 % against
Table S1's 2.1-2.3 / 3.4-4.0 / 5.3-5.4), below the shipped arm's own under-reading on all three.

## 9. One guard that passes on a thread, and one claim of this document that was wrong

* The growth arm's SM70 leucine enrichment is **3.5011 %** against Table S2's 3.5 % ceiling. The
  claim "the growth arm alone is over Rollero's ceiling on all three musts" is now carried at
  SM70 by **one thousandth of a point**. It is true and it is one re-derivation from false.
* **A first draft of this document claimed the joint arm's gate-closure fraction had moved to
  1.13 points from the blend's at SM425, and widened that guard's bound from 1.0 to 1.5 to
  admit it. Both were wrong and are reverted.** The 1.13 came from comparing the joint arm
  measured here against D-257's *recorded* blend number rather than against the blend arm
  measured in the same run. The real gap is **0.107 / 0.110 / 0.114** points across SM70 /
  SM250 / SM425 — flat, and unchanged by the composition. The guard is back at 1.0 and now
  carries a second assert pinning the gap under 0.25, four times tighter, so a sink change that
  starts steering production fails rather than passing inside a point. The bound was widened by
  a mis-measurement, and the suite caught it: the failing test is what produced this paragraph.

## 10. Vacuity

The pre-registration's vacuity check was "splits move by less than the closure control's own 2 %
tolerance". The smallest movement anywhere is threonine's high edge, 50.91 -> 51.26 % — 0.35
points on a weight that moved only 1.014x, which is the *narrowest* corner of the re-pricing.
The headline movements are 6-11 points. The beat is not vacuous.
