# D-266 findings — the joint fusel repair, four arms, every sourced observable

Numbers are from `probe_joint.py` → `findings.json` (55 s, Python 3.13, this repo at the
commit that adds this folder). Arms as in `PREREGISTER.md`. Every quadrature closure in every
arm reads 1.0000 (the D-259 control, asserted in the fixture in the shipped tests).

## Predictions, scored

| | prediction | verdict |
|---|---|---|
| P1 | harness reproduces D-257 (36.3/48.1/55.9 %, 172.3 mg/L, 72.9 % Phe) and D-260 (27.58 % / 5.900 %) | **CONFIRMED** on every number: blend 36.5 / 48.4 / 56.9, 172.37, 72.86; growth 27.58 / 5.900; shipped 81.5 / 1.507 / 2123.3 µM |
| P2 | joint clears D-257's blocker: Phe left ≤ 10 % (D-104 fixture) and ≤ 5 % (Crépin) | **CONFIRMED** — 0.00006 % and 0.000004 %, from 72.9 % under the blend alone |
| P3 | joint leucine split ≥ 25 points above growth's 27.6 %, landing 55–80 % | **PARTIAL** — +26.9 points, but 54.50 % (lo/mid/hi 46.0 / 54.5 / 61.5), half a point under the window |
| P4 | joint leucine tracer on Crépin's must falls BELOW Rollero's 3.4 % floor | **REFUTED** — 3.617 %, inside the 3.4–8.2 % band (lo/hi edges 4.29 / 3.06 %) |
| P5 | joint isoamyl on Crépin's must 1400–2000 µM (from 2123) | **REFUTED** — 2176.6 µM, +2.5 %: the denominator did not move |
| P6 | joint NT fractions within 3 points of blend's | **CONFIRMED** — 36.56 / 48.48 / 57.02 vs 36.49 / 48.39 / 56.92 |
| P7 | joint nitrogen response ≤ 1.4× | **CONFIRMED** — 1.2996× (Rollero 0.74–0.77×; shipped 2.27×) |
| P8 | joint leucine enrichment on Rollero's musts below the Table S2 floor on all three | **REFUTED on two of three** — 1.78 / 5.44 / 8.33 % against 3.4–3.5 / 4.2–4.7 / 6.8–8.2 %: under, OVER, at the top edge |
| P9 | level anchor holds under joint | **CONFIRMED** — 172.37 / 172.47 mg/L at YAN 250 / 300 |
| P10 | C and N close, net dS/dt ≤ 0, joint C refund < 1 | **CONFIRMED** — 3e-13 / 2e-11 relative, 0.0, 0.236× |

P4 and P5 are the beat. The expected verdict — "the joint arm pays for the split on the tracer" —
was wrong because the split does not reach Crépin: it stops at 54.5 %, and on D-260's line that
leaves the tracer inside Rollero's band. P5 says why the line is the same line.

## The scorecard

Crépin's own must (D-259/D-260 fixture, other precursor consumers off), mid composition:

| | shipped | blend | growth | joint | measured |
|---|---|---|---|---|---|
| leucine → lump, % | 81.5 (imposed) | 81.5 (imposed) | 27.6 | **54.5** (46.0–61.5) | Crépin 77–86 |
| isoleucine | 51.0 (imposed) | 51.0 | 38.4 | **66.2** | 51 |
| valine | 62.0 (imposed) | 62.0 | 37.2 | **65.2** | 41 |
| threonine | 82.0 (imposed) | 82.0 | 57.9 | **81.2** | 38 |
| phenylalanine | 97.5 (imposed) | 97.5 | 98.9 | **99.6** | 97.5 (sourced) |
| leucine tracer, % of isoamyl | 1.507 | 1.472 | 5.900 | **3.617** | Rollero 3.4–8.2 |
| isoamyl, µM | 2123 | 2174 | 2123 | **2177** | (joint ceiling 1170) |
| phenylalanine left, % | 1.35 | 56.2 | 0.000 | **0.000** | ~0 |
| propanol de novo | 0.878 | 0.880 | **0.714** | 0.875 | floor 0.80 |
| isoamyl de novo | 0.966 | 0.967 | 0.911 | 0.948 | — |

Rollero's own musts (D-255 fixture, full Process set):

| | shipped | blend | growth | joint | measured |
|---|---|---|---|---|---|
| isoamyl at NT / EF, % — SM70 | 100.6 | 36.5 | 100.6 | **36.6** | 51.1–53.9 |
| SM250 | 100.6 | 48.4 | 100.6 | **48.5** | 44.5–51.9 |
| SM425 | 100.6 | 56.9 | 100.6 | **57.0** | 42.3–51.3 |
| isoamyl EF, µM — SM70 | 1205 | 1577 | 1205 | **1579** | 1066–1337 |
| SM250 | 2080 | 1845 | 2081 | **1848** | 1314–1365 |
| SM425 | 2736 | 2047 | 2740 | **2052** | 793–1034 |
| response SM425/SM70 | 2.27× | 1.30× | 2.27× | **1.30×** | 0.74–0.77× |
| leucine enrichment, % — SM70 | 0.99 | 0.76 | 3.79 | **1.78** | 3.4–3.5 |
| SM250 | 2.05 | 2.31 | 7.84 | **5.44** | 4.2–4.7 |
| SM425 | 2.65 | 3.54 | 10.13 | **8.33** | 6.8–8.2 |
| leucine-labelled isoamyl, µM — SM250 | 42.6 | 42.6 | 163 | **101** | 55.2–64.2 |
| SM425 | 72.4 | 72.4 | 277 | **171** | 65.0–70.3 |
| valine enrichment, % — SM70 | 1.23 | 0.94 | 2.04 | **0.87** | 2.1–2.3 |
| SM250 | 2.44 | 2.75 | 3.64 | **2.15** | 3.4–4.0 |
| SM425 | 3.09 | 4.13 | 4.28 | **2.85** | 5.3–5.4 |

Anchors and ledgers: undosed isoamyl 170.5 (shipped/growth) and 172.4 mg/L (blend/joint) at
YAN 250, band 140–205. D-104 fixture Phe left: 20.3 / 72.9 / 0.000 / 0.000 %. Worst joint carbon
refund 0.584 / 0.263 / 0.236 / 0.236 × growth's draw. Net dS/dt ≤ 0 everywhere in every arm.

## What the numbers say

1. **The two refusals were not each other's defect, and the pre-registered mechanism for that
   was wrong.** The joint arm clears D-257's blocker (P2) and keeps its temporal gains (P6, P7)
   — but it does not pay on the Crépin-must tracer (P4 refuted) because the split stops at
   54.5 % and D-260's line then leaves 3.6 % for the tracer. It moves the model's point along
   the line; it does not move the line (P5): isoamyl on this must is 2123 → 2177 µM, the
   leucine pool still empties, so the joint-satisfaction ceiling (≤ 1170 µM at Crépin's 77 %)
   is exactly where D-260 left it. The blend is level-preserved at the Wang anchor by
   construction, and Crépin's must at 28 °C lands at the same total.
2. **What the growth sink alone breaks, the joint arm restores: the sourced propanol floor.**
   Growth-anchoring lets threonine's Ehrlich draw run un-capped and propanol reads 71.4 % de
   novo against the 80 % floor D-244/D-248 hold; the blend's cut of the Ehrlich draw brings it
   back to 87.5 %. Neither single repair could have shown this.
3. **What the joint arm breaks that the shipped model gets exactly: the three other splits.**
   The shipped sink imposes Crépin's isoleucine/valine/threonine shares (51/41/38) by
   construction; the growth sink under-shoots two of them and the joint arm over-shoots all
   three (66/65/81). "Closer to reality" is per-axis here, not overall.
4. **On Rollero's own musts the joint arm's leucine tracer brackets the measurement in
   enrichment and over-shoots it in amount.** 1.78 / 5.44 / 8.33 % against 3.4 / 4.2–4.7 /
   6.8–8.2 % reads under / over / at-the-top-edge — a nitrogen response of the enrichment
   (4.7× across his range) steeper than his (2.2×). In labelled micromoles it is 101 and 171
   against 55–64 and 65–70: the SM425 "in band" reading is an over-attributed numerator over an
   inflated denominator (2052 vs 793–1034 µM). Score numerator and denominator apart.
5. **Valine's tracer moves the wrong way under the joint arm** (0.87 / 2.15 / 2.85 % against
   2.1–2.3 / 3.4–4.0 / 5.3–5.4): the blend makes most isoamyl de novo, so the valine-labelled
   share falls. The growth sink alone is the best of the four on valine (2.04 / 3.64 / 4.28).

## Mutation arms (run after the guards were written; see the D-266 record)

Recorded in the record's §7 with the test names each arm reddens.
