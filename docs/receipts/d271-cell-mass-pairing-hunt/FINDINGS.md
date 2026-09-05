# FINDINGS — the count-plus-weight pairing hunt (D-271, 2026-09-05)

**Result: EMPTY.** The negative now has a mechanism instead of a count of files searched, and
two guards' worth of consequence.

## 1. The two-sided negative, four named instances

Each source reports the ONE currency its own endpoint needs, and never both:

| source | counts? | weighs? | why it fails |
|---|---|---|---|
| Crepin 2017 Data Set S1 (on disk since D-246) | **no** | **yes** — filtered, washed, 105 C / 48 h to constant weight | a real gravimetric weight with nothing to divide it by |
| Varela 2004 (AEM 70:3392) | inoculum only (Neubauer, 1e6 cells/mL) | yes (85 C to constant weight) | count is the INOCULUM, weight is the crop — two timepoints |
| Foster 2022 (Front. Microbiol. 13:747546) | yes (haemocytometer) | only as a trehalose normaliser (mg/g DCW) | no DCW at a fermentation timepoint |
| Tyrell 2013 (BrewingScience 66:76-83) | yes (Fig. 4 panel) | no — pitched by count, MEBAK III 10.4.3 | no weight anywhere |

**Two were already on record** at D-219 sec 4 (Foster, Tyrell) and were re-derived here before
the archive was re-read. The new instance is **Crepin** — closest anyone gets.

## 2. The chemostat ring, checked

Where the pairing IS routine, because cells-per-gram is an instrument calibration there. Does
not deliver: the Delft school works in **C-mol and g DW and does not count cells**. Neither
van Gulik & Heijnen 1995 nor Lange & Heijnen 2001 (already read at D-267) carries a conversion.

## 3. What the literature offers instead — never a weighing divided by a count

* **Klis, de Koster & Brul 2014** (Eukaryot. Cell 13(1):2-9, PMC3910951): **16.5 pg** haploid
  (44 fL), **31.2 pg** diploid (83 fL), exponential phase, rich medium, 30 C. Computed —
  "multiplying the volume with the density (1.11) ... with the dry weight fraction (0.34)".
  The route reproduces its own printed masses to 0.7 % (16.61 / 31.32 pg).
* **BNID 101795: 60 pg**, Physical Biology of the Cell Table 1.1 — a stated rule of thumb.

## 4. The headline: EVERY figure in reach lands below branch 1's 70.9-91.9 pg

Klis haploid 16.5, Klis diploid 31.2, rule of thumb 60, D-219's settled 40, D-219's band high
edge 50, D-270 sec 7's re-priced elemental estimate 47.71 and 62.12, and the sourced
cross-check's high edge 56.61 — **eight figures, all below**. A direction, not a settlement:
none is a pairing, none is an ale cell in wort.

## 5. The size route is FRAME-BROKEN

Branch 1's demand through the sourced constants: **188-244 fL**, equivalent sphere
**7.11-7.75 um**. Whether that is ordinary depends on which size source is opened:

| | Klis (microscopy volume) | Okada 2023 (FSC apparent diameter) | ratio |
|---|---|---|---|
| haploid | 44 fL (4.38 um equiv) | 7.3 um => 203.7 fL | **4.63x** |
| diploid | 83 fL | 9.4 um => 434.9 fL | **5.24x** |

The demand is **2.26-2.93x a diploid** in Klis's frame and **0.92-1.20x a haploid** in Okada's.
A route that answers both ways cannot adjudicate. Okada's brewing diploid (K7A, 12.6 um) is not
an escape — it lives in the frame that makes everything large.

## 6. The one repair that ships

D-219's cross-check used rho 1.11 and an **unsourced** "~30 % dry matter". Its printed 30-57 pg
edges imply dry fractions **0.270 and 0.342**, so the "~30 %" was doing duty as a range. Klis
sources both constants; at 0.34 the cross-check is **37.74-56.61 pg** — the value D-219's own
upper edge was already computed at. A narrowing, not a move. The settled 28-50 pg band is
**untouched**: it comes from the elemental route, not from this check.

## 7. An error corrected against itself

This file's first draft read D-270 sec 7 as re-pricing **branch 1's demand** to 47.71-62.12 pg.
It re-priced the **engine-side elemental estimate**; the demand stayed at D-230's 70.9-91.9 and
what narrowed was the GAP (2.03-2.64x -> 1.14-1.93x). Now pinned by a guard.

## 8. Scope

The engine's 4e-11 g is D-219's DEFINITION of the gram Coleman's fit is counted in. Nothing
here moves it and nothing here was scored as if it could.
