# D-273 — beer's three flux-linked acid courses, re-framed

Receipts for the record. `probe.py` prints every number below and writes `findings.json`;
`PREREGISTER.md` holds the seven predictions, written before the model was run.

Nothing here is a new transcription. `TYRELL_ACID_COURSE_PPM` (Figs 9/10/14) and
`TYRELL_FLUX_FRACTION` (Fig. 4's extract panel) were both read at D-215, off the same four
ferments of the same 12 °P wort. What changes is which of them gets compared to which.

## 1. The frame D-215 §3 used, and what it contains

D-215 §3 scored each acid's **measured** course against the **modelled** one, by day, and read
the per-acid sign off the difference: succinic late by 25 points, malic early by 25, lactic
nearly right. It concluded the errors OPPOSE and that *"there is no single timing correction
that helps all three"*.

That difference is two defects added together. The modelled column carries the model's own
day-2 flux deficit — the engine ferments this wort ~3× too slowly at day 2 (D-215 §4), which is
still open and was **deliberately parked** at D-223, because the uptake rate that closes it
loses both of beer's speed anchors. So the comparison asks "is the acid mistimed?" while
holding a stopwatch that is itself wrong.

This is the same error D-215 §7 caught itself making one section earlier, in mirror image: there
the *measured* acid was paired with the *measured* flux while the argument needed the model's;
here the *measured* acid is paired with the *modelled* flux while the question needs the
source's. [[feedback-pair-the-arm-with-its-baseline]]

## 2. Frame B — the source against itself. No model quantity at all

How far each acid still has to go, in points of its own day-0 → day-7 rise, at the moment
Tyrell's wort has fermented the stated fraction:

| day | flux done | succinic | lactic | malic |
|---|---|---|---|---|
| 1 | 15.0 % | −3.9 | +5.2 | +24.3 |
| **2** | **59.4 %** | **+13.5** | **+44.8** | **+63.5** |
| **5** | **100.3 %** | **+11.1** | **+21.9** | **+38.4** |

**All three trail the sugar, in one direction, in a stable order: succinic < lactic < malic.**
There is no opposition anywhere in the source.

Read floors, from the ±2 ppm figure tolerance as a share of each acid's whole rise: lactic 2.8,
succinic 4.3, malic 8.2 points. Every day-2 and day-5 lag clears its own floor several times
over. **Day 1 does not** — malic's +24.3 is barely three read errors on a 24.25 ppm rise, and
succinic's −3.9 is inside its floor — so day 1 is excluded from the finding and from the guard.

## 3. Why the calendar frame changes sign — and why this is not arithmetic

On the calendar, each acid's residual is *its real lag minus the model's flux deficit*. The
deficit is one number for all three; the lags are three different numbers. Where the deficit
falls relative to them decides the sign pattern:

| day | model's flux deficit | span of the measured lags | inside? |
|---|---|---|---|
| 1 | +8.8 | [−3.9, +24.3] | yes |
| 2 | +40.7 | [+13.5, +63.5] | yes |
| 5 | +24.5 | [+11.1, +38.4] | yes |

It lands **inside** the span at every day scored. That is what turns succinic's residual
negative and leaves malic's positive off one shared defect.

**The decomposition itself is a tautology and is not offered as evidence** — any three numbers
satisfy `a − c = (a − b) + (b − c)`. What is not a tautology, and is the load-bearing
measurement, is the *inclusion*: had the deficit exceeded every lag, all three acids would read
early; had it fallen below every lag, all three would read late. Either outcome was available
and neither happened. `test_the_models_speed_deficit_lies_inside_the_span_of_the_measured_lags`
asserts the inclusion, so the explanation goes red if beer's speed ever leaves the span.

## 4. Frame C — the model against its own clock

The producer is `Y·Σr`, so each acid's progress is the fermentation's progress re-scaled:

| day | acid % of rise | ferm % of wort | residual |
|---|---|---|---|
| 1 | 5.74 | 6.42 | −0.68 |
| 2 | 17.42 | 19.51 | **−2.09** |
| 3 | 36.47 | 38.27 | −1.80 |
| 4 | 58.88 | 60.04 | −1.17 |
| 5 | 78.48 | 79.09 | −0.61 |
| 7 | 100.00 | 100.00 | 0.00 |

Both normalised on day 7. The residual is small, one-signed, and has a cause: carbon routed
into the acid, fusel and ester pools leaves `S` without passing through
`fermentative_uptake_rates`, so total sugar decline slightly exceeds the flux the producer is
paid on.

The three modelled acids differ from one another by **1.1e-13 points** — they are one curve.
This is why D-215 §3's modelled column was a single number (20.5 %) for all three, and it means
a per-acid modelled *shape* is not something this engine can currently produce at all.

## 5. The sharper finding: a shape the law cannot express at any speed

Tyrell's extract panel reads **1.003 at day 5** — the wort is fermented out. His acids are not:

| acid | after day 5 | after day 4 |
|---|---|---|
| **lactic** | **21.6 %** | **22.6 %** |
| malic | 38.1 % | 16.5 % |
| succinic | 10.8 % | 4.3 % |

`organic_acid_rates` returns `[]` the instant the flux stops, so the model can place only the
residual sugar's worth there: **0.989 %** of the rise past 99 % attenuation, identical across
the three.

This is **not** the parked speed defect. It is scored on the model's own attenuation clock, so
the engine may take as long as it likes to ferment the wort out. A 20-fold gap survives.

**Only lactic's number is robust.** Which day counts as "fermented out" is a figure read, and
malic's share more than halves and succinic's falls 2.5× between the two available choices —
both dip at day 5 by more than the ±2 ppm tolerance. Lactic moves 21.6 → 22.6 %. So lactic
carries the claim, the other two corroborate its direction, and **38.1 % must not be quoted as
the size of this effect**. [[feedback-a-summary-statistic-is-not-the-curve]]

## 6. Corroboration already in the repo

D-183 retired acetic's flux-linked rate law because Tyrell's Fig. 13 put **86 % of its rise
inside the first 15 % of the flux** — acetic runs *ahead* of the sugar, and got a growth-linked
producer of its own. These three run *behind* it. Flux-linkage is the wrong shape for the whole
set, in both directions, and both halves of that statement now come from the same paper.

## 7. What is NOT built, and why

A lagging producer was not built, for three reasons:

1. The day-2 lags span roughly 5× (13.5 / 44.8 / 63.5). One shape would need three fitted
   magnitudes — a compromise rather than a mechanism, which is what D-215 §3 forbade about its
   own frame and remains right about in this one.
2. Any time-domain repair would be scored on a model whose speed is deliberately parked, so
   there is no clock to score it against.
3. The corpus's one candidate arrives spoiled. `Y_lactic_sugar_beer`'s own note already says
   *"one dataset cannot separate a late excretion from an autolytic release"*, and D-215 §1
   found autolysis reported with **opposite signs** by two of the five beer texts, on a
   weeks-scale against this 7-day window.

The post-attenuation gap ships as a strict xfail instead: it states what is true of the source
and false of the model, and it turns green when a source not paid per gram of sugar exists.
