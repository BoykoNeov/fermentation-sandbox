# D-273 pre-registration — beer's three flux-linked acid courses, re-framed

Written BEFORE any model run. Owner picked this beat from the open ledger
(`test_the_three_flux_linked_acid_courses_are_mistimed`, strict xfail, D-215 §3).

Only arithmetic on the already-transcribed `TYRELL_ACID_COURSE_PPM` /
`TYRELL_FLUX_FRACTION` has been done at the time of writing (probe0, below).

## The hypothesis

D-215 §3 scored each acid's *measured* course against the *model's* course on the
CALENDAR, and found the three timing errors have OPPOSITE signs. That comparison
carries the model's own ~2.8x day-2 flux deficit (D-215 §4, still open and
DELIBERATELY parked at D-223) inside it. Scored in the source's own frame —
measured acid rise against measured flux — the opposition should disappear.

## Predictions

* **P1** In Tyrell's own frame all three acids TRAIL his flux at day 2 and day 5;
  no acid leads. (Day 1 excluded: malic's rise is 24.25 ppm at a +-2 ppm read
  tolerance, so a day-1 lag of ~24 points is inside noise.)
* **P2** The model's acid-rise fraction equals the model's OWN fermented fraction
  to within 2 points at every whole day 1-7. If it does not, the flux-frame guard
  this beat proposes has a premise that must be restated. **BLOCKING.**
* **P3** The lag ordering succinic < lactic < malic is stable across days 2 and 5
  (not a single-day artefact).
* **P4** The model's flux deficit at day 2 lies strictly INSIDE the range of the
  three measured lags, i.e. between succinic's and malic's. This is what makes the
  calendar-frame residual change sign across the three, and it is NOT a tautology.
* **P5** Tyrell's wort is fermented out by day 5 (1.003) yet all three acids gain
  a further share of their rise between day 5 and day 7; the shipped rate law
  returns `[]` once flux stops, so it cannot express ANY of it, at any speed.
* **P6** Lactic's post-attenuation share is robust to which endpoint is taken as
  "fermented out" (day 4 vs day 5); malic's is NOT and roughly halves.
* **P7** No single lag shape fits all three: the day-2 lags span ~5x, so a shared
  correction would need a per-acid magnitude. A mechanism is therefore NOT
  buildable from this dataset, and this beat ships measurement + guards only.

## Abandonment rule

If P2 fails, the flux-frame framing is withdrawn as written rather than patched,
and the beat reports why. If P1 fails for any acid at day 2 or day 5, the
"single direction" claim is dropped and D-215 §3 stands.

---

## Outcome (appended after the probes, per the rule above)

| | prediction | outcome |
|---|---|---|
| P1 | all three trail at days 2 and 5 | **CONFIRMED** — +13.5/+44.8/+63.5 and +11.1/+21.9/+38.4 |
| P2 | model acid tracks its own flux to 2 points (**blocking**) | **CONFIRMED** — worst residual 2.09 points, at day 2 |
| P3 | ordering succinic < lactic < malic stable | **CONFIRMED** at days 2 and 5 |
| P4 | the flux deficit falls INSIDE the lag span | **CONFIRMED** at days 1, 2 and 5 (+8.8 / +40.7 / +24.5) |
| P5 | acids rise after attenuation; the law cannot express it | **CONFIRMED** — 0.989 % modelled against 21.6 % measured |
| P6 | lactic robust to the endpoint, malic not | **CONFIRMED** — 21.6/22.6 % vs 38.1/16.5 % |
| P7 | no mechanism buildable from this dataset | **HELD** — measurement + guards only, nothing in `src/` |

Nothing was withdrawn: the abandonment rule did not fire. Two numbers came out larger than the
D-215-era figures they update, both because the scenario moved under them (D-222's counted pitch,
D-223's re-anchored uptake rate): the day-2 flux deficit reads +40.7 points here against the
~38.9 implied by D-215 §4's own table, and the model's day-2 acid fraction reads 17.42 % against
the 20.5 % D-215 quoted.
