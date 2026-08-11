---
name: oxidation-cascade-and-dosing
description: "The Fe(II)+O2 cascade set and SO2 dosing (D-141 to D-148) - cascade is built and non-default, must-dosing is refused"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — oxidation cascade and SO2 dosing.** Detail split out of
`.claude/memory/project-fermentation-sandbox.md` at D-185; that file's ledger points here by
path. Read it when working on this subject. Every bullet is *what it forbids* + the D-record
to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue
past it from this file.** **Never evict an old prohibition to buy a line.**

**Cascade — BUILT and NON-DEFAULT (D-141). Do not re-build; do not flip it silently.**
- `core/kinetics/oxidative_cascade.py`, 8 Processes, `quinone` the one new slot, OFF the ledger (`h2o2` QSS is
  **`1/k` NOT `ln2/k`**). **REPLACEMENT** via `_OXIDATIVE_SETS` (`get_medium` /
  `compile_scenario(..., oxidative=)`). **`"direct"` is default and stays default.**
- **Do NOT close the fate gap by tuning** — fates move up to **25×** purely from re-homing the rate law. Both
  settling sources are **paywalled (nine hosts)**; abstract K values are **never a `source:` field**
  [[feedback-conceded-caveats-are-not-coverage]].
- **Never re-derive the activation floor as a fit** (name retired at D-172, above) — the O2 budget
  agreeing is **supply limitation (D-136), NOT a rate-law check**. Activation **reads the reductant
  pools** — never lump. The 31 D-140 guards stand — **never re-derive their pins**; edit only the two
  seams. Beer's O2 is a **floor `>=5 mg/L`**, never 5.71. Pin rtol **1e-4**
  [[feedback-pin-tolerance-vs-solver-tolerance]]; `quinone == 0.0` under direct is **exact**.
- **Benchmark EXISTS and is ACTIVE** (`tests/benchmarks/test_validation_danilewicz_so2_o2.py`, 15 tests, no
  xfail/skip — **not open work**). **Never pin 1.7** — one dataset's mode, above the other's whole range; assert
  the limits 1/2 + bands. The falsifier is the **traverse** (D-141's "structurally cannot produce" was wrong). Quinone branching **NOT settled**.
- **Operating point is load-bearing — enforce the >10 free-SO₂ floor, never state it in prose.** It is
  `SIM_CURVATURE_FLOOR_MGL`, **NOT Miao's criterion — never re-encode his**; keep the value + excluded-wines
  table. **Never report a single-dose verdict on this ratio** (direct unconditional, cascade **straddles**). The
  **"4–7× quinone shortfall" + its `xfail` are WITHDRAWN**.
- **"The sim under-binds SO₂" is WITHDRAWN (D-143) — never re-assert**; never compare an addition-method secant
  with an oxidation-path slope (concretely `oxidation_path_slope` vs `MIAO_BUFFERING_BAND`); **never re-run a
  pool sweep**. Acetaldehyde 0.0000 mM at D-142's operating point is a **scenario artefact**, not a schema gap. **A green suite proves nothing about a quoted decimal — every assertion is a band.**
- **The four D-142 artefacts are FIXED (D-144)** — Miao's **Table 4 intercepts ARE free SO₂ at O₂ exhaustion**:
  **read the intercept, never recompute from an assumed dose**. D-144's test **NAMES NO CAUSE — keep it that
  way**. Locus guard is **pure algebra** (a state route re-entangles `ph_of_state`); masses from `core.chemistry`
  (**M_SO2 = 64.058**); `_SO2_BINDERS` reads all **four**.

**Dosing schedule — REFUSED (D-145). Never re-propose must-dosing; benchmark unchanged.**
- **D-143's "five for five" is three for five.** **Import statistics from the shipped benchmark module** —
  re-derived helpers carried both errors. **The reservoir CANNOT move the Table 3 secant** — saturated both
  ends ⇒ moves the **INTERCEPT**, not the slope; `free0*` is **floor-CONTINGENT, never state it bare**.
  **Depletion gap RETIRED, not reframed** — D-143 omitted `M_SO2/M_O2`; only the secant is invariant, in band.
- **Four traps that each produce a plausible green:** `ProcessSet.disable()` is silently undone by `begin_aging`;
  `param_values` is a **property returning a fresh dict**; mol/L ×1000 is mmol/L not mg/L; **`trajectory.y` is
  `(n_states, n_times)` — `y[-1]` is the LAST SLOT's series**, not the final state (D-147).
- **§2.4 CLOSED (D-148) — never re-open as "the Brett/quench over-draw": no quench draw exists.** Live pair was
  Brett/POF; **summing is NOT the mechanism**. **Never build the shared depletion gate** — `depletion_gate`'s 8
  sites buy per-draw first-orderness `_decarboxylation_branch` already has. Out-of-band numbers are **BDF step
  artefacts**; an undershot pool **freezes negative**. Aging Process **names** are enumerated by `test_fusel_*`.
