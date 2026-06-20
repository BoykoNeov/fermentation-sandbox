# Milestone 1 — context (key files & decisions)

## Where things live (built in Milestone 0)

- **Tiers:** `src/fermentation/core/tiers.py` — `Tier`, `combine`.
- **State:** `src/fermentation/core/state.py` — `StateSchema`, `VarSpec`,
  `StateVector`. Sugar is a vector (D-4).
- **Process:** `src/fermentation/core/process.py` — subclass `Process`
  (`name`/`tier`/`touches` + `derivatives`); `ProcessSet` sums + derives tiers.
- **Runtime:** `src/fermentation/runtime/integrate.py` — `simulate`, `Trajectory`.
- **Parameters:** `src/fermentation/parameters/` — schema + `load_parameters`;
  data in `data/wine_generic.yaml` (currently **placeholders**, tagged speculative).
- **Units:** `src/fermentation/units/convert.py` — Brix/SG/Plato/ABV/°C/days.
- **Validation:** `src/fermentation/validation/` — conservation checks + benchmark
  specs (`BENCHMARKS`) + `compare_series`.
- **Benchmarks (to unskip):** `tests/benchmarks/test_milestone1.py`.
- **Toy reference model:** `tests/conftest.py` (`MassConservingFermentation`) —
  pattern for a conservation-respecting Process; replace with real kinetics.

## Canonical units (D-3) — get these right in every formula

concentration **g/L**, temperature **K**, time **hours**. μ_max etc. in `1/h`.
Benchmarks are quoted in days → convert at the assertion boundary.

## Parameters to source (each needs full provenance, handoff §2.3)

μ_max, K_s (sugar), K_n (nitrogen), ethanol-inhibition constant(s) + tolerance,
biomass yield, ethanol yield (have ~0.47), Arrhenius A + E_a per rate,
maintenance/decay rate, max population. Reconcile against:
- Coleman et al. 2007 — temperature-dependent N-limited wine model (10.1128/AEM.00845-07)
- Cramer et al. 2002; Malherbe et al. 2004 — wine kinetics
- Gee & Ramirez 1988; de Andrés-Toro et al. 1998 — dynamic beer models

When a placeholder is replaced with a sourced value, re-tag its tier (typically
`plausible`; `validated` once it reproduces a benchmark).

## Beer specifics to handle (handoff §2.1)

Sequential sugar uptake (the `S` vector), apparent-vs-real attenuation, and the
diacetyl rest (the byproduct itself is Tier-2, but the uptake order matters now).

## Open question to resolve at M1 start

- Growth form: logistic vs Monod-on-N vs both. Lean Monod-on-nitrogen since the
  defining mechanism is N-limitation; validate against the lag→exp→stationary
  shape.
