---
name: keto-acid-excretion-shape
description: "The excreted keto-acid pools' rate SHAPE (D-49 option B), measured and refused for alpha-ketobutyrate at D-189"
metadata: 
  node_type: memory
  type: project
  originSessionId: 53cf9b2e-d59e-43b1-8918-019200d03f50
  modified: 2026-08-12T10:03:55.517Z
---

**Live prohibitions — the excreted keto-acid pools' rate shape.** Split out at D-185's pattern;
the status ledger points here by path. Read it when working on this subject. Every bullet is
*what it forbids* + the D-record to read for *why*. If a prohibition looks unconvincing, **go
read its D-record — do not argue past it from this file.** **Never evict an old prohibition to
buy a line.**

**Option B is MEASURED and REFUSED for alpha-ketobutyrate (D-189). Do not re-propose it there
as untested — `keto_acids.py` said "untested, and is the next step here" for 82 records, and
that sentence is now retired.**
- **A growth-linked SOURCE with the shipped flux-linked SINK drains the pool it feeds.** Growth
  ends day 0.98, dryness is day 2.04, the sink clears `exp(-0.1 × 106.4)` ≈ 2.4e-5 across that
  window: residual **2.000 → 0.000 mg/L**, sotolon **1.599 → 0.027 µg/L**. The arithmetic is
  **structural and general** — all three sinks share `k = 0.1` and that flux integral — so
  **option B cannot be adopted for the source ALONE in pyruvate or α-KG either.**
- **Both rescues were BUILT as probes and both REFUSED.** Growth-linking the *sink* too is
  **output-identical** once re-anchored (share 11.70 %, residual 2.000, sotolon 1.597 vs 1.599)
  — refused because it asserts **yeast stop re-assimilating when they stop DIVIDING**, which
  nothing sources and which contradicts the co-metabolic reasoning in `keto_acids.yaml`.
  Re-anchoring the sink instead needs **0.0196, 1.5× below its own declared band's low edge**.
- **Never quote the shipped residual's temperature-invariance as EXACT.** It is 5.65e-7
  relative across 15–28 °C, which is the BDF solver's floor at `rtol=1e-6`. A pin at `rel=1e-6`
  would have **0.57× headroom** — it would pass on the tolerance, not the physics. Both rescues
  cost that invariance outright (**2.03×** and **16.33×**).
- **NO GUARD WAS ADDED AND NONE IS OWED.** The shipped suite already forbids option B at **five
  asserts across three files** (`test_aging_scenario.py:1004/1039/1122/1175`,
  `test_closure_ingress.py:492`) — sotolon's perceptibility, the premox threshold, the D-100
  tripwire, and the closure ordering, which collapses to **`leaky == tight`** because sotolon
  loses its substrate entirely, un-lifting the D-108 limitation D-136/D-187 spent two beats
  lifting. Do not "add coverage" here.
- **D-107's diagnosis is CORRECTED, not merely extended.** "The diagnosis is the excretion
  SHAPE" is incomplete: shape is worth **+8.4 points** (2.62 → 11.03 %), **D-104's anabolic
  threonine sink is worth +12.3** (2.62 → 14.92 %), neither reaches Crépin's **19 %** and they
  **do not sum to it** — a third factor is unidentified. The gate is a **competition for one
  molecule**. The sink is **NOT a candidate fix** (its 77–86 % protein share is sourced, and
  Crépin's yeast ran the same competition). **The gate never becomes a fitted fraction**: the
  shape alone cannot reach 19 %, so no version was ever one constant from a match.
- **`EthanolToleranceDeath` moves this share by EXACTLY NOTHING** (2.620 % with and without) —
  a killed hypothesis, do not re-run it.
- **PYRUVATE's option B is CLOSED — REFUSED at D-195** (the "stays open" sentence is retired).
  **Never argue it as invisible-gain: option B WORKS**, peaking **122.5** vs the shipped monotone
  **30.0 mg/L**. Refused on COST — that residual is D-51's SO₂ binder, so draining it takes free
  SO₂ **33.1 → 42.8 (+29 %)**, an error in the FLATTERING direction. **37 asserts across 7 files
  forbid it; no guard owed.** Decisive: **Miao's secant 1.1005 vs measured 1.2526-1.9882** — and
  a magnitude-only control pins that on **the RESIDUAL, not the shape** (1.100510 vs 1.100501),
  so **no rescue variant with kinder timing exists** and the residual is **load-bearing for
  agreement with real wine**. **α-KG INFERRED, never measured** — same sink `k`, same D-51.
- **The share has ONE measurement route** — re-evaluating the rate law along the *unperturbed*
  trajectory. Disabling the other threonine consumers to isolate `d(threonine)` holds open the
  very gate the number is about and reports **87.8 %**. Never measure it that way.
- **Nothing pins 2.62 % and nothing should.** D-107's **1.7 %** was NOT reproduced and the
  difference is **not attributable** (the defining files are unchanged since D-122; D-107's
  record never states its probe's must). **A literal whose provenance you cannot establish is
  not pinnable.** The share is scenario-bound anyway: **3.63 % / 3.15 % / 2.62 %** at
  15/20/28 °C.
- **A growth-linked driver would need naming in BOTH** `ArrheniusTemperature.for_growth`'s
  extra targets **and** `BiomassCarryingCapacity.modifies` (`biomass_growth_rate` returns the
  **base** rate — the D-32 coupling in D-183's form), which would take the pool off the
  temperature-flat stance its two siblings keep. Recorded so a future attempt does not
  rediscover it as a bug.

Receipts: `M:\claud_projects\temp\ferment\d189-growth-linked-excretion\` — `PREREGISTER.md`,
`RESULTS.md`, `probe1..6.py`, the `d189_optionb` pytest plugin, and both suite arms.
Related: [[feedback-a-shape-change-is-a-change-to-the-pair]],
[[feedback-mutate-the-premise-before-building-the-guard]],
[[feedback-pair-the-arm-with-its-baseline]].
