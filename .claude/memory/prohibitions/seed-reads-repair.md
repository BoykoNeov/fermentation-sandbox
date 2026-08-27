---
name: seed-reads-repair
description: "D-241 — the D-45 fallback seeds are DRAWN via CompiledScenario.seed_reads, derived from the y0 rules; six repaired, two subsumed, tier half deliberately untouched"
metadata: 
  node_type: memory
  type: project
  originSessionId: ea3f3770-ac9d-4299-8476-1248aca8b862
  modified: 2026-08-27T08:47:25.403Z
---

**Live prohibitions — the seed-drawability repair (D-241).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about making a compile-time seed drawable, `seed_reads`,
`y0_for_member`, "the ensemble does not vary the initial conditions", or the two Coleman biomass
coefficients. Every bullet is *what it forbids* + the record to read for *why*. **If a prohibition
looks unconvincing, go read D-241 — do not argue past it from here.** Sibling:
`banded-undrawn-census.md` (D-240, the census this repaired six rows of).

- **The repair is BUILT — never re-propose "make the seeds drawable" as open work.**
  `CompiledScenario.seed_reads` is a `tuple[str, ...]` threaded into `simulate_ensemble` by
  `run_ensemble` and unioned in `_resolve_sample_names`. Six names drawn:
  `dms_potential_initial`, `bound_h2s_initial`, `bound_methanethiol_initial`,
  `must_fermentable_fraction`, `o2_wort_aeration_beer`, `burst_antioxidant_initial`.
- **`seed_reads` is DERIVED from `_member_seed_rules()` and must stay derived.** A name is
  drawable **iff** a rule re-seeds it, so D-240 §10's "drawing a name that cannot reach `y0`,
  which is worse than the gap" is *unrepresentable*, not merely guarded. **Never introduce a
  second, hand-written list** — that reintroduces the exact failure the design excludes.
- **The union goes into the DEFAULT branch only, and BEFORE `exclude`.** An explicit `only=`
  still means exactly those names, or every one-parameter sweep in the repo silently becomes a
  joint one [[feedback-a-one-parameter-sweep-is-not-the-band]]; unioning after `exclude` would
  make a seed the one sampled name with no way to pin it. Both are guarded — do not "simplify".
- **It carries NO TIER CLAIM, on purpose, and that gap is still OPEN.** `reads` has two masters
  (D-160): sampler scope and tier propagation. This channel takes the first alone. What a seed's
  provenance tier should do to a *state slot* is unmeasured, and closing it needs its own beat.
  **Never cite D-241 as having decided it**, and never quietly add tier propagation here.
- **Each rule fires only while the compiled slot still holds the parameter's value**, and that one
  equality does three jobs: refuses a scenario-stated level (D-24), stops silently if the seam
  changes, and reproduces **D-147's burst-wiring gate for free**. Do not replace it with an
  explicit wiring check — the table deliberately does not know which oxidative set compiled.
- **`burst_antioxidant_initial` is drawn under `direct_burst` ONLY, and that is CORRECT.** Under
  `direct`/`cascade` D-147 zeroes the slot, so across its whole 50× band the **entire `y0` is
  bit-identical** — nothing to repair there. It keeps a WIRING-CONDITIONAL/INERT verdict in the
  sibling census. **Never file the single-wiring scope as a half-repair.** Consequence: D-234's
  census is now wiring-dependent for exactly one name, which
  `test_the_census_itself_is_the_same_under_every_oxidative_wiring` predicted in its own docstring.
- **Beer is the NULL CONTROL and must stay one.** Its reported band is unchanged to **1.000** on
  `X`/`E`/`o2`/`acetic` while the rule fires (~1e-8 relative). The guard asserts BOTH halves — a
  zero alone is also what a rule that never ran produces
  [[feedback-a-control-needs-mechanical-reach]]. Never drop the reach assertion.
- **Quote the PAIRED numbers, and both statistics.** Against a shipped-before arm carrying rules
  1–3 with rule 4 off (**not** `y0_for_member=None`, which would credit D-233/236/238 here): band
  ratios `dms` **2.83×**, `methanethiol` **2.07×**, `bound_h2s` 1.36×, `E` 1.04×,
  `burst_antioxidant` **6.97×**. Per-member: median `dms` shift **49.6 %**, worst `E` **3.01 g/L**.
  **The band ratio alone misdescribes it** — ethanol's band barely moves while members move by
  grams [[feedback-a-summary-statistic-is-not-the-curve]].
- **The per-member conservation guard reads the MEMBER's `y0`, not the compiled one.** It read the
  compiled array — wrong since D-233, green until a drawn `must_fermentable_fraction` moved the
  carbon denominator, then 99.043 vs 97.834 g C/L. **Never "fix" that shape by widening a
  tolerance** [[feedback-a-half-pinned-read-is-green-until-the-quantity-moves]].
