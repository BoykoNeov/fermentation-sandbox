---
name: compile-sampled-census
description: "D-234/235/236 - the compile-read AND sampled census is 32 names, enumerated and classified, and after both repairs NO row is LIVE"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1a6ea7a3-7bda-4049-b556-1015cdbe0204
  modified: 2026-08-26T12:24:03.623Z
---

**Live prohibitions — the compile-read-AND-sampled census (D-234).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about "parameters that are pinned at compile", `set_ph` under an
ensemble, `copper_typical`, the `must_aa_fraction_*` split, or a scenario override knob. Every
bullet is *what it forbids* + the record to read for *why*. **If a prohibition looks
unconvincing, go read D-234 — do not argue past it from this file.** **Both of D-234's LIVE rows
were REPAIRED on 2026-08-26 (D-235 pH, D-236 copper); the census now has NO live row.**

- **The census is RUN and CLASSIFIED — never re-propose it as unenumerated.** D-233 §10's parked
  item is CLOSED. **32** names over the battery (wine 29, beer 20, wine-with-overrides 30), each
  carrying a verdict in `tests/test_compile_sampled_census.py`'s `CENSUS`. The predicate is
  compile-time READ ∩ the sampler's own `_resolve_sample_names` — **not** D-153/156/157/159's
  drawability surface, which asks the opposite question and is **disjoint by construction**
  (pinned). **Do not re-run either audit against the other's headline.**
- **NEVER size this set with a grep.** `parameters["…"]` finds **26** and misses **21** — the 19
  `pKa_*` and both `nitrogen_uptake_charge_*` are read off `parameters.resolve()`. I
  pre-registered **12-20** from that grep and the answer was **32**
  [[feedback-size-a-set-with-the-instrument-that-can-see-it]]. The instrument is a class-level
  recorder with its own positive control, because a blinded one reports an EMPTY census and every
  membership test then passes vacuously.
- **The pH pair is REPAIRED (D-235) — never re-propose it, and never repair ONE anchor.**
  `StateMutation` is now `(schema, state, params)` and mutations get the **RUNNING** map, so
  `set_ph` re-anchors per member: worst miss **0.07896 → 1.9836e-11**, spread **0.13202 →
  2.1265e-11** (24 members), the t=0 anchor's own class. **Score an anchor AT the breakpoint**,
  not 0.05 h past it — that offset is ferment drift (7.04e-07), not anchor error. The same
  widening exposed `add_dap` reading z̄ off the compile map (**0.0108** pH), also repaired, so
  D-234's "only `set_ph` would change behaviour" is WRONG. Exactly **2** of 13 verbs read the
  argument, pinned by a recording map. D-234's blast radius was priced at definitions and missed
  **28** `event.mutate(...)` CALL sites
  [[feedback-a-blast-radius-counted-at-definitions-misses-call-sites]].
- **`copper_typical` is REPAIRED (D-236) — and its old band share was 100 % ARTEFACT.**
  `reanchor_for_member` is now `CompiledScenario.y0_for_member`, a list of per-slot RULES; rule 2
  re-seeds `copper` from the member's draw, so D-134's `f(Cu) == 1` holds for every member and
  aged `A420` went **16.65 % spread → bit-identical** (12 members, `only=["copper_typical"]`).
  **Never quote the old sampler figure as copper uncertainty** — the coherent channel measured
  **exactly 0.000000 %**. The rule is CONDITIONAL on the scenario not naming `copper_gpl`: there
  the two are genuinely independent (**14.34 %**, honest, must survive) and re-seeding would
  breach D-24. Both halves of the branch are guarded (Arms E and F).
- **A defect pin that stays GREEN through its repair is about the PATH it drives.** The copper pin
  did exactly that — its arms drive `simulate_scheduled` by hand and never reach `y0_for_member`.
  It was **re-scoped, not deleted**, and the repair got its own ensemble-frame guard. Do not read
  such a green as "the repair was inert"
  [[feedback-a-defect-pin-can-outlive-its-defect-by-driving-another-path]].
- **The eight `must_aa_fraction_*` are CLOSED at D-206 — do not re-measure them.** All eight were
  enumerated from the registry there, not argued from methionine. Verified current here, to the
  digit: every gate **0.888889** at nominal, methionine **0.952381→0.833333** (widest),
  phenylalanine **0.892430→0.885375** (narrowest). The honest channel was **0.06 %**.
  (The nominal gate moves with the DOSE — 0.909091 at `amino_acids_gpl` 1.0. That is not drift.)
- **`biomass_carrying_capacity` / `k_autolysis` under a scenario override are BY DESIGN — never
  file them as a defect.** The override is the **MODE** of `triangular(low, value, high)`, not a
  value the sampler discards (measured: drawn mean **3.8296** at an override of 4.0 against a 2.5
  reference). That is what D-164's in-band bound exists to keep well-formed. **D-24's surviving
  exclusion — scenario inputs are never sampled — is NOT breached**, and D-233 §9's warning
  against reading its correction as licence still stands.
- **NOT covered, deliberately:** the cascade / `direct_burst` wirings are still an unmeasured
  third battery member — but D-236's repair reaches them mechanically (it moves the SEED, so
  `copper − copper_typical` is zero for every member under any wiring centred on the same
  reference); that is an argument, not a measurement. `pKa_peptide_buffer`'s band re-derivation is
  still untouched, and `peptide_buffer_capacity_beer` stays back-solved offline at the nominal.
