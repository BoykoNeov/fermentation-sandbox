---
name: compile-sampled-census
description: "D-234 - the compile-read AND sampled census is 32 names, enumerated and classified; two rows stay LIVE and priced, the rest are closed"
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
unconvincing, go read D-234 — do not argue past it from this file.**

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
- **`set_ph` under an ensemble is LIVE, and D-186 forbade the state we are in.** Its event closes
  over the compile-time resolved map, so members land up to **0.07896** pH from the target and
  span **0.13202**, against **2.03e-11** at t=0 since D-233 repaired that anchor. D-186's own
  docstring says *"do not fix this one anchor alone"* — **D-233 did**. The repair is a
  `StateMutation` signature widening (13 `def mutate` + 13 sites in `src`, 1 in `tests`, 3 type
  refs; only `set_ph` changes behaviour) and it is **priced, not shipped — owner's call**.
  The guard PINS THE DEFECT: a RED means it was repaired, **delete it, never revert**.
- **`copper_typical` is LIVE and its band share is 100 % ARTEFACT.** It seeds the `copper` slot at
  compile and mean-centres `PhenolicBrowning`'s `f(Cu)`, so D-134's `f == 1` invariant holds only
  at the nominal draw. Control arm (recompile, both roles) moves aged `A420` by **exactly
  0.000000 %** — bit-identical; sampler arm **−14.04 %** at the band top, and a wine drawn
  *higher* in typical copper browns *less*. **Never quote the sampler figure as copper
  uncertainty.** The fix must be CONDITIONAL on the scenario not naming `copper_gpl` (then the
  two are genuinely independent), which is why it needs the `y0` generalisation D-233 declined.
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
- **NOT covered, deliberately:** the cascade / `direct_burst` wirings are an unmeasured third
  battery member (their sinks declare `copper_typical` too — expected to WIDEN §5, not flip it,
  and that is an expectation); `pKa_peptide_buffer`'s band re-derivation is still untouched.
