# Milestone 2 — Tier-2 (plausible mechanisms)

> Status: **scoping**. Milestone 1 (the Tier-1 validated core) is closed: all three
> §2.2 benchmarks green, D-8 → D-17. This plan opens Tier-2. The two opening calls
> are recorded in **DECISIONS D-18**: pH is a full **charge-balance solver**
> (derived-algebraic, couplings emerge — not a scripted tracked-pH), and the
> **byproducts/temperature beat is built first**, ahead of pH (a deliberate inversion
> of the handoff §6 order, justified below).
>
> Goal: add Tier-2 mechanisms **one Process at a time, each behind its own tests**,
> tagged `plausible`/`speculative`, while the Tier-1 suite stays green and the
> validated core stays togglable-isolable (prime directive #3).

## Build order (dependency-ordered; handoff §6 step 4, re-sequenced per D-18)

```
  byproducts + temperature axis   ← FIRST beat (this milestone's active work)
        │  (T + N only; no pH dependency)
  pH / acid charge-balance solver ← keystone infrastructure (D-18); unblocks ↓
        ├── SO₂ (free/bound/molecular; pKa-driven speciation)
        └── MLF — Oenococcus oeni (2nd organism; malic→lactic; pH-sensitive)
                  └── mixed cultures / Brett / sour consortium
  stochastic ensemble wrapper      ← physics-free runtime layer; parallel, any time
```

Why byproducts before pH (full rationale in D-18): it **closes the last skipped
benchmark**, **exercises the Arrhenius temperature machinery** that M1 built but never
tested (inert at the isothermal `T_ref`), and is the **most self-contained** Tier-2
physics (esters/fusels are T- and N-driven, pH-independent) — so it costs the pH chain
nothing and defers the heavy charge-balance commitment until its design is locked.

---

## Active beat: temperature-/metabolism-driven byproducts (handoff §3.2)

### Definition of done

1. **Unskip and pass `test_lower_temperature_is_slower_but_cleaner`** in
   `tests/benchmarks/test_milestone1.py` — the directional spec: lower fermentation
   temperature ⇒ longer time-to-dryness **and** fewer fusel/ester byproducts.
   - The *"slower"* half already holds today (a constant sub-`T_ref` run activates the
     dormant `ArrheniusTemperature` modifiers); this beat must keep it true and add a
     real multi-temperature comparison.
   - The *"cleaner"* half is the new physics: ester + fusel pools that fall with
     temperature.
2. Each new Process ships its own unit tests (monotone in T, produced-only/non-negative,
   **togglable-off ⇒ the validated core is byte-for-byte unchanged**).
3. Carbon conservation still holds to tolerance (see "carbon accounting" below for how
   trace byproducts are handled). The Tier-1 §2.2 trio stays green.
4. `pytest` / `ruff` / `mypy` green; new parameters carry real provenance and honest
   (mostly `speculative`, directional-only — handoff §3.5) tiers.

### Scope of this beat — esters and fusels only

The benchmark names exactly *"fusel/ester byproducts"*, so this beat models those two
and **stubs** the rest of §3.2 (diacetyl, acetaldehyde, H₂S — see stubs below). Keeps
the first Tier-2 Process tight and benchmark-aligned.

- **Higher / fusel alcohols** (Ehrlich pathway from amino acids) — rise with temperature;
  the N-relationship is **non-monotonic**, so model the pathway, *do not fit a slope*
  (handoff §3.2). First-pass form ties fusel production to the amino-acid/growth flux
  scaled by a temperature factor, with the non-monotonic-N behaviour flagged as a known
  simplification to revisit.
- **Esters** (isoamyl acetate ≈ banana, ethyl hexanoate, ethyl acetate) — favoured by
  warmth and strain, coupled to nitrogen dynamics. First-pass form ties ester synthesis
  to the fermentative flux scaled by a temperature factor.

These are **additive** Processes (they *produce* a compound), so they fit `ProcessSet`'s
additive model directly — no new mechanism is needed (contrast the multiplicative
`RateModifier` hook the M1 inhibition/Arrhenius work introduced).

### State & schema

New **produced-only pools** (0 at pitch, accumulate during ferment): `esters`, `fusels`.
Use `VarSpec.default = 0.0` — the D-16 pattern that let `Gly`/`Byp`/`X_dead` land without
touching the ~37 initial-condition call sites. Both wine and beer carry them (esters/
fusels are produced in both); beer-specific markers (diacetyl) come in a later beat.

### Carbon accounting (a build sub-decision to settle, flagged here)

Esters/fusels carry carbon but are **trace** (mg/L–low-hundreds-mg/L) beside the g/L
ethanol flux — the same order as the M1 biomass diversion (~1–2 %, D-8). Two options,
to decide at build time:
- **(a) Route their carbon from sugar** into the tracked pools and scale the
  ethanol/CO₂ split, exactly as D-16 does for `Gly`/`Byp` ⟹ `total_carbon` closes to
  machine precision for any yields. **Recommended** — keeps the project's carbon-first
  discipline; the bookkeeping is the D-16 algebra reused.
- (b) Treat as trace and carbon-unaccounted (magnitude immaterial), tier `speculative`.
  Simpler but leaves a (tiny) carbon leak the conservation tests would have to exempt.
Note fusels derive from amino acids (Ehrlich), so their carbon skeleton is arguably
N-pool-linked rather than sugar-linked — the build must pick one source cleanly.

### Parameters to source (provenance, like the D-12 sweep)

Ester-synthesis and fusel (Ehrlich) rate constants + their temperature sensitivities.
Reading list (reconcile, don't transcribe): de Andrés-Toro et al. 1998 (dynamic beer
model — includes byproduct terms), Malherbe et al. 2004 (wine), plus ester/higher-alcohol
kinetic literature. Tiers: `plausible` only where a source measures *our* form; otherwise
`speculative`. The benchmark is a *directional* check, so validation is qualitative
(handoff §3.5) — these do not earn promotion past `plausible`.

### Approach (test-driven, mirrors M1)

1. Add `esters`/`fusels` to the wine & beer schemas as defaulted produced-only pools.
2. Implement the fusel (Ehrlich) Process + unit test; then the ester Process + unit test.
   Each: additive, produced-only, monotone-increasing in T, togglable-off.
3. Settle the carbon-accounting sub-decision (recommend route-from-source, D-16 style);
   extend `total_carbon` if pools are carbon-routed.
4. Source + reconcile the rate parameters; replace placeholders; tag tiers honestly.
5. Add the multi-temperature comparison and unskip
   `test_lower_temperature_is_slower_but_cleaner`. Confirm the §2.2 trio + conservation
   stay green.
6. Record outcomes in a DECISIONS entry (D-19) and update this plan + ARCHITECTURE.

---

## Scoped but not yet designed (stubs, in dependency order)

- **pH / acid charge-balance solver** — the keystone (DECISIONS D-18). Designed deep at
  its own beat. Track tartaric/malic/lactic/acetic (± carbonic) as carbon-accounted
  state; solve `Σ charge = 0` for `[H⁺]` each RHS call; `pH` is a derived pure function.
  Must resolve the three D-18 couplings (evolved-vs-dissolved CO₂; acid carbon vs the
  D-16 `Byp`=succinic sink; pKa(T)).
- **SO₂** — free/bound/molecular equilibrium, **pH-dependent** (molecular fraction via
  pKa ≈ 1.81 is the antimicrobial one); binds acetaldehyde. Needs pH first.
- **Malolactic fermentation (*Oenococcus oeni*)** — a second-organism Process activated
  by a "pitch MLF" event: L-malic → L-lactic + CO₂, deacidifies (pH rises). Its growth
  is pH/ethanol/SO₂/T-sensitive. Needs pH first; first consumer of the multi-organism
  competition extension.
- **Mixed cultures / Brett / sour consortium** — resource competition (extended Monod /
  Lotka–Volterra). After MLF.
- **Stochastic ensemble wrapper** — physics-free runtime layer sampling each parameter
  within its provenance `Uncertainty` band and running ensembles → confidence bands on
  every output. Lives in `runtime` over `simulate`, outside the deterministic core
  (handoff §1.6). Buildable in parallel with the byproducts beat; API is an engineering
  choice (no scoping gate). See `milestone-2-context.md`.

## Out of scope for Milestone 2 (Tier-3)

Aging (the years axis), oak extraction, and the heuristic sensory/OAV layer — all Tier-3,
isolated and `speculative`, after Tier-2 settles (handoff §4, §6 step 5).

## Risks

- **Validation gets scarce here** (handoff §3.5): good Tier-2 datasets are often
  proprietary. Where curves are unavailable, validate **directionally** (warmer ⇒ more
  esters; MLF ⇒ malic falls / pH rises) and tag `plausible`, never `validated`.
- **Parameter sourcing** again dominates effort — byproduct kinetics are scattered and
  strain-specific (handoff §2.3 risk, recurring).
- **pH solver robustness** (later beat): the per-RHS `[H⁺]` root-find must stay smooth
  and fast for the stiff BDF loop; a kink there chatters the solver (the same smoothness
  concern as D-9's catabolite repression and D-10's C¹ inhibition touchdown).
