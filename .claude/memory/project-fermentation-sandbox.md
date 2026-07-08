---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation
engine in Python (uv, scipy/numpy/pydantic). Repo:
https://github.com/BoykoNeov/fermentation-sandbox (default branch `main`).

**This file is session-boot context, not a changelog.** The full per-decision
narrative lives in the repo and is the single source of truth — do NOT copy it
here (see [[feedback-batch-end-ritual]]):

- `docs/DECISIONS.md` — every engineering decision D-1 … D-50, with the fork,
  the choice, and the reasoning. **The canonical archive.**
- `docs/ARCHITECTURE.md` — layering, package map, the core/runtime/scenario seams.
- `docs/plans/milestone-*.md` — task checklists per milestone.
- `CLAUDE.md` — prime directives (tiers, provenance, one-directional deps).

**Status (2026-07-07):**
- **Milestone 0 & 1 COMPLETE.** M1 = single-strain isothermal nitrogen-limited
  primary fermentation passing both §2.2 benchmarks (wine ~24°Brix→dry 8–14d,
  beer ~1.048→~1.010 in ~5.5d). Tier sweep (D-17): VALIDATED reserved for
  independent measured data; the benchmark pass earns PLAUSIBLE.
- **Milestone 2 — ALL physics beats COMPLETE**, plus post-M2 refinement work
  through **D-55** (POF v2, closes the last D-40 pt4 deferrals): **661 passed**,
  ruff + mypy strict clean. Brett
  arc (D-40) DONE: pt1 phenol
  pathway + pt2 BrettGrowth + pt3 BrettDeath (SO₂ kill) + **pt4 POF+ yeast
  decarboxylase + emergent reservoir**. Built so far: pH charge-balance solver
  (D-18), aroma byproducts esters/fusels (D-19–21), SO₂ speciation + free/bound
  split (D-22/28), MLF conversion + O.oeni diacetyl (D-23/31), stochastic
  ensemble + spread attribution + LHS/Sobol (D-24/25/37), diacetyl/acetaldehyde/
  H₂S (D-26/27/29), residual-N carrying cap (D-30), amino-acid ledger + fusel
  re-route + autolysis (D-32/33/34), event loop + temperature ramp + discrete
  intervention verbs `add_dap`/`add_so2`/`rack`/`pitch_mlf` (D-35/36),
  MLF-growth: dynamic `X_mlf` bacterial biomass, autocatalytic conversion (D-38),
  **MLF death + rack-removes-`X_mlf`: the MLF arc (D-23→D-31→D-38→D-39) is
  COMPLETE (D-39)**.

**D-39 (MLF death, the crux worth remembering):** `MalolacticDeath` moves viable
`X_mlf` → new `X_mlf_dead` pool, driven by **molecular SO₂ ALONE**
(`1 − g_SO₂`), Arrhenius temperature (warm kills, cold preserves — NOT the
cardinal γ(T)). A probe **killed** the first `1 − toxicity` draft: the ethanol
Luong wall already drives `1 − toxicity` ≈ 0.92 at post-AF ethanol, so ethanol
alone wiped bacteria in ~1 week (O.oeni actually persists weeks–months; SO₂/rack
is the deliberate kill). SO₂-only **decoupled** `k_death_mlf` from the unsulfited
conversion, so k re-tuned 0.02→0.05/h. **v1 tradeoff (owned):** no SO₂ ⇒ no death
(bacteria persist); slow ethanol/age decline deferred v2. Commit 2: `rack`
removes both `X_mlf` pools off the lees (asymmetry: viable *yeast* stays — it's
in suspension; O.oeni is lees-associated). Diacetyl lock-in demo is confounded
(kill removes diacetyl source AND sink) ⇒ validated on `X_mlf` directly.

**D-41 (MLF v2 benign senescence, landed — closes the MLF arc's last deferral):**
new `MalolacticSenescence` gives *O. oeni* a small always-on **baseline mortality**
`r_sen = k_senescence_mlf·X_mlf·arrhenius(T)` into the same `X_mlf_dead` pool, so a
pitched untreated wine slowly declines over **weeks-to-months** (lifts D-39's "no
SO₂ ⇒ never die"). **Environment-free** (no pH/ethanol/SO₂ — the D-39 Luong-wall
wipeout crux reused: any ethanol term ⇒ ~1-week kill) and **Arrhenius not γ(T)**
(warm accelerates, cold preserves; γ(T) would peak at the growth optimum —
backwards). **Separate isolable Process** ⇒ D-39 SO₂ kill stays byte-for-byte;
total mortality now `r_sen + r_death`. New speculative `k_senescence_mlf` 5e-4/h
(t½ ~58 d, ~100× below the SO₂ kill), reuses `E_a_death_mlf`/`T_ref`; no `brentq`.
`X_mlf_dead` a **terminal sink** (autolysis reads only yeast `X_dead` — no recycling
loop, advisor-checked). §2.2 + 0.1813 deacidification control-diff unmoved. Growth-
isolation tests difference senescence out (mid-run `pitch_mlf` re-enables the whole
gated set). **MLF arc D-23→D-31→D-38→D-39→D-41 complete.**

**D-38 open fork for owner (still reversible):** MLF-growth conservation is
**N-anchored, carbon-from-sugar** (advisor's higher-fidelity pick); rejected
branch was C-anchored / arginine-deiminase. Contained swap if owner prefers ADI
— see DECISIONS D-38 "anchoring fork".

**D-40 pt2 (BrettGrowth, the numeric crux worth remembering):** X_brett is now
dynamic — nitrogen-anchored on `amino_acids`, carbon shortfall drawn from
**ETHANOL** (owner fork), so Brett grows in a *dry* finished wine → phenols
accelerate autocatalytically. Brett has NO self-arrest (no sugar Monod / ethanol
wall, unlike MLF), so an intrinsic logistic carrying-capacity brake
`(1 − X_brett/K)` is the only ceiling. **Advisor-caught BDF blow-up** (X_brett→23,
aa→−4.5 under default BDF, while RK45/LSODA were correct): the `E≤0` guard had no
smooth shadow, so BDF's finite-difference Jacobian straddled an on/off step at E=0
during primary AF. Fix = ethanol Monod `E/(K_E_brett+E)` (the missing shadow, and
physically right — Brett grows on ethanol). **Invariant learned: every hard guard
must be shadowed by a smooth factor that reaches zero first.** Regression pinned
under BDF + a BDF/RK45/LSODA agreement test.

**D-40 pt3 (BrettDeath, landed):** `X_brett → X_brett_dead` on **molecular SO₂
alone** (`1 − g_SO₂`), Arrhenius temperature (cold preserves, not cardinal γ(T))
— the MalolacticDeath (D-39) form. Cleaner than MLF: Brett's gate has NO
ethanol/pH term to strip out (Brett is ethanol/acid-tolerant), so "SO₂ alone
kills" is *directly* correct physics, not a confounder-fix. Neutral transfer
(both pools weighted since pt2). Pitch-gated (dies whether or not it grew);
racking already removes it (pt1). New speculative `k_death_brett` (0.03/h, below
`k_death_mlf` — Brett more SO₂-tolerant) + `E_a_death_brett`. Headline: SO₂
crashes a *growing* population (`X_brett_dead` accumulates + `X_brett` < dose,
vs still-growing control).

**D-40 pt4 (POF+ yeast, landed — closes the Brett arc + M2 physics):**
`YeastPOFDecarboxylation` = a POF+ *S. cerevisiae* cinnamate decarboxylase (the
*same* reaction/routing as `BrettDecarboxylation`, hydroxycinnamics→vinylphenols
+CO2, 9=8+1) but with **no reductase**, so it fills the shared `vinylphenols`
reservoir it cannot drain → with no Brett the vinylphenols **strand** (ethylphenols
stays 0 + VALIDATED); a later Brett gets a **head start**. **Forks (owner):**
separate opt-in Process (not strain-flag); pure-enable key **`pof_positive`**
(binary trait; rate stays YAML), **independent of `brett_pitch_gpl`**; carbon from
hydroxycinnamics (forced). **Flux-coupled** (catalyst = viable yeast `X` via
`fermentative_flux_shape`, NOT X_brett — runs during AF, stops at dryness);
**temperature-flat** (no E_a_pof, the `AcetolactateExcretion` precedent).
**Advisor test-design crux:** stranding is the PRIMARY headline (timing-independent);
head-start is an **early-time/time-to-threshold** claim — conservation forces
*asymptotic* ethylphenols **equal** (same total precursor), so never assert higher
*final* ep. New speculative `k_pof_decarb`=2.5e-6 (~49% conversion during AF).

**D-42 (H₂S CO₂-stripping sink, landed — closes the §3.2 aroma beat + M2 physics):**
`HydrogenSulfideVolatilization` sweeps volatile dissolved `h2s` into a new carbon-free
`h2s_gas` headspace pool on the CO₂ flux (`-k_h2s_volatil·flux·f_gas(E_a_uptake)·
f_part(dH_h2s_volatil)·h2s`), so `h2s` is now the µg/L **residual** and `h2s + h2s_gas`
is cumulative produced. The **exact ester D-20/D-21 precedent but carbon-free** ⇒ both
pools on no ledger, transfer neutral by construction (**no conservation.py change**;
produced-total invariant replaces the ester carbon-closure test). **Flux cancels in
the residual** (`h2s_ss = k_h2s·gate/(k_h2s_volatil·f_gas·f_part)`) ⇒ rises as N
depletes, freezes at dryness. New speculative `k_h2s_volatil`=1.0 (~99.7% stripped) +
**sourced** `dH_h2s_volatil`=17.5 kJ/mol (Sander −d ln kH/d(1/T)≈2100 K, exothermic ⇒
+sign). **Honest artifact flagged:** T-flat production + T-rising strip ⇒ residual
falls with warmer ferment (unbenchmarked). Always-on both media; new `h2s_gas` slot
(wine 32→33, beer 19→20); flipped run-level h2s reads → h2s+h2s_gas across 3 test files.

**D-44 (autolytic H₂S source + copper fining, landed — the two D-42 deferred items):**
`AutolyticHydrogenSulfide` = a **yield on the D-34 autolysis flux** (`y_h2s_autolysis·
k_autolysis·f_T·X_dead`, the D-33 recompute-the-producer idiom so `autolysis_rate_per_h`
moves peptide + sulfide release on ONE clock). **Non-flux-linked** is the crux: the D-42
CO₂ sink gates off at dryness, so this sulfide **accumulates un-stripped as residual**
post-dryness — the reductive fault (verified: rises day15→40, >5× default; default freezes).
Opt-in + wine-only, disabled *with* `YeastAutolysis` (both now in `_AUTOLYSIS_PROCESSES`);
undosed run byte-for-byte core (5 benchmarks pass, RUN not inferred). New speculative
`y_h2s_autolysis`=2e-5 g/g (~1% of biomass-S ceiling). `add_copper` verb precipitates H₂S as
CuS (1:1, `copper_h2s_binding`=0.536 g/g, plausible), removes `min(present, capacity)`,
ledger-neutral (carbon-free, `add_so2` precedent). **Two honesty caveats:** (1) emergent
post-dryness, NOT an AF/post-AF switch (X_dead exists during AF, where the sink still strips);
(2) `h2s_gas` "produced total" broadens to sulfate-reduction + stripped-autolytic once opted in.
**Mercaptan pool DEFERRED to owner** (carbon-bearing thiols — must be a real total_carbon
species, no carbon-free lump; formation murky) — see the open question below.

**D-45 (mercaptan pool + copper mercaptide, landed — closes the reductive-sulfur beat):**
owner chose **Option A** (over B): a carbon-bearing `mercaptans` pool (methanethiol, C1,
N-free; wine 33→34) filled by `AutolyticMercaptan` as a yield on the shared autolysis flux,
drawing the thiol carbon from `amino_acids` and **deaminating** the nitrogen to `N` (the D-33
FuselAminoAcidReroute idiom). Carbon+N close by construction (C into mercaptans = C out of
amino_acids; arg N → N pool). **Caveats worth remembering:** (1) the carbon is drawn from the
*arginine* `amino_acids` lump, NOT literal methionine (exact on ledger, approximate on
provenance; ~0.66 vs ~1 mol N/mol MeSH); (2) it's the **first autolysis-gated N-writer**, so
autolysis-on drops structural `tier_of("N")` PLAUSIBLE→SPECULATIVE even amino-dose-off (D-27 E
parallel) — verified + pinned. `add_copper` extended: binds **H₂S first** (CuS Ksp≫mercaptide),
leftover Cu binds mercaptans (Cu(SR)₂ 1:2, `copper_mercaptan_binding`=1.514 g/g, banded HARD —
copper is incomplete on thiols, useless on disulfides). **Mercaptan removal books a carbon
external flow** (thiols carry carbon), so its conservation test uses the flow-identity `final ==
initial + Σ flows`, NOT assert_conserved. Shared `autolysis_flux(y,schema,params)` helper now
backs all 3 autolysis Processes (D-34/D-44/D-45). New speculative `y_mercaptan`=1e-5 g/g (below
the H₂S yield). Opt-in with autolysis; 5 benchmarks byte-for-byte (run, not inferred).

**D-46 (solve_ph totality fix — D-45 shipped main RED, this greens it):** D-45's
extra `mercaptans` slot (wine 33→34) shifted **BDF's** adaptive stepping and exposed
a **latent** fragility (NOT a D-44/D-45 derivative bug — no-op'ing both new Processes
still repro'd it). BDF's `num_jac` Jacobian probe perturbs the `cation_charge` slot far
outside its physical ~0.03 mol/L range (to **3.81**); there `charge_residual` is positive
across the whole [0,14] bracket, so `brentq` in `solve_ph` threw "f(a) and f(b) must have
different signs" → 3 `test_brett.py` integration tests red. **RK45/LSODA green, cation
constant 0.0254** — pins it as a BDF-num_jac-only artifact. **Fix:** make `solve_ph`
**total** — `charge_residual` is strictly monotone-decreasing in pH, so both-ends-positive
⇒ return 14, both-negative ⇒ return 0, else brentq (physiological cation falls through ⇒
**bit-for-bit pH, byte-for-byte curves**; clamp tames only the Jacobian probe, never the
exact RHS). +3 direct `test_acidbase.py` unit tests pin totality at the FUNCTION level (the
Brett tests catch it only incidentally via a 120-day run). **Invariant learned: num_jac
probes every state var outside its physical range — any core helper with a bracket/`log`/
`sqrt` that assumes a physical domain must be total, not throw.** Sibling to the D-40-pt2
"every hard guard needs a smooth shadow" BDF lesson. Commit e087a3d.

**D-47 (SO₂-bound acetaldehyde protected from ADH — lands the D-28 deferred RHS
coupling):** `AcetaldehydeReduction` now reduces only **free** (unbound) acetaldehyde
(`acidbase.free_acetaldehyde`), so dosed SO₂ **locks in** acetaldehyde — a sulfited wine
strands a residual (~27 mg/L at a 50 mg/L pitch dose, ~0.78 mol/mol SO₂) and free SO₂
stays depressed. **Literature-grounded** (Han et al. 2020: "bound acetaldehyde could not be
metabolized … only free acetaldehyde"; K=1.5e-6 = the literature value; the ~0.76×
degradation slowdown reproduced at sub-stoichiometric field doses). Owner chose **bake-in
default-on** (MLF SO₂-gate precedent) — this **intentionally RETIRES the D-22/D-28 "SO₂
readout-only" invariant** for sulfited runs. Preserved exactly where it still holds: undosed
= byte-for-byte D-27 (exact `so2_total>0` guard, no per-RHS brentq), **no §2.2 benchmark
doses SO₂**, carbon still closes, pH still not a charge actor (SO₂ couples only via
charge-free acetaldehyde), core-ferment footprint ≤1e-3 (E→viability ripple only).
**Emergent downstream (faithful, re-pinned):** SO₂ dosed *during* AF is now only a
**partial** MLF brake (stranded acetaldehyde sequesters the antimicrobial pool — "bound SO₂
isn't antimicrobial"), and molecular SO₂ nets *down* over the run — the counterintuitive
flip. CAVEAT (speculative): bound treated inert-to-ADH ⇒ upper bound on persistence; pitch
dose = maximal stranding, post-AF dose strands ~nothing. **No new params.** `touches`/`reads`
unchanged (SO₂/pH params read inside the helper; output already speculative). Commit next.

**D-48 (SO₂-induced acetaldehyde over-production — the peak half D-47 scoped out;
commit f4b6c0b):** the induced-over-production term (glyceropyruvic redox pull, Han
2020) added to `AcetaldehydeProduction`: `+k_acet_so2_induced·flux·so2_total`, a
carbon-exact **borrow from E** (D-27 invariant), exact `so2_total>0` guard (undosed =
byte-for-byte core). **The crux worth remembering — the task premise was wrong, data
reshaped it:** D-47's caveat said protection captured only "half" and a production
"half" was missing from the finished-wine level. FALSE — **D-47 protection ALONE
already delivers 1.3–1.5× the field 0.39 mg/mg end-state slope** (25.7/56.1/119 mg/L
at 50/100/200 vs 19.5/39/78). End state is **capped by the SO₂-binding equilibrium
(D-28), not production** — an over-produced slice is reduced away; only the bound pool
survives; both drivers leave end-state ~unchanged (25.7→25.8). So no additive end-state
term fits; **owner chose Option 3 = scope D-48 to the transient PEAK** (real, distinct,
no end-state anchor). **Driver = TOTAL SO₂, reversing the owner's earlier free-SO₂
pick:** free SO₂ collapses to ~0 at the peak (all sulfite bound), so it's **empirically
inert** on the very observable (+0.1–1.4 mg/L) — a stability premise the data refuted;
total SO₂ is open-loop but flux-gated (no runaway, verified to 200 mg/L). **Magnitude
sizing = a cross-process reality constraint (owner asked "what is closer to reality"):**
the raised peak sequesters SO₂ and weakens the molecular-SO₂ MLF brake (real: "bound SO₂
isn't antimicrobial"); at the initial 1e-2 that brake dropped to 44% malic retained,
crossing below the literature "partial brake, >half retained" regime. **k set to the
LARGEST value keeping that regime: 4e-3** (malic 2.09/~48% converted; ceiling ~5e-3 at
the 2.0 floor) — represents BOTH real effects at max self-consistent magnitude, and the
MLF test stays in its **original** `2.0<malic<3.0` band (no re-pin). Peak lift +3.8 mg/L
@50, ~+15 @200. Speculative tier. Advisor caught the MLF scope-leak (surfaced to owner,
not silently re-pinned). **Method beat worth remembering: two advisor→owner escalations
turned a "just add a bump" task into a driver reversal + a reality-anchored magnitude —
the empirical check (end-state increment table) refuted the spec before code shipped.**

**D-49 (excreted overflow-pyruvate pool — part 1 of the D-47/D-28 overshoot fix):**
the 1.3–1.5× finished-wine acetaldehyde overshoot is a **real missing mechanism, not
a mis-calibration** — the model routes 100% of bound SO₂ onto acetaldehyde, but real
wine shares SO₂ with competing carbonyls (pyruvate, α-KG; Jackowetz & Mira de Orduña
2013). Owner chose to **build the competitor pools**; D-49 = pyruvate (D-50 α-KG,
D-51 the coupled multi-carbonyl SO₂ equilibrium that reads them). **Crux #1 — side
pool, NOT on-pathway precursor:** the "max-fidelity" rework (route acetaldehyde's
carbon *through* pyruvate) was designed and **rejected** (advisor retracted its own
rec) — it conflates the intracellular flux intermediate (never persists/binds SO₂)
with the extracellular excreted residual (persists/binds), its persistence needs a
fake mechanism ("SO₂ shields pyruvate from PDC" — PDC is intracellular), and it would
make dosed SO₂ **suppress** acetaldehyde (opposite of D-48). So **acetaldehyde/D-27/
D-47/D-48 stay untouched.** Model = D-19/D-26 side-pool idiom: `PyruvateExcretion`
(flux-linked, T-flat) draws C3 from `S`; `PyruvateReassimilation` returns to **E+CO2**
(C3→C2+C1, one mole each, carbon-closing like malic→lactic+CO2) NOT `S` (post-dryness
S=0 ⇒ refund-to-sugar destroys carbon). Wine-only, both speculative. **Crux #2 —
flux-link the reassimilation (mid-build mechanism fix, advisor-confirmed):** first
build borrowed ADH's **no-flux viable-X gate** → residual came out **0.0 mg/L**,
because a clean ferment ends with yeast **still viable** (~0.4 g/L) so no-flux
draining clears the pool over the long tail. Overflow-pyruvate reassimilation is
**co-metabolic** (tracks active ferment), the *opposite* of ADH. Fix: flux-link it
(share excretion's `X·S/(K+S)`) ⇒ at dryness **both terms freeze** the pool → residual
pegged to *end-of-ferment*, **crash- AND duration-independent** (30.0 mg/L at both 21
and 40 days). v1 simplification (documented): monotonic rise to plateau
`k_exc/k_reassim`, drops the real peak-then-decline (nothing reads the peak; D-51 reads
only the residual; growth-coupled excretion = option B, deferred). Sized by the
**ratio** 3e-3/1e-1 = 30 mg/L. Isolability: not byte-for-byte (routes a trace of sugar
carbon to ethanol) but ABV/CO2 endpoint delta **rel ~4.4e-5** ≪0.1%; carbon closes to
machine precision. New `test_keto_acids.py` (19); 6 existing `full_params`/schema tests
updated for the new shared `keto_acids.yaml` + `pyruvate` slot (wine 34→35).

**D-50 (excreted α-ketoglutarate pool — the crux worth remembering):** same excreted-
side-pool structure as D-49's pyruvate (`AlphaKetoglutarateExcretion`/
`AlphaKetoglutarateReassimilation`), with **one fork, advisor-caught before coding:**
pyruvate's `C3→C2(ethanol)+C1(CO2)` reassimilation is mole-for-mole *only because* 3
carbons happens to equal one Gay-Lussac fermentation unit (2 carbon→ethanol : 1
carbon→CO2) — that's *why* the detour is nearly isolable (rel ~4e-5 delta), not
because "return to E" is inherently safe. α-KG's C5 doesn't divide 1:1, so naively
copying mole-for-mole would divert reassimilation **throughput** (not just the
residual, ~10–20× larger) away from ethanol, threatening §2.2 ABV/CO2. Fix: return
carbon at the **same 2:1 ratio** instead (5/3 mol ethanol + 5/3 mol CO2 per mole).
The tempting "more faithful" alternative — route to succinate/`Byp` via the real
α-KG-dehydrogenase reaction — was also rejected: that enzyme is repressed under the
same anaerobic conditions that make α-KG overflow, and the real reassimilation fate
is N-coupled glutamate synthesis (unmodelled), so neither destination is "more true"
— fidelity wasn't the tiebreaker, isolability was. **Measured, not just threshold-
checked:** α-KG residual lands at exactly 20.0 mg/L (pyruvate unchanged 30.0),
combined ABV/CO2 isolability delta rel ~7.3e-5 (~2× pyruvate-alone's ~4e-5, as
expected from two detours). **CALIBRATION-PENDING flag for D-51:** both residuals
are order-of-magnitude estimates, not fits — D-51 must re-derive the pyruvate/α-KG
ratio against the field slope, and must work in **moles** (SO₂ binds molar
concentration): α-KG's higher molar mass (146.1 vs pyruvate's 88.06 g/mol) means
20 mg/L is only ~40% of pyruvate's molar contribution despite being 67% of it by
mass. New `total_carbon` weighting term needed in `conservation.py` (caught
immediately by the carbon test). 17 new tests (36 total in `test_keto_acids.py`).

**D-51 (coupled multi-carbonyl SO₂ equilibrium, worked in moles — the crux worth
remembering):** generalized `bound_so2_molar` from D-28's single-carbonyl closed-
form quadratic to N competing carbonyls sharing one bisulfite pool via a
**competitive-Langmuir partition reduced to one shared `brentq` root** (`h` =
reactive bisulfite; each carbonyl's bound share `Aᵢ·h/(Kᵢ+h)`) — proven (20-trial
numeric + regression-anchor test) to reduce **exactly** to the D-28 quadratic at
n=1. Now reads acetaldehyde + pyruvate + α-KG together, natively in **moles**
(resolving the D-50 mole-vs-mass flag). New `K_pyruvate_so2`/`K_alpha_kg_so2`
(Burroughs & Sparks 1973 — same paper whose acetaldehyde Kd matches the existing
`K_acetaldehyde_so2` exactly, a direct cross-check). **The honest finding (re-
derived, not inherited, per the D-50 calibration-pending flag):** measured against
the field 0.39 mg/mg slope at 50/100/200 mg/L SO₂, D-51 narrows the D-48 overshoot
(1.32/1.44/1.53× → **1.15/1.32/1.45×**) but does **not** close it. A sensitivity
check pushed both keto-acid residuals to the top of their literature-sourced bands
(pyruvate 100 mg/L, α-KG 70 mg/L — frozen state verified to land exactly there) and
got 0.86/1.10/1.29×: the 50 mg/L point now *undershoots* while 200 mg/L still
misses — the crossover signature of finite-capacity Langmuir competitors saturating
against a field regression that stays linear across the dose range. **This is
structural, not a value not yet found** — no in-band pool resizing closes it
uniformly. Per the owner's guardrail ("do not force-fit beyond the literature-
sourced pool ranges") and advisor concurrence, **shipped D-49/D-50 residuals
(30/20 mg/L) unchanged**; the remaining gap is a recorded open deferred item (needs
a different structure — e.g. a dose-scaling binder — not more pool mass).
**Reshapes its own task premise** ("D-51 is the actual fix") the same way D-48 did
— real, correct, load-bearing progress, not closure. Genuine side effect fixed
honestly: the always-on keto-acid pools now also compete for bisulfite in the
SO₂-dosed MLF integration test (`test_malolactic.py`), so free/molecular SO₂ ends
lower (~21%→~15% of an 80 mg/L dose) and MLF edges just past halfway converted
(~51%, was ~48%) — measured and re-pinned, not loosened blindly. 650 passed
(646+4 new D-51 tests, incl. 5 benchmarks), ruff+mypy clean.

**D-52 (MLF v2 — bounded ethanol/starvation stress multiplier on `MalolacticSenescence`, landed
2026-07-07):** with M2 physics complete through D-51, owner asked to close "whichever [MLF v2
deferred item] is closer to reality." **Advisor reversed the first-pass pick before any code was
written:** the initial read favoured a `BrettSenescence` twin, but DECISIONS D-40 pt3 already states
Brett's *indefinite persistence without SO₂* is an intentional fidelity choice ("an honest reflection
of how tenacious a barrel Brett infection is") — so a senescence twin there would be a fidelity
*downgrade*. Verified `malolactic.py`'s `r_sen` is a *tiny* rate (~100× below `k_death_mlf`) with no
multiplier, unlike D-39's crux (a *large* full-kill-calibrated rate × unbounded `1−toxicity`≈0.92) —
so a *bounded* stress factor could lift D-41's "environment-free" deferral safely. Model:
`stress = 1 + k_senescence_ethanol_scale·[E/(E+ethanol_tolerance_mlf)] +
k_senescence_starvation_scale·[K_aa_mlf/(K_aa_mlf+amino_acids)]`, both terms smooth Monod-type factors
in [0,1) (not a Luong-wall near-binary shape), hard-capped at 2.5× (worst-case half-life ~23 d,
verified empirically — never near the ~1-week wipeout regime). Reuses existing concentration scales
(`ethanol_tolerance_mlf`/`K_aa_mlf`); only 2 new dimensionless ceiling params. Genuine side effect
measured + banded honestly (not loosened blindly, the D-51 discipline): MLF-diacetyl clearing test's
final/peak ratio rose to ~0.861 (faster X_mlf decline late in a 30-d run ⇒ less bacterial reductase).
**654 green** (650+4 net new tests) + 5 benchmark, ruff+mypy clean. **Method beat worth remembering:**
THREE advisor() passes shaped this arc, each catching what the previous one missed. Pass 1
(pre-build) caught a wrong delegated-judgment pick by checking the codebase's own prior decisions —
"closer to reality" required reading D-40's Brett characterization, not just reasoning from first
principles. Pass 2 (post-commit "done?" check) caught that the wipeout-guard test implicitly fixed
T=20°C while its name claimed "worst case" — split into a T_ref-scoped test + a temperature-invariant
kill-ratio test (commit 3e47830). Pass 3 (post-research) caught that the D-53 fix (below) needed a
test *assertion flip*, not a re-band, and surfaced the "D-52 is now inert" honesty point before it
could be silently absorbed.

**D-53 (2026-07-07, correction — k_senescence_mlf was wrong by ~50×):** D-52 left an open calibration
question ("typical wine now declines in ~29 d, not D-41's ~2mo target — re-anchor k, or leave it?").
Owner asked for deep research instead of picking a number — **this overturned the premise entirely.**
Real, unsulfited finished wine shows NO detectable spontaneous O. oeni decline for 3–5 months
(Windholtz et al. 2025, OENO One, doi:10.20870/oeno-one.2025.59.3.9346; Millet 2001 thesis: population
"maintained at ~10^6 CFU/mL" at 0 mg/L SO2 over 3 months). The steep decline D-41's citations implied
is actually SO2-DRIVEN (Kioroglou et al. 2020, doi:10.3389/fmicb.2020.562560) — already
`MalolacticDeath`'s (D-39) territory. D-41's citations (Ribereau-Gayon; Bartowsky & Henschke 2004)
supported general "SO2 controls spoilage LAB" practice, not a specific weeks-to-months spontaneous-
decline claim — a misread that propagated uncaught into D-52. **Fix:** `k_senescence_mlf` 5e-4 → 1e-5
(round, upper-bound-consistent — no source measures decline past 5 months). D-52's stress-multiplier
MECHANISM unchanged, only the baseline magnitude. **Honest consequence surfaced to owner, not
buried:** at this magnitude D-52's multiplier is empirically inert on every simulated timescale
(worst-case combined stress still gives a multi-year half-life) — owner chose to KEEP the structure
(least churn) over stripping it back to D-41's flat form. Test consequence was an assertion flip:
`test_so2_crashes_bacteria_over_the_slow_senescence_baseline` asserted a *measurable* decline the new
evidence contradicts — renamed to `..._near_stable_...` and flipped to assert near-stability (~0.990
at day 21, was ~0.608); diacetyl-clearing test rationale corrected too (~0.742, was ~0.861). 654 green
(same count, reassigned not added), ruff+mypy clean. **Deep-research infra note:** the workflow's
`resumeFromRunId` path silently dropped `args` after a session restart (errored "No research question
provided") — had to relaunch fresh with `name`+`args` instead of `scriptPath`+`resumeFromRunId`.

**D-54 (POF v2 pt1, `E_a_pof` — direction-checked before calibrated):** `YeastPOFDecarboxylation`
was D-40-pt4 temperature-flat; owner picked "POF v2" as the next work. **Advisor caught a
wrong-by-default move before any code:** cloning a nearby `E_a_*` value (e.g. `E_a_decarb`) would
have picked the wrong DIRECTION, because POF's rate is flux-coupled to sugar uptake — this
codebase's own D-19 "ordering constraint" (`E_a_esters`/`E_a_fusels`) says a flux-coupled
byproduct's NET (time-integrated-to-dryness) total scales as
`exp(-((E_a_byproduct-E_a_uptake)/R)(1/T-1/T_ref))`, so the sign is set by `E_a_pof` **relative to**
`E_a_uptake` (55,100 J/mol), not `E_a_pof` alone. Researched (not guessed) the real direction:
brewing practice on this exact enzyme (Pad1/Fdc1, wheat-beer/Weizen fermentation) is unambiguous —
**cooler ferments retain MORE clove/4-vinylguaiacol character**, opposite the esters/fusels
direction. So `E_a_pof`=25000 J/mol, set BELOW `E_a_uptake` (not above, as naive analogy suggests).
Two new tests split the two effects: raw-rate direction (`E_a_pof>0`, isolated via fixed
flux/precursor) vs. net finished-wine direction (full 12°C-vs-28°C scenario runs, empirically
confirms the algebra). v1's implicit `E_a_pof=0` already had the right direction by luck; v2 makes
it a genuine sourced term instead of an accidental placeholder.

**D-55 (POF v2 pt2 — the ferulic-acid split, 3 commits):** owner said "full precursor split, break
it into smaller pieces, your call." **The scope-collapsing fact:** `hydroxycinnamics` is booked as
p-coumaric acid LITERALLY (9 carbons, `M_P_COUMARIC` used throughout); ferulic acid is a DIFFERENT
10-carbon molecule (`10C → 9C vinylguaiacol + 1C CO2`) — so any fixed-ratio split of the existing
pool's OUTPUT breaks carbon closure by construction (a 9C precursor can't yield a 9C product + CO2).
This collapsed 3 candidate designs to one real option (genuine second precursor pool) vs. "document
as a limit" (the D-51 precedent) — not three co-equal choices. Built as: (1) scaffolding — new
species `M_FERULIC`/`M_VINYLGUAIACOL`/`M_ETHYLGUAIACOL`, 3 new wine-only slots (schema 36→39), their
`total_carbon` weighting; (2) decarboxylation branch — both `BrettDecarboxylation` and
`YeastPOFDecarboxylation` gained a ferulic branch via a shared `_decarboxylation_branch` helper;
(3) reduction branch + scenario wiring — `BrettVinylphenolReduction` gained the
vinylguaiacol→ethylguaiacol branch via `_reduction_branch`, plus `ferulic_acid_gpl` wired into the
compiler. **Relative kinetics sourced, not cloned:** Edlin et al. 1998 gives PAIRED Vmax/Km for both
substrates on a homologous decarboxylase in the SAME assay (~0.606× rate, ~0.742× half-saturation) —
real ratio-derived params (`k_brett_decarb_ferulic`, `K_hydroxycinnamic_ferulic`,
`k_pof_decarb_ferulic`), not independently re-estimated. **One honest gap left unclosed:** Tchobanov
et al. 2008 confirms the SAME reductase acts on both vinylguaiacol/vinylphenol (upgrades an
assumption to a fact) but gives no paired rate to derive a ratio, so `k_brett_reduction` is reused
unchanged for both branches — documented as a simplification, not silently assumed. 661 passed
(44 in `test_brett.py`, 8 net new across D-54+D-55), ruff+mypy clean, undosed default runs stay
byte-for-byte the validated core (new slots default 0, no benchmark doses them). Both D-40 pt4
POF-v2 deferrals (E_a_pof + the split) are now closed.

**D-56 (2026-07-07, first independent-data validation attempt — Varela et al. 2004): model runs
2-4x too fast, diagnosed not fixed.** Owner picked "validation against real data" as the next
direction post-M2. Found (via 104-agent deep-research sweep) two datasets genuinely independent of
Coleman 2007 (the paper the wine parameters are fit to): Varela 2004 (Chile, same EC1118 strain,
Table 1 endpoint values w/ replicate uncertainty — no digitization needed) and Palma 2012 (both
figure-only for real time series, contrary to the research report's initial claim — verified by
direct WebFetch). Ran the model at Varela's exact conditions (28°C, 240 g/L sugar, 300 vs 50 mg
N/L): model is **~2x too fast** even in-range (N=300, traces to the ALREADY-documented
`q_sugar_max` total-vs-active-biomass caveat) and **~4x too fast** at severe N-deficiency (N=50,
below Coleman's 70-350 mg N/L fit floor — an ADDITIONAL ~2x gap, isolated via a biomass-hours
integral argument to the nitrogen-limited growth kinetics, consistent with Bisson's hexose-
transporter-turnover literature). A single-term ethanol-driven N-gated fix was prototyped
(monkeypatched, no core files touched) and **structurally disproved** — narrowing the term to fit
the N-differentiation starves the N=300 correction, widening it to fix N=300 collapses the
differentiation; at least two mechanisms are needed. Stopped there deliberately: Varela is the
project's only independent dataset, and tuning ≥2 free params against 2 data points would burn it
as a validation set (the firewall: never fit to the one thing you're validating against). Landed as
`tests/benchmarks/test_validation_varela2004.py` — a real-data regression benchmark that
characterizes the CURRENT gap (catches it narrowing OR widening), not a physics fix. No tier
promotion (matching an aggregate endpoint doesn't license per-parameter VALIDATED bumps — non-
identifiability; also `ProcessSet.tier_of`'s honest path already floors wine S/X at SPECULATIVE via
`K_s`/`K_repression`/`Y_byproduct_sugar`, independent of fit quality). **Method beat: 3 advisor()
catches in one arc** — mechanical tier-promotion reality-check before coding, a confounded
20°C-vs-28°C comparison caught before it misattributed the gap to Arrhenius (contradicting D-14's
own line-for-line Coleman reproduction), and the validation/calibration firewall catch mid-sweep.
664 passed (661+3), ruff+mypy clean. Full diagnosis in DECISIONS.md D-56.

**D-57 (2026-07-07, correction + fix — D-56 finding 1 was a stale-note misdiagnosis; the real bug was
`k_prime_d`'s missing quadratic temperature scaling):** picked up D-56's "build a Bisson N-gated
transporter mechanism" as the next task; checking the premise against current code (not the D-56
record) before building anything found the premise half-wrong. **Mechanism 1 doesn't exist to fix**:
the "`q_sugar_max` applies to total biomass, no active/inactive split" note was written in D-12,
*before* D-13 added exactly that split (`X`/`X_dead`, verified byte-for-byte equivalent to Coleman's
own eqs 1-2 via `test_coleman_reconstruction.py`) — stale documentation, not a live gap. **The real,
sourced bug:** `k_prime_d` (Coleman's death-rate constant, his one QUADRATIC-in-T parameter) shipped
with NO temperature modifier at all ("M1 is isothermal" — true when written, stale once M2 added
temperature ramps). Every non-20°C run since has driven growth/uptake with Arrhenius scaling while
death stayed frozen at the 20°C rate — inert on short/high-N runs, compounding badly on long/
low-N ones, exactly why D-56 misread the gap as "worse at low N" (a nitrogen-transporter symptom that
was actually a T-scaling bug). Fixed: `ColemanQuadraticDeathTemperature`, a new `RateModifier`
implementing Coleman's own regression directly (not an Arrhenius approximation — the note itself
says a single E_a can't reproduce the curvature), 3 new sourced params
(`k_prime_d_a1`/`a2`/`t_floor`), wired into both wine and beer. Verified decisive: Coleman's own
reference model at Varela's 28°C inputs matches the engine almost exactly, both pre- and post-fix.
**Measured result:** N=50 duration gap to Varela narrowed ~4x→~2.2x; the N50/N300 ratio shortfall
against Varela's real 4.12x fell from ~1.94x-too-small to ~1.17x-too-small — real progress, not
closure. A second, independent correction surfaced while finishing this properly: the benchmark
compared Varela's biomass to viable `X` alone; WebFetch of the primary source confirmed Varela
measures TOTAL dry cell weight (gravimetric, dead+viable) — corrected to `X + X_dead`, which also
reproduces D-56 finding 3's already-known Y_X/N cross-study gap (~42% low at N=300, ~7% at N=50)
almost exactly. **Flagged, not fixed:** at N=50 the model's viable/dead split implies ~94-98% dead by
dryness while Varela reports >97% viable (LIVE/DEAD staining) — read as a vitality-vs-viability
category mismatch (X_dead = lost catalytic capacity, not membrane integrity), not evidence k_prime_d's
magnitude is wrong (Coleman's own model shows the identical crash); changing k_prime_d would break the
Coleman line-for-line reconstruction, out of scope. **Recommendation, not yet decided by owner:** the
residual ~1.17x N-specific gap is small enough that a Bisson-sourced transporter mechanism (D-56's
original ask) is no longer clearly worth the calibration/validation-firewall risk — owner's call
whether to pursue further or accept as a documented model limit. 664 passed (unchanged count — a fix,
not new tests), ruff+mypy clean. Two advisor() passes, both premise corrections (see D-57 for the
full method-beat writeup).

**D-58 (2026-07-08): MLF v2 sub-items resolved — `BrettSenescence` re-confirmed declined,
`BrettEthanolToxicity` built.** Picked up D-52's two remaining MLF-v2 sub-items (`BrettSenescence`
twin; separate `molecular_so2_death_scale`). Two parallel Opus research agents (opposite angles,
mirroring the D-53 method) **re-confirmed D-52: no `BrettSenescence` twin** — no literature source shows
Brett declining from elapsed time alone; every observed decline traces to SO₂, ethanol toxicity, or
substrate exhaustion. D-40/D-52's "persists indefinitely" wording softened: read as "no positive evidence
for spontaneous decline without SO₂," not literal immortality. **Side-effect finding, verified against
`brett.py` directly (not just agent claims), then BUILT:** Barata et al. 2008 shows Brett grows at 8%
v/v ethanol, dies above ~14% (no SO₂ needed) — the model's old `BrettGrowth` (logistic plateau brake,
never declines) + SO₂-only `BrettDeath` couldn't reproduce that bloom-then-death dynamic. Owner chose
to build it: a shared `brett_ethanol_survival_factor` helper (a **threshold** form, deliberately NOT
the whole-range MLF Luong wall — that would spuriously suppress Brett at ordinary wine strength,
verified via a 22-Brix probe topping out at E≈106.6 g/L, below the 110 g/L onset) feeds BOTH a new
`BrettGrowth` upper wall (arrests growth near the ~118 g/L ceiling) and a new sibling `BrettEthanolToxicity`
Process (kills X_brett→X_brett_dead above the onset, reusing `k_death_brett`/`E_a_death_brett`/`T_ref`
rather than sourcing new magnitude params — no independent activation energy exists since Barata ran at
one fixed 25°C). Scope-limited on purpose: does NOT attempt Barata's confounded 12%-ethanol/50-day
starvation-plus-ethanol result (below the onset by design). 12 new tests incl. a headline integration
test (26-Brix/~13%ABV unsulfited wine crashes; 22-Brix/~11%ABV control keeps growing). **676 passed
(664+12), ruff+mypy clean, zero existing test needed a re-band** (onset sits above every existing
scenario's finished-wine ethanol). Full design/build record in DECISIONS.md D-58.

**D-59 (2026-07-08, research only, no code changed): 6-agent validation-direction sweep.** Owner picked
"validation" (of validation/UX/new-physics) as the next direction, then asked to research all 3 open
validation sub-threads before building anything: N-gap (D-56/D-57), SO2/acetaldehyde overshoot (D-51),
broadening coverage. Ran 6 parallel Opus agents (2 opposite-angle per thread, the D-53/D-58 method).
**Headline finding: Coleman (fit) and Varela (validation) are the SAME strain lineage (Prise de
Mousse)** — the project's "independent" validation has been monostrain. **SO2 overshoot: RESOLVED
(accepted, not closed)** — affinity arithmetic proves no binder pool can help (acetaldehyde ~99% bound,
100-370x tighter than competitors) and the field anchor is a cross-sectional survey, not a titration;
documented as a structural limit, no new pool. **N-gap: not resolved, sharper framing** — the real
miss is that the model never reproduces fermentation *arrest* (Varela N50 got stuck at 16 g/L residual;
model always finishes dry), not the 1.17x ratio; recommended next step is a free internal diagnostic
(model N50 viable-biomass vs Varela's measured cells), not yet run. **Coverage: Palma 2012** (different
strain, n=3 + SD error bars) is a strong validation candidate, not yet digitized; **beer has no
accessible independent in-regime dataset** (best option is an off-regime lager cross-check, Speers 2003).
**Advisor caught a firewall conflict:** Palma 2012 was proposed as both a validation dataset AND a
mechanism-parameter source by two different agents — resolved to validation-only (a per-cell N
mechanism, if ever built, must source from Salmon 1989 instead). Full findings + citations: DECISIONS.md
D-59. Nothing built yet — owner has not chosen which recommended next step (N50 diagnostic / Palma
digitization / D-51 citation fix / beer lager check) to start.

**Open / next candidates** (details in DECISIONS "Deferred" + latest D-record):
- **D-56/D-57/D-59 N-specific gap:** real next step is the free N50 viable-biomass-vs-Varela diagnostic
  (not yet run) — resolves whether this is a biomass/death-rate issue (recalibrate existing terms) or
  a genuine per-cell transporter gap (build from Salmon 1989, not Palma — see D-59 firewall note).
- **Palma 2012 digitization** (D-59): different-strain validation dataset, glucose+ethanol only,
  scoped and ready to build whenever owner picks it up.
- **Beer independent validation** (D-59): no in-regime dataset exists publicly; decide whether to
  build the off-regime Speers 2003 lager cross-check or defer.
- ~~`BrettSenescence` twin~~ — **RE-CONFIRMED DECLINED in D-58 (2026-07-08),** now literature-checked
  (two independent research agents), not just folk wisdom: no source shows Brett declining from
  elapsed time alone. Revisit only if a source shows genuine age-based (not ethanol/SO₂/starvation)
  decay.
- A separate `molecular_so2_death_scale` for `MalolacticDeath` — pure parameter-architecture split,
  ruled out in D-52 as zero fidelity gain, still mechanically available if ever wanted (deprioritized
  again in D-58).
- ~~The residual D-51 overshoot gap (1.15–1.45× field, worst at high SO₂ dose)~~ — **RESOLVED
  (accepted, not closed) in D-59 (2026-07-08).** Affinity arithmetic + field-anchor category-mismatch
  both independently rule out a new binder pool; documented as a structural model limit. See D-59.
- **Milestone 3 / next direction** — **validation was picked (2026-07-08)** and researched in D-59;
  **UX and new-physics-scope remain untouched.** Next: owner picks among D-59's concrete next steps
  (N50 diagnostic, Palma digitization, D-51 citation fix, beer lager check) or pivots to UX/physics.
- ~~MLF v2 further refinements~~ — **ALL RESOLVED.** Benign baseline mortality: DONE in D-41.
  Ethanol/starvation modulation of the senescence baseline: DONE in D-52 (later found empirically
  inert at real-world magnitudes, D-53). `BrettSenescence` twin: re-confirmed declined in D-58.
  `molecular_so2_death_scale`: deprioritized (zero fidelity gain, D-52/D-58).
- ~~Mercaptan pool~~ — **CLOSED in D-45 (2026-07-06).** Owner chose Option A: real
  carbon-bearing `mercaptans` pool (methanethiol, C1, N-free), filled by `AutolyticMercaptan`
  off the shared autolysis flux; `add_copper` extended to bind H₂S first, leftover Cu to
  mercaptans (Cu(SR)₂ 1:2). Reductive-sulfur beat (H₂S + mercaptans) fully complete —
  see DECISIONS.md D-45. (Superseded the D-44 "OPEN" flag, which went stale after D-45 landed.)
- ~~Default-on residual-N model~~ — **DECLINED in D-43 (2026-07-06, spike, no source
  change).** A mass-balance argument proved default-on residual *assimilable* N is
  **Coleman-incompatible regardless of mechanism** (two-pool/cell-quota/satiation):
  Coleman builds biomass by ~day 1.3, pinning external assimilable N to ~0 by then for
  every dose, so nothing biomass-preserving can widen the H₂S lever or leave a
  late-window residual without breaking the Coleman sugar curve. Keep the D-30 opt-in
  cap as-is. If the H₂S cross-must lever is ever wanted default-on, re-point the *H₂S
  gate* onto a dose-correlated proxy (an H₂S-model change), NOT the N model. The
  negative result **is** the deliverable — closes the thread open since D-23/D-29/D-30.
- ~~POF v2~~ — **CLOSED in D-54/D-55 (2026-07-07).** `E_a_pof` temperature dependence
  (direction sourced from brewing practice, set BELOW `E_a_uptake` per the D-19 ordering
  constraint) + the vinylguaiacol/vinylphenol split into a genuine second ferulic-acid
  precursor branch (ratio-sourced kinetics from Edlin et al. 1998). See DECISIONS.md D-54/D-55.
- **Next milestone/direction after POF v2** — M2 physics + its refinements are now all done
  through D-55; next direction (validation/UX/new physics) is owner's call. Ask before picking.
