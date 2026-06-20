# Fermentation Research Sandbox — Engineering Handoff

> **Audience:** Claude Code (and any human collaborators) building a simulation engine for
> wine and beer fermentation, grounded in published science where possible and clearly
> labelled as speculative where not.
>
> **Status:** Greenfield. This document defines the architecture, the modelling tiers, the
> validation discipline, and a build sequence. It does **not** prescribe code; treat the
> equation references as a reading list, not a spec to transcribe blindly.

---

## 0. Project intent and prime directives

This is a **research sandbox**, not a game and not a homebrew calculator. The bar is
correspondence with reality, not fun or convenience. Three principles override everything
else and should be enforced in code and tests, not just honored in spirit:

1. **Fidelity is tiered and explicit.** Every modelled quantity belongs to one of three
   confidence tiers — *validated*, *plausible*, *speculative*. The tier is metadata that
   travels with the value all the way to any output or plot. The engine must never silently
   blend a validated concentration with a speculative one and present the result as equally
   trustworthy.

2. **Parameters are data with provenance, never magic numbers in code.** Every kinetic
   constant carries its source (DOI/citation), the conditions it was measured under, units,
   an uncertainty range, and a confidence tier. A parameter we guessed is tagged as a guess.
   This is what makes the difference between a sandbox a researcher can trust and a toy.

3. **The validated core is built first and protected.** We start with the narrow, well-
   established case (single-strain isothermal sugar fermentation), validate it against known
   curves, and only then expand outward. Speculative layers are physically isolated so they
   cannot contaminate the core's numerics or its tests.

The intended growth path is **validated → plausible → speculative**, and the architecture
exists to make each expansion an *addition* rather than a rewrite.

---

## 1. Architecture

### 1.1 Layering

Four layers, strictly one-directional dependencies (lower layers know nothing of higher):

```
  scenario / analysis  (declarative recipes, parameter sweeps, plotting, Monte Carlo)
          │  consumes time-series, owns no physics
  ────────┼─────────────────────────────────────────
  runtime  (time-stepping, event queue, phase switching, stochastic wrapper)
          │  integrates the core, knows nothing of UI
  ────────┼─────────────────────────────────────────
  domain core  (state vector + Process objects that contribute rates)
          │  pure, deterministic, no I/O, heavily unit-tested
  ────────┼─────────────────────────────────────────
  parameter store  (versioned data files: values + provenance + tier)
```

The **domain core has no UI, no file I/O, no global state, and no randomness.** Given a state
and a parameter set it returns derivatives. That purity is precisely what makes it testable
against benchmark curves, which is the whole game for a research tool.

### 1.2 The state vector

A single typed structure holding everything that evolves in time. Start minimal; it grows
as tiers are added. Initial contents:

- viable biomass `X`
- fermentable sugar(s) `S` — scalar for wine v1; a vector (glucose / maltose / maltotriose)
  for beer, because uptake is sequential
- ethanol `E`
- yeast assimilable nitrogen (YAN) `N`
- temperature `T`
- evolved CO₂ (the experimentally measurable proxy — fermentation is tracked by weight loss
  / CO₂ release, so this is a primary validation channel, not an afterthought)

Each scalar should be able to carry its confidence tier and, optionally, an uncertainty band.

### 1.3 The Process abstraction

The core's central idea: a **Process** is anything that contributes to the time derivative of
the state. Primary fermentation is a process. Malolactic fermentation is a process. Oxidation
is a process. Oak extraction is a process. Each Process:

- reads the current state and parameters,
- returns its contribution to `d(state)/dt`,
- declares which state variables it touches and which tier it belongs to,
- can be individually enabled, disabled, or swapped.

The total derivative is the sum of active Processes' contributions. This compositionality is
the key design bet: additive × organism × temperature × vessel × time is a combinatorial space
far too large to validate exhaustively, so we model **mechanisms** and let combinations *emerge*
rather than scripting outcomes. It also means a speculative Process can be toggled off and the
validated core still runs and still passes its tests.

### 1.4 Runtime: events + continuous integration

Between interventions the state evolves continuously (ODE integration). Interventions —
*add SO₂, pitch an MLF culture, rack to barrel, step the temperature, add DAP, dose oxygen* —
are **timed events** that either mutate the state or change which Processes are active. So the
runtime is an event-driven loop: integrate to the next event, apply it, re-evaluate the active
process set, continue. Phase transitions (e.g. primary fermentation completing) are themselves
events, possibly triggered by state conditions rather than wall-clock time.

### 1.5 Scenarios as data

A run is fully described by a declarative document (YAML/JSON or a small DSL): initial
composition, organism(s) and strain, temperature schedule, vessel geometry/material, and a
timeline of interventions. **No physics lives here.** This makes parameter sweeps, Monte Carlo,
and scenario sharing trivial, and keeps the engine reusable across wine, beer, cider, and mead
without code changes.

### 1.6 Numerical strategy (read before integrating anything)

- **Stiffness is guaranteed.** Fast transients (CO₂ dissolution, acetaldehyde spikes) coexist
  with slow drift (months of aging). Use an implicit adaptive solver — BDF or Radau
  (`scipy.integrate.solve_ivp` is the obvious first home). Do not reach for explicit RK and
  fight instability.
- **Multi-scale time spans minutes → years.** Don't integrate aging at fermentation
  resolution. Use **phase-based integration**: each phase activates a relevant Process set and
  an appropriate step-size regime. Consider quasi-steady-state approximations for fast
  variables during slow phases.
- **Conservation as an invariant.** Carbon, nitrogen, and mass balances should hold to
  tolerance. Encode these as runtime assertions and as tests — a model that quietly creates
  carbon is broken regardless of how good the curves look.
- **Determinism by default, stochasticity as a wrapper.** The core is deterministic. Realism
  and replicate variation come from a runtime layer that samples parameters within their
  provenance-declared uncertainty and runs ensembles. Keep this *outside* the core so single
  runs stay reproducible and debuggable.

### 1.7 Language

Prototype and validate the kinetics in **Python** (scipy/numpy ecosystem is the fastest path
to a validated core, and the literature's reference implementations live there). Keep the
core pure so that, if performance or shipping later demands it, it can be ported to Rust/TS
behind the same scenario interface. Validate first; optimize never-prematurely.

---

## 2. Tier 1 — Solidly modelable (build this first)

This is real, published, validatable science. The first milestone lives entirely here.

### 2.1 Primary alcoholic fermentation

A coupled ODE system over `{X, S, E, N, T}`. The standard, well-supported ingredients:

- **Growth** — logistic or Monod-type, biomass tied to substrate and especially nitrogen
  availability. Wine fermentations are typically **nitrogen-limited**: cells grow while YAN
  lasts, then stop dividing but continue fermenting at a declining rate. This mechanism is the
  textbook cause of *sluggish* and *stuck* fermentations and must be in v1.
- **Ethanol inhibition** — growth and fermentation rate decline as ethanol rises; viability
  collapses past a strain-specific tolerance (commonly ~14–16% ABV). Various functional forms
  exist (linear, exponential); pick one with literature support and tag the constant.
- **Temperature dependence** — Arrhenius on the rate constants; temperature also shifts max
  population and ethanol tolerance. Warmer ferments faster but stresses yeast earlier.
- **Sugar → ethanol stoichiometry** — theoretical Gay-Lussac yield ≈ 0.51 g ethanol / g
  glucose; realised ≈ 0.46–0.48 once biomass and glycerol are accounted for.
- **Early oxygen / sterol requirement** — a brief aerobic phase supports sterol and
  unsaturated-fatty-acid synthesis that underpins later ethanol tolerance. Can be simplified
  or stubbed in v1, but leave the hook.

**Beer differs in three concrete ways** worth handling early because they're well-characterized:
sequential sugar uptake (glucose, then maltose, then maltotriose) which makes `S` a vector;
the apparent-vs-real attenuation distinction; and the **diacetyl rest** (yeast produces then
reabsorbs diacetyl — a defining lager parameter, see Tier 2).

**Reading list (kinetic models, not gospel — reconcile, don't transcribe):**
- Coleman et al. (2007), temperature-dependent nitrogen-limited wine fermentation model.
- Cramer et al. (2002); Malherbe et al. (2004) — wine fermentation kinetics.
- Gee & Ramirez (1988); de Andrés-Toro et al. (1998) — dynamic beer fermentation models.

### 2.2 Validation benchmarks (write these as tests before tuning the model)

- A ~24 °Brix must (≈240 g/L sugar) at 20 °C ferments to dryness in roughly **10–14 days**,
  with a visible lag → exponential → stationary biomass trajectory.
- A ~1.048 OG ale wort at 20 °C attenuates to roughly 1.010 over **~5–7 days**.
- CO₂ evolution rate rises to a peak, then tails off — its integral tracks sugar consumed.
- Lower fermentation temperature ⇒ slower but generally cleaner ferment.

These are the acceptance criteria for the first milestone. **Test-drive the model:** encode the
benchmark, then iterate parameters and functional forms until it passes — this is exactly where
Claude Code is strongest.

### 2.3 Parameters required (each needs full provenance)

μmax, half-saturation constants (Ks for sugar and nitrogen), ethanol-inhibition constant(s) and
tolerance threshold, biomass and ethanol yield coefficients, Arrhenius pre-exponential and
activation energy per rate, maintenance/decay rate, maximum population. Expect to spend more
time sourcing and reconciling these (strain- and condition-specific, scattered, sometimes
contradictory) than writing solver code. Budget for it.

---

## 3. Tier 2 — Modelable with effort (plausible)

Sound mechanisms with literature support, but more parameters, more coupling, and harder
validation (good datasets are often proprietary). Add these one Process at a time, each behind
its own tests.

### 3.1 Microbial ecology (multi-organism competition)

Extend single-organism kinetics to resource competition (extended Monod; Lotka–Volterra-style
interaction for shared substrates). Each added organism multiplies the parameter count, so add
deliberately:

- **Malolactic fermentation** — *Oenococcus oeni* converts L-malic → L-lactic acid + CO₂;
  deacidifies and nudges pH up ~0.1–0.3. Its own growth is sensitive to pH, ethanol, SO₂, and
  temperature. Naturally modelled as a second organism Process activated by a "pitch MLF" event.
- **Brettanomyces** — slow grower producing volatile phenols (4-ethylphenol, 4-ethylguaiacol)
  from hydroxycinnamic acids. The canonical "funk" mechanism.
- **Sour-beer consortium** — *Lactobacillus*/*Pediococcus* (lactic acid), *Brettanomyces*, and
  *Acetobacter* (acetic acid, oxygen-dependent). Genuine microbial succession; a good stress
  test of the competition model.

### 3.2 Temperature- and metabolism-driven byproducts

Model the **major** flavor-active compounds (not the full volatome) as outputs of yeast
metabolism modulated by temperature and nitrogen status:

- **Esters** (isoamyl acetate ≈ banana, ethyl hexanoate, ethyl acetate) — favored by warmth and
  strain; coupled to nitrogen dynamics.
- **Higher / fusel alcohols** (Ehrlich pathway from amino acids) — rise with temperature; the
  relationship to nitrogen is non-monotonic, so model the pathway, don't fit a slope.
- **Diacetyl (VDK)** — valine pathway; produced then reabsorbed (the diacetyl rest). Time-and-
  temperature dependent; central to lager quality.
- **Acetaldehyde** — transient intermediate with a characteristic early peak.
- **Glycerol** — ~5–8 g/L byproduct; mouthfeel contribution; useful carbon-balance check.
- **H₂S** — linked to nitrogen/sulfur deficiency; a useful early-warning signal to model.

### 3.3 Additives with clear mechanisms

- **SO₂** — free/bound/molecular equilibrium, strongly **pH-dependent** (molecular SO₂ is the
  antimicrobial fraction; governed by pKa ≈ 1.81). Binds acetaldehyde and other carbonyls.
  Antimicrobial action above a molecular-SO₂ threshold. This couples directly to the pH model.
- **Nutrient additions (DAP, complex nutrients)** — feed straight into `N`; a lever on stuck
  fermentations.
- **Hop bittering** — alpha-acid → iso-alpha-acid isomerization (utilization as a function of
  boil time, gravity, temperature) → IBU. Pre-fermentation/wort-side but belongs in the engine.
- **Acid/sugar adjustments** — tartaric acid additions, chaptalization — simple state mutations
  via events.

### 3.4 pH / acid–base system (harder than it looks; prioritize it)

A proper pH requires a proton/charge balance across the weak-acid system (tartaric, malic,
lactic, acetic, carbonic) with titratable acidity tracked alongside. pH then feeds back into
*almost everything*: yeast and bacterial growth, SO₂ speciation, and reaction rates. Treat the
pH solver as core infrastructure, not a byproduct — many Tier-2 mechanisms are wrong without it.

### 3.5 Validation note

Exact-curve validation gets scarce here. Where you can't validate quantitatively, validate
**qualitatively and directionally** (warmer ⇒ more esters; SO₂ ⇒ suppressed spoilage; MLF ⇒
malic falls, pH rises) and tag outputs *plausible*, not *validated*.

---

## 4. Tier 3 — Frontier (speculative; isolate and label loudly)

Real chemistry underlies all of this, but integrating it into a trustworthy prediction is **not
solved science.** Build it as clearly-fenced, swappable Processes tagged *speculative*, never
wired so tightly that they can perturb the validated core or its tests.

### 4.1 Aging — the "years" axis

Individual reactions are modelable; their *sum over years into a verdict* is not. Candidates:

- ester formation/hydrolysis equilibria shifting over time,
- oxidation (phenolics, acetaldehyde generation, browning, Strecker degradation),
- tannin–anthocyanin polymerization (red-wine color and astringency evolution),
- **oak extraction** — diffusion-limited release of vanillin, whiskey lactones, ellagitannins,
  and toast compounds, as a function of surface-to-volume ratio, toast level, and barrel age,
- micro-oxygenation / reductive–oxidative balance and sulfide evolution,
- long-aging Maillard chemistry, sotolon in oxidative styles.

Each is a defensible mechanism; the honest framing is "tracking a handful of key compounds with
literature-based kinetics," **not** "predicting what five years in barrel tastes like."

### 4.2 Sensory / flavor prediction

This is explicitly **heuristic, not physics.** A reasonable approach: map compound
concentrations to **odor activity values** (OAV = concentration ÷ perception threshold), then
project OAVs onto a descriptor space. Implement it as a **separate, swappable sensory model**
tagged *speculative*, with a clean seam so it can later be replaced by an ML model trained on
real sensory-panel data should such data become available. The cardinal rule: the sensory layer
consumes the chemistry; the chemistry never depends on the sensory layer.

### 4.3 Keeping it honest

The temptation in a "based on real science" tool is to let plausible-looking speculation borrow
the credibility of the validated core. The tier metadata (§0.1) exists to prevent exactly that.
Any plot, export, or report must surface the tier of what it's showing.

---

## 5. Cross-cutting requirements (apply to every tier)

- **Provenance store.** Parameters live in versioned data files; each entry: value, units,
  source citation/DOI, measurement conditions, uncertainty range, confidence tier. No exceptions.
- **Confidence propagation.** Tier metadata travels with values into outputs. Outputs declare
  their lowest constituent tier.
- **Test-driven model development.** Benchmarks and conservation laws are tests. New Processes
  ship with their own qualitative or quantitative checks. The Tier-1 suite must stay green as
  later tiers are added.
- **Stochastic ensembles** as a runtime wrapper over the deterministic core (parameter sampling
  within provenance bounds; Monte Carlo for replicate spread).
- **Compositionality over scripting.** Model mechanisms; never hardcode the outcome of a
  specific additive/organism/temperature combination.

---

## 6. Suggested build sequence

1. **Skeleton + parameter store + validation harness.** Layering, the Process interface, the
   provenance-backed parameter store, and a test harness that can assert against benchmark
   curves and conservation laws — *before* any real kinetics. Get the architecture honest first.
2. **Tier-1 minimal core.** Single strain, isothermal, single sugar, nitrogen-limited
   fermentation. Pass the §2.2 wine and beer benchmarks. This is Milestone 1 and the project's
   foundation.
3. **Runtime maturation.** Event queue, temperature schedules, phase switching, stochastic
   wrapper. Re-validate.
4. **Tier-2, one Process at a time,** each behind its own tests: pH/acid system first (it
   unblocks the rest), then SO₂, then byproducts, then second organisms (MLF, then mixed
   cultures). Tag everything *plausible*.
5. **Tier-3 as isolated, clearly-labelled modules:** aging chemistry, then a swappable sensory
   layer. Tag *speculative*; keep the Tier-1 suite green throughout.

---

## 7. Open decisions for the human (please resolve before/early in build)

- **Wine first, beer first, or shared core from day one?** A shared core is the design intent,
  but the first validated benchmark must target one. Recommend wine for Tier-1 (cleaner single-
  sugar nitrogen-limited literature), beer immediately after (it exercises the sugar-vector and
  diacetyl machinery).
- **Do we have, or can we obtain, any real fermentation datasets** (time-series of
  Brix/gravity, cell counts, temperature)? The credibility of the "research" claim scales with
  this. Proprietary winery/brewery/yeast-supplier data, if reachable, is gold.
- **How rich must the pH/acid model be at first** — full proton balance, or a tracked-pH
  approximation that we upgrade later? This affects how soon Tier-2 mechanisms become trustworthy.
- **Units and conventions** to standardize up front: SI internally vs. industry units
  (°Brix, SG, °P, ABV, IBU) at the boundaries; a single canonical internal representation with
  conversion only at I/O edges.
