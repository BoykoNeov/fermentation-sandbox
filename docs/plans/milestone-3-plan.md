# Milestone 3 — Tier-3 (speculative frontier): sensory/OAV + aging chemistry

> Status: **beat 1a built (D-67); aging axis opened (D-68); first aging Process built (D-69);
> aging-phase scenario wiring landed (D-70); oxidative aging axis opened (D-71); SO₂ scavenging
> — the first O₂ sink — built (D-72); O₂ sub-axis reworked for always-on sinks (D-73);
> `PhenolicBrowning` — the first always-on sink, and the first VISIBLE oxidative product (the
> `A420` browning index) — built (D-74); `StreckerDegradation` — the Strecker aldehydes
> (methional/phenylacetaldehyde), a wine-only substrate-gated O₂ sink — built (D-75).**
> The OAV sensory readout ships (`fermentation.sensory` + `sensory.yaml` + `tests/test_sensory_oav.py`).
> **D-69 built `EsterHydrolysis` (`core/kinetics/aging.py` + shared `aging.yaml` + `tests/test_aging.py`):
> net decay of the `esters` pool toward a lower equilibrium floor, Arrhenius warmer-ages-faster, released
> carbon split 5:2 → `fusels` + `Byp` (advisor-settled crux; carbon closes to machine precision).
> **D-70 wired it into the scenario pipeline (753 tests green): `EsterHydrolysis` into both media
> (disabled at compile — aging is inherently post-ferment), a `begin_aging` intervention verb (the
> `pitch_mlf` reconfigure pattern MINUS the state mutation) that enables it over a post-fermentation
> segment whose span is set by `duration_days`, and `aging.yaml` into `compile.py`'s `shared_files`. The
> advisor reframe: the load-bearing call was the deferred "what runs during aging" — settled Stance A
> (leave the ferment set on) because every producer of `esters`/`fusels`/`Byp` is fermentative-flux-gated
> and quiescent at dryness, so the aging signal is unconfounded. §7 slow-phase = the segment restart (no
> new machinery). An un-aged run is byte-for-byte the pre-aging core (isolable).**
> **D-71 opened the OXIDATIVE aging axis: `OxidativeAcetaldehyde` (dissolved O₂ oxidises ethanol →
> acetaldehyde, the 'sherry'/oxidised note the D-67 lens already reads) on a new carbon-free `o2` state
> pool + an `add_oxygen` dosing verb + 3 `aging.yaml` params (768 tests). Advisor crux: O₂, not ethanol,
> is the rate-limiter — first-order in the finite `o2` pool ⇒ SATURATING acetaldehyde, not the unbounded
> ethanol-first alternative; the owner chose the dissolved-O₂ pool (Approach B) over the smaller
> ethanol-first form, so `o2` is the shared substrate the whole future oxidative sub-axis (browning /
> Strecker / SO₂ consumption) will draw down. Carbon closes machine-precision (E→acetaldehyde, the D-27
> reduction reversed; O₂ off every ledger); reductive aging (`begin_aging` with no `add_oxygen`) is
> byte-for-byte the ester-only case.** D-72 built `SulfiteOxidation` (the first O₂ sink, wine-only,
> substrate-gated: "SO₂ protects until exhausted"); D-73 reworked the O₂ sub-axis so `k_ethanol_oxidation`
> is a *share* not the total, letting *always-on* sinks compose; D-74 built `PhenolicBrowning` — the first
> always-on sink and the first VISIBLE oxidative product: medium-agnostic, dominant O₂ share (diverts O₂
> from — and suppresses — oxidative acetaldehyde), accumulating a new `A420` browning-index state slot (an
> optical absorbance, off every ledger; the `iso_alpha` pattern, not the D-67 post-hoc OAV). **D-75 built
> `StreckerDegradation` — the O₂/amino-acid Strecker aldehydes methional (cooked-potato off-note) +
> phenylacetaldehyde (honey), the first aging Process to add aroma pools the D-67 lens did not read (two new
> wine slots + thresholds). CRUX (owner fork, superseding the D-71→D-74 forward-guess): gating the O₂ draw on
> `amino_acids` makes it DOUBLY substrate-gated (o2 AND amino_acids), like `SulfiteOxidation`, so it ADDS ON
> TOP with NO re-baseline of `k_ethanol_oxidation` — the "reduce k_ethanol again" note is retired. Carbon +
> nitrogen close (the D-45 mercaptan draw+deaminate idiom + a CO₂ decarb term); wine-only, silent without both
> substrates. Owner also chose TWO pools over one lump (opposite sensory valence).**
> **D-76 closed the emergent SUR-LIE → Strecker pathway with NO new physics: lees autolysis (`YeastAutolysis`,
> D-34) and `StreckerDegradation` (D-75) simply COMPOSE — opting into `autolysis_rate_per_h` + `add_oxygen` +
> `begin_aging` with no `amino_acids_gpl` dose lets dead lees self-digest, refilling `amino_acids` from the
> wine's own dead yeast, which the Strecker route degrades (both aldehydes emerge from the physically-real
> nitrogen source — the DIRECTIONAL result; the absolute level runs ~8× the D-75 dosed literature anchor, an
> order-of-magnitude figure not a prediction, a recorded open item in D-76). Owner fork (chose A: verify +
> document + test) over B (re-gate autolysis to
> the aging phase); a discriminating measurement decided it — the pre-dryness active-ferment release is ~15 mg/L
> (bounded-small), the ~385 mg/L breakpoint pool is legit post-dryness sur-lie, so autolysis-from-t0 needs no
> re-gating. Carbon + nitrogen close on the new triple-draw compose (autolysis refills while Strecker + mercaptan
> draw). Test-only change (helper kwarg + 3 scenario tests); zero production code.
>
> **D-77 built `OakExtraction` — the barrel/chip aroma-extractive axis, the FIRST NON-oxidative aging Process and a
> SEPARATE axis (draws no O₂). Four wood extractives — `whiskey_lactone` (coconut, light-toast), `vanillin` (vanilla,
> medium peak), `guaiacol` (smoky, heavy — the OAK guaiacol, distinct from Brett 4-EG), `eugenol` (clove, heavy) —
> diffuse in and rise toward a per-compound saturation ceiling: `d(C_i)/dt = k_oak·f(T)·max(0, ceiling_i − C_i)`, the
> inverse of EsterHydrolysis's decay-to-floor. `E_a_oak_extraction` deliberately WEAK (20 kJ/mol — diffusion, not
> reaction). The ceilings are SET-AND-HOLD off-ledger STATE slots the new `add_oak {oak_gpl, toast}` verb writes
> (`oak_gpl × toast-specific yield` from `oak.yaml`), enabled by `begin_aging`. KEY DESIGN TURN: the advisor's
> first-pass "mint the ceilings as provenance-backed params in the verb" recipe was overturned by primary-source
> evidence (verbs can't inject into the compiled ParameterSet; `param_update` is plain floats absent from
> `param_tiers`, and `begin_aging`-before-`add_oak` would KeyError mid-integrate) — a 2nd advisor pass agreed and
> dropped it; the state-slot dose (the `cation_charge`/`add_oxygen` idiom) meets every goal (provenance in the
> yields, D-1 moot since all oak pools floor SPECULATIVE, no KeyError window). OFF EVERY LEDGER (exogenous
> wood-derived, the `iso_alpha` precedent — cleanest aging Process, no `chemistry.py` change). Owner forks: 4
> compounds (+eugenol) and WINE-ONLY wired (like Strecker — medium-agnostic physics ≠ agnostic wiring). Explicit
> `ceiling ≤ 0` undershoot guard (the floor is 0, so `max()` alone fabricates extract). 8 wine-only slots, 838
> tests.**
> **D-78 (latest) built `EllagitanninOxidation` — the oak-tannin O₂-scavenging sink, the BRIDGE from the oak
> extractive axis to the O₂ sub-axis D-77 kept separate (854 tests, +16). Ellagitannin is oak's hydrolysable tannin:
> `OakExtraction` extracts a 5th pool (identical diffusion), and a new `EllagitanninOxidation` draws its share of the
> shared `o2` budget (bilinear `[o2]·[ellagitannin]`, the `SulfiteOxidation` form) and CONSUMES it. THE SPINE is oak
> PROTECTION: an oaked+oxygenated wine browns LESS (lower A420) and makes LESS oxidative acetaldehyde than an un-oaked
> wine at the same O₂ dose (the D-72 SO₂-protection threshold with a RENEWABLE buffer — the wood re-supplies tannin
> below the ceiling, so protection is sustained for months-to-years, unlike SO₂'s finite pool). SUBSTRATE-GATED on
> `[ellagitannin]` ⇒ adds ON TOP with NO re-baseline (the `k_ethanol+k_browning=5e-4` anchor untouched) even though
> it's a DOMINANT sink when present — proving the substrate-gated/always-on distinction, not the magnitude, is
> load-bearing. MASS-based consumption yield `y_ellag_per_o2` (g/g, no fake molar mass for the lumped macromolecule);
> off every ledger (both slots), so nothing conserved moves. ASTRINGENCY is a TASTE readout `analysis.astringency_series`
> (mg/L ellagitannin, IBU-exact, reads no threshold — the iso_alpha exclusion), NOT the OAV odor lens and NOT an
> A420-style slot (it tracks the current pool). Softening = ONE contributor (oxidative consumption); tannin–anthocyanin
> polymerization is the separate deferred beat (the `tannin` namespace left free). Toast ordering: ellagitannin DECLINES
> with toast (thermolabile). 1 advisor pass (endorsed all 4 forks). Next: beat 1b (descriptor projection), non-oxidative
> Maillard Strecker, tannin–anthocyanin polymerization (now unblocked), barrel fill-number, or barrel-beer oak — all
> deferred.**
>
> **PROGRESS (D-79, 2026-07-13): `TanninAnthocyaninCondensation` built — the DOMINANT red-wine astringency-softening +
> colour-stabilization mechanism, the EIGHTH aging Process, the SECOND non-oxidative one, and a THIRD separate axis on
> GRAPE pools (868 tests, +14). Free monomeric `anthocyanin` + condensed grape `tannin` (two GRAPE must inputs,
> `anthocyanin_gpl`/`tannin_gpl`, default 0 ⇒ white wine) condense (bilinear `k·f(T)·[anthocyanin]·[tannin]`,
> reaction-scale E_a) into a soft, SO₂/pH-STABLE polymeric pigment — the young-purple → aged-brick-red evolution.
> OAK- AND O₂-INDEPENDENT (the correctness crux): grape condensed `tannin` ≠ oak hydrolysable `ellagitannin`, and the
> Process draws NO o2 — a steel-tank red still polymerizes (reusing ellagitannin would wrongly require `add_oak`); it's
> the grape `tannin` the D-78 note left the namespace free for. DOUBLY substrate-gated ⇒ white/no-tannin wine
> byte-for-byte inert; off every ledger (grape-derived, no `chemistry.py` change). The polymeric pigment is a POST-HOC
> readout (`polymeric_pigment_series` = anthocyanin₀ − anthocyanin), NOT a slot — the A420 discriminator applied
> (anthocyanin's single fate ⇒ reconstructible), keeping v1 to 2 slots. Readouts: `astringency_series` now =
> (tannin + ellagitannin)×1000 (softens as tannin condenses; polymeric excluded = soft); new `color_series` counts
> free anthocyanin + polymeric pigment (colour RETAINED — the monomeric→polymeric shift, not vanishing; bleaching
> deferred). ACETALDEHYDE BRIDGE DEFERRED — the advisor caught it as a conservation trap (acetaldehyde is ON the carbon
> ledger; an off-ledger pigment consuming it fails `assert_conserved`) — it's the explicit named next beat. 1 advisor
> pass (before writing; adjusted 2 leanings, 1 conservation-breaking). Next: the acetaldehyde-bridged route, tannin
> self-polymerization, SO₂/pH bleaching (promotes pigment to a slot), beat 1b, Maillard Strecker, or barrel fill-number.**
>
> **PROGRESS (D-80, 2026-07-13): `AcetaldehydeBridgedCondensation` built — the acetaldehyde-bridged (ethylidene) route,
> the NINTH aging Process, the THIRD non-oxidative one, the SPLIT-LEDGER beat, and the first link from the oxidative
> sub-axis to red-wine colour (887 tests). Trilinear `k·f(T)·[free acetaldehyde]·[anthocyanin]·[tannin]`; a new ON-ledger
> `ethyl_bridge` slot captures the acetaldehyde carbon (weighted at cf(ethylidene)) so carbon closes non-trivially.
> HONEST-FRAMING: v1 delivered the carbon + O₂→pigment MECHANISM, NOT a colour behaviour change — `color_series` was
> O₂-invariant (superseded at D-81). 1 pre-work + 1 done-call advisor pass.**
>
> **PROGRESS (D-81, 2026-07-13): the SO₂/pH anthocyanin-BLEACHING beat — `AnthocyaninFading` built + polymeric pigment
> PROMOTED to a slot, so `color_series` now GENUINELY DECLINES (899 tests, two commits). User chose "Both (C)" at the
> design fork = reversible masking readout (A) + irreversible fade sink (B), a two-beat split; D-81 delivers B, D-82 (the
> masking readout) is the COMMITTED second half, still owed. The fade is O₂-COUPLED (bilinear `k·f(T)·[o2]·[anthocyanin]`
> on the shared o2 pool, anthocyanin → colourless `faded_anthocyanin`), so SO₂ colour-protection is EMERGENT (SO₂
> scavenges o2 via D-72, leaving less to fade) — advisor caught that the reversible-masking alternative would make colour
> RISE, not decline. HONEST O₂-GATING: colour declines UNDER O₂ exposure; an anaerobic red still holds flat via
> condensation. Three-slot identity anthocyanin + polymeric_pigment + faded_anthocyanin ≡ anthocyanin₀ (by construction).
> SUPERSEDES D-80's O₂-invariance framing. 1 pre-work fork-resolution + 1 done-call advisor pass. Next: D-82 (masking
> readout, owed), O₂-independent thermal fade, tannin self-polymerization, beat 1b, Maillard Strecker, barrel fill-number.**
>
> **PROGRESS (D-82, 2026-07-13): the masking readout — `observed_color_series` built, DELIVERING D-81's COMMITTED beat A
> (906 tests). A pure READOUT (no state slot, no fate): the Somers reversible SO₂/pH bleaching. Free monomeric anthocyanin
> is masked by a coloured fraction `χ = 1/(1 + K_h/h + K·[HSO₃⁻])` (`acidbase.anthocyanin_coloured_fraction`), polymeric
> pigment counted FULL (SO₂/pH-resistant); `observed = χ·anthocyanin·1000 + pigment·1000`. Advisor's load-bearing catch:
> the COMPETITIVE single denominator, NOT a product of two fractions — carbinol + bisulfite adduct are parallel drains of
> the flavylium pool, so a product form carries a spurious cross-term (bisulfite bleaching the colourless carbinol, ~4×
> error at pH 3.4/20 mg/L). Reads FREE bisulfite (`bisulfite_so2_at_ph`, after carbonyl binding), so reversibility is
> emergent (SO₂ bound/oxidised → mask lifts). OPPOSITE SO₂-sign to D-81 (here SO₂ MASKS; D-81 SO₂ PROTECTS) — different
> series, both real, comment cross-refs guard it. Two params (`pKa_flavylium_hydration` plausible, `K_anthocyanin_bisulfite`
> speculative). 1 pre-work + 1 done-call advisor pass. THE COLOUR AXIS'S "Both (C)" REQUEST IS NOW COMPLETE (B=D-81 fade +
> A=D-82 mask). Next: O₂-independent thermal fade, tannin self-polymerization / tannin-ethyl-tannin, beat 1b, Maillard
> Strecker, barrel fill-number, barrel-beer oak.**
> Milestone 1 (Tier-1 validated core) and Milestone 2 (Tier-2
> plausible mechanisms) are closed — the §2.2 benchmark trio is green and §3.3
> "additives with clear mechanisms" completed at D-65 (717 tests). This plan opens
> **Tier-3**, the handoff's §4 frontier: "real chemistry, but integrating it into a
> trustworthy prediction is *not solved science*." Everything here is `speculative`,
> isolated, and clearly labelled — it must never perturb the validated core or its
> tests (prime directive #3).
>
> The two opening calls are recorded in **DECISIONS D-66**: (1) build the
> **sensory/OAV readout layer first**, aging chemistry second — *inverting* the
> handoff §6-step-5 order ("aging then sensory"); (2) handle the lumped aroma pools
> with a **representative-compound threshold per lump** (owner call, over the
> single-compound-only alternative).

## Build order (dependency-ordered; handoff §6 step 5, re-sequenced per D-66)

```
  sensory / OAV readout layer      ← FIRST beat (this milestone's active work)
        │  (pure readout over compounds already tracked; zero core risk)
        │  1a. OAV ratio (sourced thresholds)      ← the honest, sourced part
        └─ 1b. descriptor-space projection         ← deferred: a further heuristic leap
  aging chemistry (the "years" axis) ← subsequent beats, one Process at a time
        ├── ester formation/hydrolysis equilibria over time
        ├── oxidation (acetaldehyde/phenolic browning, Strecker)
        ├── oak extraction (vanillin, whiskey lactones — D-77; ellagitannins — D-78)
        ├── tannin–anthocyanin polymerization (red colour / astringency)   ← D-79 (grape pools, oak-/O₂-independent)
        └── micro-oxygenation / reductive–oxidative sulfide evolution; Maillard/sotolon
```

**Why sensory before aging** (full rationale in D-66): the sensory layer is a **pure
readout** over aroma-active compounds the model *already* tracks (esters, fusels,
diacetyl, acetaldehyde, H₂S, 4-ethylphenol, 4-ethylguaiacol, mercaptans) — so it adds
**no new ODE physics and zero risk to the validated core**. It ships first because it
then becomes the **acceptance lens for aging**: once OAVs exist, every aging Process's
effect on the aroma profile is visible immediately. Aging chemistry is heavier — new
speculative RHS Processes on a years-scale phase (phase-based integration, handoff §7
multi-scale) with scattered parameter sourcing — so it comes second, one Process at a
time behind its own tests. The handoff's "aging then sensory" order is reference, not
gospel (CLAUDE.md); the owner's own framing put sensory first too.

---

## Active beat: the sensory / OAV readout layer (handoff §4.2)

### Placement & the isolation firewall (§4.2 cardinal rule)

A **new top-layer package `fermentation.sensory`**, a sibling of `fermentation.analysis`
in the dependency graph:

```
scenario / validation  →  runtime  →  core  →  parameters / units
                    sensory  ─┘  (consumes Trajectory + thresholds; imported by NOTHING lower)
```

- It consumes a `runtime.Trajectory` (state series) plus a threshold table and returns
  OAV series. It imports the core/runtime *downward* only; **nothing in core/runtime/
  scenario imports it back** — the handoff §4.2 cardinal rule: *the sensory layer
  consumes the chemistry; the chemistry never depends on the sensory layer.*
- **Thresholds load directly into the sensory layer, NOT through the compile seam.**
  Unlike `acidbase.yaml` / `vicinal_diketones.yaml` (merged into every compiled scenario
  at `compile.py`'s `shared_files` because a *Process* reads them), **no RHS reads a
  perception threshold** — so a new `sensory.yaml` is loaded by the sensory module on its
  own, never merged into `CompiledScenario.param_values`. The chemistry never even sees
  the sensory params. This is a stronger isolation than any Tier-2 readout.
- **Tier floor.** Every OAV output tier is `Tier.combine(chemistry_input_tier,
  SPECULATIVE)` → **speculative**, *even when the input chemistry is validated*. The
  sensory mapping is itself speculative (§0.1 / Tier docstring names "sensory mapping"
  as the canonical speculative case), so it caps everything it touches. Read the input
  compound's tier via `ProcessSet.tier_of(pool, ...)` and combine with speculative.

### Definition of done (beat 1a — OAV ratio only)

1. `OAV_i(t) = concentration_i(t) / threshold_i` for each aroma-active compound, mapped
   over a trajectory (mirroring `analysis.ibu_series` / `molecular_so2_series`): an
   `oav_series(traj, thresholds, compound)` per compound plus a finished-profile view over
   the medium's active compounds at a chosen time. Dimensionless. **The aggregate reports
   per-compound OAVs and flags which exceed 1 (above-threshold) — NOT a single summed
   scalar**: summing OAVs assumes perceptual additivity, which is contested, so a summed
   number would over-claim (settle the exact aggregate shape in D-67).
2. A `sensory.yaml` provenance file with **real, sourced perception thresholds** (see
   sourcing below), each carrying value, unit, `source`, `conditions` (**the matrix** —
   see below), `uncertainty`, `tier: speculative`.
3. Unit tests: OAV monotone-increasing in its pool; **identically 0 when the pool is 0**
   (an unhopped/unspoiled/clean run has no false aroma); **tier is the speculative floor**
   even for a validated input; the reported matrix matches the medium. Plus a golden
   sanity check that a known concentration → the literature OAV (e.g. diacetyl at ~2×
   its lager threshold reads OAV ≈ 2).
4. **The Tier-1 suite and all conservation tests are byte-for-byte untouched** — the
   readout adds **no state slot, no RHS, no ledger entry** (prime directive #3, trivially:
   nothing in core changes). `pytest` / `ruff` / `mypy` green.

### Compounds and their thresholds

Aroma-active pools already in state (`core/media.py`), classified by medium and by how
cleanly OAV applies. **The set is medium-specific** (mirroring the beer-only `iso_alpha`):
`_common_specs` (both media) carries only `esters`, `fusels`, `diacetyl`, `acetaldehyde`,
`h2s`; `ethylphenols`, `ethylguaiacols`, and `mercaptans` are appended in `wine_schema`
**only**. So the **beer OAV set = the 5 common compounds**; the **wine OAV set = those 5
+ 4-EP + 4-EG + mercaptans**. A beer `oav_series` must not reach for the three wine-only
slots (they do not exist in the beer schema).

| Pool | Medium | Representative compound (threshold key) | Descriptor | Note |
|------|--------|-----------------------------------------|-----------|------|
| `diacetyl` | wine + beer | 2,3-butanedione | buttery | single molecule — clean OAV |
| `acetaldehyde` | wine + beer | acetaldehyde | green apple / bruised | single molecule — clean OAV |
| `h2s` | wine + beer | hydrogen sulfide | rotten egg | single molecule — clean OAV |
| `esters` | wine + beer | **isoamyl acetate** (stand-in) | banana / fruity | **lumped** → representative threshold |
| `fusels` | wine + beer | **isoamyl alcohol** (3-methylbutan-1-ol, stand-in) | solventy / fusel | **lumped** → representative threshold |
| `ethylphenols` | **wine only** | 4-ethylphenol | horse-sweat / barnyard | single molecule — clean OAV |
| `ethylguaiacols` | **wine only** | 4-ethylguaiacol | clove / smoky | single molecule — clean OAV |
| `mercaptans` | **wine only** | **methanethiol** (stand-in, already the pool's named stand-in) | reductive / drains | **lumped** → representative threshold |

- **The lumped-pool call (D-66, owner-chosen).** `esters`/`fusels`/`mercaptans` are single
  g/L pools that really mix several molecules whose thresholds span ~3 orders of magnitude.
  We assign each lump the threshold of one **named representative compound** — the stand-in
  its `VarSpec` description *already* names (fusels = Ehrlich higher alcohols → isoamyl
  alcohol; mercaptans = "methanethiol stand-in") — compute OAV uniformly, and carry
  **"assumes fixed lump composition"** loudly in that threshold's provenance `notes`. This
  keeps the dominant young-product aromas (esters, fusels) in the numeric readout; the
  honesty cost is the fixed-composition assumption, flagged at the source.
- **`iso_alpha` / IBU is excluded** — it is a **taste** (bitterness), not an odor, and is
  already a direct mg/L→IBU readout (`analysis.ibu_series`, D-64). Do not shoehorn a
  bitterness into an odor-threshold OAV. (A future taste-intensity readout is separate.)

### Units & matrix (both load-bearing)

- **Units:** state is g/L; literature odor thresholds are µg/L–mg/L. Convert at the sensory
  boundary via `fermentation.units` (add g/L↔µg/L helpers if absent). OAV itself is
  dimensionless.
- **Matrix-specificity (a provenance requirement, not optional).** Ethanol and the wine/beer
  matrix shift most odor thresholds substantially — a wine-matrix threshold ≠ beer ≠ water/
  model solution. Every threshold's `conditions` **must record the matrix it was measured
  in**; the wine profile reads wine-matrix thresholds where they exist, beer reads beer-
  matrix, and any water/model-solution fallback is flagged as a matrix gap in `notes`.

### Parameters to source (provenance, like the D-12 sweep)

Perception (odor) thresholds, matrix-specific. Reading list (reconcile, don't transcribe):
Guth 1997 (wine aroma thresholds), Francis & Newton 2005 (wine flavour compounds review),
Meilgaard 1975 / Meilgaard et al. (beer flavour thresholds), Ferreira et al. 2000 (wine
odor-activity), plus diacetyl (~0.1 mg/L lager) and 4-EP/4-EG (~425 / ~110 µg/L red-wine
sensory) from the spoilage literature already cited in `vicinal_diketones.yaml` / the Brett
beats. Tiers: **all `speculative`** — thresholds are panel- and reference-dependent and swing
by orders of magnitude, and the OAV *mapping* is speculative regardless; carry wide
uncertainty bands.

### Approach (test-driven, mirrors the D-64 IBU readout)

1. Create `src/fermentation/sensory/` (package) with an `oav` module; add `sensory.yaml` to
   `parameters/data/` and load it standalone (a `load_thresholds()` helper, NOT via
   `compile.py`'s `shared_files`).
2. Implement `oav_series` per compound + a finished-profile aggregate; each: dimensionless,
   monotone in its pool, 0 when the pool is 0, tier = `combine(input_tier, SPECULATIVE)`.
3. Source + reconcile the thresholds; replace any placeholders; record the matrix in
   `conditions`; flag the lump-composition assumption in the three lumped thresholds.
4. Unit tests (per the DoD) + a golden OAV sanity check. Confirm the §2.2 trio +
   conservation + full suite stay byte-for-byte green (readout touches no core state).
5. Record outcomes in a DECISIONS entry (D-67) and update this plan + ARCHITECTURE.

---

## Deferred / later beats (in order)

- **1b. Descriptor-space projection** — map the OAV vector onto a descriptor vocabulary
  ("fruity/buttery/barnyard/reductive/…"). This is a *further* heuristic leap beyond the
  sourced OAV ratio ("OAVs → this smells like leather and banana"), so it is fenced as a
  **separate, even-more-speculative swappable sub-model** (handoff §4.2: a clean seam so an
  ML model trained on sensory-panel data could later replace it). Keeping it out of beat 1a
  keeps the sourced-ratio layer honest.
- **Aging chemistry (§4.1), one Process at a time on a slow/years phase** — ester
  formation/hydrolysis equilibria; oxidation (acetaldehyde generation, phenolic browning,
  Strecker degradation); **oak extraction** (diffusion-limited vanillin / whiskey lactones /
  ellagitannins as a function of surface-to-volume ratio, toast level, barrel age — new
  extracted pools, akin to dosed inputs); tannin–anthocyanin polymerization (red-wine colour
  and astringency evolution); micro-oxygenation / reductive–oxidative sulfide balance;
  long-aging Maillard / sotolon in oxidative styles. Each: `speculative`, isolable/togglable,
  **phase-based integration** (handoff §7 — do *not* integrate years at fermentation
  resolution), and validated only by the sensory lens built above. The honest framing stays
  "tracking a handful of key compounds with literature-based kinetics," **not** "predicting
  what five years in barrel tastes like."

## Risks

- **Validation is essentially absent at Tier-3** (handoff §4.3). Odor thresholds vary by
  matrix, panel, and reference by orders of magnitude; aging verdicts are not solved science.
  Carry wide uncertainty bands, tag `speculative`, **never** `validated`, and lean on
  directional checks (more diacetyl ⇒ higher buttery OAV) not magnitudes.
- **The §4.3 credibility firewall.** The Tier-3 temptation is to let plausible-looking
  speculation borrow the validated core's credibility. Every OAV / descriptor / aging output
  **must surface its speculative tier** in any plot, export, or report — the tier-floor rule
  above enforces this at the API.
- **Lump-composition assumption** (the accepted D-66 (a) call): a lump's OAV is only as
  meaningful as its assumed fixed composition. Flagged in every lumped threshold's provenance;
  revisit if/when the esters/fusels pools are ever speciated (a *chemistry*-layer change, which
  would be motivated on its own merits, never to serve the sensory layer — §4.2).
- **Aging multi-scale stiffness** (handoff §7): the years phase must use phase-based Process
  activation and an appropriate step regime; do not integrate aging at ferment resolution.
- **Thresholds sit outside the D-24 ensemble sweep** — a deliberate consequence of loading
  `sensory.yaml` standalone rather than through the compile seam: `simulate_ensemble` samples
  only compiled-scenario params, so it will *not* propagate threshold uncertainty into the OAV
  band. Defensible for a speculative readout (the OAV floor is already speculative), but state
  it explicitly in D-67 so it never later reads as an oversight, not a choice.
