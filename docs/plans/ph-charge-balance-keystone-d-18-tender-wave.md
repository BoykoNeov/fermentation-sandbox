# pH charge-balance keystone (D-18)

## Context

Milestone 2's keystone. Today the engine has **no pH** — no pH/pKa/charge symbol
anywhere in `src/`. DECISIONS **D-18** settled *how* pH must be built: a **full
proton/charge-balance solver**, not a tracked-pH approximation. Each weak acid is
carried as state; at a given state we solve `Σ charge = 0` for `[H⁺]` and read
`pH = −log₁₀[H⁺]`. The discriminator is **compositionality**: the two couplings
Tier-2 needs — MLF deacidification (pH rises ~0.1–0.3 as malic→lactic) and SO₂
speciation (molecular fraction set by pKa) — must *emerge* from the balance, not be
scripted per event. This keystone unblocks **SO₂ → MLF → mixed cultures**.

pH is the **derived algebraic pure function** of state, exactly like
`total_carbon`/ABV: there is no `dpH/dt`, the core stays pure. For D-18 there is
**no RHS consumer yet** (SO₂/MLF wire pH into rates in later beats), so the
deliverable is the *solver + a post-hoc pH/TA readout*. We still write the solver
with the in-loop signature `(y, schema, params)` and keep it medium-agnostic so the
future `RateModifier` hook and beer's phosphate set reuse it for free.

### Decisions baked in (owner-confirmed)

1. **Wine-only acid state.** Beer pH is a phosphate-buffered different acid system
   with no sourced data yet — explicitly deferred. All D-18 acids are wine acids.
2. **A strong-cation term is mandatory, not optional.** Weak acids alone give
   pH ≈ 2.3 at must tartaric levels (~33 mM, pKa1≈3.04); real must is ~3.3 — K⁺ as
   bitartrate supplies the counter-charge. Without it the solver is *qualitatively*
   wrong.
3. **Anchoring = inverse (now).** Scenario gives acid concentrations + a measured
   `initial_ph`; compile **back-solves the strong-cation charge** to reproduce it,
   then stores it as a constant state slot; pH evolves emergently as acids change.
   Honest claim: **D-18 predicts pH *changes*, not absolute initial pH** (initial
   pH is an input). This absorbs activity-coefficient & cation uncertainty into one
   fitted term (how Boulton's wine-pH model is anchored). *Forward-from-cation is a
   documented future option* — the core solver is anchoring-agnostic and the cation
   stays a state slot, so adding a forward `cation_meq_l` input later is additive.
4. **Byp = include-by-reading (coupling #2 closed).** The charge balance reads the
   *existing* `Byp` pool as its succinic-equivalent — **zero new carbon**, so
   `total_carbon` is unchanged and the double-count is closed, not deferred. Caveat:
   `Byp` lumps neutral 2,3-butanediol, slightly overstating acid charge (~1–1.5 mM
   vs ~20 mM buffer — minor).

### Scope caveats (documented, with numbers — *justified* scope, not hand-waves)

- **Carbonic omitted (coupling #1).** At pH 3.3, bicarbonate charge ~0.03 mM vs a
  ~20 mM buffer (~0.1%). Correct to omit below pH ~4; `CO2` state stays the evolved
  proxy (D-15). Revisit threshold: deacidified/low-acid musts above pH ~4.
- **Constant pKa (coupling #3).** Carboxylic ΔH_ionization ≈ 0; the pKa shift over
  10–30 °C is <0.05 units, inside the pKa uncertainty. (We omit carbonic, which is
  the one acid with real T-dependence.)
- **Ionic strength / activity (fourth caveat).** Wine I ≈ 0.05–0.1 M; rigorous
  balance needs Davies/Debye-Hückel. Concentration-based *apparent* pKa is the
  standard plausible-tier simplification; inverse anchoring folds the activity error
  into the fitted cation at t=0, leaving it to affect only the *slope* (buffer
  capacity), where we claim only directional fidelity.
- **pH tier = `plausible`.** The CRC pKa measurement is validated, but applying
  25 °C / I=0 constants to wine is extrapolation. pH's tier is computed **explicitly**
  (see step 6) — it must *not* inherit the `VALIDATED` default that `tier_of` returns
  for the inert acid slots no Process touches.

## Implementation

One-directional deps respected: pure solver in **core**; trajectory-level series in a
new top-layer **analysis** module; pKa via the **parameters** store; scenario inputs
at the **compile** seam.

### 1. State: add wine-only acid + cation slots — `src/fermentation/core/media.py`

`_common_specs` is shared by both media (media.py:67), so do **not** touch it. Append
the new slots in `wine_schema()` (media.py:113) only, so `beer_schema()` is untouched:

```python
def wine_schema() -> StateSchema:
    specs = _common_specs(VarSpec("S", "g/L", description="fermentable sugar"))
    specs += [
        VarSpec("tartaric", "g/L", default=0.0, description="tartaric acid (must input; diprotic)"),
        VarSpec("malic",    "g/L", default=0.0, description="L-malic acid (must input; diprotic; MLF substrate)"),
        VarSpec("lactic",   "g/L", default=0.0, description="L-lactic acid (produced-only; MLF product)"),
        VarSpec("cation_charge", "mol/L", default=0.0,
                description="net strong-cation charge (K+-dominant), constant; back-solved from initial_ph (D-18)"),
    ]
    return StateSchema(specs)
```

- **`default=0.0` is load-bearing**: existing wine scenarios/tests that don't name
  acids still compile (all four → 0). With acids and cation at 0 the slots are inert
  and contribute 0 to every conservation sum, so the validated core and its tests are
  untouched (prime directive #3). pH is simply not meaningful for a no-acid scenario
  and is only *computed* when requested.
- No Process touches these slots in D-18 → their derivatives are 0 → constant
  trajectory. They exist so MLF can later deplete `malic`/grow `lactic`, and so the
  charge balance and `total_carbon` can read them.
- `cation_charge` is a charge density (mol⁺/L), not a mass concentration — state is
  already heterogeneous (`T` in K), so this is consistent. Keeping it in state (not a
  param) is right: it's scenario-specific like sugar-from-Brix, and a future
  cold-stabilization/racking Process could deplete it (KHT precipitation).

### 2. The solver — new `src/fermentation/core/acidbase.py`

Pure, no I/O, no state. `scipy.optimize.brentq` (scipy already a dep via `solve_ivp`).
Concrete surface:

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class AcidSpec:
    molar_mass: float            # g/mol, from chemistry.py
    pka_param_names: tuple[str, ...]   # 1 entry monoprotic, 2 diprotic

# Which state slots are charge-active acids, + how Byp is read.
ACID_STATE: dict[str, AcidSpec]   # {"tartaric": .., "malic": .., "lactic": ..}
BYP_AS_SUCCINIC: AcidSpec          # read the existing Byp slot as succinic-equiv

KW = 1.0e-14                       # in-code physical constant (per D-3; like GAS_CONSTANT)

def mean_charge(h: float, pkas: tuple[float, ...]) -> float:
    """Mean anion charge magnitude per mole at proton conc h (Henderson-Hasselbalch:
    monoprotic α1 = Ka/(Ka+h); diprotic via D = h^2 + Ka1 h + Ka1 Ka2)."""

def charge_residual(ph: float, totals_molar, cation, byp_succinic_molar, pka_map) -> float:
    """Net charge as a function of pH: (cation + [H+]) − ([OH-] + Σ acid-anion charge).
    Monotonically decreasing in pH ⇒ a single smooth root. Solve in pH-space (smooth
    variable) not [H+]-space."""

def solve_ph(totals_molar, cation, byp_succinic_molar, pka_map) -> float:
    """brentq(charge_residual, 0.0, 14.0) — wide bracket guarantees a sign change."""

def solve_cation_charge(totals_molar, byp_succinic_molar, pka_map, target_ph) -> float:
    """Inverse anchoring — CLOSED FORM (no root-find): at target pH,
    cation = Σ acid-anion charge + Kw/[H+] − [H+]. Raises ValueError if < 0
    (initial_ph inconsistent with the acid load — unphysical)."""

def ph_of_state(y, schema, params: Mapping[str, float]) -> float:
    """In-loop signature, pure — `params` is the RESOLVED `{name: float}` map that
    `Process.derivatives`/`RateModifier.factor` receive (NOT a ParameterSet), so the
    future in-loop hook is genuinely free. Reads acid g/L from y (only the ACID_STATE
    slots that exist in schema), converts to mol/L via chemistry molar masses, reads
    cation from the cation_charge slot and Byp as succinic-equiv, pulls pKa from params,
    returns pH."""

def titratable_acidity(y, schema, params: Mapping[str, float]) -> float:
    """TA as g/L tartaric-equivalent — a SECOND derived pure function of the same acid
    state (no new state). Σ C_a·(z_max − mean_charge(sample pH)) in eq/L, × M_tartaric/2
    (≈75.04 g/eq). Approximates the pH-8.2 endpoint as full carboxylic deprotonation and
    omits the free-[H⁺] term (~0.4 mM, ~0.4% at must pH) — note it; conventional TA-to-8.2
    includes it. Fine for the plausible tier and the 6–9 g/L band."""

def ph_tier(params: ParameterSet) -> Tier:
    """combine(all pKa param tiers, Tier.PLAUSIBLE) — computed EXPLICITLY, takes the full
    ParameterSet (kept separate from the resolved-map hot-loop signature above). Does NOT
    read tier_of the inert acid slots (that returns VALIDATED via empty combine)."""
```

- **Smoothness:** `charge_residual` is a sum of C∞ Henderson-Hasselbalch terms →
  smooth root in the acid totals. The C¹ check (step 7) guards the *future* in-loop
  BDF use; there is no consumer in D-18. If `brentq`-per-RHS ever bites later, the
  closed-form `d(charge)/d[H⁺]` → Newton in 3–4 iters is the perf escape hatch — note
  it, don't build it.
- **Medium-agnostic:** the solver iterates only the `ACID_STATE` slots present in the
  passed schema, so beer's future phosphate set drops in by extending `ACID_STATE`.

### 3. Acid stoichiometry — `src/fermentation/core/chemistry.py`

Add pure formula constants (per D-3 these live in code, like `M_SUCCINIC`), then
register in `MOLAR_MASS` and `CARBON_ATOMS` (chemistry.py:71-96):

```python
M_TARTARIC = 4*_M_C + 6*_M_H + 6*_M_O   # C4H6O6 ≈ 150.086
M_MALIC    = 4*_M_C + 6*_M_H + 5*_M_O   # C4H6O5 ≈ 134.087
M_LACTIC   = 3*_M_C + 6*_M_H + 3*_M_O   # C3H6O3 ≈  90.078
# MOLAR_MASS += {"tartaric_acid": M_TARTARIC, "malic_acid": M_MALIC, "lactic_acid": M_LACTIC}
# CARBON_ATOMS += {"tartaric_acid": 4, "malic_acid": 4, "lactic_acid": 3}
```

**MLF-ready carbon:** malic (C₄) → lactic (C₃) + CO₂ (C₁) balances (4 = 3+1), so these
weights make a future MLF Process carbon-closing. `succinic_acid` already exists for
the Byp read.

### 4. Carbon weights — `src/fermentation/validation/conservation.py`

In `total_carbon` (conservation.py:40), mirror the existing `Gly`/`Byp` block (:69-72):

```python
if "tartaric" in schema: w[schema.slice("tartaric")] = carbon_mass_fraction("tartaric_acid")
if "malic"    in schema: w[schema.slice("malic")]    = carbon_mass_fraction("malic_acid")
if "lactic"   in schema: w[schema.slice("lactic")]   = carbon_mass_fraction("lactic_acid")
```

`cation_charge` gets weight 0 automatically (`schema.zeros()` + not a carbon species);
`total_nitrogen`/`total_mass` need no change (cation/acids are neither N nor in
{S,E,CO2}). `Byp` weighting is **unchanged** — we only *read* it for charge, add no
carbon. In D-18 the acids are constant ⇒ drift 0 ⇒ existing conservation tests pass;
the weights are there for MLF.

### 5. Parameters: pKa set + compile-seam inputs

**New `src/fermentation/parameters/data/acidbase.yaml`** — shared, medium-agnostic,
loaded alongside the medium file (the store merges paths; names don't collide). pKa
set (sanity anchors): `pKa_tartaric_1` 3.04 / `pKa_tartaric_2` 4.37; `pKa_malic_1`
3.40 / `pKa_malic_2` 5.11; `pKa_lactic` 3.86; `pKa_succinic_1` 4.21 /
`pKa_succinic_2` 5.64. Each a full `Parameter` (value, unit `"pH units"`, tier
`plausible`, uncertainty `{low, high, note}`, provenance `{source: CRC/textbook,
conditions: "25 °C, I≈0", doi?, notes}`). `KW` stays an in-code constant. (Acetic pKa
4.76 is folded into the succinic-equivalent Byp read for now; add a separate `acetic`
species only if Byp is later speciated.)

**Load seam — `_load_parameters` (compile.py:195):** on the default-lookup branch,
merge the shared file:

```python
base = Path(data_dir) if data_dir is not None else default_data_dir()
medium_path = base / f"{scenario.medium}_{scenario.strain}.yaml"
shared = base / "acidbase.yaml"
return load_parameters(medium_path, shared)   # merge; pKa names are collision-free
```

The explicit `parameter_paths` override branch (compile.py:200) stays caller-controlled
— document that callers wanting pH must include `acidbase.yaml` in their paths.

**Compile inputs — `_ALLOWED_KEYS["wine"]` (compile.py:86)** gains optional keys:
`tartaric_gpl`, `malic_gpl`, `initial_ph` (`lactic` is produced-only → not an input;
cation is back-solved). All optional so existing scenarios still validate.

**`_wine_initial` (compile.py:109)** packs the acid slots and back-solves the cation
(only when acids + `initial_ph` are supplied; otherwise leave all four at 0):

```python
tartaric = _optional(values, "tartaric_gpl", 0.0)
malic    = _optional(values, "malic_gpl", 0.0)
# ... pack "tartaric"/"malic"/"lactic": 0.0 ...
if "initial_ph" in values:
    cation = acidbase.solve_cation_charge(  # mol/L; converts g/L→mol/L internally
        totals_molar=..., byp_succinic_molar=0.0,  # Byp=0 at pitch
        pka_map=_pka_map(parameters), target_ph=float(values["initial_ph"]),
    )  # raises ValueError (caught/re-raised as a clear compile error) if unphysical (<0)
    out["cation_charge"] = cation
```

`solve_cation_charge`'s negative-cation guard surfaces an `initial_ph` inconsistent
with the acid load at compile, rather than packing a nonsense state.

### 6. Derived pH/TA exposure — new `src/fermentation/analysis.py` (top layer)

`ph_of_state`/`titratable_acidity` are pure and live in core. The *series* helpers need
`Trajectory` (a runtime type), so they go one layer up — a thin new module mirroring how
`units` provides scalar conversions and how benchmarks map ABV over `traj.series("E")`:

```python
from fermentation.core import acidbase
from fermentation.runtime.integrate import Trajectory

def ph_series(traj: Trajectory, params) -> FloatArray: ...        # ph_of_state over columns
def titratable_acidity_series(traj, params) -> FloatArray: ...
# tier reported via acidbase.ph_tier(params)
```

In D-18 the **new** acid slots (tartaric/malic/lactic) are constant, but the pH series
is **not flat** — and that is the point. `Byp` is core (the D-16 realised-yield diversion
in `SugarUptakeToEthanolCO2` grows it 0 → ~1–4 g/L over the ferment, even with the aroma
tuple off), and the include-by-reading decision makes that growing succinate anion charge
count. With the cation frozen at its pitch value (back-solved at Byp=0), pH therefore
**drifts mildly down** (~0.05–0.1 units for ~2 g/L succinic-equiv at must pH) as
fermentation proceeds. This is a realistic, *emergent* second demonstration — beyond the
malic→lactic headline — that the solver responds to acid dynamics with no scripting.
(Verify the magnitude against the actual realised-yield numbers during implementation.)

### 7. Tests — `tests/` (`test_acidbase.py` + extend conservation/compile tests)

Ranked, headline first:

1. **HEADLINE — malic→lactic ΔpH ∈ [0.1, 0.3]** (proves MLF-enablement *without* MLF
   built): use a **malic-rich** composition (the case where MLF actually matters — ΔpH
   at fixed cation is composition-dependent, so a tartaric-heavy/low-malic must can land
   below 0.1) + anchored cation; substitute X mM malic → X mM lactic 1:1; recompute pH;
   assert the rise lands in [0.1, 0.3]. D-18's proof-of-purpose. **If the first numbers
   give e.g. 0.08, fix by choosing a more malic-rich must — NOT by widening the band**
   (CLAUDE.md forbids weakening benchmark tests).
2. **Charge residual ≈ 0** at the solved `[H⁺]` (the balance actually balances).
3. **Monotonicity** — more tartaric ⇒ lower pH; more cation ⇒ higher pH.
4. **Smoothness / C¹** — finite-difference `dpH/d(acid)` continuous (guards future BDF).
5. **Round-trip** — a wine scenario with `initial_ph` compiles and `ph_of_state(y0)`
   reproduces it to tolerance (inverse anchoring is exact at t=0). *Note this test is
   tautological w.r.t. unit conversion* — `solve_cation_charge`/`solve_ph` are inverses
   applying the same g/L→mol/L factor, so a conversion bug cancels. It does NOT guard
   absolute correctness; test 6 does.
6. **Back-solved cation is physical (the units guard)** — for a textbook must the
   back-solved `cation_charge` lands in a physical K⁺ range **~25–50 meq/L** (K⁺ ~1–2 g/L
   ÷ 39.1 g/mol). This is the clean catch for a g/L↔mol/L factor error that the round-trip
   cannot see and TA-in-band only partly does.
7. **Compile guard** — `initial_ph` inconsistent with acid amounts (→ negative cation)
   raises at compile.
8. **Carbon conservation unchanged** — run an existing wine conservation path with the
   new acid slots (zero, then nonzero-constant); drift stays at tolerance. (Byp still
   grows and is still weighted as before — unaffected.)
9. **TA consistency** — computed TA for a representative must lands ~6–9 g/L
   tartaric-equiv.
10. **Tier** — `ph_tier(params)` / TA tier resolves to `plausible`, never `validated`.

### 8. Docs

- **`docs/DECISIONS.md` D-18**: flip status open→resolved; record the cation term
  (with the 2.3-vs-3.3 number), inverse anchoring + the "predicts pH *changes*" claim,
  the four resolved caveats with their numbers, and the Byp include-by-reading
  resolution of coupling #2. Note forward-anchoring as a future option.
- **`docs/plans/milestone-2-tasks.md`**: check off the pH-beat boxes (:118-125);
  flag the malic→lactic ΔpH acceptance test as the gate.
- **`docs/ARCHITECTURE.md`**: add a short pH/acid-state subsection (pH as a derived
  pure function alongside `total_carbon`; wine acid + cation slots; the new
  `analysis.py` observable layer).

## Verification

```bash
uv run pytest -q                         # all green incl. new tests/test_acidbase.py
uv run pytest tests/test_acidbase.py -q  # the new solver suite, headline test first
uv run ruff check . && uv run ruff format --check .
uv run mypy                              # strict on src — new acidbase.py/analysis.py typed
```

End-to-end smoke (manual): compile a wine scenario with `tartaric_gpl`, `malic_gpl`,
`initial_ph: 3.4`; `simulate`; map `analysis.ph_series` over the trajectory and confirm
it starts at 3.4 (anchoring exact at t=0) and **drifts mildly down** (~0.05–0.1) as `Byp`
accumulates — *not* flat (see step 6); confirm `titratable_acidity` lands in 6–9 g/L. Then
in a unit test do the malic→lactic substitution on a malic-rich must and confirm pH rises
0.1–0.3 — the keystone demonstrably enables what it exists for.

**"Done" boundary:** solver + post-hoc pH/TA readout, no RHS consumer. SO₂ and MLF
(which wire pH into rates) are the next beats.
