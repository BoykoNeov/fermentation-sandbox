"""Wine pH as a charge-balance solution — the Milestone-2 keystone (decision D-18).

pH is a **derived algebraic pure function of state**, exactly like ``total_carbon``
or ABV: there is no ``dpH/dt`` and the core stays pure. Given the charge-active acid
concentrations carried in the state vector, we solve the electroneutrality condition

    Σ charge = 0   ⇔   (cation + [H⁺]) − ([OH⁻] + Σ acid-anion charge) = 0

for ``[H⁺]`` and report ``pH = −log₁₀[H⁺]``. Building pH as a *full proton balance*
(rather than a tracked-pH approximation) is the decision's discriminator: the two
couplings Tier-2 needs — MLF deacidification (pH rises as diprotic malic → monoprotic
lactic) and SO₂ speciation (molecular fraction set by pKa) — *emerge* from the balance
rather than being scripted per event. This module is the solver + post-hoc readout;
there is no RHS consumer in D-18 (SO₂/MLF wire pH into rates in later beats).

**SO₂ speciation (decisions D-22, D-28)** is the first such readout to land. In D-22 the
dosed slot was *free* SO₂ and :func:`molecular_so2` partitioned it into its molecular
(antimicrobial) fraction *at the solved pH*, delivering the keystone's promised "dose SO₂
→ speciation falls out of the current pH" coupling. **D-28** completes it once acetaldehyde
became real state (D-27): the slot is reinterpreted as **total** SO₂ (conserved, inert),
and :func:`speciate_so2` derives the free/bound split by the acetaldehyde-bisulfite binding
equilibrium — ``free = total − bound``, ``molecular = free × fraction(pH)`` — so the early
acetaldehyde peak *emergently* sequesters SO₂ and depresses free/molecular, which recovers
as acetaldehyde is reduced. It remains **readout-only** — total SO₂ is a state slot but is
*not* in the charge balance (its minor bisulfite charge would be absorbed by the inverse-
anchored cation at t=0 regardless; see D-22), so the D-18 solver signatures are untouched,
and the split feeds back into the acetaldehyde reduction (**decision D-47**): bound acetaldehyde
is protected from alcohol dehydrogenase, so :func:`free_acetaldehyde` gates the reduction onto the
unbound share and SO₂ **locks in** acetaldehyde — the D-22/D-28 "SO₂ is readout-only" invariant is
intentionally retired for sulfited runs (an *un*-dosed run is still byte-for-byte the D-27 core, and
no §2.2 benchmark doses SO₂). At acetaldehyde = 0 the split collapses to D-22 exactly (``free ==
total``). The one live consumer is the pre-existing MLF
antimicrobial gate, which now reads the *derived* free-molecular SO₂ via
:func:`molecular_so2_at_ph` (bound SO₂ is not antimicrobial — a correct consequence, D-28).

**Multi-carbonyl competition (decision D-51).** Acetaldehyde is not the only carbonyl that binds
SO₂ in real wine: the excreted keto-acid overflow pools pyruvate (D-49) and α-ketoglutarate (D-50)
persist in finished wine and compete for the *same* shared bisulfite pool (Burroughs & Sparks 1973;
Jackowetz & Mira de Orduña 2013). :func:`bound_so2_molar` generalizes the D-28 single-carbonyl
quadratic to a competitive-Langmuir partition over all three carbonyls at once — one shared
"reactive bisulfite" root, solved by ``brentq`` — which reduces *exactly* to the original D-28
quadratic when only one carbonyl is present (byte-for-byte on a keto-acid-pool-off run). The
practical effect: acetaldehyde's bound (and hence ADH-protected) share is now *smaller* than the
D-28/D-47 single-carbonyl form, because pyruvate/α-KG soak up part of the dosed SO₂ — the mechanism
D-48's finished-wine overshoot flagged as missing.

**Acid speciation (Henderson-Hasselbalch).** Each weak acid's mean anion charge per
mole at proton concentration ``h`` is a smooth (C∞) function of ``h``:

  * monoprotic (one pKa): ``α₁ = Ka/(Ka + h)`` → mean charge ``α₁`` (one −1 site);
  * diprotic (two pKas): with ``D = h² + Ka₁·h + Ka₁·Ka₂`` the species fractions are
    ``[H₂A] = h²/D``, ``[HA⁻] = Ka₁·h/D``, ``[A²⁻] = Ka₁·Ka₂/D``; mean charge magnitude
    ``(1·Ka₁·h + 2·Ka₁·Ka₂)/D``.

Because every term is a smooth, monotone function of pH and the cation/H⁺ terms are
too, :func:`charge_residual` is monotonically decreasing in pH ⇒ a single smooth root;
we solve in pH-space (the well-conditioned variable) with ``brentq`` over ``[0, 14]``,
where the wide bracket guarantees a sign change.

**The strong-cation term is mandatory (D-18).** Weak acids alone give pH ≈ 2.3 at must
tartaric levels (~33 mM, pKa₁ ≈ 3.04); real must is ~3.3. Potassium (as bitartrate)
supplies the counter-charge. We carry it as a net strong-cation charge density
(``cation_charge`` state slot, mol⁺/L) and **back-solve it from a measured initial pH**
at compile (:func:`solve_cation_charge`): the honest claim is that D-18 predicts pH
*changes*, not absolute initial pH (initial pH is an input). Inverse anchoring absorbs
activity-coefficient and cation uncertainty into one fitted term — how Boulton's
wine-pH model is anchored — and leaves the concentration-based *apparent* pKa
simplification to affect only the slope (buffer capacity), where we claim only
directional fidelity. Tier: **plausible** (CRC pKa values are measured, but applying
25 °C / I=0 constants to wine is extrapolation).

**Medium-agnostic.** The solver iterates only the :data:`ACID_STATE` slots present in
the passed schema, so beer's future phosphate set drops in by extending that mapping;
the in-loop signature ``(y, schema, params)`` matches ``Process.derivatives`` /
``RateModifier.factor`` so a future in-loop pH→rate hook is free.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from scipy.optimize import brentq

from fermentation.core.chemistry import (
    M_ACETALDEHYDE,
    M_ALPHA_KETOGLUTARATE,
    M_LACTIC,
    M_MALIC,
    M_PYRUVATE,
    M_SO2,
    M_SUCCINIC,
    M_TARTARIC,
)
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier, combine

#: Ionic product of water [mol²/L²] at 25 °C. A universal physical constant, so it
#: lives in code with its citation (like ``GAS_CONSTANT`` in the Arrhenius modifier,
#: decision D-3) rather than the empirical-parameter store. At wine pH the OH⁻ term is
#: negligible (~1e-11 M) but it is carried so the balance is exact across all pH and
#: ``solve_ph`` brackets cleanly up to pH 14.
KW = 1.0e-14


@dataclass(frozen=True)
class AcidSpec:
    """How one charge-active acid maps onto chemistry constants and pKa parameters.

    ``molar_mass`` converts the state slot's g/L to mol/L; ``pka_param_names`` lists
    the parameter-store names of its pKa(s) — one entry for a monoprotic acid, two
    (in ascending order pKa₁ < pKa₂) for a diprotic one. ``protons`` (= number of
    pKas) is the maximum anion charge, used by :func:`titratable_acidity`.
    """

    molar_mass: float
    pka_param_names: tuple[str, ...]

    @property
    def protons(self) -> int:
        return len(self.pka_param_names)


#: The charge-active acids carried as wine state slots (decision D-18, wine-only).
#: Beer's phosphate-buffered system is a different acid set with no sourced data yet
#: and is explicitly deferred — extend this mapping when it lands. ``lactic`` is
#: produced-only (the MLF product), but it is charge-active the moment it exists.
ACID_STATE: dict[str, AcidSpec] = {
    "tartaric": AcidSpec(M_TARTARIC, ("pKa_tartaric_1", "pKa_tartaric_2")),
    "malic": AcidSpec(M_MALIC, ("pKa_malic_1", "pKa_malic_2")),
    "lactic": AcidSpec(M_LACTIC, ("pKa_lactic",)),
}

#: The existing ``Byp`` minor-byproduct pool is *read* as a succinic-equivalent diprotic
#: acid for the charge balance (decision D-18, coupling #2 "include-by-reading"): zero
#: new carbon, so ``total_carbon`` is unchanged and the former double-count is closed.
#: Caveat: ``Byp`` lumps neutral 2,3-butanediol too, slightly overstating acid charge
#: (~1–1.5 mM vs a ~20 mM buffer — minor).
BYP_KEY = "Byp"
BYP_AS_SUCCINIC = AcidSpec(M_SUCCINIC, ("pKa_succinic_1", "pKa_succinic_2"))

#: Every pKa parameter name the solver may read — the single list :func:`ph_tier`
#: and the compile/analysis ``build_pka_map`` adapters iterate, so the sites cannot drift.
PKA_PARAM_NAMES: tuple[str, ...] = tuple(
    name for spec in (*ACID_STATE.values(), BYP_AS_SUCCINIC) for name in spec.pka_param_names
)

#: **Total** SO₂ — a wine-only state slot (g/L of SO₂-equivalent), dosed at pitch and
#: **inert** (no Process touches it ⇒ conserved, exactly like the D-18 acids). Read ONLY
#: by the SO₂-speciation readouts below, never by the pH solver or :func:`titratable_acidity`
#: (the readout-only design, decisions D-22/D-28). In D-22 this slot *was* free SO₂; once
#: acetaldehyde became real state (D-27) it was reinterpreted as total, and free/bound are
#: derived by the binding equilibrium (decision D-28) — at acetaldehyde = 0, free == total,
#: so D-22's molecular-SO₂ curve is recovered exactly. NB "total" here is really "free +
#: acetaldehyde-bound": other carbonyl binders (pyruvate, α-ketoglutarate, sugars) are not
#: modelled, so this slightly *under*-binds — a documented v1 overclaim on the name (D-28).
SO2_STATE_KEY = "so2_total"

#: The acetaldehyde state slot (g/L) the free/bound SO₂ split reads as the SO₂ binder
#: (decision D-28). Built as real state in D-27; absent ⇒ no binding ⇒ free == total.
ACETALDEHYDE_KEY = "acetaldehyde"

#: The excreted keto-acid state slots (g/L) that compete with acetaldehyde for the shared
#: bisulfite pool (decision D-51). Built as real state in D-49/D-50; wine-only and absent
#: whenever the ``_KETO_ACID_PROCESSES`` tuple is off (or on the beer schema), in which case
#: they contribute zero binding and the equilibrium collapses to the D-28 acetaldehyde-only form.
PYRUVATE_KEY = "pyruvate"
ALPHA_KG_KEY = "alpha_ketoglutarate"

#: The acetaldehyde-bisulfite adduct dissociation constant (mol/L), referenced to bisulfite
#: HSO₃⁻ (decision D-28). Read only by the binding solver, never by the pH balance.
SO2_BINDING_PARAM = "K_acetaldehyde_so2"

#: The pyruvate- and α-ketoglutarate-bisulfite adduct dissociation constants (mol/L), the two
#: competing carbonyl binders D-51 adds to the same equilibrium (Burroughs & Sparks 1973).
PYRUVATE_SO2_BINDING_PARAM = "K_pyruvate_so2"
ALPHA_KG_SO2_BINDING_PARAM = "K_alpha_kg_so2"

#: The sulfurous-acid pKa parameter names the molecular-SO₂ readout reads — DELIBERATELY
#: kept out of :data:`PKA_PARAM_NAMES`. SO₂ does not enter the charge balance in D-22, so
#: ``build_pka_map`` / :func:`charge_residual` never see these; pH is solved from the
#: organic acids and the dosed free SO₂ is then partitioned at that pH (readout-only).
SO2_PKA_PARAM_NAMES: tuple[str, ...] = ("pKa_sulfurous_1", "pKa_sulfurous_2")


# -- speciation ---------------------------------------------------------------


def mean_charge(h: float, pkas: tuple[float, ...]) -> float:
    """Mean anion charge magnitude per mole of an acid at proton concentration ``h``.

    Henderson-Hasselbalch: ``pkas`` are the pKa value(s) in ascending order — one for a
    monoprotic acid (``α₁ = Ka/(Ka + h)``), two for a diprotic one (via the partition
    ``D = h² + Ka₁·h + Ka₁·Ka₂``). Returns a value in ``[0, len(pkas)]`` that rises
    smoothly as ``h`` falls (higher pH ⇒ more dissociated ⇒ more negative charge).
    """
    if len(pkas) == 1:
        ka = 10.0 ** (-pkas[0])
        return float(ka / (ka + h))
    if len(pkas) == 2:
        ka1 = 10.0 ** (-pkas[0])
        ka2 = 10.0 ** (-pkas[1])
        denom = h * h + ka1 * h + ka1 * ka2
        return float((ka1 * h + 2.0 * ka1 * ka2) / denom)
    raise ValueError(f"only mono- and diprotic acids supported, got {len(pkas)} pKa(s)")


def neutral_fraction(h: float, pkas: tuple[float, ...]) -> float:
    """Fraction of an acid present as its fully-protonated (neutral) species at ``h``.

    The complement of :func:`mean_charge`'s dissociation: the undissociated H₂A (or HA)
    share. Monoprotic ``h/(h + Ka)``; diprotic ``h²/D`` with the same partition
    ``D = h² + Ka₁·h + Ka₁·Ka₂``. Used by the molecular-SO₂ readout (decision D-22): the
    *molecular* (antimicrobial) fraction of free SO₂ is exactly the undissociated
    SO₂·H₂O share. Rises smoothly toward 1 as ``h`` rises (lower pH ⇒ more protonated).
    """
    if len(pkas) == 1:
        ka = 10.0 ** (-pkas[0])
        return float(h / (h + ka))
    if len(pkas) == 2:
        ka1 = 10.0 ** (-pkas[0])
        ka2 = 10.0 ** (-pkas[1])
        denom = h * h + ka1 * h + ka1 * ka2
        return float(h * h / denom)
    raise ValueError(f"only mono- and diprotic acids supported, got {len(pkas)} pKa(s)")


def bisulfite_fraction(h: float, pkas: tuple[float, ...]) -> float:
    """Fraction of an acid present as its **singly-ionized** species (HA⁻) at ``h``.

    The middle species between :func:`neutral_fraction` (H₂A) and the fully-dissociated
    form: monoprotic ``Ka/(Ka + h)`` (= its only anion), diprotic ``Ka₁·h/D`` with the
    same partition ``D = h² + Ka₁·h + Ka₁·Ka₂``. Used by the free/bound-SO₂ split (decision
    D-28): the bisulfite HSO₃⁻ share of free SO₂ is the reactive nucleophile that binds
    acetaldehyde, so ``[HSO₃⁻] = free_SO₂ · bisulfite_fraction(pH)``. At wine pH (well above
    the sulfurous pKa₁ ≈ 1.81 and below pKa₂ ≈ 7.2) this is ~0.94–0.99 — nearly all free SO₂
    is bisulfite — so it varies only mildly across the wine range.
    """
    if len(pkas) == 1:
        ka = 10.0 ** (-pkas[0])
        return float(ka / (ka + h))
    if len(pkas) == 2:
        ka1 = 10.0 ** (-pkas[0])
        ka2 = 10.0 ** (-pkas[1])
        denom = h * h + ka1 * h + ka1 * ka2
        return float(ka1 * h / denom)
    raise ValueError(f"only mono- and diprotic acids supported, got {len(pkas)} pKa(s)")


# -- the charge balance -------------------------------------------------------


def charge_residual(
    ph: float,
    totals_molar: Mapping[str, float],
    cation: float,
    byp_succinic_molar: float,
    pka_map: Mapping[str, tuple[float, ...]],
) -> float:
    """Net charge [eq/L] as a function of pH — zero at electroneutrality.

    ``(cation + [H⁺]) − ([OH⁻] + Σ acid-anion charge)``. ``totals_molar`` maps each acid
    slot name to its mol/L; ``pka_map`` maps the same names (plus :data:`BYP_KEY`) to
    pKa tuples; ``byp_succinic_molar`` is the ``Byp`` pool read as succinic-equivalent.
    Monotonically decreasing in pH (cation/H⁺ fall, anion charge rises) ⇒ a single
    smooth root.
    """
    h = 10.0 ** (-ph)
    oh = KW / h
    anion = byp_succinic_molar * mean_charge(h, pka_map[BYP_KEY])
    for name, conc in totals_molar.items():
        anion += conc * mean_charge(h, pka_map[name])
    return float((cation + h) - (oh + anion))


def solve_ph(
    totals_molar: Mapping[str, float],
    cation: float,
    byp_succinic_molar: float,
    pka_map: Mapping[str, tuple[float, ...]],
) -> float:
    """Solve ``charge_residual = 0`` for pH, clamped to the physical window ``[0, 14]``.

    ``charge_residual`` is strictly monotone *decreasing* in pH (cation/[H⁺] fall, anion
    charge rises), so residual(0) is its maximum and residual(14) its minimum. For a
    physiological cation the curve crosses zero exactly once inside the bracket and
    ``brentq`` finds it. This function is *also* called with a NON-physiological cation:
    BDF's ``num_jac`` perturbs the ``cation_charge`` state slot far outside its ~0.03 mol/L
    range while probing the Jacobian, which can push the whole curve positive (cation
    swamps all buffering ⇒ the electroneutral pH lies *above* 14) or negative (an acid load
    with no counter-cation ⇒ *below* 0). Returning the boundary is then the exact answer —
    the root lies outside the physical window — and keeps ``ph_of_state`` a TOTAL, bounded
    function of state so the Jacobian probe cannot raise (D-46). The real trajectory never
    leaves the bracket (RK45/LSODA hold ``cation_charge`` constant), so every physiological
    call falls through to the identical ``brentq`` — bit-for-bit pH, byte-for-byte curves.
    """
    args = (totals_molar, cation, byp_succinic_molar, pka_map)
    if charge_residual(0.0, *args) <= 0.0:
        return 0.0  # net-negative even fully protonated ⇒ electroneutral pH ≤ 0
    if charge_residual(14.0, *args) >= 0.0:
        return 14.0  # still net-positive at pH 14 ⇒ electroneutral pH ≥ 14
    return float(brentq(charge_residual, 0.0, 14.0, args=args, xtol=1e-10))


def solve_cation_charge(
    totals_molar: Mapping[str, float],
    byp_succinic_molar: float,
    pka_map: Mapping[str, tuple[float, ...]],
    target_ph: float,
) -> float:
    """Inverse anchoring — the strong-cation charge that reproduces ``target_ph``.

    Closed form (no root-find): electroneutrality at the target pH rearranges to
    ``cation = [OH⁻] + Σ acid-anion charge − [H⁺]``. Raises ``ValueError`` if the result
    is negative — that means ``target_ph`` is *below* what the acid load alone produces
    (an unphysical negative strong-cation charge), surfaced at compile rather than
    packed into a nonsense state.
    """
    h = 10.0 ** (-target_ph)
    oh = KW / h
    anion = byp_succinic_molar * mean_charge(h, pka_map[BYP_KEY])
    for name, conc in totals_molar.items():
        anion += conc * mean_charge(h, pka_map[name])
    cation = oh + anion - h
    if cation < 0.0:
        raise ValueError(
            f"initial_ph={target_ph} is below the acid load's intrinsic pH: it implies a "
            f"negative strong-cation charge ({cation:.4e} mol/L). Lower the acids or raise "
            "initial_ph."
        )
    return float(cation)


# -- in-loop / readout helpers (pure, signature matches Process.derivatives) ---


def build_pka_map(params: Mapping[str, float]) -> dict[str, tuple[float, ...]]:
    """Adapter: resolved ``{name: float}`` params → ``{acid: (pKa, …)}`` for the solver.

    Single source of truth reused by ``ph_of_state``, the analysis series helpers and the
    compile back-solve, so the three call sites cannot drift in how they read pKa. Keyed
    by acid slot name plus :data:`BYP_KEY`.
    """
    out: dict[str, tuple[float, ...]] = {
        name: tuple(params[n] for n in spec.pka_param_names) for name, spec in ACID_STATE.items()
    }
    out[BYP_KEY] = tuple(params[n] for n in BYP_AS_SUCCINIC.pka_param_names)
    return out


def _totals_molar(y: FloatArray, schema: StateSchema) -> dict[str, float]:
    """Acid slot concentrations present in ``schema``, converted g/L → mol/L."""
    return {
        name: float(y[schema.slice(name)][0]) / spec.molar_mass
        for name, spec in ACID_STATE.items()
        if name in schema
    }


def _byp_succinic_molar(y: FloatArray, schema: StateSchema) -> float:
    """The ``Byp`` pool (g/L) read as succinic-equivalent mol/L, or 0 if absent."""
    if BYP_KEY not in schema:
        return 0.0
    return float(y[schema.slice(BYP_KEY)][0]) / BYP_AS_SUCCINIC.molar_mass


def _cation(y: FloatArray, schema: StateSchema) -> float:
    """The net strong-cation charge (mol⁺/L) from the state slot, or 0 if absent."""
    if "cation_charge" not in schema:
        return 0.0
    return float(y[schema.slice("cation_charge")][0])


def ph_of_state(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """pH of a single state vector — pure, in-loop signature ``(y, schema, params)``.

    ``params`` is the RESOLVED ``{name: float}`` map (what ``Process.derivatives`` /
    ``RateModifier.factor`` receive), so the future in-loop pH→rate hook is genuinely
    free. Reads the ACID_STATE slots present in ``schema`` (g/L → mol/L via chemistry
    molar masses), the ``cation_charge`` slot, and ``Byp`` as succinic-equivalent.
    """
    return solve_ph(
        _totals_molar(y, schema),
        _cation(y, schema),
        _byp_succinic_molar(y, schema),
        build_pka_map(params),
    )


def titratable_acidity(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Titratable acidity [g/L tartaric-equivalent] — a second derived pure function.

    Approximates the pH-8.2 titration endpoint as *full* carboxylic deprotonation: each
    acid contributes ``C_a·(z_max − mean_charge(current pH))`` eq/L of still-titratable
    protons, summed and weighted by the tartaric equivalent weight (M_tartaric/2 ≈ 75.04
    g/eq). Omits the free-[H⁺] term (~0.4 mM, ~0.4 % at must pH) that conventional
    TA-to-8.2 includes — fine for the plausible tier and the 6–9 g/L must band. Uses the
    same acid state as ``ph_of_state`` (no new state).

    CAVEAT — trust the *must* (t=0) value, not the end-of-ferment series. The whole ``Byp``
    pool is read as a fully-titratable diprotic succinic acid, so as ``Byp`` accumulates
    over a ferment (the D-16/D-19 realised-yield diversion, ~3 g/L) the computed TA *rises*
    ~3–4 g/L. Real wine TA is flat-to-*declining* during fermentation (tartrate
    precipitation, malic metabolism), so the end-of-ferment TA here is an **over-estimate,
    not a fidelity-grade readout**. The cause is upstream pool sizing/booking (``Byp``
    lumps neutral 2,3-butanediol yet is booked diprotic; the pool itself exceeds real
    succinic 0.5–1.5 g/L), not this function — which is exact given its inputs. Bounded for
    *pH* as minor (~1–1.5 mM vs ~20 mM buffer) by D-18; the *TA* impact is direct and
    larger. Use the t=0 must TA as the band check; treat the series as directional only.
    """
    pka_map = build_pka_map(params)
    totals = _totals_molar(y, schema)
    cation = _cation(y, schema)
    byp = _byp_succinic_molar(y, schema)
    ph = solve_ph(totals, cation, byp, pka_map)
    h = 10.0 ** (-ph)

    eq_per_l = byp * (BYP_AS_SUCCINIC.protons - mean_charge(h, pka_map[BYP_KEY]))
    for name, conc in totals.items():
        eq_per_l += conc * (ACID_STATE[name].protons - mean_charge(h, pka_map[name]))
    return eq_per_l * (M_TARTARIC / 2.0)


def ph_tier(params_tier_of: Mapping[str, Tier]) -> Tier:
    """The tier of a derived pH/TA value — computed EXPLICITLY (decision D-18).

    The lowest of the pKa parameter tiers floored at ``PLAUSIBLE`` (applying 25 °C / I=0
    constants to wine is an extrapolation, so pH is never ``VALIDATED`` however good the
    pKa source). Takes a ``{name: Tier}`` map (e.g. ``ParameterSet.tier_map()``) — kept
    separate from the resolved-float hot-loop signature above. It must NOT read the tier
    of the inert acid *state* slots (no Process touches them, so ``ProcessSet.tier_of``
    returns ``VALIDATED`` for them — which would over-report pH's confidence).
    """
    tiers = [params_tier_of[n] for n in PKA_PARAM_NAMES if n in params_tier_of]
    return combine([*tiers, Tier.PLAUSIBLE])


# -- SO₂ speciation readout (pH-coupled molecular fraction, decision D-22) ------


def molecular_so2_fraction(ph: float, pkas: tuple[float, ...]) -> float:
    """Fraction of free SO₂ present as molecular (antimicrobial) SO₂·H₂O at ``ph``.

    The molecular species is the undissociated acid, so this is just
    :func:`neutral_fraction` evaluated at ``h = 10⁻ᵖᴴ``. With the sulfurous pKa₁ ≈ 1.81
    (and pKa₂ ≈ 7.2 negligible at wine pH) it reduces to the textbook
    ``1/(1 + 10^(pH − pKa₁))`` and falls ~3× per 0.5 pH unit — the pH coupling the D-18
    charge-balance keystone exists to make emerge ("dose SO₂ → speciation falls out of
    the current pH"). Bisulfite/sulfite are the remaining ``mean_charge`` share.
    """
    return neutral_fraction(10.0 ** (-ph), pkas)


def _so2_total(y: FloatArray, schema: StateSchema) -> float:
    """The total-SO₂ pool (g/L), or 0 if the schema has no :data:`SO2_STATE_KEY` slot."""
    if SO2_STATE_KEY not in schema:
        return 0.0
    return float(y[schema.slice(SO2_STATE_KEY)][0])


def _acetaldehyde_molar(y: FloatArray, schema: StateSchema) -> float:
    """The acetaldehyde pool as mol/L, or 0 if the :data:`ACETALDEHYDE_KEY` slot is absent.

    The SO₂ binder (decision D-28). Clamped ≥ 0 against solver undershoot, exactly as the
    D-27 acetaldehyde reduction clamps it.
    """
    if ACETALDEHYDE_KEY not in schema:
        return 0.0
    return max(float(y[schema.slice(ACETALDEHYDE_KEY)][0]), 0.0) / M_ACETALDEHYDE


def _pyruvate_molar(y: FloatArray, schema: StateSchema) -> float:
    """The excreted-overflow pyruvate pool as mol/L, or 0 if the :data:`PYRUVATE_KEY` slot
    is absent (decision D-51: the second SO₂ binder, built as real state in D-49)."""
    if PYRUVATE_KEY not in schema:
        return 0.0
    return max(float(y[schema.slice(PYRUVATE_KEY)][0]), 0.0) / M_PYRUVATE


def _alpha_kg_molar(y: FloatArray, schema: StateSchema) -> float:
    """The excreted-overflow α-ketoglutarate pool as mol/L, or 0 if the :data:`ALPHA_KG_KEY`
    slot is absent (decision D-51: the third SO₂ binder, built as real state in D-50)."""
    if ALPHA_KG_KEY not in schema:
        return 0.0
    return max(float(y[schema.slice(ALPHA_KG_KEY)][0]), 0.0) / M_ALPHA_KETOGLUTARATE


def bound_so2_molar(
    total_molar: float,
    carbonyls: tuple[tuple[float, float], ...],
    bisulfite_frac: float,
) -> tuple[float, ...]:
    """Bound SO₂ [mol/L] for each of N carbonyls competing for one shared bisulfite pool.

    Generalizes the D-28 single-carbonyl equilibrium to a competitive-Langmuir partition
    (decision D-51): every carbonyl ``i`` forms a 1:1 adduct with bisulfite,
    ``K_i = (A_i − x_i)·[HSO₃⁻] / x_i``, and all of them draw on the SAME unbound-SO₂ pool
    ``u = C − Σx_i`` (``C`` = total SO₂, mol/L), whose reactive bisulfite share is
    ``h = u·β``. Solving each adduct for its own bound amount gives the Langmuir form
    ``x_i = A_i·h / (K_i + h)``, and substituting into the shared-pool definition of ``h``
    collapses the whole system to ONE monotone equation in ``h``:

        ``β·C − β·Σᵢ(A_i·h/(K_i+h)) − h = 0``

    (strictly decreasing: the sum term rises with ``h``, ``−h`` falls), so a single
    ``brentq`` root-find over ``[0, β·C]`` — where the sign change is guaranteed — replaces
    the old per-carbonyl quadratic. With one carbonyl this reduces *exactly* to the D-28
    quadratic root (test-verified), so a keto-acid-pool-off run is byte-for-byte the D-28/
    D-47 form. ``carbonyls`` is an ordered ``(species_molar, k_diss)`` sequence (e.g.
    ``(acetaldehyde_molar, K_acetaldehyde_so2)``); the return tuple matches its order and
    length. Each entry is clamped to ``[0, min(A_i, C)]``. Zero for every entry when
    ``total_molar``/``bisulfite_frac`` is zero, or (per-entry) when that carbonyl's own molar
    amount or ``K`` is zero. Pure algebra (no state) so the equilibrium is unit-testable in
    isolation.
    """
    n = len(carbonyls)
    if total_molar <= 0.0 or bisulfite_frac <= 0.0:
        return (0.0,) * n
    active = [(a, k) for a, k in carbonyls if a > 0.0 and k > 0.0]
    if not active:
        return (0.0,) * n

    def shared_bisulfite_residual(h: float) -> float:
        bound_sum = sum(a * h / (k + h) for a, k in active)
        return bisulfite_frac * total_molar - bisulfite_frac * bound_sum - h

    h_max = bisulfite_frac * total_molar
    h = brentq(shared_bisulfite_residual, 0.0, h_max, xtol=1e-14)
    return tuple(
        min(max(a * h / (k + h), 0.0), min(a, total_molar)) if a > 0.0 and k > 0.0 else 0.0
        for a, k in carbonyls
    )


@dataclass(frozen=True)
class So2Speciation:
    """The full free/bound/molecular SO₂ split of one state vector (decisions D-28, D-51).

    All concentrations are g/L expressed *as SO₂* (mass-preserving): ``total`` is the dosed
    conserved slot, ``bound`` the combined acetaldehyde/pyruvate/α-ketoglutarate-hydroxysulphonate
    adducts (decision D-51 — the three carbonyls share ONE bisulfite pool, see
    :func:`bound_so2_molar`), ``free = total − bound`` the analytically-measured free SO₂, and
    ``molecular`` the antimicrobial (undissociated SO₂·H₂O) share of *free* at ``ph``.
    ``molecular_fraction`` is that share of free (the D-22 curve). CAVEAT (D-51, narrower than
    the pre-D-51 D-28 caveat): ``bound`` still omits SUGARS, the remaining known wine SO₂
    binder — so ``bound`` still under-estimates and ``free``/``molecular`` still slightly
    over-estimate the protective pool, just by less than before; the ``total`` slot is really
    "free + acetaldehyde/pyruvate/α-KG-bound".
    """

    ph: float
    total: float
    bound: float
    free: float
    molecular: float
    molecular_fraction: float


def _bound_molar_split(
    total_so2_molar: float,
    acetaldehyde_molar: float,
    pyruvate_molar: float,
    alpha_kg_molar: float,
    bisulfite_frac: float,
    params: Mapping[str, float],
) -> tuple[float, float, float]:
    """The 3-carbonyl coupled bind (decision D-51) — the one place all three ``K``s meet.

    Shared by :func:`_speciate_at_ph` (which only needs the *sum*, for ``bound``) and
    :func:`free_acetaldehyde` (which needs acetaldehyde's *own* share of that sum) so the two
    call sites can never disagree about how the SO₂ is partitioned across carbonyls.
    """
    acet, pyr, akg = bound_so2_molar(
        total_so2_molar,
        (
            (acetaldehyde_molar, params[SO2_BINDING_PARAM]),
            (pyruvate_molar, params[PYRUVATE_SO2_BINDING_PARAM]),
            (alpha_kg_molar, params[ALPHA_KG_SO2_BINDING_PARAM]),
        ),
        bisulfite_frac,
    )
    return acet, pyr, akg


def _speciate_at_ph(
    total: float,
    acetaldehyde_molar: float,
    pyruvate_molar: float,
    alpha_kg_molar: float,
    ph: float,
    params: Mapping[str, float],
) -> So2Speciation:
    """Split a known total SO₂ (g/L) at an already-solved ``ph`` — the shared inner kernel.

    Factored out so the pH is solved exactly once per call site: the public readouts solve
    it via :func:`ph_of_state`, while an in-loop consumer that has *already* solved pH (the
    MLF antimicrobial gate) passes it straight in. SO₂ is expressed as SO₂, so g/L ⇄ mol/L
    is a plain ``M_SO2`` scale and the split is mass-preserving.
    """
    pkas = tuple(params[n] for n in SO2_PKA_PARAM_NAMES)
    h = 10.0 ** (-ph)
    beta = bisulfite_fraction(h, pkas)
    acet, pyr, akg = _bound_molar_split(
        total / M_SO2, acetaldehyde_molar, pyruvate_molar, alpha_kg_molar, beta, params
    )
    bound = (acet + pyr + akg) * M_SO2
    free = max(total - bound, 0.0)
    frac = molecular_so2_fraction(ph, pkas)
    return So2Speciation(
        ph=ph, total=total, bound=bound, free=free, molecular=free * frac, molecular_fraction=frac
    )


def speciate_so2(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> So2Speciation:
    """Full free/bound/molecular SO₂ split of a state vector — a derived pure function.

    The headline D-28 readout, generalized by D-51 to three competing carbonyls: pH is solved
    from the organic acids (:func:`ph_of_state` — SO₂ is **not** in the charge balance,
    decisions D-22/D-28), the dosed total SO₂ is split into acetaldehyde/pyruvate/α-KG-bound
    vs free at that pH (the coupled binding equilibrium), and the free pool partitioned into
    its molecular antimicrobial share. In-loop signature ``(y, schema, params)`` like
    :func:`ph_of_state`. Returns an all-zero split (at the solved pH) when no SO₂ is dosed or
    the slot is absent. At acetaldehyde = pyruvate = α-KG = 0 the bound term vanishes and
    ``free == total`` — so D-22's molecular-SO₂ curve is reproduced exactly (the regression
    anchor), and on a keto-acid-pool-off run (pyruvate/α-KG slots absent from the schema) this
    is byte-for-byte the pre-D-51 acetaldehyde-only split. Report the conventional mg/L via
    :func:`fermentation.units.convert.gpl_to_mgl`.
    """
    total = _so2_total(y, schema)
    ph = ph_of_state(y, schema, params)
    if total <= 0.0:
        frac = molecular_so2_fraction(ph, tuple(params[n] for n in SO2_PKA_PARAM_NAMES))
        return So2Speciation(
            ph=ph, total=0.0, bound=0.0, free=0.0, molecular=0.0, molecular_fraction=frac
        )
    return _speciate_at_ph(
        total, _acetaldehyde_molar(y, schema), _pyruvate_molar(y, schema),
        _alpha_kg_molar(y, schema), ph, params,
    )  # fmt: skip


def bound_so2(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Combined acetaldehyde/pyruvate/α-KG-bound SO₂ [g/L] (decisions D-28, D-51) — 0 when
    undosed."""
    return speciate_so2(y, schema, params).bound


def free_so2(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Free SO₂ [g/L] of a state vector — total minus the combined carbonyl-bound share
    (decisions D-28, D-51)."""
    return speciate_so2(y, schema, params).free


def molecular_so2_at_ph(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], ph: float
) -> float:
    """Molecular (antimicrobial) SO₂ [g/L] at an **already-solved** ``ph`` (decisions D-28, D-51).

    For in-loop consumers that have solved pH themselves (the MLF antimicrobial gate reuses
    its single brentq solve): free = total − (acetaldehyde+pyruvate+α-KG)-bound at ``ph``, then
    the molecular share of free. Returns 0 when no SO₂ is dosed. The bound term reads all three
    carbonyls from the same state, so as any of them peaks it sequesters SO₂ and the
    antimicrobial free pool (correctly) drops — bound SO₂ is not antimicrobial.
    """
    total = _so2_total(y, schema)
    if total <= 0.0:
        return 0.0
    return _speciate_at_ph(
        total, _acetaldehyde_molar(y, schema), _pyruvate_molar(y, schema),
        _alpha_kg_molar(y, schema), ph, params,
    ).molecular  # fmt: skip


def free_acetaldehyde(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float], ph: float
) -> float:
    """Acetaldehyde NOT bound to SO₂ [g/L] at an already-solved ``ph`` — the ADH-reducible pool.

    The RHS side of the D-28 free/bound split (decision D-47), reading the D-51 multi-carbonyl
    equilibrium. The acetaldehyde-bisulfite adduct (1-hydroxyethanesulphonate) is **protected
    from alcohol dehydrogenase** — the literature is explicit that "acetaldehyde bound to SO₂
    could not be metabolized by yeast during fermentation; only free acetaldehyde could impact
    metabolism" — so the enzymatic reduction
    (:class:`fermentation.core.kinetics.acetaldehyde.AcetaldehydeReduction`) acts on this free
    share alone. ``free = total_acetaldehyde − bound``, with ``bound`` acetaldehyde's OWN share
    of the same :func:`_bound_molar_split` the SO₂ readout uses — pyruvate and α-KG now compete
    for the shared bisulfite pool too (decision D-51), so this ``bound`` is *smaller* than the
    pre-D-51 acetaldehyde-only form whenever a keto-acid pool is present, leaving MORE free
    acetaldehyde for ADH (the mechanism behind D-51's smaller finished-wine stranding).

    Returns the **total** acetaldehyde (g/L) whenever no SO₂ is dosed or no acetaldehyde is
    present — the binding vanishes and the reduction is byte-for-byte the D-27 core (unaffected
    by whether pyruvate/α-KG pools exist, since with no acetaldehyde there is nothing here to
    protect). In-loop signature ``(y, schema, params, ph)`` like :func:`molecular_so2_at_ph`;
    the caller supplies the already-solved pH so this never triggers a second ``brentq``.
    Clamped ≥ 0.
    """
    acet_molar = _acetaldehyde_molar(y, schema)
    total_so2 = _so2_total(y, schema)
    if acet_molar <= 0.0 or total_so2 <= 0.0:
        return acet_molar * M_ACETALDEHYDE  # no binding ⇒ all acetaldehyde is reducible (D-27)
    pkas = tuple(params[n] for n in SO2_PKA_PARAM_NAMES)
    beta = bisulfite_fraction(10.0 ** (-ph), pkas)
    bound_acet, _, _ = _bound_molar_split(
        total_so2 / M_SO2, acet_molar, _pyruvate_molar(y, schema), _alpha_kg_molar(y, schema),
        beta, params,
    )  # fmt: skip
    return max(acet_molar - bound_acet, 0.0) * M_ACETALDEHYDE


def molecular_so2(y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> float:
    """Molecular (antimicrobial) SO₂ [g/L] of a state vector — a derived pure function.

    The D-22 headline, now reading the D-51 multi-carbonyl free/bound split: pH is solved from
    the organic acids, the dosed *total* SO₂ split into acetaldehyde/pyruvate/α-KG-bound vs
    free, and the molecular share of **free** returned (bound SO₂ is not antimicrobial).
    Convenience wrapper over :func:`speciate_so2`. Returns 0 when no SO₂ is dosed (or the slot
    is absent); at acetaldehyde = pyruvate = α-KG = 0 it collapses to the exact D-22
    ``free × fraction(pH)``.
    """
    return speciate_so2(y, schema, params).molecular


def molecular_so2_tier(params_tier_of: Mapping[str, Tier]) -> Tier:
    """Tier of the derived SO₂-speciation values — computed EXPLICITLY (decisions D-22/D-28/D-51).

    The readout solves pH (so it inherits every pH-solver pKa tier, exactly as
    :func:`ph_tier`), partitions by the sulfurous pKa(s), *and* (D-28/D-51) splits free/bound by
    all three carbonyl-bisulfite binding constants — so its tier is the lowest of **all these**
    parameter sets (pH pKas, sulfurous pKas, the 3 binding ``K``s), floored at ``PLAUSIBLE``
    (apparent constants applied to wine are extrapolation; SO₂ speciation is never ``VALIDATED``
    however good the source). Covers ``molecular``/``free``/``bound`` alike (they share these
    inputs). Like :func:`ph_tier` it must NOT read the inert SO₂/acetaldehyde/keto-acid/acid
    *state*-slot tiers (which ``tier_of`` reports ``VALIDATED`` for, since no Process touches
    them directly here).
    """
    names = (
        *PKA_PARAM_NAMES,
        *SO2_PKA_PARAM_NAMES,
        SO2_BINDING_PARAM,
        PYRUVATE_SO2_BINDING_PARAM,
        ALPHA_KG_SO2_BINDING_PARAM,
    )
    tiers = [params_tier_of[n] for n in names if n in params_tier_of]
    return combine([*tiers, Tier.PLAUSIBLE])
