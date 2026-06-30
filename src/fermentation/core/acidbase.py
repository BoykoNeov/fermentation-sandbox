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

from fermentation.core.chemistry import M_LACTIC, M_MALIC, M_SUCCINIC, M_TARTARIC
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
    """Solve ``charge_residual = 0`` for pH over the bracket ``[0, 14]``.

    The wide bracket guarantees a sign change (positive net charge at pH 0, negative at
    pH 14), so ``brentq`` converges on the single root without an initial guess.
    """
    return float(
        brentq(
            charge_residual,
            0.0,
            14.0,
            args=(totals_molar, cation, byp_succinic_molar, pka_map),
            xtol=1e-10,
        )
    )


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
