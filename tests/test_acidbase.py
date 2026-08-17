"""The wine pH charge-balance solver and its derived pH/TA readout (decision D-18).

Ranked headline-first. The keystone's proof-of-purpose is
``test_headline_malic_to_lactic_raises_ph``: the malic→lactic deacidification (the
chemistry MLF performs) raises pH by 0.1–0.3, demonstrated *without* an MLF Process
built — the solver responds to acid dynamics on its own. The rest pin the balance
(residual ≈ 0, monotonicity, smoothness), the inverse anchoring (round-trip, physical
back-solved cation, the unphysical-initial_ph compile guard), that the new acid slots
leave carbon conservation intact, that TA lands in the textbook band, and that the
derived pH tier is computed explicitly as ``plausible`` (never ``validated``).
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.analysis import ph_series, titratable_acidity_series
from fermentation.core import acidbase
from fermentation.core.chemistry import M_LACTIC, M_MALIC, M_NITROGEN, M_TARTARIC
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime.integrate import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import assert_conserved, total_carbon

#: Potassium molar mass [g/mol] — used only to state test 6's expected cation range
#: from first principles (K⁺ 1–2 g/L is the physical wine range), NOT read by the solver.
M_K = 39.0983


@pytest.fixture
def pset():
    """Real wine kinetic params + the shared pKa set, as a ParameterSet."""
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml")


@pytest.fixture
def params(pset):
    """Resolved ``{name: float}`` map the solver hot-loop signature consumes."""
    return pset.resolve()


@pytest.fixture
def pka(params):
    return acidbase.build_pka_map(params)


@pytest.fixture
def beer_params():
    """Beer's parameter set — needed for D-209's beer nitrogen-charge derivation."""
    data = default_data_dir()
    return load_parameters(
        data / "beer_generic.yaml", data / "acidbase.yaml", data / "beer_acids.yaml"
    )


#: The ``N`` these hand-built states carry. Named rather than inlined because since D-209 it is
#: no longer pH-irrelevant: the assimilable-nitrogen pool is on the cation side of the balance,
#: so :func:`_anchor_cation` has to subtract its charge to hit a target pH.
_WINE_STATE_NITROGEN_GPL = 0.2


def _wine_state(schema: StateSchema, **acids: float) -> FloatArray:
    """A wine state vector with arbitrary bulk values + given acids.

    ``N`` is NOT pH-irrelevant since D-209 — see :data:`_WINE_STATE_NITROGEN_GPL`.
    """
    base: dict[str, float | list[float]] = {
        "X": 0.5,
        "S": [240.0],
        "E": 0.0,
        "N": _WINE_STATE_NITROGEN_GPL,
        "T": 293.15,
        "CO2": 0.0,
    }
    base.update(acids)
    return schema.pack(base)


def _anchor_cation(
    pka: Mapping[str, tuple[float, ...]],
    tartaric_gpl: float,
    malic_gpl: float,
    target_ph: float,
    params: Mapping[str, float] | None = None,
) -> float:
    """The ``cation_charge`` SLOT value that anchors a :func:`_wine_state` at ``target_ph``.

    ``params`` is optional only so the pKa-map-only callers below keep working; pass it whenever
    the anchored state is going to be read back through :func:`acidbase.ph_of_state`. Since
    D-209 the assimilable-nitrogen pool is itself on the cation side, and ``_wine_state`` seeds
    ``N`` at 0.2 g/L, so the slot is the solved TOTAL minus what that nitrogen supplies —
    exactly the subtraction ``cation_charge_for_ph`` and both compile-seam anchors make. Without
    it these states read back ~0.12 pH ABOVE their target.
    """
    totals = {
        "tartaric": tartaric_gpl / M_TARTARIC,
        "malic": malic_gpl / M_MALIC,
        "lactic": 0.0,
    }
    total = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, target_ph)
    if params is None:
        return total
    return total - acidbase.nitrogen_charge_from_gpl(_WINE_STATE_NITROGEN_GPL, "wine", params)


# -- 1. HEADLINE: malic→lactic deacidification raises pH 0.1–0.3 ---------------


def test_headline_malic_to_lactic_raises_ph(params, pka):
    # A malic-rich must (the case where MLF matters): tartaric 4 / malic 4 g/L, anchored
    # to a measured pH 3.4. Full MLF converts all malic → lactic mole-for-mole; pH must
    # rise into the deacidification band. Settled empirically — a tartaric-heavy must
    # would land below 0.1; the fix for an out-of-band number is a more malic-rich must,
    # NOT widening the band (CLAUDE.md forbids weakening benchmark tests).
    schema = wine_schema()
    cation = _anchor_cation(pka, 4.0, 4.0, 3.4, params)
    y0 = _wine_state(schema, tartaric=4.0, malic=4.0, lactic=0.0, cation_charge=cation)
    ph0 = acidbase.ph_of_state(y0, schema, params)
    assert ph0 == pytest.approx(3.4, abs=1e-3)  # anchoring exact at t=0

    # all malic → lactic, conserving moles (g/L scales by the molar-mass ratio)
    lactic_gpl = (4.0 / M_MALIC) * M_LACTIC
    y1 = _wine_state(schema, tartaric=4.0, malic=0.0, lactic=lactic_gpl, cation_charge=cation)
    ph1 = acidbase.ph_of_state(y1, schema, params)

    delta = ph1 - ph0
    assert 0.1 <= delta <= 0.3, f"malic→lactic ΔpH {delta:.3f} outside [0.1, 0.3]"


# -- 2. the balance actually balances -----------------------------------------


def test_charge_residual_zero_at_solved_ph(pka):
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    assert acidbase.charge_residual(ph, totals, cation, 0.0, 0.0, pka) == pytest.approx(
        0.0, abs=1e-9
    )


# -- 3. monotonicity ----------------------------------------------------------


def test_more_acid_lowers_ph_more_cation_raises_ph(pka):
    base = {"tartaric": 4.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(base, 0.0, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(base, cation, 0.0, 0.0, pka)

    more_tartaric = {**base, "tartaric": 6.0 / M_TARTARIC}
    assert acidbase.solve_ph(more_tartaric, cation, 0.0, 0.0, pka) < ph  # more acid → lower pH

    assert acidbase.solve_ph(base, cation * 1.2, 0.0, 0.0, pka) > ph  # more cation → higher pH


# -- 4. smoothness / C¹ (guards a future in-loop BDF consumer) ----------------


def test_ph_is_smooth_in_acid(pka):
    # No ``params`` argument here on purpose: this test feeds the cation straight to
    # ``solve_ph`` with a hand-built totals map and never reads a STATE back, so there is no
    # nitrogen slot in play and the D-209 subtraction would be wrong rather than merely unneeded.
    cation = _anchor_cation(pka, 5.0, 3.0, 3.4)
    tartaric = np.linspace(3.0, 7.0, 41) / M_TARTARIC
    ph = np.array(
        [
            acidbase.solve_ph(
                {"tartaric": t, "malic": 3.0 / M_MALIC, "lactic": 0.0}, cation, 0.0, 0.0, pka
            )
            for t in tartaric
        ]
    )
    d1 = np.diff(ph)  # first difference (∝ dpH/d tartaric)
    d2 = np.diff(d1)  # second difference — small & sign-stable ⇒ no kink
    assert np.all(d1 < 0.0)  # strictly monotone
    assert np.max(np.abs(d2)) < 1e-2  # no derivative jump on a fine grid


# -- 4b. totality: solve_ph clamps a non-physiological probe cation, never raises (D-46) --
# BDF's num_jac perturbs the ``cation_charge`` state slot far outside its ~0.03 mol/L range
# while building the Jacobian, which can push ``charge_residual`` positive (or negative) across
# the whole [0, 14] bracket. ``solve_ph`` must stay a TOTAL, bounded function and clamp to the
# window rather than let ``brentq`` throw "f(a) and f(b) must have different signs". The three
# Brett integration tests only catch this incidentally (a 120-day run happens to drive the
# probe there); these pin it at the function level so a refactor cannot silently un-total it.


def test_solve_ph_clamps_huge_probe_cation_to_14(pka):
    # A probe cation two orders of magnitude above the physical ~0.03 mol/L: no acid load can
    # neutralise it, so the electroneutral pH lies above the window ⇒ clamp to 14, not raise.
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    assert acidbase.solve_ph(totals, 3.81, 0.0, 0.0, pka) == 14.0


def test_solve_ph_clamps_negative_probe_cation_to_0(pka):
    # The mirror probe: a large acid load with a strongly NEGATIVE strong-cation charge is
    # net-negative even fully protonated ⇒ electroneutral pH below the window ⇒ clamp to 0.
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    assert acidbase.solve_ph(totals, -2.0, 0.0, 0.0, pka) == 0.0


def test_solve_ph_physiological_cation_falls_through_to_brentq(pka):
    # The untouched path: a physiological cation still returns an interior root, unclamped and
    # bit-for-bit identical to the brentq result (the clamp branches are never taken).
    totals = {"tartaric": 6.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    ph = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    assert 0.0 < ph < 14.0
    assert ph == pytest.approx(3.4, abs=1e-6)  # inverts solve_cation_charge exactly


# -- 5. round-trip: a compiled scenario reproduces its measured initial_ph -----
# NB tautological w.r.t. the g/L→mol/L factor (solve_cation_charge / solve_ph are
# inverses applying the same conversion, so a unit bug cancels). Test 6 is the guard.


def _wine_scenario(**initial_extra: float) -> Scenario:
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.5}
    initial.update(initial_extra)
    return Scenario(
        name="wine-ph",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=14.0,
    )


def test_compiled_scenario_reproduces_initial_ph():
    compiled = compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4))
    ph0 = acidbase.ph_of_state(compiled.y0, compiled.schema, compiled.param_values)
    assert ph0 == pytest.approx(3.4, abs=1e-3)


# -- 6. THE UNITS GUARD: back-solved cation is physical (K⁺ 1–2 g/L) -----------


def test_back_solved_cation_is_physical():
    # Independent of the solver's own arithmetic: a textbook must (TA ~6–9 g/L, pH 3.4)
    # is counter-charged by K⁺ ~1–2 g/L, i.e. ~25.6–51.2 meq/L (1–2 g/L ÷ 39.0983).
    # A g/L↔mol/L factor error (which the round-trip can't see) lands the cation orders
    # of magnitude outside this band, so this is the clean catch.
    compiled = compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4))
    cation = compiled.schema.get(compiled.y0, "cation_charge")  # mol⁺/L
    meq_per_l = cation * 1000.0
    lo, hi = 1.0 / M_K * 1000.0, 2.0 / M_K * 1000.0  # ~25.6 .. 51.2 meq/L
    assert lo <= meq_per_l <= hi, (
        f"cation {meq_per_l:.1f} meq/L outside physical K⁺ {lo:.1f}–{hi:.1f}"
    )


# -- 7. compile guard: initial_ph below the acid load's intrinsic pH raises ----


def test_unphysical_initial_ph_raises_at_compile():
    # pH 2.0 with 6 g/L tartaric + 3 g/L malic needs a NEGATIVE strong cation — unphysical.
    with pytest.raises(ValueError, match="unphysical"):
        compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=2.0))


def test_initial_ph_without_pka_params_raises_clearly():
    # The explicit parameter_paths override is caller-controlled: a caller asking for
    # initial_ph but NOT including acidbase.yaml must get a clear, actionable error
    # (the missing-pKa KeyError is caught and re-raised), not a raw KeyError.
    wine_only = default_data_dir() / "wine_generic.yaml"  # deliberately omit acidbase.yaml
    with pytest.raises(ValueError, match="acidbase.yaml"):
        compile_scenario(
            _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4),
            parameter_paths=[wine_only],
        )


# -- 8. carbon conservation unchanged by the new acid slots -------------------


def test_carbon_conserved_with_constant_acids():
    compiled = compile_scenario(
        _wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4), strict=True
    )
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    carbon = total_carbon(
        compiled.schema, biomass_carbon_fraction=compiled.parameters["biomass_C_fraction"].value
    )
    # Acids are inert (no Process touches them) — a constant offset that drifts 0; the
    # rest of the carbon ledger still closes to machine precision with them present.
    assert_conserved(traj, carbon, rtol=1e-6, atol=1e-9, label="total carbon (with acids)")


# -- 9. TA lands in the textbook band -----------------------------------------


def test_titratable_acidity_in_band(params):
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4, params)
    y = _wine_state(schema, tartaric=6.0, malic=3.0, lactic=0.0, cation_charge=cation)
    ta = acidbase.titratable_acidity(y, schema, params)
    assert 6.0 <= ta <= 9.0, f"TA {ta:.2f} g/L tartaric-equiv outside the 6–9 band"


# -- 10. tier is computed explicitly as plausible, never validated ------------


def test_ph_tier_is_plausible(pset):
    # Computed explicitly (not inherited): the lowest pKa tier floored at plausible. The
    # pKa params are all plausible, and pH is never validated however good the pKa source.
    #
    # SCOPED TO WINE since D-179. The schema argument is not decoration: beer's lumped
    # peptide buffer is honestly speculative, so an unscoped call combines it in and reports
    # SPECULATIVE — for an acid this wine does not carry. Asking about wine's pH now means
    # saying "wine".
    assert acidbase.ph_tier(pset.tier_map(), wine_schema()) is Tier.PLAUSIBLE


def test_ph_tier_unscoped_is_conservative_not_wine(pset):
    """The no-schema default reports the WORST tier across all media, and that is deliberate.

    A caller with no schema does not know which medium it holds. Defaulting to wine's set
    would be back-compatible and would quietly over-report a beer pH as plausible; defaulting
    to the union can only ever under-claim. Prime directive #1 says an output's tier is the
    LOWEST of its inputs, so the conservative default is the correct one — but it is a real
    behaviour change from D-178, so it is pinned rather than left implicit.
    """
    assert acidbase.ph_tier(pset.tier_map()) is Tier.SPECULATIVE
    assert acidbase.ph_tier(pset.tier_map(), beer_schema()) is Tier.SPECULATIVE


# -- analysis series + the emergent Byp pH drift (the second demonstration) ----


def test_ph_series_drifts_down_from_byp_and_nitrogen_uptake():
    """A wine's pH falls over a ferment for TWO reasons, and this decomposes them exactly.

    The acid slots are constant, yet pH is not flat. Until D-209 the only driver was ``Byp``
    (the core realised-yield diversion), which grows over the ferment and which the charge
    balance reads by inclusion — worth a mild ~0.06 pH. D-209 adds the larger one: the
    assimilable-nitrogen pool carries net positive charge, so the yeast consuming it removes
    cation charge and the pH falls further. Total measured drift ~0.19 pH.

    The decomposition needs no parameter surgery and no second run. Re-reading the END state
    with ``N`` put back to its t=0 value makes the nitrogen term contribute exactly what it
    contributed at the anchor, which is the frozen-constant case — i.e. bit-for-bit the
    pre-D-209 balance — so the difference between that pH and the real end pH is the nitrogen
    term's whole contribution, and the rest is ``Byp``.
    """
    compiled = compile_scenario(_wine_scenario(tartaric_gpl=6.0, malic_gpl=3.0, initial_ph=3.4))
    traj = simulate(compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h)
    ph = ph_series(traj, compiled.param_values)
    ta = titratable_acidity_series(traj, compiled.param_values)
    assert ph[0] == pytest.approx(3.4, abs=1e-3)  # anchored at pitch
    drift = ph[0] - ph[-1]
    assert 0.10 <= drift <= 0.30, f"pH drift {drift:.3f} outside the expected mild fall"

    schema, params_ = compiled.schema, compiled.param_values
    end = np.asarray(traj.y[:, -1], dtype=float).copy()
    end[schema.slice("N")] = traj.y[schema.slice("N"), 0]
    byp_only_drift = ph[0] - acidbase.ph_of_state(end, schema, params_)
    nitrogen_share = drift - byp_only_drift
    assert 0.02 <= byp_only_drift <= 0.09, (
        f"the Byp-only component is {byp_only_drift:.4f} pH; D-18 through D-208 measured this "
        "whole drift at 0.02-0.15 and it was all Byp, so this component must not have moved"
    )
    assert nitrogen_share > byp_only_drift, (
        f"nitrogen uptake contributes {nitrogen_share:.4f} pH against Byp's "
        f"{byp_only_drift:.4f}; since D-209 it is the LARGER of the two drivers of a wine's "
        "falling pH, and a reversal here means the term has been scoped away"
    )

    # The MUST (t=0) TA is the fidelity-grade value, in the 6-9 g/L band. The TA SERIES
    # then RISES as Byp accumulates (whole pool read as titratable diprotic succinic) —
    # an ACKNOWLEDGED upstream artifact (D-16/D-19 pool booking), backwards to real wine
    # (TA flat-to-declining during ferment). Pinned here as known/directional, NOT
    # fidelity: see acidbase.titratable_acidity caveat. Don't "fix" it by changing D-18.
    assert ta.shape == traj.t.shape and np.all(ta > 0.0)
    assert 6.0 <= ta[0] <= 9.0  # must value is trustworthy
    assert ta[-1] > ta[0]  # documents the artifact rise (do not assert end-of-ferment band)


# ======================================================================================
# D-182 — dissolved CO2 as carbonic acid in the charge balance
# ======================================================================================


def test_carbonic_is_not_a_member_of_any_acid_registry():
    """The structural claim, asserted because getting it wrong FAILS SILENTLY.

    Dissolved CO2 has no state slot: it is derived from ``CO2`` and ``T``. The acid registries
    are lists of SLOTS, and ``_totals_molar`` skips any name not in the schema — so putting
    ``carbonic`` in ``WINE_ACIDS``/``BEER_ACIDS`` would not raise, it would resolve to nothing
    and the term would quietly vanish from the balance while looking present in the source.
    That is why it travels as its own argument, the way ``Byp`` does.
    """
    for registry in (acidbase.WINE_ACIDS, acidbase.BEER_ACIDS, acidbase.ALL_ACIDS):
        assert acidbase.CARBONIC_KEY not in registry
    assert acidbase.CARBONIC_KEY not in wine_schema()
    assert acidbase.CARBONIC_KEY not in beer_schema()
    # ...but it IS a key in the pKa map, exactly like Byp, or charge_residual would KeyError.
    data = default_data_dir()
    pmap = acidbase.build_pka_map(
        load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml").resolve()
    )
    assert acidbase.CARBONIC_KEY in pmap and acidbase.BYP_KEY in pmap


def test_the_carbonic_parameters_are_all_in_the_sampled_read_set():
    """The D-160 guard, in the one place the derivation could not reach.

    ``PKA_PARAM_NAMES`` derives itself from the acid registries so a new SLOT cannot be
    forgotten. Carbonic is not a registry member, so the derivation would have skipped it and
    ``pKa_carbonic_1`` would sit outside ``PH_SYSTEM_READS`` — undeclared, therefore never
    sampled, therefore narrowing the reported spread of every pH-dependent output below what
    its own provenance justifies. The three solubility parameters are not pKas at all and
    could not have arrived by that route under any derivation.
    """
    for name in ("pKa_carbonic_1", *acidbase.CO2_SOLUBILITY_PARAMS):
        assert name in acidbase.PH_SYSTEM_READS, f"{name} would never be sampled"
    # Still disjoint from the SO2 set — a pH-only caller must not declare the binding Ks.
    assert not set(acidbase.PH_SYSTEM_READS) & set(acidbase.SO2_BINDING_READS)


def test_carbonic_is_monoprotic_and_never_reaches_the_polyprotic_branch(params):
    """pKa2 (10.3) is omitted, not deferred: at any beverage pH the carbonate fraction is
    below 1e-9, so a second proton would be arithmetic on nothing. Asserted rather than
    asserted-in-prose, because "add the second pKa for completeness" is the obvious edit."""
    assert acidbase.CARBONIC_AS_CO2.protons == 1
    # What the omission is worth, computed with a hand-built diprotic set rather than quoted,
    # at the HIGHEST pH any medium here starts from (a 5.65 wort, taken at 5.7 to be generous
    # to the second proton). The number that matters is the RATIO to the first proton, not the
    # absolute: an absolute threshold would silently pass on a pool that had grown.
    diprotic = (params["pKa_carbonic_1"], 10.3)
    h = 10.0 ** (-5.7)
    mono = acidbase.mean_charge(h, (params["pKa_carbonic_1"],))
    di = acidbase.mean_charge(h, diprotic)
    assert (di - mono) / mono < 1e-4, (
        f"the second proton adds {(di - mono) / mono:.2e} of the first proton's charge "
        f"({di - mono:.2e} per mole against {mono:.3f}); D-182 measured 4.6e-05"
    )


def test_saturation_is_the_printed_constant_at_the_reference_and_falls_with_temperature(params):
    """Henry's law with a van 't Hoff transfer — exact at T_ref, monotone decreasing above it.

    At the reference temperature the transfer is exactly 1, so the model uses the PRINTED
    in-beer constant unscaled (the ``T_ref`` idiom of the Arrhenius modifiers). Away from it
    solubility falls as temperature rises, which is the direction every source states.
    """
    ref = params["T_ref_co2_solubility"]
    assert acidbase.co2_saturation_gpl(ref, params) == pytest.approx(
        params["H_co2_beverage"], rel=1e-12
    ), "the transfer must be exactly 1 at the temperature the constant is quoted at"
    temps = [275.0, 283.15, 288.15, 293.15, 303.15]
    sats = [acidbase.co2_saturation_gpl(t, params) for t in temps]
    assert all(a > b for a, b in zip(sats[:-1], sats[1:], strict=True)), (
        "CO2 solubility must FALL as temperature rises"
    )
    # A lager at 10 C holds meaningfully more than an ale at 20 C — the reason this ships
    # temperature-dependent rather than as one number.
    assert sats[1] / sats[3] > 1.25


def test_saturation_is_total_under_a_jacobian_probe(params):
    """Bounded like ``solve_ph``'s bracket (D-46): ``num_jac`` perturbs the ``T`` slot, and a
    non-positive temperature would make 1/T blow up. The probe must not be able to raise."""
    assert acidbase.co2_saturation_gpl(0.0, params) == 0.0
    assert acidbase.co2_saturation_gpl(-5.0, params) == 0.0
    assert acidbase.co2_saturation_gpl(1e-6, params) >= 0.0


def test_the_evolved_slot_is_capped_at_saturation_and_is_zero_before_fermentation(params):
    """``min(evolved, saturation)`` — the two ends of the term's whole behaviour.

    The ``CO2`` slot is a cumulative evolved-gas integral reaching ~40 g/L in a beer and
    ~100 g/L in a wine, against a saturation near 2. Reading it directly would put twenty to
    fifty times too much carbonic acid in the balance, which is the single most likely
    misreading of this module.
    """
    schema = wine_schema()
    m_co2 = acidbase.CARBONIC_AS_CO2.molar_mass
    sat = acidbase.co2_saturation_gpl(293.15, params)

    unfermented = _wine_state(schema, CO2=0.0, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(unfermented, schema, params) == 0.0

    finished = _wine_state(schema, CO2=100.0, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(finished, schema, params) * m_co2 == pytest.approx(sat)

    # Below saturation the liquid holds ALL of it — the mass-conservative half of the min.
    early = _wine_state(schema, CO2=0.5, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(early, schema, params) * m_co2 == pytest.approx(0.5)

    # Solver undershoot on the evolved pool cannot make the term negative.
    undershot = _wine_state(schema, CO2=-1e-9, tartaric=5.0)
    assert acidbase.dissolved_co2_molar(undershot, schema, params) == 0.0


def test_the_inverse_anchor_cannot_absorb_the_carbonic_term(pka, params):
    """THE STRUCTURAL REASON THIS TERM MOVES pH AT ALL, and the contrast that proves it.

    D-178 rejected malt phosphate on the ground that a species of CONSTANT charge is absorbed
    outright by the back-solved cation and is a near no-op. Dissolved CO2 escapes that only
    because it is ~0 in the must/wort where the anchor is taken and saturated afterwards.

    So this asserts BOTH halves: anchoring CO2-free and then finishing saturated moves the
    finished pH, while anchoring a state that ALREADY carries the same CO2 gives it back —
    the phosphate result, reproduced on purpose as the counterfactual.
    """
    totals = {"tartaric": 5.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    sat_molar = acidbase.co2_saturation_gpl(293.15, params) / acidbase.CARBONIC_AS_CO2.molar_mass

    # Real ordering: anchor with no CO2, finish with it.
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    ph_free = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka)
    ph_saturated = acidbase.solve_ph(totals, cation, 0.0, sat_molar, pka)
    assert ph_free == pytest.approx(3.4, abs=1e-9)
    assert ph_saturated < ph_free, "dissolved CO2 must acidify"

    # Counterfactual: the same CO2 present when the anchor is taken. The cation absorbs it and
    # the finished pH is the anchor again — a no-op, exactly as D-178 measured for phosphate.
    absorbed_cation = acidbase.solve_cation_charge(totals, 0.0, sat_molar, pka, 3.4)
    ph_absorbed = acidbase.solve_ph(totals, absorbed_cation, 0.0, sat_molar, pka)
    assert ph_absorbed == pytest.approx(3.4, abs=1e-9)
    assert absorbed_cation > cation, "absorbing the anion charge takes MORE counter-cation"


def test_titratable_acidity_models_a_degassed_sample(params):
    """A titration is run on a DEGASSED sample, so carbonic must not count — and the check has
    to be that TA is BITWISE unchanged by the CO2 in the state, not merely close.

    OIV titratable acidity explicitly excludes carbonic (and sulfurous) acid; brewers and
    winemakers degas before titrating for exactly this reason. The asymmetry with
    ``ph_of_state``, which DOES read it, is the measurement convention rather than an
    oversight — and it is the kind of asymmetry a later reader tidies away.
    """
    schema = wine_schema()
    flat = _wine_state(schema, CO2=0.0, tartaric=5.0, malic=3.0, cation_charge=0.05)
    carbonated = _wine_state(schema, CO2=100.0, tartaric=5.0, malic=3.0, cation_charge=0.05)

    assert acidbase.titratable_acidity(flat, schema, params) == acidbase.titratable_acidity(
        carbonated, schema, params
    ), "TA must be bit-for-bit identical — it models a degassed sample"
    # ...while pH is NOT, which is what makes the exclusion a choice rather than a no-op.
    assert acidbase.ph_of_state(carbonated, schema, params) < acidbase.ph_of_state(
        flat, schema, params
    )


def test_the_model_sits_below_the_printed_volumes_band_and_the_head_free_scope_is_why(params):
    """The one cross-check that is independent of every shipped number.

    "The Chemistry of Beer" sec 6.7.1 prints, as an OUTCOME rather than a constant, that beer
    contains 1.2-1.7 volumes of CO2 after a normal nonpressurized fermentation (1 volume =
    0.196 % w/w). The model gives ~1.02 volumes at 15 C. It lands BELOW that band, and the
    reason is the scope stated in the YAML: this term is fixed at 1 atm CO2 partial pressure,
    while a real tank carries metres of hydrostatic head. Pinned as DIRECTIONAL — model at or
    under the printed low edge — never as a target to hit, because the two are not the same
    quantity. Reproducing 1.7 volumes at 1 atm would need a fermentation at about -2 C.
    """
    gpl_per_volume = 1.96  # 1 volume = 0.196 % by weight, printed in the same sentence
    for celsius in (10.0, 15.0, 20.0):
        volumes = acidbase.co2_saturation_gpl(273.15 + celsius, params) / gpl_per_volume
        assert 0.8 < volumes <= 1.2, (
            f"at {celsius:g} C the model holds {volumes:.2f} volumes; it must stay at or below "
            "the printed 1.2-1.7 band's low edge, which the head-free 1 atm scope explains. "
            "Above it would mean the model had acquired vessel pressure from somewhere."
        )


def test_the_apparent_pka_is_the_right_one_and_the_true_one_would_overstate_the_charge(params):
    """APPARENT, and the word is load-bearing: the constant is referenced to TOTAL dissolved
    CO2, of which under 1 % is actually H2CO3. Our pool IS total dissolved CO2. Using the TRUE
    carbonic pKa (3.7) against that pool would overstate the anion charge by hundreds of
    times — the single most plausible wrong edit to this parameter, so its size is measured."""
    h = 10.0 ** (-4.9)  # a finished beer
    apparent = acidbase.mean_charge(h, (params["pKa_carbonic_1"],))
    true_h2co3 = acidbase.mean_charge(h, (3.7,))
    assert true_h2co3 / apparent == pytest.approx(32.7, rel=0.05), (
        f"substituting the TRUE carbonic pKa would multiply the anion charge by "
        f"{true_h2co3 / apparent:.1f}x at beer pH; D-182 measured 32.7x. The two constants "
        "differ by 2.7 pKa units because under 1 % of dissolved CO2 is actually H2CO3 — if "
        "this ratio has collapsed toward 1, the shipped value is no longer the apparent one."
    )


def test_a_wine_moves_but_barely_and_the_geometry_is_why(params, pka):
    """The measured cost of shipping this MEDIUM-AGNOSTIC rather than beer-scoped.

    Wine's pH sits ~3 units BELOW the apparent pKa, so under 0.1 % of its dissolved CO2 is
    dissociated, against ~4.5 % at beer's pH. That is the identical geometric argument D-178
    used to REJECT malt phosphate for beer, pointing the other way — which is what makes
    "medium-agnostic" a physics statement rather than a convenience.
    """
    totals = {"tartaric": 5.0 / M_TARTARIC, "malic": 3.0 / M_MALIC, "lactic": 0.0}
    sat_molar = acidbase.co2_saturation_gpl(293.15, params) / acidbase.CARBONIC_AS_CO2.molar_mass
    cation = acidbase.solve_cation_charge(totals, 0.0, 0.0, pka, 3.4)
    shift = acidbase.solve_ph(totals, cation, 0.0, 0.0, pka) - acidbase.solve_ph(
        totals, cation, 0.0, sat_molar, pka
    )
    assert 0.0 < shift < 0.005, (
        f"the carbonic term moved a wine by {shift:.4f} pH; D-182 measured ~0.0007. It is "
        "small BY GEOMETRY, and a large value here means either the pKa or the saturation "
        "has moved somewhere it should not have."
    )
    # The dissociated fractions that produce that asymmetry, stated as the check itself.
    frac_wine = acidbase.mean_charge(10.0 ** (-3.33), (params["pKa_carbonic_1"],))
    frac_beer = acidbase.mean_charge(10.0 ** (-4.96), (params["pKa_carbonic_1"],))
    assert frac_wine < 0.001 < 0.02 < frac_beer < 0.10


# -- 12. the assimilable-nitrogen pool's charge (decision D-209) ---------------
#
# The pool is not electrically neutral: it is ammonium plus amino acids whose side chains are
# charged at fermentation pH, so it sits on the CATION side and the yeast consuming it acidifies
# the liquid. These tests pin the arithmetic, the re-allocation that keeps t=0 anchoring exact,
# the opt-in gate, and — the one that matters most — the DERIVATION, so the provenance in
# ``acidbase.yaml`` is executable rather than merely written down.

#: Peyer 2017 Table 16 control column, mg/L, the 18 free amino acids of a malt wort. Transcribed
#: here so ``nitrogen_uptake_charge_beer`` can be re-derived from its own cited source rather than
#: trusted as a literal [[feedback-transcribe-tables-not-prose]]. Proline is absent from that
#: table, which is exactly right: it is Jones & Pierce Group D and brewing yeast does not
#: assimilate it, so the ``N`` pool (assimilable nitrogen by definition) excludes it.
_PEYER_WORT_AMINO_ACIDS_MGL = {
    "alanine": 36.9, "arginine": 47.6, "asparagine": 32.0, "aspartic": 27.5,
    "glutamic": 22.2, "glutamine": 41.5, "glycine": 11.3, "histidine": 22.0,
    "isoleucine": 23.4, "leucine": 50.7, "lysine": 30.2, "methionine": 10.2,
    "phenylalanine": 41.7, "serine": 23.6, "threonine": 20.1, "tryptophan": 14.0,
    "tyrosine": 30.8, "valine": 42.3,
}
#: ``(molar mass, nitrogen atoms, pKa_COOH, pKa_NH3, (side-chain pKa, sign) | None)``.
#: The nitrogen COUNT is the load-bearing column: ``zbar``'s denominator is ELEMENTAL nitrogen
#: because that is what the ``N`` slot holds, so arginine's +1 spreads over FOUR nitrogens and
#: contributes +0.25 per mole N. Getting that convention backwards inflates the cationic half
#: roughly fourfold, which is why it is asserted separately below.
_AMINO_ACID_CHEMISTRY = {
    "alanine": (89.09, 1, 2.34, 9.69, None),
    "arginine": (174.20, 4, 2.17, 9.04, (12.48, +1)),
    "asparagine": (132.12, 2, 2.02, 8.80, None),
    "aspartic": (133.10, 1, 1.99, 9.90, (3.90, -1)),
    "glutamic": (147.13, 1, 2.10, 9.47, (4.07, -1)),
    "glutamine": (146.15, 2, 2.17, 9.13, None),
    "glycine": (75.07, 1, 2.34, 9.60, None),
    "histidine": (155.16, 3, 1.82, 9.17, (6.04, +1)),
    "isoleucine": (131.17, 1, 2.36, 9.68, None),
    "leucine": (131.17, 1, 2.36, 9.60, None),
    "lysine": (146.19, 2, 2.18, 8.95, (10.53, +1)),
    "methionine": (149.21, 1, 2.28, 9.21, None),
    "phenylalanine": (165.19, 1, 1.83, 9.13, None),
    "serine": (105.09, 1, 2.21, 9.15, None),
    "threonine": (119.12, 1, 2.09, 9.10, None),
    "tryptophan": (204.23, 2, 2.83, 9.39, None),
    "tyrosine": (181.19, 1, 2.20, 9.11, (10.07, -1)),
    "valine": (117.15, 1, 2.32, 9.62, None),
}
_WORT_PH = 5.65
_NH4_PKA = 9.25
_PEYER_WORT_AMMONIUM_MG_N_PER_L = (25.0, 30.0)  # Peyer Table 2, read as mg N/L
_PEYER_DILUTION = 2.0  # Table 16's wort is CW0.5, diluted 50:50 with water


def _fraction_protonated(ph: float, pka: float) -> float:
    return float(1.0 / (1.0 + 10.0 ** (ph - pka)))


def _amino_acid_charge(ph: float, name: str) -> float:
    _, _, pka_cooh, pka_nh3, side = _AMINO_ACID_CHEMISTRY[name]
    charge = -(1.0 - _fraction_protonated(ph, pka_cooh))  # alpha-COOH: -1 once deprotonated
    charge += _fraction_protonated(ph, pka_nh3)  # alpha-NH3+: +1 while protonated
    if side is not None:
        pka, sign = side
        if sign > 0:
            charge += _fraction_protonated(ph, pka)
        else:
            charge -= 1.0 - _fraction_protonated(ph, pka)
    return charge


def _wort_amino_acid_pool(ph: float) -> tuple[float, float]:
    """``(mmol elemental N, mmol charge)`` per litre of FULL-STRENGTH Peyer wort."""
    nitrogen = charge = 0.0
    for name, mgl in _PEYER_WORT_AMINO_ACIDS_MGL.items():
        molar_mass, n_atoms = _AMINO_ACID_CHEMISTRY[name][0], _AMINO_ACID_CHEMISTRY[name][1]
        mmol = mgl / molar_mass * _PEYER_DILUTION
        nitrogen += mmol * n_atoms
        charge += mmol * _amino_acid_charge(ph, name)
    return nitrogen, charge


def test_the_beer_nitrogen_charge_is_reproduced_from_its_cited_composition(beer_params):
    """THE PROVENANCE, RE-DERIVED — the shipped number must fall out of its own sources.

    ``nitrogen_uptake_charge_beer`` is a DERIVED parameter, which is the kind most easily
    corrupted by a later edit: nothing about a lone float says which convention produced it.
    So this recomputes it from Peyer's Table 16 amino acids plus Table 2's ammonium and checks
    the shipped value and BOTH band edges land where the derivation puts them. It is also the
    guard on the arginine convention: swap the elemental-nitrogen denominator for an
    alpha-amino one and this fails by roughly the factor that error is worth.
    """
    nitrogen, charge = _wort_amino_acid_pool(_WORT_PH)
    # Units cross-check, and it is the reason Table 2's ammonium row is read as mg N/L: the
    # amino acids summed as ELEMENTAL nitrogen land inside Table 2's own printed 150-230 mg/L
    # for the free-amino-acid row, whereas as compound mass they would be ~1037 and far outside.
    assert 150.0 <= nitrogen * M_NITROGEN <= 230.0, (
        f"the transcribed composition gives {nitrogen * M_NITROGEN:.1f} mg N/L, outside Peyer "
        "Table 2's printed 150-230 for free amino acids — the two tables must be consistent or "
        "the ammonium row's units are not established either"
    )

    param = beer_params["nitrogen_uptake_charge_beer"]
    edges = []
    for ammonium_mg_n in _PEYER_WORT_AMMONIUM_MG_N_PER_L:
        n_ammonium = ammonium_mg_n / M_NITROGEN
        total_n = nitrogen + n_ammonium
        total_charge = charge + n_ammonium * _fraction_protonated(_WORT_PH, _NH4_PKA)
        edges.append(total_charge / total_n)
    assert min(edges) == pytest.approx(param.uncertainty.low, abs=5e-4)
    assert max(edges) == pytest.approx(param.uncertainty.high, abs=5e-4)
    assert param.value == pytest.approx(sum(edges) / 2.0, abs=5e-4), (
        "the nominal must be the derivation's own midpoint, not a fitted value"
    )

    # THE NEAR-CANCELLATION, asserted because it is what makes this parameter essentially a
    # measurement of wort AMMONIUM. The cationic amino acids (arginine, lysine, histidine) and
    # the anionic ones (aspartate, glutamate) very nearly cancel, so the amino-acid half is a
    # small residue and ~78 % of the shipped value is the ammonium term.
    assert charge / nitrogen == pytest.approx(0.0395, abs=5e-4), (
        f"the amino-acid-only charge is {charge / nitrogen:.4f} per mole N; D-209 measures "
        "+0.0395. If this has grown, check the arginine convention first: it carries FOUR "
        "nitrogens, so its +1 is worth +0.25 per mole N and not +1"
    )
    assert (charge / nitrogen) / param.value < 0.25, (
        "the amino-acid half must stay the minor term — the ammonium share is what this "
        "parameter mostly measures, and that is the single-source exposure the record names"
    )


def test_arginine_carries_four_nitrogens_per_unit_charge():
    """The convention, isolated — a one-line guard on the easiest fourfold error here.

    Named separately from the derivation above because a reader changing the denominator would
    otherwise see only a band-edge mismatch, with nothing saying which convention is right.
    """
    molar_mass, n_atoms = _AMINO_ACID_CHEMISTRY["arginine"][:2]
    assert n_atoms == 4
    charge = _amino_acid_charge(_WORT_PH, "arginine")
    assert charge == pytest.approx(1.0, abs=0.01), "arginine is +1 at wort pH (guanidinium 12.48)"
    assert charge / n_atoms == pytest.approx(0.25, abs=0.01), (
        "per mole of ELEMENTAL nitrogen — which is what the N slot holds — arginine is +0.25"
    )
    assert molar_mass == pytest.approx(174.2, abs=0.1)


@pytest.mark.parametrize(
    ("medium", "initial"),
    [
        (
            "wine",
            {"brix": 22.0, "yan_mgl": 200.0, "pitch_gpl": 0.25,
             "tartaric_gpl": 6.0, "malic_gpl": 3.0, "initial_ph": 3.4},
        ),
        (
            "beer",
            {"glucose_gpl": 12.0, "maltose_gpl": 57.0, "maltotriose_gpl": 13.0,
             "yan_mgl": 200.0, "pitch_gpl": 1.0, "initial_ph": 5.65},
        ),
    ],
)
def test_the_nitrogen_term_is_a_reallocation_so_the_anchor_is_untouched(medium, initial):
    """t=0 IS UNCHANGED BY D-209, and this is the assert that would catch it if it were not.

    ``cation_charge`` is back-solved from ``initial_ph``, so before D-209 it already contained
    the nitrogen pool's charge — lumped and frozen as a constant. Making the pool explicit moves
    that share out of the slot and onto the ``N`` state. Two consequences, both checked here for
    BOTH media because the subtraction is duplicated at three sites (each compile-seam anchor and
    ``cation_charge_for_ph``) and a sign error at any of them would silently anchor every must
    and wort to the wrong pH rather than raise:

    * the slot FALLS by exactly the nitrogen charge, so the cation SIDE is unchanged;
    * the anchored pH still reproduces ``initial_ph`` to the closed form's own precision.

    ``titratable_acidity`` is checked too: it solves its own starting pH off the same cation, so
    a must's t=0 TA — the fidelity-grade value, per that function's caveat — must not move.
    """
    compiled = compile_scenario(
        Scenario(
            name=f"d209-anchor-{medium}",
            medium=medium,
            initial=dict(initial),
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=1.0,
        )
    )
    schema, params = compiled.schema, compiled.param_values
    y0 = np.asarray(compiled.y0, dtype=float)

    nitrogen_charge = acidbase.nitrogen_charge_molar(y0, schema, params)
    assert nitrogen_charge > 0.0, "a nitrogen-bearing anchored state must carry positive charge"

    assert acidbase.ph_of_state(y0, schema, params) == pytest.approx(
        initial["initial_ph"], abs=1e-9
    ), "the anchor must still reproduce initial_ph — the D-209 subtraction is what makes it so"

    # The cation SIDE, reconstructed: slot + nitrogen must equal what the acids alone demand.
    slot = float(y0[schema.slice("cation_charge")][0])
    demanded = acidbase.solve_cation_charge(
        acidbase._totals_molar(y0, schema),
        acidbase._byp_succinic_molar(y0, schema),
        acidbase.dissolved_co2_molar(y0, schema, params),
        acidbase.build_pka_map(params),
        float(initial["initial_ph"]),
    )
    assert slot + nitrogen_charge == pytest.approx(demanded, rel=1e-12), (
        "the slot plus the nitrogen charge must be the total the target pH demands; if these "
        "diverge the term has been double-counted or dropped at one of its three sites"
    )
    assert slot < demanded, "the slot must have FALLEN — that is what re-allocation means"

    # And the round trip through the public state-level inverse, which must return the SLOT.
    assert acidbase.cation_charge_for_ph(
        y0, schema, params, float(initial["initial_ph"])
    ) == pytest.approx(slot, rel=1e-12)


def test_an_unanchored_beer_gets_no_nitrogen_charge(params):
    """THE GATE, and it is load-bearing rather than defensive.

    An un-anchored beer's charge balance is empty by construction (``_beer_acids`` seeds every
    acid slot from ``initial_ph`` or not at all, D-179) but its ``N`` slot still holds ~200 mg/L.
    Ungated, the nitrogen term would hand that empty balance ~0.0025 mol/L of cation charge with
    no acid to meet it, and ``solve_ph`` would answer around 11 where such a state answered 7.0 —
    a strong-base artefact in exactly the place D-179's gate exists to protect, and one that
    reaches the aging trajectory through ``EsterHydrolysis``. So the term rides the same opt-in.
    """
    schema = beer_schema()
    y = schema.zeros()
    y[schema.slice("N")] = 0.2
    y[schema.slice("T")] = 293.15
    assert not acidbase.charge_balance_is_populated(y, schema)
    assert acidbase.nitrogen_charge_molar(y, schema, params) == 0.0
    assert acidbase.ph_of_state(y, schema, params) == pytest.approx(7.0, abs=1e-6), (
        "an empty balance must still read 7.0 — nitrogen alone is not pH information"
    )

    # ...and with the balance populated by a single acid, the term switches on.
    y[schema.slice("lactic")] = 0.3
    assert acidbase.charge_balance_is_populated(y, schema)
    assert acidbase.nitrogen_charge_molar(y, schema, params) > 0.0


def test_a_nitrogen_free_state_is_bitwise_the_pre_d209_balance(params):
    """Inertness where it can be claimed exactly: no nitrogen, no term.

    The honest inertness check for this parameter, because a state with ``N`` = 0 makes the term
    identically zero by construction rather than merely small — and unlike a band-edge screen
    that is a real algebraic zero [[feedback-nominal-on-a-band-edge-is-not-inertness]].
    """
    schema = wine_schema()
    cation = _anchor_cation(acidbase.build_pka_map(params), 6.0, 3.0, 3.4)
    y = _wine_state(schema, tartaric=6.0, malic=3.0, lactic=0.0, cation_charge=cation, N=0.0)
    assert acidbase.nitrogen_charge_molar(y, schema, params) == 0.0
    assert acidbase.ph_of_state(y, schema, params) == pytest.approx(3.4, abs=1e-9), (
        "with no nitrogen the balance is the pre-D-209 one, so the un-subtracted anchor is right"
    )
