"""The JOINT fusel repair, measured for the first time (decision D-266).

D-257 built the temporal fusel repair — the Ehrlich rate blended with a de-novo term,
``flux · [f·N/(K_n+N) + (1−f)] · arrhenius``, ``f = 0.79``, the five ``k`` re-anchored ×0.4033,
the amino-acid attribution scaled by the Ehrlich share of the rate — and REVERTED it, because the
shipped precursor sink rides the alcohol rate and slowing production stranded 72.9 % of the
must's phenylalanine. D-259/D-260 measured the growth-anchored sink that clears exactly that
blocker and REFUSED it, because under the shipped rate law Crépin's leucine split and Rollero's
leucine tracer trade one-for-one along ``tracer = (1 − split) · consumed / isoamyl``.

Each refusal's ground is the other repair's defect, and nobody had run the pair. This module
does: four arms — ``shipped``, ``blend`` (D-257 alone), ``growth`` (D-259/D-260 alone), ``joint``
— on every sourced fusel observable the suite already transcribes. **Nothing in ``src/``
changes**; every arm is a test-side counterfactual, and the ``blend`` and ``growth`` arms must
first reproduce D-257's and D-260's own published numbers before any ``joint`` number is read
(:func:`test_the_harness_reproduces_d257_and_d260_before_it_measures_anything`).

Receipts, in the repo this time: ``docs/receipts/d266-joint-fusel-repair/`` (the
pre-registration written before the first number, the probe, its JSON, the findings).
"""

from __future__ import annotations

import contextlib
from collections.abc import Iterator, Mapping
from typing import Any

import numpy as np
import pytest

from fermentation.core.chemistry import (
    M_ISOAMYL_OH,
    MOLAR_MASS,
    carbon_mass_fraction,
    sugar_species,
)
from fermentation.core.kinetics import byproducts as byproducts_module
from fermentation.core.kinetics import precursor_fates as precursor_fates_module
from fermentation.core.kinetics.amino_acid_pools import SPEC_BY_SPECIES, depletion_gate
from fermentation.core.kinetics.amino_acids import AminoAcidAssimilation
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.kinetics.carbon_routing import (
    FUSEL_SPECS,
    fermentative_flux_shape,
    non_ehrlich_fraction_param,
)
from fermentation.core.kinetics.growth import (
    GrowthNitrogenLimited,
    assimilable_nitrogen_pools,
    biomass_growth_rate,
)
from fermentation.core.kinetics.precursor_fates import PrecursorNonEhrlichFates
from fermentation.core.process import RateModifier
from fermentation.core.state import FloatArray, StateSchema
from fermentation.runtime import simulate_scheduled
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.validation import total_carbon, total_nitrogen
from tests.test_defined_media import (
    ROLLERO_SUGAR_GPL,
    ROLLERO_TEMP_C,
    SOURCED_PITCH_GPL,
    _brix_for,
    commensurate_pools,
    commensurate_scenario,
    rollero_must_mm,
)
from tests.test_fusel_keto_acid_node import _OTHER_PRECURSOR_CONSUMERS
from tests.test_fusel_reroute import (
    _ROLLERO_EF_ISOAMYL_UM,
    _ROLLERO_EF_LEUCINE_ENRICHMENT,
    _ROLLERO_NT_ISOAMYL_UM,
)
from tests.test_precursor_fates import (
    _CREPIN_PROTEIN_PCT,
    _D257_LEVEL_PRESERVING_K_FACTOR,
    _D259_EDGES,
    _ROLLERO_LEUCINE_SHARE_OF_ISOAMYL_PCT,
    _d259_weights,
    _GrowthAnchoredFates,
)

#: D-257 §1: the catabolic share, fitted on Rollero's SM250 — the ONE value D-257 measured.
_D257_CATABOLIC_SHARE = 0.79
#: The sourced propanol de-novo floor D-244 asserted and D-248 cleared (Crépin 81 % newly
#: synthesised 2-KB; Rollero > 90 % CCM). Read from ``test_fusel_catabolic_shape``'s own band.
_PROPANOL_DE_NOVO_FLOOR = 0.80
#: Rollero Table S1's isoamyl-alcohol valine enrichment at EF, ``{yan: (2 mg/l, 8 mg/l)}`` —
#: the same transcription ``test_fusel_reroute._ROLLERO_EF_ENRICHMENT`` carries, as percent.
_ROLLERO_EF_VALINE_ENRICHMENT_PCT = {70.0: (2.1, 2.3), 250.0: (4.0, 3.4), 425.0: (5.3, 5.4)}
#: Rollero Table S2's leucine-labelled isoamyl at EF, µM, ``{yan: (2 mg/l, 8 mg/l)}`` — the
#: NUMERATOR of the enrichment, transcribed at D-256 §3, so an enrichment can be scored apart
#: from the total it is divided by.
_ROLLERO_EF_LEUCINE_LABELLED_UM = {70.0: (37.3, 45.5), 250.0: (55.2, 64.2), 425.0: (65.0, 70.3)}
_YANS = (70.0, 250.0, 425.0)

_ORIG_SHAPE = byproducts_module.fusel_rate_shape
_ORIG_DRAWS = byproducts_module.ehrlich_draws
_DRAWS_BINDING = "ehrlich_draws"  # precursor_fates imports it by value; patched by name
_LABEL_BINDING = "LABELLED_PRECURSOR"  # not re-exported; reached by name, as D-256 did
_ORIG_LABEL = getattr(byproducts_module, _LABEL_BINDING)


# ---------------------------------------------------------------------------------------------
# The blended rate law, as a patch on the bindings the code actually reads.
# ---------------------------------------------------------------------------------------------
def _terms(
    y: FloatArray, schema: StateSchema, params: Mapping[str, float]
) -> tuple[float, float, float]:
    flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
    if flux <= 0.0:
        return 0.0, 0.0, 0.0
    n = sum(assimilable_nitrogen_pools(y, schema))
    gate = n / (params["K_n"] + n) if n > 0.0 else 0.0
    temp = float(y[schema.slice("T")][0])
    return flux, gate, arrhenius_factor(temp, params["E_a_fusels"], params["T_ref"])


class _Blend:
    """D-257's rate law and its attribution scaling, installed on the bindings the code reads.

    ``byproducts.fusel_rate_shape`` is what the producer's pool rates and its sugar draw read
    (through ``fusel_carbon_draw_by_species``), so it becomes the blend. ``ehrlich_draws`` is
    what the re-route and the D-104 sink read — in TWO modules, because ``precursor_fates``
    imports it by value ([[feedback-patch-the-binding-the-code-reads]]) — and it becomes the
    original helper evaluated under the Ehrlich-only shape, which is D-257 §4's "attribution
    scaled by the Ehrlich share of the current rate": an identity at ``f = 1``.
    """

    def __init__(self, f: float) -> None:
        self.f = f

    def blended(
        self,
        y: FloatArray,
        schema: StateSchema,
        params: Mapping[str, float],
        *,
        growth_coupled: bool = False,
    ) -> float:
        if growth_coupled:
            return _ORIG_SHAPE(y, schema, params, growth_coupled=True)
        flux, gate, f_t = _terms(y, schema, params)
        if flux <= 0.0:
            return 0.0
        return float(flux * (self.f * gate + (1.0 - self.f)) * f_t)

    def ehrlich_only(
        self,
        y: FloatArray,
        schema: StateSchema,
        params: Mapping[str, float],
        *,
        growth_coupled: bool = False,
    ) -> float:
        if growth_coupled:
            return _ORIG_SHAPE(y, schema, params, growth_coupled=True)
        flux, gate, f_t = _terms(y, schema, params)
        if flux <= 0.0 or gate <= 0.0:
            return 0.0
        return float(flux * self.f * gate * f_t)

    def scaled_draws(self, y: FloatArray, schema: StateSchema, params: Mapping[str, float]) -> Any:
        byproducts_module.fusel_rate_shape = self.ehrlich_only
        try:
            return _ORIG_DRAWS(y, schema, params)
        finally:
            byproducts_module.fusel_rate_shape = self.blended


@contextlib.contextmanager
def _rate_law(blend: bool, f: float = _D257_CATABOLIC_SHARE) -> Iterator[None]:
    if not blend:
        yield
        return
    patch = _Blend(f)
    byproducts_module.fusel_rate_shape = patch.blended
    byproducts_module.ehrlich_draws = patch.scaled_draws
    setattr(precursor_fates_module, _DRAWS_BINDING, patch.scaled_draws)
    try:
        yield
    finally:
        byproducts_module.fusel_rate_shape = _ORIG_SHAPE
        byproducts_module.ehrlich_draws = _ORIG_DRAWS
        setattr(precursor_fates_module, _DRAWS_BINDING, _ORIG_DRAWS)


def _rescale_k(params: dict[str, float], blend: bool, k_factor: float) -> None:
    if blend:
        for spec in FUSEL_SPECS:
            params[spec.k_param] = params[spec.k_param] * k_factor


def _install_growth_sink(cs: Any, edge: str) -> tuple[_GrowthAnchoredFates, list[RateModifier]]:
    """D-259's counterfactual, with growth's own modifiers attached (D-260 §5's correction)."""
    sink = _GrowthAnchoredFates(_d259_weights(edge))
    cs.process_set._processes[sink.name] = sink
    attached: list[RateModifier] = []
    for modifier in cs.process_set.active_modifiers:
        if GrowthNitrogenLimited.name in modifier.modifies:
            modifier.modifies = (*modifier.modifies, sink.name)
            attached.append(modifier)
    assert attached, "no modifier scales growth: the D-32 attachment measures nothing"
    return sink, attached


def _integrate(cs: Any, params: dict[str, float], points: int) -> Any:
    traj = simulate_scheduled(
        cs.process_set,
        params,
        cs.y0,
        cs.t_span_h,
        events=cs.events,
        param_tiers=cs.parameters.tier_map(),
        t_eval=np.linspace(0.0, cs.t_span_h[1], points),
    )
    assert traj.success, traj.message
    return traj


_ARMS: dict[str, tuple[bool, bool]] = {  # name -> (blend, growth sink)
    "shipped": (False, False),
    "blend": (True, False),
    "growth": (False, True),
    "joint": (True, True),
}


# ---------------------------------------------------------------------------------------------
# Crépin's must — the D-259/D-260 fixture, all five precursors, both draws reconstructed.
# ---------------------------------------------------------------------------------------------
def _crepin_arm(
    blend: bool,
    growth: bool,
    *,
    edge: str = "mid",
    f: float = _D257_CATABOLIC_SHARE,
    k_factor: float = _D257_LEVEL_PRESERVING_K_FACTOR,
    points: int = 4001,
) -> dict[str, Any]:
    """Splits, the leucine tracer, isoamyl, what is left of each pool, de-novo shares.

    Both draws are reconstructed from their own formulae along the trajectory and tied to the
    depletion the solver realised by a per-species closure, which the fixture asserts — D-259
    §7's discipline, inherited unchanged.
    """
    with _rate_law(blend, f):
        cs = compile_scenario(commensurate_scenario("crepin", days=14.0))
        for name in _OTHER_PRECURSOR_CONSUMERS:
            cs.process_set.disable(name)
        params = cs.param_values
        _rescale_k(params, blend, k_factor)
        sink: _GrowthAnchoredFates | None = None
        attached: list[RateModifier] = []
        if growth:
            sink, attached = _install_growth_sink(cs, edge)
        traj = _integrate(cs, params, points)
        y, t, schema = np.asarray(traj.y, float), np.asarray(traj.t, float), cs.schema
        species = ("leucine", "isoleucine", "valine", "threonine", "phenylalanine")
        lump = {s: np.zeros_like(t) for s in species}
        ehrlich = {s: np.zeros_like(t) for s in species}
        sourced = {spec.pool: np.zeros_like(t) for spec in FUSEL_SPECS}
        draws_fn = byproducts_module.ehrlich_draws  # the binding the run used
        for i in range(t.size):
            column = y[:, i]
            for draw in draws_fn(column, schema, params):
                sp = draw.precursor.species
                if sp in ehrlich:
                    ehrlich[sp][i] += draw.precursor_carbon / carbon_mass_fraction(sp)
                sourced[draw.alcohol.pool][i] += draw.alcohol_carbon
            if sink is not None:
                base_dx = max(biomass_growth_rate(column, schema, params), 0.0)
                factor = 1.0
                for m in attached:
                    factor *= m.factor(float(t[i]), column, schema, params)
                for s in species:
                    gate = depletion_gate(column, schema, params, (SPEC_BY_SPECIES[s],))
                    lump[s][i] = sink.weights[s] * base_dx * gate * factor
            else:
                for s in species:
                    share = params[non_ehrlich_fraction_param(s)]
                    lump[s][i] = ehrlich[s][i] * share / (1.0 - share)
    out: dict[str, Any] = {"split_pct": {}, "closure": {}, "left_pct": {}, "de_novo": {}}
    for s in species:
        mw = MOLAR_MASS[s]
        lump_um = float(np.trapezoid(lump[s], t)) / mw * 1e6
        ehrlich_um = float(np.trapezoid(ehrlich[s], t)) / mw * 1e6
        pool = y[schema.slice(s), :][0]
        consumed_um = (float(pool[0]) - float(pool[-1])) / mw * 1e6
        assert consumed_um > 0.0, f"vacuous: the must's {s} was not consumed"
        out["split_pct"][s] = 100.0 * lump_um / (lump_um + ehrlich_um)
        out["closure"][s] = (lump_um + ehrlich_um) / consumed_um
        out["left_pct"][s] = 100.0 * float(pool[-1]) / float(pool[0])
        if s == "leucine":
            out["leucine_supply_um"] = float(pool[0]) / mw * 1e6
            out["leucine_consumed_um"] = consumed_um
            out["leucine_ehrlich_um"] = ehrlich_um
    isoamyl_um = float(y[schema.slice("isoamyl_alcohol"), -1][0]) / M_ISOAMYL_OH * 1e6
    assert isoamyl_um > 0.0, "vacuous: no isoamyl alcohol was made"
    out["isoamyl_um"] = isoamyl_um
    out["leucine_tracer_pct"] = 100.0 * out["leucine_ehrlich_um"] / isoamyl_um
    for spec in FUSEL_SPECS:
        made = float(y[schema.slice(spec.pool), -1][0]) * carbon_mass_fraction(spec.species)
        assert made > 0.0, f"vacuous: no {spec.pool} was made"
        out["de_novo"][spec.pool] = 1.0 - float(np.trapezoid(sourced[spec.pool], t)) / made
    out["biomass_peak_gpl"] = float(y[schema.slice("X"), :].max())
    return out


@pytest.fixture(scope="module")
def crepin():
    """The four arms plus the joint arm's composition bracket — computed ONCE, closure asserted
    HERE so every consumer inherits the check that ties a split to the run that produced it."""
    arms = {name: _crepin_arm(*cfg) for name, cfg in _ARMS.items()}
    for edge in _D259_EDGES:
        arms[f"joint_{edge}"] = _crepin_arm(True, True, edge=edge)
    for name, arm in arms.items():
        for s, closure in arm["closure"].items():
            assert closure == pytest.approx(1.0, abs=0.02), (
                f"quadrature control failed for {s} in arm {name!r}: the two integrated draws "
                f"account for {closure:.4f} of the depletion the solver realised. Either the "
                "Processes that ran are not the ones the arm names, or the grid is too coarse "
                "(D-103's trapezoid trap). Nothing read off this arm is trustworthy until it "
                "closes"
            )
    return arms


# ---------------------------------------------------------------------------------------------
# Rollero's musts — the D-255 fixture, full Process set, both tracer labels.
# ---------------------------------------------------------------------------------------------
def _rollero_scenario(yan: float) -> Scenario:
    pools, model_yan = commensurate_pools(rollero_must_mm(yan))
    initial = {
        "brix": _brix_for(ROLLERO_SUGAR_GPL),
        "yan_mgl": model_yan,
        "pitch_gpl": SOURCED_PITCH_GPL,
        "amino_acids_gpl": 1.0,
    } | pools
    return Scenario(
        name=f"d266-rollero-{yan:.0f}",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=ROLLERO_TEMP_C)],
        duration_days=14.0,
    )


def _gate_closure_index(traj: Any, schema: StateSchema, params: Mapping[str, float]) -> int:
    """``test_fusel_reroute._fusel_nitrogen_gate_closure``'s landmark, on the gate's own
    quantity (extracellular N plus the store, no amino-acid pool)."""
    for i in range(traj.y.shape[1]):
        n = sum(assimilable_nitrogen_pools(traj.y[:, i], schema))
        if n / (params["K_n"] + n) < 0.01:
            return i
    raise AssertionError("vacuous: the fusel nitrogen gate never closes over this run")


def _rollero_arm(blend: bool, growth: bool, yan: float, label: str) -> dict[str, Any]:
    with _rate_law(blend):
        setattr(byproducts_module, _LABEL_BINDING, label)
        try:
            cs = compile_scenario(_rollero_scenario(yan), strict=True)
            params = cs.param_values
            _rescale_k(params, blend, _D257_LEVEL_PRESERVING_K_FACTOR)
            if growth:
                _install_growth_sink(cs, "mid")
            traj = _integrate(cs, params, int(cs.t_span_h[1]) + 1)
        finally:
            setattr(byproducts_module, _LABEL_BINDING, _ORIG_LABEL)
    schema = cs.schema
    i = _gate_closure_index(traj, schema, params)
    iso = traj.y[schema.slice("isoamyl_alcohol")][0]
    assert iso[-1] > 0.0, "vacuous: no isoamyl alcohol was made"
    labelled = float(traj.y[schema.slice("isoamyl_alcohol_valine"), -1][0])
    return {
        "nt_fraction_pct": 100.0 * float(iso[i] / iso[-1]),
        "isoamyl_ef_um": float(iso[-1]) / M_ISOAMYL_OH * 1e6,
        "enrichment_pct": 100.0 * labelled / float(iso[-1]),
        "labelled_um": labelled / M_ISOAMYL_OH * 1e6,
    }


@pytest.fixture(scope="module")
def rollero():
    """``rollero[arm][yan]`` with the valine tracer's numbers, plus the leucine arm's
    enrichment reached the D-256 way (re-pointing the module global) — with the arms-differ
    check that catches a patch that never reached the sourcing layer."""
    out: dict[str, dict[float, dict[str, Any]]] = {}
    for name, (blend, growth) in _ARMS.items():
        out[name] = {}
        for yan in _YANS:
            valine = _rollero_arm(blend, growth, yan, "valine")
            leucine = _rollero_arm(blend, growth, yan, "leucine")
            assert abs(leucine["enrichment_pct"] - valine["enrichment_pct"]) > (
                0.05 * valine["enrichment_pct"]
            ), (
                f"{name} SM{yan:.0f}: the leucine and valine label arms read the same number, "
                "so the label patch did not reach the sourcing layer (D-256's anti-vacuity)"
            )
            valine["leucine_enrichment_pct"] = leucine["enrichment_pct"]
            valine["leucine_labelled_um"] = leucine["labelled_um"]
            out[name][yan] = valine
    return out


# ---------------------------------------------------------------------------------------------
# The D-104 fixture (full set) and the D-112 anchor.
# ---------------------------------------------------------------------------------------------
def _d104_arm(blend: bool, growth: bool) -> dict[str, Any]:
    """Phenylalanine, the D-117 refund guard, net dS/dt, and both ledgers, on D-257's own
    fixture — brix 24, YAN 250, pitch 0.25, 1 g/L amino acids, 20 °C, 14 d."""
    with _rate_law(blend):
        scenario = Scenario(
            name="d266-d104-fixture",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25, "amino_acids_gpl": 1.0},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=14.0,
        )
        cs = compile_scenario(scenario, strict=True)
        params = cs.param_values
        _rescale_k(params, blend, _D257_LEVEL_PRESERVING_K_FACTOR)
        sink: _GrowthAnchoredFates | None = None
        if growth:
            sink, _ = _install_growth_sink(cs, "mid")
        traj = _integrate(cs, params, 337)
        schema, ps = cs.schema, cs.process_set
        sink_proc = sink if sink is not None else ps._processes[PrecursorNonEhrlichFates.name]
        swap = ps._processes[AminoAcidAssimilation.name]
        s_slice = schema.slice("S")
        worst_refund, worst_ds = 0.0, -np.inf
        for i in range(traj.y.shape[1]):
            yi, ti = traj.y[:, i], float(traj.t[i])
            worst_ds = max(worst_ds, float(ps.total_derivatives(ti, yi, params)[s_slice].sum()))
            base_dx = biomass_growth_rate(yi, schema, params)
            if base_dx <= 1e-6:
                continue
            carbon = 0.0
            for proc in (swap, sink_proc):
                if not ps.is_enabled(proc.name):
                    continue
                d = proc.derivatives(ti, yi, schema, params)
                for offset, sp in enumerate(sugar_species(schema)):
                    carbon += float(d[s_slice.start + offset]) * carbon_mass_fraction(sp)
            worst_refund = max(worst_refund, carbon / (params["biomass_C_fraction"] * base_dx))
    c_fn = total_carbon(schema, biomass_carbon_fraction=params["biomass_C_fraction"])
    n_fn = total_nitrogen(schema, biomass_nitrogen_fraction=params["biomass_N_fraction"])
    c = np.array([c_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])
    n = np.array([n_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])
    phe = traj.y[schema.slice("phenylalanine")][0]
    assert float(phe[0]) > 0.0, "vacuous: the must carries no phenylalanine"
    return {
        "phe_left_pct": 100.0 * float(phe[-1] / phe[0]),
        "worst_joint_c_refund": worst_refund,
        "max_net_ds_dt": worst_ds,
        "carbon_rel_drift": float(np.max(np.abs(c - c[0])) / c[0]),
        "nitrogen_rel_drift": float(np.max(np.abs(n - n[0])) / n[0]),
        "biomass_final_gpl": float(traj.y[schema.slice("X"), -1][0]),
    }


@pytest.fixture(scope="module")
def d104():
    return {name: _d104_arm(*cfg) for name, cfg in _ARMS.items()}


def _anchor_isoamyl_mgl(blend: bool, yan: float) -> float:
    """``test_fusel_catabolic_shape._isoamyl_no_dose`` — D-112's undosed anchor, 20 °C."""
    with _rate_law(blend):
        scenario = Scenario(
            name="d112-anchor",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": yan, "pitch_gpl": 0.25, "amino_acids_gpl": 0.0},
            temperature_schedule=[
                TemperaturePoint(day=0.0, celsius=20.0),
                TemperaturePoint(day=14.0, celsius=20.0),
            ],
            duration_days=14.0,
        )
        cs = compile_scenario(scenario)
        params = cs.param_values
        _rescale_k(params, blend, _D257_LEVEL_PRESERVING_K_FACTOR)
        traj = _integrate(cs, params, 337)
    return float(traj.y[cs.schema.slice("isoamyl_alcohol"), -1][0]) * 1e3


# ---------------------------------------------------------------------------------------------
# The guards.
# ---------------------------------------------------------------------------------------------
def test_the_harness_reproduces_d257_and_d260_before_it_measures_anything(crepin, rollero, d104):
    """P1, and it is the licence for every other number in this module.

    The ``blend`` arm must read D-257's own numbers — 36.3 / 48.1 / 55.9 % of the isoamyl
    present at nitrogen-gate closure on Rollero's three musts, 172.3 mg/L at the D-112 anchor,
    72.9 % of phenylalanine stranded on the D-104 fixture — and the ``growth`` arm must read
    D-260's: leucine split 27.58 %, tracer 5.900 % on Crépin's must. The ``shipped`` arm must
    read what D-256/D-260 recorded for it. If any of these has moved, the patch bindings or the
    counterfactual are not the ones the archive measured, and nothing below is attributable.
    ([[feedback-reproduce-a-published-number-before-trusting-the-new-column]])
    """
    for yan, recorded in zip(_YANS, (36.3, 48.1, 55.9), strict=True):
        got = rollero["blend"][yan]["nt_fraction_pct"]
        assert got == pytest.approx(recorded, abs=1.5), (
            f"blend SM{yan:.0f}: {got:.1f} % of isoamyl at gate closure against D-257's "
            f"{recorded} %. The blended rate law here is not the one D-257 measured"
        )
        assert rollero["shipped"][yan]["nt_fraction_pct"] > 99.0, (
            "the shipped arm no longer makes all its isoamyl before the gate shuts (D-256)"
        )
    assert _anchor_isoamyl_mgl(True, 250.0) == pytest.approx(172.3, abs=1.0), (
        "the blend does not hold the Wang anchor at 172.3 mg/L; D-257's ×0.4033 rescale is "
        "no longer level-preserving on this tree"
    )
    assert d104["blend"]["phe_left_pct"] == pytest.approx(72.9, abs=2.0), (
        f"blend strands {d104['blend']['phe_left_pct']:.1f} % of phenylalanine against "
        "D-257's 72.9 % — the blocker this module measures the clearing of has moved"
    )
    assert d104["shipped"]["phe_left_pct"] == pytest.approx(20.3, abs=2.0)

    growth, shipped = crepin["growth"], crepin["shipped"]
    assert growth["split_pct"]["leucine"] == pytest.approx(27.58, abs=0.3), (
        f"growth-anchored leucine split {growth['split_pct']['leucine']:.2f} % against "
        "D-260's 27.58 %: the counterfactual is not D-260's (are growth's modifiers attached?)"
    )
    assert growth["leucine_tracer_pct"] == pytest.approx(5.900, abs=0.05)
    assert shipped["split_pct"]["leucine"] == pytest.approx(81.5, abs=0.1)
    assert shipped["leucine_tracer_pct"] == pytest.approx(1.507, abs=0.02)
    assert shipped["isoamyl_um"] == pytest.approx(2123.3, abs=5.0)

    # The four arms are four different models: a patch that silently did not apply would make
    # two of them read alike. The shipped and blend arms share the leucine split BY
    # CONSTRUCTION (the shipped sink imposes `f`), so they are told apart on the axis the blend
    # moves — how much phenylalanine it strands — and the sink swap on the split itself.
    assert crepin["blend"]["left_pct"]["phenylalanine"] > (
        crepin["shipped"]["left_pct"]["phenylalanine"] + 40.0
    ), "the blend arm reads like the shipped arm: the rate-law patch did not reach the code"
    assert crepin["joint"]["split_pct"]["leucine"] > crepin["growth"]["split_pct"]["leucine"] + 20
    assert crepin["shipped"]["split_pct"]["leucine"] > crepin["growth"]["split_pct"]["leucine"] + 40


def test_the_blend_at_f_1_and_k_1_IS_the_shipped_rate_law(crepin):
    """The positive control on the patch itself: D-257's form is isolable by a parameter.

    At ``f = 1`` the blended shape is ``flux · gate · arrhenius`` and the Ehrlich-only shape is
    the same expression, so with the five ``k`` unscaled the ``joint`` machinery must reproduce
    the ``growth`` arm — the same counterfactual under the ORIGINAL rate law — to solver noise.
    This is what proves the patch reaches the code AND does nothing it was not asked to.
    """
    identity = _crepin_arm(True, True, f=1.0, k_factor=1.0)
    growth = crepin["growth"]
    for s in growth["split_pct"]:
        assert identity["split_pct"][s] == pytest.approx(growth["split_pct"][s], rel=1e-6), (
            f"{s}: the blend at f=1, k=1 reads {identity['split_pct'][s]:.6f} % against the "
            f"unpatched {growth['split_pct'][s]:.6f} % — the patch changes something at the "
            "identity, so every blended number carries an artefact"
        )
    assert identity["isoamyl_um"] == pytest.approx(growth["isoamyl_um"], rel=1e-6)
    assert identity["leucine_tracer_pct"] == pytest.approx(growth["leucine_tracer_pct"], rel=1e-6)


def test_the_joint_arm_CLEARS_the_blocker_d257_refused_the_blend_on(crepin, d104):
    """P2. The ``blend`` arm is the control: it is the arm that HAS the blocker.

    D-257 reverted the repair because the shipped sink consumes precursor at ``f/(1−f)`` times
    the Ehrlich draw, so slowing production strands the must's phenylalanine (72.9 % left).
    Growth-anchoring the sink removes that coupling, and the joint arm leaves essentially none —
    on D-257's own fixture and on Crépin's must. The anti-vacuity arm: the rest of the ferment
    does not move with it (biomass within 2 % across all four arms), so the change is the
    sink's and not the whole ferment's.
    """
    assert d104["blend"]["phe_left_pct"] > 55.0, "the control arm has lost the blocker"
    assert d104["joint"]["phe_left_pct"] < 0.05, (
        f"the joint arm leaves {d104['joint']['phe_left_pct']:.3f} % of the D-104 fixture's "
        "phenylalanine — D-257's blocker is NOT cleared by the pair, contrary to D-266"
    )
    assert crepin["joint"]["left_pct"]["phenylalanine"] < 0.05
    assert crepin["blend"]["left_pct"]["phenylalanine"] > 40.0, (
        "the blend-alone arm no longer strands phenylalanine on Crépin's must; the control "
        "this attribution rests on has moved"
    )
    biomass = [d104[name]["biomass_final_gpl"] for name in _ARMS]
    assert max(biomass) - min(biomass) < 0.02 * min(biomass), (
        f"final biomass spans {min(biomass):.4f}-{max(biomass):.4f} g/L across the arms; the "
        "phenylalanine change is then not attributable to the sink alone"
    )


def test_the_joint_arm_KEEPS_the_blends_temporal_gains(rollero):
    """P6 and P7: the sink swap does not touch production.

    Within 1 point of the blend on Rollero's three musts, and the residual over-response is
    pinned as what it is: SM250 sits inside his 44.5-51.9 %, SM70 below his 51-54 % and SM425
    above his 42-51 % — a model spread of 36 → 57 % across his nitrogen range where he measures
    a near-flat 42-54 %, the same shape D-257 §1 recorded for the blend alone. The isoamyl
    response across the range is 1.30× against his 0.74-0.77× (the shipped 2.27×).
    """
    measured = {
        yan: [
            nt / ef
            for nt, ef in zip(_ROLLERO_NT_ISOAMYL_UM[yan], _ROLLERO_EF_ISOAMYL_UM[yan], strict=True)
        ]
        for yan in _YANS
    }
    for yan in _YANS:
        joint, blend = (
            rollero["joint"][yan]["nt_fraction_pct"],
            rollero["blend"][yan]["nt_fraction_pct"],
        )
        assert joint == pytest.approx(blend, abs=1.0), (
            f"SM{yan:.0f}: the joint arm's {joint:.1f} % at gate closure is more than a point "
            f"from the blend's {blend:.1f} %; the growth sink is moving production"
        )
    lo250, hi250 = 100.0 * min(measured[250.0]), 100.0 * max(measured[250.0])
    assert lo250 < rollero["joint"][250.0]["nt_fraction_pct"] < hi250
    assert rollero["joint"][70.0]["nt_fraction_pct"] < 100.0 * min(measured[70.0])
    assert rollero["joint"][425.0]["nt_fraction_pct"] > 100.0 * max(measured[425.0])
    response = rollero["joint"][425.0]["isoamyl_ef_um"] / rollero["joint"][70.0]["isoamyl_ef_um"]
    assert response == pytest.approx(1.30, abs=0.05), (
        f"joint isoamyl response across Rollero's range is {response:.2f}×; D-257 measured "
        "1.30× for the blend and the sink swap is not supposed to move it"
    )
    shipped = rollero["shipped"][425.0]["isoamyl_ef_um"] / rollero["shipped"][70.0]["isoamyl_ef_um"]
    assert shipped > 2.0, "the shipped over-response (D-256's 2.27×) has closed; re-measure"


def test_the_blend_moves_the_point_along_d260s_line_and_does_not_move_the_line(crepin):
    """P3, P4 and P5 together — the beat's finding, and it reverses the pre-registered one.

    D-260's identity ``tracer = (1 − split) · consumed / isoamyl`` binds because the leucine
    pool empties whatever draws it and isoamyl is inert under a sink change. The blend cuts the
    Ehrlich draw (×0.79 gated, ×0.4033 on ``k``), so under the growth sink the lump absorbs the
    residue and the leucine split rises 27.6 → 54.5 % (46.0-61.5 across the composition
    bracket). It was pre-registered that the tracer would then fall below Rollero's floor. It
    does not: 3.62 % (3.06-4.29), inside 3.4-8.2 %, because the split stops well short of
    Crépin's 77 %. And isoamyl on this must is 2123 → 2177 µM (+2.5 %) — the blend is
    level-preserved at the Wang anchor by construction — so the LINE is the same line and the
    joint-satisfaction ceiling (≤ 1170 µM at 77 %) is exactly where D-260 left it. **The joint
    arm is a different point on D-260's line, not an escape from it.**
    """
    shipped = crepin["shipped"]
    crepin_lo, _ = _CREPIN_PROTEIN_PCT["leucine"]
    rollero_lo, rollero_hi = _ROLLERO_LEUCINE_SHARE_OF_ISOAMYL_PCT
    for name in _ARMS:
        arm = crepin[name]
        assert arm["leucine_consumed_um"] == pytest.approx(arm["leucine_supply_um"], rel=0.01), (
            f"{name}: leucine is no longer fully consumed; the identity below rests on that"
        )
        assert arm["isoamyl_um"] == pytest.approx(shipped["isoamyl_um"], rel=0.03), (
            f"{name}: isoamyl {arm['isoamyl_um']:.0f} µM against the shipped "
            f"{shipped['isoamyl_um']:.0f}. The denominator is supposed to be inert across all "
            "four arms; if it moves, the line moves and D-260's collision is re-openable"
        )
        predicted = (1.0 - arm["split_pct"]["leucine"] / 100.0) * arm["leucine_consumed_um"]
        predicted = 100.0 * predicted / arm["isoamyl_um"]
        assert predicted == pytest.approx(arm["leucine_tracer_pct"], rel=1e-4), (
            f"{name}: split and tracer are not two readings of one integral"
        )
        ceiling = arm["leucine_supply_um"] * (1.0 - crepin_lo / 100.0) / (rollero_lo / 100.0)
        assert arm["isoamyl_um"] > ceiling, (
            f"{name}: isoamyl {arm['isoamyl_um']:.0f} µM is under the joint-satisfaction "
            f"ceiling {ceiling:.0f}; both targets are reachable and the collision has gone"
        )

    joint = crepin["joint"]
    assert joint["split_pct"]["leucine"] == pytest.approx(54.5, abs=1.5), (
        f"joint leucine split {joint['split_pct']['leucine']:.1f} % (recorded 54.5)"
    )
    assert joint["split_pct"]["leucine"] > crepin["growth"]["split_pct"]["leucine"] + 20.0
    assert joint["split_pct"]["leucine"] < crepin_lo, (
        "the joint arm reaches Crépin's band — D-260's collision is relieved and this record "
        "is out of date"
    )
    assert rollero_lo < joint["leucine_tracer_pct"] < rollero_hi, (
        f"joint tracer {joint['leucine_tracer_pct']:.3f} % is outside Rollero's "
        f"{rollero_lo}-{rollero_hi} %. It was PREDICTED to fall below the floor and did not; "
        "if it now does, the pre-registered mechanism was right after all — re-record"
    )
    assert joint["leucine_tracer_pct"] == pytest.approx(3.62, abs=0.15)
    assert crepin["joint_lo"]["split_pct"]["leucine"] == pytest.approx(46.0, abs=1.5)
    assert crepin["joint_hi"]["split_pct"]["leucine"] == pytest.approx(61.5, abs=1.5)
    assert crepin["joint_hi"]["split_pct"]["leucine"] < crepin_lo


def test_the_joint_arm_OVER_shoots_crepin_on_the_three_precursors_the_shipped_sink_imposes(crepin):
    """The other side of the split, which the shipped model gets exactly BY CONSTRUCTION.

    ``f_non_ehrlich_isoleucine/valine/threonine`` are Crépin's own 51 / 41 / 38 imposed; the
    growth sink alone under-shoots two of them and over-shoots threonine; the joint arm
    over-shoots all three (66 / 65 / 81 at the mid composition, and above her at every edge).
    The end inversion D-259 found survives — leucine is still the lowest and threonine the
    highest — while isoleucine and valine, which D-259 corrected to the measured order, now
    read within a point of each other. Phenylalanine reads 99.6 % against its sourced 0.975:
    just above D-259's 96-99.5 window, on the same side as the growth arm's 98.9.
    """
    for edge in ("lo", "mid", "hi"):
        arm = crepin[f"joint_{edge}"]["split_pct"]
        for species in ("isoleucine", "valine", "threonine"):
            hi = _CREPIN_PROTEIN_PCT[species][1]
            assert arm[species] > hi, (
                f"joint ({edge}) sends {arm[species]:.1f} % of consumed {species} to the lump "
                f"against Crépin's {hi:.0f} %; it was measured to over-shoot at every edge"
            )
        assert arm["leucine"] < arm["threonine"], "the end inversion has closed; re-open D-259"
    mid = crepin["joint"]["split_pct"]
    assert mid["isoleucine"] == pytest.approx(66.2, abs=1.5)
    assert mid["valine"] == pytest.approx(65.2, abs=1.5)
    assert mid["threonine"] == pytest.approx(81.2, abs=1.5)
    assert abs(mid["isoleucine"] - mid["valine"]) < 2.0
    assert 99.0 < mid["phenylalanine"] < 100.0
    growth = crepin["growth"]["split_pct"]
    assert growth["isoleucine"] < 51.0 and growth["valine"] < 41.0 and growth["threonine"] > 38.0


def test_rolleros_leucine_tracer_brackets_the_enrichment_and_overshoots_the_amount(rollero):
    """P8, refuted on two of three — and the amount says why the SM425 "hit" is not one.

    Enrichment on Rollero's own musts: 1.78 / 5.44 / 8.33 % against Table S2's 3.4-3.5 /
    4.2-4.7 / 6.8-8.2 % — under at SM70, over at SM250, at the top edge at SM425. The shipped
    arm is under on all three and the growth arm over on all three; the joint arm sits between
    and its enrichment responds 4.7× across his nitrogen range where his responds ~2.2×. Scored
    apart from the total ([[feedback-a-near-constant-ratio-can-be-two-errors-growing-together]]),
    the labelled AMOUNT is over at both higher levels — 101 µM against 55-64, 171 against 65-70
    — so the SM425 enrichment lands inside the band only because the total it is divided by is
    2× his (2052 against 793-1034 µM). Two errors, cancelling.
    """
    joint = rollero["joint"]
    lo70 = min(_ROLLERO_EF_LEUCINE_ENRICHMENT[70.0]) * 100.0
    lo250, hi250 = (x * 100.0 for x in sorted(_ROLLERO_EF_LEUCINE_ENRICHMENT[250.0]))
    lo425, hi425 = (x * 100.0 for x in sorted(_ROLLERO_EF_LEUCINE_ENRICHMENT[425.0]))
    assert joint[70.0]["leucine_enrichment_pct"] < lo70
    assert joint[250.0]["leucine_enrichment_pct"] > hi250, (
        f"SM250: joint leucine enrichment {joint[250.0]['leucine_enrichment_pct']:.2f} % was "
        f"measured ABOVE Rollero's {lo250}-{hi250} %; the pre-registration expected below"
    )
    assert joint[425.0]["leucine_enrichment_pct"] == pytest.approx(hi425, abs=0.5)
    assert joint[425.0]["leucine_enrichment_pct"] > lo425
    for yan in _YANS:
        assert (
            rollero["shipped"][yan]["leucine_enrichment_pct"]
            < min(_ROLLERO_EF_LEUCINE_ENRICHMENT[yan]) * 100.0
        ), f"SM{yan:.0f}: the shipped arm is no longer under Table S2's floor (D-256)"
        assert (
            rollero["growth"][yan]["leucine_enrichment_pct"]
            > max(_ROLLERO_EF_LEUCINE_ENRICHMENT[yan]) * 100.0
        ), f"SM{yan:.0f}: the growth arm alone is no longer over Table S2's ceiling"
    response = joint[425.0]["leucine_enrichment_pct"] / joint[70.0]["leucine_enrichment_pct"]
    assert response > 4.0, f"joint enrichment response {response:.2f}× (recorded 4.7×)"
    for yan in (250.0, 425.0):
        measured_hi = max(_ROLLERO_EF_LEUCINE_LABELLED_UM[yan])
        assert joint[yan]["leucine_labelled_um"] > 1.5 * measured_hi, (
            f"SM{yan:.0f}: joint leucine-labelled isoamyl {joint[yan]['leucine_labelled_um']:.0f} "
            f"µM against Rollero's ≤ {measured_hi} µM — the amount was measured 1.6-2.5× over"
        )
    assert joint[425.0]["isoamyl_ef_um"] > 1.9 * max(_ROLLERO_EF_ISOAMYL_UM[425.0]), (
        "the SM425 total is no longer ~2× Rollero's; the cancellation argument above is void"
    )


def test_valines_tracer_moves_the_WRONG_way_under_the_joint_arm(rollero):
    """The blend makes most isoamyl de novo, so the valine-labelled share falls: 0.87 / 2.15 /
    2.85 % against Table S1's 2.1-2.3 / 3.4-4.0 / 5.3-5.4 %, below the shipped arm's own
    under-reading on all three musts. The growth arm alone is the best of the four here. Pinned
    so that "the joint arm is closer to reality" cannot be quoted without this axis.
    """
    for yan in _YANS:
        floor = min(_ROLLERO_EF_VALINE_ENRICHMENT_PCT[yan])
        joint, shipped = rollero["joint"][yan], rollero["shipped"][yan]
        assert joint["enrichment_pct"] < shipped["enrichment_pct"] < floor, (
            f"SM{yan:.0f}: valine enrichment joint {joint['enrichment_pct']:.2f} / shipped "
            f"{shipped['enrichment_pct']:.2f} / floor {floor}. The joint arm was measured to "
            "sit BELOW the shipped arm, which is itself under the band"
        )
        assert rollero["growth"][yan]["enrichment_pct"] > joint["enrichment_pct"]


def test_the_growth_sink_alone_breaks_the_propanol_floor_and_the_joint_arm_restores_it(crepin):
    """A finding neither single repair could show. Growth-anchoring leaves threonine's Ehrlich
    draw un-capped by ``f``, and propanol reads 71.4 % de novo on Crépin's must against the
    sourced 80 % floor (D-244, cleared at D-248); the blend's cut of the Ehrlich draw brings
    the pair back to 87.5 %, beside the shipped 87.8 %."""
    growth = crepin["growth"]["de_novo"]["propanol"]
    joint = crepin["joint"]["de_novo"]["propanol"]
    assert growth < _PROPANOL_DE_NOVO_FLOOR, (
        f"growth-anchored propanol reads {growth:.3f} de novo — measured 0.714, under the floor"
    )
    assert joint > _PROPANOL_DE_NOVO_FLOOR, (
        f"joint propanol reads {joint:.3f} de novo, under the sourced {_PROPANOL_DE_NOVO_FLOOR}"
    )
    assert joint == pytest.approx(crepin["shipped"]["de_novo"]["propanol"], abs=0.02)
    assert crepin["joint"]["de_novo"]["isoamyl_alcohol"] > 0.9


def test_the_joint_arm_closes_both_ledgers_and_never_creates_sugar(d104):
    """P10. The growth sink refunds carbon to ``S`` and nitrogen to ``N`` for the precursor it
    draws, so the D-117 sparing-credit ceiling is the guard that has to hold: 0.236× growth's
    draw, against the shipped 0.584× and the hard 1.0. Carbon and nitrogen close to solver
    precision and the summed right-hand side never makes sugar appear."""
    joint = d104["joint"]
    assert joint["carbon_rel_drift"] < 1e-9 and joint["nitrogen_rel_drift"] < 1e-9, (
        f"ledger drift C {joint['carbon_rel_drift']:.1e} / N {joint['nitrogen_rel_drift']:.1e}"
    )
    assert joint["max_net_ds_dt"] <= 0.0, f"net dS/dt reached {joint['max_net_ds_dt']:+.3e}"
    assert joint["worst_joint_c_refund"] < 1.0
    assert joint["worst_joint_c_refund"] == pytest.approx(0.236, abs=0.03)
    assert d104["shipped"]["worst_joint_c_refund"] == pytest.approx(0.584, abs=0.03), (
        "the shipped joint refund (D-118's 0.584×) has moved; re-derive before comparing"
    )
