"""D-266 probe: four arms x every sourced fusel observable. No src/ change; test-side patches.

Run from the repo root:  uv run --group ui --python 3.13 python docs/receipts/d266-joint-fusel-repair/probe_joint.py
Writes findings.json beside itself.
"""

from __future__ import annotations

import contextlib
import json
import pathlib
import sys
import time
from collections.abc import Iterator, Mapping

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from fermentation.core.chemistry import (  # noqa: E402
    M_ISOAMYL_OH,
    MOLAR_MASS,
    carbon_mass_fraction,
    sugar_species,
)
from fermentation.core.kinetics import byproducts as bp  # noqa: E402
from fermentation.core.kinetics import precursor_fates as pf  # noqa: E402
from fermentation.core.kinetics.amino_acid_pools import SPEC_BY_SPECIES, depletion_gate  # noqa: E402
from fermentation.core.kinetics.arrhenius import arrhenius_factor  # noqa: E402
from fermentation.core.kinetics.carbon_routing import (  # noqa: E402
    FUSEL_SPECS,
    ISOAMYL_ALCOHOL,
    fermentative_flux_shape,
    non_ehrlich_fraction_param,
)
from fermentation.core.kinetics.growth import (  # noqa: E402
    GrowthNitrogenLimited,
    assimilable_nitrogen_pools,
    biomass_growth_rate,
)
from fermentation.core.kinetics.precursor_fates import PrecursorNonEhrlichFates  # noqa: E402
from fermentation.core.process import RateModifier  # noqa: E402
from fermentation.runtime import simulate_scheduled  # noqa: E402
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario  # noqa: E402
from fermentation.validation import total_carbon, total_nitrogen  # noqa: E402
from tests.test_defined_media import (  # noqa: E402
    ROLLERO_SUGAR_GPL,
    ROLLERO_TEMP_C,
    SOURCED_PITCH_GPL,
    _brix_for,
    commensurate_pools,
    commensurate_scenario,
    rollero_must_mm,
)
from tests.test_fusel_keto_acid_node import _OTHER_PRECURSOR_CONSUMERS  # noqa: E402
from tests.test_precursor_fates import _D259_EDGES, _GrowthAnchoredFates, _d259_weights  # noqa: E402

F_CATABOLIC = 0.79  # D-257 §1, fitted on SM250
K_FACTOR = 0.4033  # D-257, the level-preserving rescale of the five k
_ORIG_SHAPE = bp.fusel_rate_shape
_ORIG_DRAWS = bp.ehrlich_draws
_ORIG_LABEL = bp.LABELLED_PRECURSOR


def _terms(y, schema, params):
    flux = fermentative_flux_shape(y, schema, params["K_sugar_uptake"])
    if flux <= 0.0:
        return 0.0, 0.0, 0.0
    n = sum(assimilable_nitrogen_pools(y, schema))
    g_n = n / (params["K_n"] + n) if n > 0.0 else 0.0
    temp = float(y[schema.slice("T")][0])
    return flux, g_n, arrhenius_factor(temp, params["E_a_fusels"], params["T_ref"])


def blended_shape(y, schema, params, *, growth_coupled: bool = False) -> float:
    """D-257's rate law: the catabolic share gated on nitrogen, the rest de novo from sugar."""
    if growth_coupled:
        return _ORIG_SHAPE(y, schema, params, growth_coupled=True)
    flux, g_n, f_t = _terms(y, schema, params)
    if flux <= 0.0:
        return 0.0
    return float(flux * (F_CATABOLIC * g_n + (1.0 - F_CATABOLIC)) * f_t)


def ehrlich_only_shape(y, schema, params, *, growth_coupled: bool = False) -> float:
    """The Ehrlich part of the blend alone — what the sourcing layer may attribute to amino acids."""
    if growth_coupled:
        return _ORIG_SHAPE(y, schema, params, growth_coupled=True)
    flux, g_n, f_t = _terms(y, schema, params)
    if flux <= 0.0 or g_n <= 0.0:
        return 0.0
    return float(flux * F_CATABOLIC * g_n * f_t)


def scaled_draws(y, schema, params):
    """D-257 §4: attribution scaled by the Ehrlich share of the rate (identity at f = 1)."""
    bp.fusel_rate_shape = ehrlich_only_shape
    try:
        return _ORIG_DRAWS(y, schema, params)
    finally:
        bp.fusel_rate_shape = blended_shape


@contextlib.contextmanager
def rate_law(blend: bool) -> Iterator[None]:
    if not blend:
        yield
        return
    bp.fusel_rate_shape = blended_shape
    bp.ehrlich_draws = scaled_draws
    pf.ehrlich_draws = scaled_draws
    try:
        yield
    finally:
        bp.fusel_rate_shape = _ORIG_SHAPE
        bp.ehrlich_draws = _ORIG_DRAWS
        pf.ehrlich_draws = _ORIG_DRAWS


def apply_k_factor(params: dict[str, float], blend: bool) -> None:
    if blend:
        for spec in FUSEL_SPECS:
            params[spec.k_param] = params[spec.k_param] * K_FACTOR


def install_growth_sink(cs, edge: str) -> tuple[_GrowthAnchoredFates, list[RateModifier]]:
    sink = _GrowthAnchoredFates(_d259_weights(edge))
    cs.process_set._processes[sink.name] = sink
    attached: list[RateModifier] = []
    for modifier in cs.process_set.active_modifiers:
        if GrowthNitrogenLimited.name in modifier.modifies:
            modifier.modifies = (*modifier.modifies, sink.name)
            attached.append(modifier)
    assert attached, "no modifier scales growth"
    return sink, attached


ARMS = {
    "shipped": dict(blend=False, growth=False),
    "blend": dict(blend=True, growth=False),
    "growth": dict(blend=False, growth=True),
    "joint": dict(blend=True, growth=True),
}


def _integrate(cs, params, points: int):
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


def _closure_drift(traj, schema, params) -> dict[str, float]:
    c_fn = total_carbon(schema, biomass_carbon_fraction=params["biomass_C_fraction"])
    n_fn = total_nitrogen(schema, biomass_nitrogen_fraction=params["biomass_N_fraction"])
    c = np.array([c_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])
    n = np.array([n_fn(traj.y[:, i]) for i in range(traj.y.shape[1])])
    return {
        "carbon_rel_drift": float(np.max(np.abs(c - c[0])) / c[0]),
        "nitrogen_rel_drift": float(np.max(np.abs(n - n[0])) / n[0]),
    }


def crepin_arm(blend: bool, growth: bool, edge: str = "mid", points: int = 4001) -> dict:
    """The D-259/D-260 fixture: Crepin's must, other precursor consumers off."""
    with rate_law(blend):
        cs = compile_scenario(commensurate_scenario("crepin", days=14.0))
        for name in _OTHER_PRECURSOR_CONSUMERS:
            cs.process_set.disable(name)
        params = cs.param_values
        apply_k_factor(params, blend)
        sink = attached = None
        if growth:
            sink, attached = install_growth_sink(cs, edge)
        traj = _integrate(cs, params, points)
        y, t, schema = np.asarray(traj.y, float), np.asarray(traj.t, float), cs.schema
        species = ("leucine", "isoleucine", "valine", "threonine", "phenylalanine")
        lump = {s: np.zeros_like(t) for s in species}
        ehr = {s: np.zeros_like(t) for s in species}
        alcohol_from = {spec.pool: np.zeros_like(t) for spec in FUSEL_SPECS}
        draws_fn = bp.ehrlich_draws  # the binding the run used
        for i in range(t.size):
            col = y[:, i]
            for draw in draws_fn(col, schema, params):
                sp = draw.precursor.species
                if sp in ehr:
                    ehr[sp][i] += draw.precursor_carbon / carbon_mass_fraction(sp)
                alcohol_from[draw.alcohol.pool][i] += draw.alcohol_carbon
            if growth:
                base_dx = max(biomass_growth_rate(col, schema, params), 0.0)
                factor = 1.0
                for m in attached:
                    factor *= m.factor(float(t[i]), col, schema, params)
                for s in species:
                    gate = depletion_gate(col, schema, params, (SPEC_BY_SPECIES[s],))
                    lump[s][i] = sink.weights[s] * base_dx * gate * factor
            else:
                for s in species:
                    f = params[non_ehrlich_fraction_param(s)]
                    lump[s][i] = ehr[s][i] * f / (1.0 - f)
        out: dict = {"splits": {}, "closure": {}, "left_pct": {}}
        for s in species:
            mw = MOLAR_MASS[s]
            lump_um = float(np.trapezoid(lump[s], t)) / mw * 1e6
            ehr_um = float(np.trapezoid(ehr[s], t)) / mw * 1e6
            pool = y[schema.slice(s), :][0]
            consumed_um = (float(pool[0]) - float(pool[-1])) / mw * 1e6
            out["splits"][s] = 100.0 * lump_um / (lump_um + ehr_um) if (lump_um + ehr_um) > 0 else float("nan")
            out["closure"][s] = (lump_um + ehr_um) / consumed_um if consumed_um > 0 else float("nan")
            out["left_pct"][s] = 100.0 * float(pool[-1]) / float(pool[0])
            if s == "leucine":
                out["leucine_supply_um"] = float(pool[0]) / mw * 1e6
                out["leucine_consumed_um"] = consumed_um
                out["leucine_ehrlich_um"] = ehr_um
                out["leucine_lump_um"] = lump_um
        iso_um = float(y[schema.slice("isoamyl_alcohol"), -1][0]) / M_ISOAMYL_OH * 1e6
        out["isoamyl_um"] = iso_um
        out["isoamyl_mgl"] = float(y[schema.slice("isoamyl_alcohol"), -1][0]) * 1e3
        out["leucine_tracer_pct"] = 100.0 * out["leucine_ehrlich_um"] / iso_um
        # de-novo shares: 1 - amino-acid-sourced alcohol carbon / alcohol carbon made
        out["de_novo_share"] = {}
        for spec in FUSEL_SPECS:
            made_c = float(y[schema.slice(spec.pool), -1][0]) * carbon_mass_fraction(spec.species)
            sourced_c = float(np.trapezoid(alcohol_from[spec.pool], t))
            out["de_novo_share"][spec.pool] = 1.0 - sourced_c / made_c if made_c > 0 else float("nan")
        out["biomass_peak_gpl"] = float(y[schema.slice("X"), :].max())
        out |= _closure_drift(traj, schema, params)
        return out


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


def _gate_closure_index(traj, schema, params) -> int:
    for i in range(traj.y.shape[1]):
        n = sum(assimilable_nitrogen_pools(traj.y[:, i], schema))
        if n / (params["K_n"] + n) < 0.01:
            return i
    raise AssertionError("gate never closes")


def rollero_arm(blend: bool, growth: bool, yan: float, label: str = "valine") -> dict:
    with rate_law(blend):
        bp.LABELLED_PRECURSOR = label
        try:
            cs = compile_scenario(_rollero_scenario(yan), strict=True)
            params = cs.param_values
            apply_k_factor(params, blend)
            if growth:
                install_growth_sink(cs, "mid")
            traj = _integrate(cs, params, int(cs.t_span_h[1]) + 1)
        finally:
            bp.LABELLED_PRECURSOR = _ORIG_LABEL
        schema = cs.schema
        i = _gate_closure_index(traj, schema, params)
        iso = traj.y[schema.slice("isoamyl_alcohol")][0]
        tracer = float(traj.y[schema.slice("isoamyl_alcohol_valine"), -1][0])
        return {
            "nt_fraction_pct": 100.0 * float(iso[i] / iso[-1]),
            "gate_closure_h": float(traj.t[i]),
            "isoamyl_ef_um": float(iso[-1]) / M_ISOAMYL_OH * 1e6,
            f"{label}_enrichment_pct": 100.0 * tracer / float(iso[-1]),
        }


def anchor_arm(blend: bool, yan: float) -> float:
    with rate_law(blend):
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
        apply_k_factor(params, blend)
        traj = _integrate(cs, params, 337)
        return float(traj.y[cs.schema.slice("isoamyl_alcohol"), -1][0]) * 1e3


def d104_arm(blend: bool, growth: bool) -> dict:
    """The D-104/D-257 fixture, FULL process set: phenylalanine, the refund guard, net dS/dt."""
    with rate_law(blend):
        scenario = Scenario(
            name="d266-d104-fixture",
            medium="wine",
            initial={"brix": 24.0, "yan_mgl": 250.0, "pitch_gpl": 0.25, "amino_acids_gpl": 1.0},
            temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
            duration_days=14.0,
        )
        cs = compile_scenario(scenario, strict=True)
        params = cs.param_values
        apply_k_factor(params, blend)
        sink = None
        if growth:
            sink, _ = install_growth_sink(cs, "mid")
        traj = _integrate(cs, params, 337)
        schema, ps = cs.schema, cs.process_set
        phe = traj.y[schema.slice("phenylalanine")][0]
        worst_c = 0.0
        worst_ds = -np.inf
        sink_proc = sink if sink is not None else ps._processes[PrecursorNonEhrlichFates.name]
        from fermentation.core.kinetics.amino_acids import AminoAcidAssimilation

        swap = ps._processes[AminoAcidAssimilation.name]
        s_slice = schema.slice("S")
        for i in range(traj.y.shape[1]):
            yi = traj.y[:, i]
            ti = float(traj.t[i])
            worst_ds = max(worst_ds, float(ps.total_derivatives(ti, yi, params)[s_slice].sum()))
            base_dx = biomass_growth_rate(yi, schema, params)
            if base_dx <= 1e-6:
                continue
            c = 0.0
            for proc in (swap, sink_proc):
                if not ps.is_enabled(proc.name):
                    continue
                d = proc.derivatives(ti, yi, schema, params)
                # a growth-anchored sink is scaled by growth's modifiers inside the run; the
                # refund ratio here is the PRE-modifier one for both (both scale identically)
                for offset, sp in enumerate(sugar_species(schema)):
                    c += float(d[s_slice.start + offset]) * carbon_mass_fraction(sp)
            worst_c = max(worst_c, c / (params["biomass_C_fraction"] * base_dx))
        out = {
            "phe_left_pct": 100.0 * float(phe[-1] / phe[0]),
            "worst_joint_c_refund": worst_c,
            "max_net_dS_dt": worst_ds,
            "isoamyl_mgl": float(traj.y[schema.slice("isoamyl_alcohol"), -1][0]) * 1e3,
            "biomass_final_gpl": float(traj.y[schema.slice("X"), -1][0]),
        }
        out |= _closure_drift(traj, schema, params)
        return out


def main() -> None:
    t0 = time.time()
    results: dict = {}
    for arm, cfg in ARMS.items():
        print(f"== {arm} ==", flush=True)
        r: dict = {}
        r["crepin"] = crepin_arm(cfg["blend"], cfg["growth"])
        print("  crepin:", json.dumps({k: v for k, v in r["crepin"].items() if k in ("splits", "leucine_tracer_pct", "isoamyl_um", "left_pct")}, default=float), flush=True)
        r["rollero"] = {}
        for yan in (70.0, 250.0, 425.0):
            v = rollero_arm(cfg["blend"], cfg["growth"], yan, "valine")
            leu = rollero_arm(cfg["blend"], cfg["growth"], yan, "leucine")
            v["leucine_enrichment_pct"] = leu["leucine_enrichment_pct"]
            v["arms_differ"] = abs(leu["leucine_enrichment_pct"] - v["valine_enrichment_pct"]) > 0.05 * v["valine_enrichment_pct"]
            r["rollero"][str(int(yan))] = v
            print(f"  rollero SM{yan:.0f}:", json.dumps(v, default=float), flush=True)
        r["anchor_mgl"] = {str(int(yan)): anchor_arm(cfg["blend"], yan) for yan in (250.0, 300.0)}
        print("  anchor:", r["anchor_mgl"], flush=True)
        r["d104"] = d104_arm(cfg["blend"], cfg["growth"])
        print("  d104:", json.dumps(r["d104"], default=float), flush=True)
        results[arm] = r
        print(f"  [{time.time() - t0:.0f} s]", flush=True)
    results["joint_edges"] = {}
    for edge in _D259_EDGES:
        results["joint_edges"][edge] = crepin_arm(True, True, edge)
        print(f"joint {edge}: splits", json.dumps(results["joint_edges"][edge]["splits"], default=float), flush=True)
    out = pathlib.Path(__file__).with_name("findings.json")
    out.write_text(json.dumps(results, indent=1, default=float), encoding="utf-8")
    print("wrote", out, f"[{time.time() - t0:.0f} s]")


if __name__ == "__main__":
    main()
