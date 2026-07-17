"""Tests for the mercaptan (thiol) beat — the carbon-bearing autolytic reductive off-aroma (D-45).

Beyond H₂S (D-44), the other "reduction" off-aroma the model carries is **methanethiol**.
:class:`AutolyticMercaptan` fills the ``methanethiol`` pool as a yield on the shared autolysis flux
— but because methanethiol carries **carbon** (unlike H₂S), it draws that carbon from the
**methionine** pool (decision D-100 — the *actual* precursor; D-45 had to draw from the lumped
arginine pool and document that arginine contains no sulfur) and **deaminates** the nitrogen to
``N`` (Option A, the D-33 fusel-reroute idiom). This suite pins the closed form, the
carbon+nitrogen closure (both by construction), the
availability gate + guards, the new ``tier_of("N")`` drop (the first autolysis-gated N-writer), the
opt-in isolability, and the emergent post-dryness accumulation.

**The pool was ``mercaptans``, a lumped plural, through D-109 — and the lump was FALSE, not merely
coarse** (decision D-110): nothing in the model produces ethanethiol or any other thiol, so the
fixed-composition caveat described a mixture the mass balance never contained. The rename to the
one molecule it holds moved no number, because the pool already *was* that molecule at every layer.
See ``test_the_pool_is_one_molecule_and_names_it`` and
``test_the_model_makes_no_thiol_but_this_one`` — the latter is the one that keeps the claim honest
as the model grows.
"""

from collections.abc import Mapping

import numpy as np
import pytest

from fermentation.core.chemistry import (
    MOLAR_MASS,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics import AutolyticMercaptan, GrowthNitrogenLimited
from fermentation.core.kinetics.arrhenius import arrhenius_factor
from fermentation.core.media import wine_schema
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate
from fermentation.scenario import Scenario, TemperaturePoint, compile_scenario
from fermentation.sensory.oav import AROMA_COMPOUNDS, load_thresholds
from fermentation.validation import (
    assert_conserved,
    assert_nonnegative,
    total_carbon,
    total_nitrogen,
)
from tests.conftest import seed_amino_acids

_MERCAPTAN_SPECIES = "methanethiol"
#: The thiol's real precursor since D-100 (was the lumped arginine stand-in).
_AA_SPECIES = "methionine"
#: Demethiolation's C4 co-product (D-107) — the keto-acid node. Its existence is what let the
#: methionine draw become the honest 1:1: the other four carbons finally have somewhere to go.
_CO_PRODUCT_SPECIES = "alpha_ketobutyrate"


@pytest.fixture
def store():
    # wine_generic carries every mercaptan read (y_methanethiol, k_autolysis, E_a_autolysis, T_ref,
    # K_amino_acids); the shared H₂S file rounds out the resolved set for tier maps.
    return load_parameters(
        default_data_dir() / "wine_generic.yaml",
        default_data_dir() / "hydrogen_sulfide.yaml",
    )


@pytest.fixture
def params(store):
    return store.resolve()


def _mercaptan_y0(
    schema: StateSchema,
    params: Mapping[str, float],
    *,
    x_dead: float = 1.5,
    amino_acids: float = 1.0,
    s: float = 200.0,
    x: float = 2.0,
    t: float = 293.15,
) -> FloatArray:
    # A post-AF-ish state: dead biomass to autolyse, an amino-acid pool (autolysis-refilled) to
    # source the mercaptan carbon, ethanol high, nitrogen exhausted. The amino acids are seeded at
    # MUST-SPECTRUM composition (D-100), the state in which every per-species gate provably equals
    # the pre-split lumped gate — so the closed form below asserts the same numbers it always did.
    y = schema.pack(
        {
            "X": x,
            "S": [s],
            "E": 100.0,
            "N": 0.0,
            "T": t,
            "CO2": 0.0,
            "X_dead": x_dead,
        }
    )
    return seed_amino_acids(y, schema, params, amino_acids)


def _expected(params: Mapping[str, float], *, x_dead: float, amino_acids: float, t: float):
    """The closed form: (r_merc, met_mass, keto_rate, n_release).

    **This helper used to encode the bug (decision D-107).** Through D-106 it read
    ``aa_mass = merc_carbon / c_methionine`` — sizing the methionine draw to the *thiol's single
    carbon*, i.e. 0.2 mol methionine per mol thiol. It matched the Process exactly and passed for
    sixty-odd decisions, because it was the same arithmetic written twice. That is the D-105 lesson
    in test form: **a closed-form check written from the code cannot find a stoichiometry error in
    the code** — only the reaction can, and the reaction is
    ``1 mol methionine → 1 mol methanethiol + 1 mol 2-oxobutyrate + NH₃``.

    So this is now written from the **reaction**, not the implementation: one mole of methionine per
    mole of thiol, its C4 to α-ketobutyrate, its nitrogen to N — and the carbon identity ``5 ==
    1+4``
    is asserted independently in :func:`test_carbon_closes_at_the_derivative_level` rather than
    assumed here.
    """
    f_t = arrhenius_factor(t, params["E_a_autolysis"], params["T_ref"])
    r_autolysis = params["k_autolysis"] * f_t * x_dead
    # At must-spectrum composition methionine's relative-depletion gate aa_i/(K·f_i + aa_i) is
    # ALGEBRAICALLY the pre-split lumped gate aa/(K + aa) (decision D-100), so this closed form is
    # unchanged by the split — which is exactly the property being asserted.
    gate = amino_acids / (params["K_amino_acids"] + amino_acids)
    r_merc = params["y_methanethiol"] * r_autolysis * gate
    # STOICHIOMETRY, not carbon-sizing: 1 mol methionine per mol thiol (D-107).
    n_merc = r_merc / MOLAR_MASS[_MERCAPTAN_SPECIES]  # mol/L/h
    met_mass = n_merc * MOLAR_MASS[_AA_SPECIES]
    keto_rate = n_merc * MOLAR_MASS[_CO_PRODUCT_SPECIES]  # the C4 that used to be discarded
    n_release = met_mass * nitrogen_mass_fraction(_AA_SPECIES)
    return r_merc, met_mass, keto_rate, n_release


# -- metadata + closed form ---------------------------------------------------


def test_metadata():
    p = AutolyticMercaptan()
    assert p.name == "autolytic_mercaptan"
    assert p.tier is Tier.SPECULATIVE
    # Draws METHIONINE (D-100), not the retired lumped arginine pool — and since D-107 books the
    # C4 co-product demethiolation really releases, which is what makes that draw 1:1.
    assert set(p.touches) == {"methanethiol", "methionine", "alpha_ketobutyrate", "N"}
    assert "amino_acids" not in p.touches
    assert set(p.reads) == {
        "y_methanethiol",
        "k_autolysis",
        "E_a_autolysis",
        "T_ref",
        "K_amino_acids",
        "must_aa_fraction_methionine",
    }


def test_mercaptan_is_carbon_bearing():
    # THE structural difference from H₂S (which is carbon-free): methanethiol carries one carbon, so
    # the pool sits on total_carbon and the Process must source that carbon from a tracked pool.
    assert carbon_mass_fraction(_MERCAPTAN_SPECIES) == pytest.approx(12.011 / 48.107, rel=1e-4)
    assert nitrogen_mass_fraction(_MERCAPTAN_SPECIES) == 0.0  # N-free ⇒ the drawn N is deaminated


def test_matches_closed_form(params):
    schema = wine_schema()
    y = _mercaptan_y0(schema, params, x_dead=1.5, amino_acids=1.0)
    d = AutolyticMercaptan().derivatives(0.0, y, schema, params)
    r_merc, met_mass, keto_rate, n_release = _expected(
        params, x_dead=1.5, amino_acids=1.0, t=293.15
    )
    assert schema.get(d, "methanethiol") == pytest.approx(r_merc)
    assert schema.get(d, "methanethiol") > 0.0
    assert schema.get(d, "methionine") == pytest.approx(-met_mass)
    assert schema.get(d, "alpha_ketobutyrate") == pytest.approx(keto_rate)
    assert schema.get(d, "N") == pytest.approx(n_release)
    # touches ONLY those four — nothing else on the state moves (not X_dead, S, E, h2s, …)
    for name in schema.names:
        if name in ("methanethiol", "methionine", "alpha_ketobutyrate", "N"):
            continue
        assert schema.get(d, name) == pytest.approx(0.0, abs=1e-18), name


def test_carbon_closes_at_the_derivative_level(params):
    """Carbon out of methionine == carbon into methanethiol + alpha-ketobutyrate (D-107).

    **The meaning of this test changed even though its form did not.** Through D-106 it compared two
    numbers that were *defined* to be equal: the draw was sized from the thiol's carbon, so closure
    was a tautology and the test could never fail — which is precisely why it never caught the 5×
    under-draw sitting next to it (D-105: "a draw defined to close the ledger can never violate
    it"). Now the draw is sized by MOLES, and closure is a real constraint that holds only because
    methionine's 5 carbons genuinely equal the thiol's 1 + the keto-acid's 4. Break the split — book
    the C4 at the wrong molar mass, or drop the pool — and this fails.
    """
    schema = wine_schema()
    d = AutolyticMercaptan().derivatives(0.0, _mercaptan_y0(schema, params), schema, params)
    c_merc = schema.get(d, "methanethiol") * carbon_mass_fraction(_MERCAPTAN_SPECIES)
    c_keto = schema.get(d, "alpha_ketobutyrate") * carbon_mass_fraction(_CO_PRODUCT_SPECIES)
    c_aa = schema.get(d, "methionine") * carbon_mass_fraction(_AA_SPECIES)
    assert c_merc + c_keto + c_aa == pytest.approx(0.0, abs=1e-18)  # gain == loss


def test_demethiolation_draws_one_mole_of_methionine_per_mole_of_thiol(params):
    """The D-105 conviction, now the D-107 fix — measured off ``dy/dt``, not declared (D-107).

    D-45 drew **0.2 mol methionine per mol methanethiol**: the draw was sized to the thiol's single
    carbon, so methionine's other four (the 2-oxobutyrate) were never charged and its nitrogen was
    released at 1/5 the real rate. Carbon closed the whole time — the draw was defined to close it —
    which is why 1140 tests and every conservation check missed it for sixty decisions.

    The mechanism the module names in prose is what convicted it, with no literature needed:
    demethiolation is ``1 mol methionine → 1 mol methanethiol + 1 mol 2-oxobutyrate + NH₃``, and a
    1:1 reaction was being run at 0.2:1 — an internal contradiction. This test pins the ratio by
    DRIVING the Process and reading the debit, so it fails on implementation drift and not merely on
    a table edit (the D-105 two-layer lesson, where layer 1 alone proved insufficient).
    """
    schema = wine_schema()
    d = AutolyticMercaptan().derivatives(0.0, _mercaptan_y0(schema, params), schema, params)
    n_thiol = schema.get(d, "methanethiol") / MOLAR_MASS[_MERCAPTAN_SPECIES]
    n_met = -schema.get(d, "methionine") / MOLAR_MASS[_AA_SPECIES]
    n_keto = schema.get(d, "alpha_ketobutyrate") / MOLAR_MASS[_CO_PRODUCT_SPECIES]
    assert n_thiol > 0.0
    assert n_met / n_thiol == pytest.approx(1.0, abs=1e-12), (
        "demethiolation must consume 1 mol methionine per mol thiol; D-45 consumed 0.2 by sizing "
        "the draw to the thiol's single carbon"
    )
    # The co-product is 1:1 too — it IS the other four carbons, so it cannot be anything else.
    assert n_keto / n_thiol == pytest.approx(1.0, abs=1e-12)


def test_nitrogen_closes_at_the_derivative_level(params):
    # The nitrogen leaving amino_acids (arginine) all lands in the N pool (methanethiol is N-free):
    # the DEAMINATION branch, so total_nitrogen is unchanged by the transfer.
    schema = wine_schema()
    d = AutolyticMercaptan().derivatives(0.0, _mercaptan_y0(schema, params), schema, params)
    n_out_of_aa = schema.get(d, "methionine") * nitrogen_mass_fraction(_AA_SPECIES)
    n_into_pool = schema.get(d, "N")  # the N pool weight is 1.0
    assert n_out_of_aa + n_into_pool == pytest.approx(0.0, abs=1e-18)


# -- flux coupling, gate, temperature, guards ---------------------------------


def test_scales_with_dead_biomass(params):
    # First-order in X_dead (via the shared autolysis flux): 2× the dead biomass ⇒ 2× the rate.
    schema = wine_schema()
    r1 = AutolyticMercaptan().derivatives(
        0.0, _mercaptan_y0(schema, params, x_dead=1.0), schema, params
    )
    r2 = AutolyticMercaptan().derivatives(
        0.0, _mercaptan_y0(schema, params, x_dead=2.0), schema, params
    )
    assert schema.get(r2, "methanethiol") == pytest.approx(2.0 * schema.get(r1, "methanethiol"))


def test_availability_gate_ramps_with_amino_acids(params):
    # The smooth gate aa/(K+aa): more amino acids ⇒ closer to the full yield; near-empty ⇒ near 0.
    schema = wine_schema()
    lo = AutolyticMercaptan().derivatives(
        0.0, _mercaptan_y0(schema, params, amino_acids=0.01), schema, params
    )
    hi = AutolyticMercaptan().derivatives(
        0.0, _mercaptan_y0(schema, params, amino_acids=2.0), schema, params
    )
    assert 0.0 < schema.get(lo, "methanethiol") < schema.get(hi, "methanethiol")


def test_is_not_flux_linked(params):
    # Like the D-44 H₂S source (and unlike the D-29 producer): reads only X_dead/amino_acids/T, so
    # it fires at S=0 and X=0 (post-fermentation) — which is why the flux-linked stripping sink
    # cannot sweep the thiols and they accumulate as residual.
    schema = wine_schema()
    dry = AutolyticMercaptan().derivatives(
        0.0, _mercaptan_y0(schema, params, s=0.0, x=0.0), schema, params
    )
    r_merc, _, _, _ = _expected(params, x_dead=1.5, amino_acids=1.0, t=293.15)
    assert schema.get(dry, "methanethiol") == pytest.approx(r_merc)
    assert schema.get(dry, "methanethiol") > 0.0


def test_rises_with_temperature(params):
    # Autolysis is enzymatic (shares E_a_autolysis via the flux), so warmer lees release faster.
    schema = wine_schema()
    cold = AutolyticMercaptan().derivatives(
        0.0, _mercaptan_y0(schema, params, t=283.15), schema, params
    )
    warm = AutolyticMercaptan().derivatives(
        0.0, _mercaptan_y0(schema, params, t=303.15), schema, params
    )
    assert schema.get(warm, "methanethiol") > schema.get(cold, "methanethiol") > 0.0


def test_zero_without_dead_biomass_or_amino_acids(params):
    # Two guards: no dead cells (autolysis_flux ≤ 0) ⇒ 0; empty amino-acid pool ⇒ 0 (no carbon
    # source, the D-33 no-op). Solver undershoots (negative) are clamped.
    schema = wine_schema()
    p = AutolyticMercaptan()
    no_dead = p.derivatives(0.0, _mercaptan_y0(schema, params, x_dead=0.0), schema, params)
    no_aa = p.derivatives(0.0, _mercaptan_y0(schema, params, amino_acids=0.0), schema, params)
    neg_dead = p.derivatives(0.0, _mercaptan_y0(schema, params, x_dead=-1e-6), schema, params)
    for d in (no_dead, no_aa, neg_dead):
        assert schema.get(d, "methanethiol") == 0.0
        assert schema.get(d, "methionine") == 0.0
        assert schema.get(d, "N") == 0.0


# -- D-110: the false lump, retired ---------------------------------------------


def test_the_pool_is_one_molecule_and_names_it(store):
    """The D-110 retire, pinned where it can fail: the slot, the weight and the threshold agree.

    The `mercaptans` lump was FALSE, not merely coarse — the pool held exactly one molecule
    under a plural name — and the rename is only legitimate because of that. So the tripwire is
    not "the flag is False" (which a careless edit could restore) but the fact underneath it:
    the pool's slot name IS its molecule, so the carbon weight, the OAV threshold and the
    Stevens exponent cannot drift apart from the compound the Process actually makes. The
    single-molecule identity `pool == representative` is what D-96 established after the
    `esters` lump's split identity produced a non-physical OAV; this pins it here.
    """
    (compound,) = [c for c in AROMA_COMPOUNDS["wine"] if c.pool == _MERCAPTAN_SPECIES]
    assert compound.lumped is False
    # The identity that makes the lump-composition question un-askable for this pool.
    assert compound.pool == compound.representative == _MERCAPTAN_SPECIES
    # The chemistry module knows it as a real molecule, so the weight is read, never assumed.
    assert MOLAR_MASS[_MERCAPTAN_SPECIES] > 0.0
    assert carbon_mass_fraction(_MERCAPTAN_SPECIES) > 0.0
    # And the sensory layer resolves it by that same name — the drift the rename forecloses.
    assert f"threshold_{_MERCAPTAN_SPECIES}_wine" in load_thresholds()


def test_the_model_makes_no_thiol_but_this_one(store):
    """WHY the lump was false, asserted rather than recited — the load-bearing half of D-110.

    The retire rests on a claim about the whole model, not about this Process: *nothing produces
    ethanethiol or any other thiol*, so a fixed-lump-composition caveat described a mixture the
    mass balance never contained. Recited in prose that claim is unfalsifiable, so it is checked
    against the state schema itself — if a future beat adds a second thiol pool, the pool really
    does become a mixture again, and this fires to say the `lumped=False` flag is now a lie.

    Deliberately checked over the SCHEMA rather than a curated list: a hardcoded list would pass
    by construction and guard nothing (D-105's `4 + 1 == 5`).
    """
    schema = wine_schema()
    # Thiols the wine literature attributes to reduction, none of which the model carries. This
    # is the set the retired `mercaptans` plural implicitly claimed to contain.
    other_thiols = {"ethanethiol", "propanethiol", "butanethiol", "benzenethiol", "furfurylthiol"}
    present = other_thiols & set(schema.names)
    assert not present, (
        f"the model now carries a second thiol {sorted(present)} — the `methanethiol` pool is a "
        "mixture again and its lumped=False flag (D-110) no longer tells the truth"
    )
    assert _MERCAPTAN_SPECIES in schema


# -- tier: the new structural drop on N (the D-27 E parallel) ------------------


def test_output_tier_is_speculative(store):
    schema = wine_schema()
    ps = ProcessSet(schema, [AutolyticMercaptan()])
    assert ps.tier_of("methanethiol", store.tier_map()) is Tier.SPECULATIVE


def test_is_the_first_autolysis_gated_n_writer(store):
    # THE tier headline (D-45, advisor-flagged): AutolyticMercaptan is the first autolysis-gated
    # Process to WRITE N (via deamination), so adding it to a set drops the structural tier_of("N")
    # from growth's PLAUSIBLE to SPECULATIVE — the D-27 E / D-26 CO2 parallel.
    schema = wine_schema()
    growth_only = ProcessSet(schema, [GrowthNitrogenLimited()])
    with_mercaptan = ProcessSet(schema, [GrowthNitrogenLimited(), AutolyticMercaptan()])
    assert growth_only.tier_of("N") is Tier.PLAUSIBLE
    assert with_mercaptan.tier_of("N") is Tier.SPECULATIVE


# -- integrated: opt-in isolability + post-dryness accumulation + conservation ------------------


def _run_autolysis(*, rate_per_h: float = 2.0e-3, days: float = 40.0):
    scenario = Scenario(
        name="wine-mercaptan",
        medium="wine",
        initial={
            "brix": 24.0,
            "yan_mgl": 80.0,
            "pitch_gpl": 0.25,
            "autolysis_rate_per_h": rate_per_h,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=days,
    )
    compiled = compile_scenario(scenario, strict=True)
    t_eval = np.linspace(0.0, compiled.t_span_h[1], int(days) + 1)
    traj = simulate(
        compiled.process_set, compiled.param_values, compiled.y0, compiled.t_span_h, t_eval=t_eval
    )
    assert traj.success, traj.message
    return traj, compiled


def test_disabled_without_opt_in():
    # Opt-in isolability (D-45/D-34): a default wine run has AutolyticMercaptan disabled, so the
    # mercaptans pool never fills and stays exactly 0 — byte-for-byte the validated core.
    scenario = Scenario(
        name="wine-default",
        medium="wine",
        initial={"brix": 24.0, "yan_mgl": 80.0, "pitch_gpl": 0.25},
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=30.0,
    )
    compiled = compile_scenario(scenario, strict=True)
    assert "autolytic_mercaptan" not in {p.name for p in compiled.process_set.active}
    traj = compiled.run()
    assert float(np.max(np.abs(traj.series("methanethiol")))) == 0.0


def test_mercaptans_accumulate_post_dryness():
    # The reductive-fault headline: with autolysis opted in, mercaptans build up and keep RISING
    # deep post-dryness (the autolysis-refilled amino_acids feed the thiol; not flux-linked, so no
    # CO2 sweeps it). Reaches the sensory scale (methanethiol threshold ~2-3 µg/L).
    traj, _ = _run_autolysis(rate_per_h=2.0e-3, days=40.0)
    merc = np.asarray(traj.series("methanethiol"))
    assert_nonnegative(traj, ("methanethiol", "methionine"), atol=1e-12)
    assert merc[-1] > 5.0e-6  # > 5 µg/L — a clear reductive signal above threshold
    i15 = int(np.argmin(np.abs(traj.t / 24.0 - 15.0)))
    assert merc[-1] > 1.5 * float(merc[i15]) > 0.0  # still rising well after dryness


def test_conserves_carbon_and_nitrogen_on_a_compiled_run():
    # Mercaptan formation is a pure transfer (amino_acids carbon → mercaptans; amino_acids nitrogen
    # → N), both weighted, so a full autolysis-on run (no copper removal) conserves BOTH elements to
    # machine precision — the closure the Option-A routing was chosen to guarantee.
    traj, compiled = _run_autolysis(rate_per_h=2.0e-3, days=40.0)
    f_c = compiled.param_values["biomass_C_fraction"]
    f_n = compiled.param_values["biomass_N_fraction"]
    assert_conserved(
        traj, total_carbon(compiled.schema, biomass_carbon_fraction=f_c), label="carbon"
    )
    assert_conserved(
        traj, total_nitrogen(compiled.schema, biomass_nitrogen_fraction=f_n), label="nitrogen"
    )
    assert float(np.asarray(traj.series("methanethiol"))[-1]) > 0.0  # it really did accumulate
