"""The uptake surplus is held INSIDE the cell, so it stops titrating the must (decision D-250).

**The defect D-248 shipped, and the invariant it broke.**
:class:`~fermentation.core.kinetics.amino_acids.AssimilableNitrogenUptake` books amino-acid
nitrogen the yeast has transported in ahead of anabolic demand. D-248 parked it in the ``N``
slot — and ``N`` is read by the acid-base charge balance (``acidbase.NITROGEN_KEY``, D-209),
at the must's mean charge per mole of nitrogen. So nitrogen already inside a cell went on
titrating the liquid around it. ``nitrogen_charge_excess``'s own docstring states the invariant
this violates: *"Constant except at dose events … only an addition of differently-charged
nitrogen moves it. No Process touches this slot"* (D-210). Uptake **is** such an addition.

**Measured before the repair**, on a wine dosed with 2 g/L amino acids: the ``N`` slot ran
300 → **436.8** mg N/L (1.456x its own starting share, never above the must's declared total —
mass was always conserved, the defect is charge), and pH ran 3.030 → **3.216**, a **+0.215**
excursion against the same run with uptake disabled. At 0.5 g/L it was +0.045, and on an undosed
must exactly 0.0000 — the Process is disabled at the compile seam there, so D-248's isolability
claim was never in question. It is a **mid-run transient**: ``N`` still reaches ~0 by day 2, so
every endpoint-scored benchmark was blind to it, and nothing in the wine suite scored pH mid-run
on a dosed must. That is why this file exists: the artefact was invisible to every shipped guard.

**The repair and its whole observable footprint.** ``stored_nitrogen`` is a wine-only slot in no
charge balance; growth's Monod and draw read it together with ``N``, split in proportion to what
each holds, and the D-32 swap refunds that draw on the **same** split. Because both pools have
the same source and the same sinks, the SUM ``N + stored_nitrogen`` follows exactly the
trajectory ``N`` alone followed before D-250 — so beer is bit-identical and wine moves nowhere
except pH. That prediction was pre-registered and is pinned below.

**What this does NOT repair, stated so a later beat does not re-find it as a bug.** MLF and Brett
stay blind to the store, deliberately: it is inside a yeast cell. Co-inoculated bacteria still
lose ~96 % of their growth increment to yeast uptake, and that is the model being **right** —
real yeast take essentially all the assimilable nitrogen (Crépin's 0.2 % residual). The gap is a
bacterial nitrogen source this model lacks (peptides, which yeast do not take), not blindness to
the yeast's own store. D-248's ``Flags: D-100`` residue closes with that reframing.
"""

from __future__ import annotations

import numpy as np
import pytest

from fermentation.core import acidbase
from fermentation.core.kinetics import AminoAcidAssimilation, GrowthNitrogenLimited
from fermentation.core.kinetics.brett import BrettGrowth
from fermentation.core.kinetics.growth import (
    STORED_NITROGEN_KEY,
    add_assimilable_nitrogen,
    assimilable_nitrogen_pools,
    biomass_growth_rate,
)
from fermentation.core.kinetics.malolactic import MalolacticGrowth
from fermentation.core.media import beer_schema, wine_schema
from fermentation.core.process import Process, ProcessSet
from fermentation.core.tiers import Tier
from fermentation.parameters.store import default_data_dir, load_parameters
from fermentation.runtime import simulate_scheduled
from fermentation.runtime.integrate import simulate
from fermentation.scenario import (
    Scenario,
    TemperaturePoint,
    amino_acid_dose_nitrogen_mgl,
    compile_scenario,
)
from fermentation.validation import assert_conserved, assert_nonnegative, total_nitrogen
from tests.test_defined_media import _assimilable_n_mgl, commensurate_scenario
from tests.test_fusel_keto_acid_node import _OTHER_PRECURSOR_CONSUMERS

UPTAKE = "assimilable_nitrogen_uptake"

#: The pH excursion D-248 shipped on the 2 g/L must, measured against the same run with uptake
#: disabled. Quoted so the guard below is a statement about a known magnitude rather than an
#: unanchored inequality.
D248_PH_EXCURSION = 0.2147

#: ``N``'s peak as a multiple of its own starting value on the same must, before the repair.
D248_AMMONIUM_INFLATION = 1.456


@pytest.fixture(scope="module")
def wine_params():
    data = default_data_dir()
    return load_parameters(data / "wine_generic.yaml", data / "acidbase.yaml").resolve()


def _dosed_must(amino_acids_gpl: float) -> dict[str, float]:
    """A wine must carrying ``amino_acids_gpl`` of dosed amino acids.

    Migrated the D-244 way — the dose's nitrogen is ADDED to the declared YAN, because since
    D-244 ``yan_mgl`` is the must's total assimilable nitrogen and the pools are carved out of
    it. Re-authoring the composition instead is the recorded error.
    """
    initial = {
        "brix": 24.0,
        "yan_mgl": 300.0,
        "pitch_gpl": 0.25,
        "tartaric_gpl": 3.0,
        "malic_gpl": 3.0,
        "initial_ph": 3.4,
        "amino_acids_gpl": amino_acids_gpl,
    }
    initial["yan_mgl"] += amino_acid_dose_nitrogen_mgl(initial)
    return initial


def _run(
    initial: dict[str, float],
    *,
    uptake: bool,
    days: float = 14.0,
    rtol: float | None = None,
    atol: float | None = None,
):
    """One arm. Compiled INSIDE the call: a reused CompiledScenario carries event state.

    ``rtol``/``atol`` are for the one caller that needs the SAME arm at two solver accuracies
    (D-253's noise-versus-over-draw discriminator); left unset the engine's defaults apply and
    every other caller integrates exactly what it always did.
    """
    scenario = Scenario(
        name="d250",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=days,
    )
    compiled = compile_scenario(scenario, strict=True)
    if not uptake:
        compiled.process_set.disable(UPTAKE)
    t_eval = np.linspace(0.0, days * 24.0, 400)
    traj = simulate(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
        t_eval=t_eval,
        # These literals MIRROR simulate()'s own defaults (integrate.py) rather than overriding
        # them; only the one caller that passes explicit values integrates anything different.
        rtol=1e-6 if rtol is None else rtol,
        atol=1e-9 if atol is None else atol,
    )
    assert traj.success, traj.message
    return compiled, traj


def _ph_course(compiled, traj) -> np.ndarray:
    return np.array(
        [
            acidbase.ph_of_state(traj.y[:, i], compiled.schema, compiled.param_values)
            for i in range(traj.y.shape[1])
        ]
    )


# -- 1. HEADLINE (fail-first): uptake no longer titrates the must ---------------------------


@pytest.mark.parametrize("dose", [0.5, 2.0])
def test_the_uptake_surplus_no_longer_raises_the_musts_ammonium(dose):
    """``N`` must not rise above its own starting value because the yeast ate amino acids.

    The fail-first control is the same must with uptake disabled: it fixes what ``N`` does when
    nothing is refunding to it, so the assertion is a contrast rather than an absolute. Before
    D-250 the uptake arm peaked at 1.456x its start on the 2 g/L must; the control peaks at
    ~1.004x, which is the D-32 swap's own refund and is bounded by growth's draw by construction.
    """
    must = _dosed_must(dose)
    _, on = _run(dict(must), uptake=True)
    _, off = _run(dict(must), uptake=False)

    n_on, n_off = on.series("N"), off.series("N")
    inflation_on = float(n_on.max() / n_on[0])
    inflation_off = float(n_off.max() / n_off[0])

    assert inflation_on < 1.02, (
        f"the must's ammonium peaks at {inflation_on:.4f}x its start with uptake on. D-248 "
        f"shipped {D248_AMMONIUM_INFLATION} on the 2 g/L must because it parked the surplus in "
        "`N`; above 1.02 the surplus is back in the medium"
    )
    assert inflation_on == pytest.approx(inflation_off, abs=0.01), (
        f"uptake moves the ammonium peak {inflation_on:.4f} vs {inflation_off:.4f} without it. "
        "The two arms should differ only in how much nitrogen the CELLS hold"
    )


def test_the_ph_excursion_d248_shipped_is_gone():
    """The pH gap between the uptake and no-uptake arms, on the must where it was 0.215.

    Not asserted to zero, and deliberately: uptake genuinely changes the run — more biomass, a
    faster ferment — so the two arms' pH curves legitimately differ in TIME. What must be gone is
    the *level* artefact, the extra ammonium. The peak pH of the two arms must now agree closely,
    which a timing difference does not disturb; D-248's peaks differed by 0.108.
    """
    must = _dosed_must(2.0)
    on_c, on = _run(dict(must), uptake=True)
    off_c, off = _run(dict(must), uptake=False)

    ph_on, ph_off = _ph_course(on_c, on), _ph_course(off_c, off)
    gap = float(np.abs(ph_on - ph_off).max())
    peak_gap = float(abs(ph_on.max() - ph_off.max()))

    assert gap < 0.5 * D248_PH_EXCURSION, (
        f"the pH gap between the arms is {gap:.4f}, not meaningfully below the "
        f"{D248_PH_EXCURSION} D-248 shipped — the surplus nitrogen is titrating the must again"
    )
    assert peak_gap < 0.01, (
        f"peak pH differs by {peak_gap:.4f} between the arms (uptake {ph_on.max():.4f} vs "
        f"{ph_off.max():.4f}); D-248's was 0.108. A LEVEL difference is the artefact"
    )


@pytest.mark.parametrize("dose", [0.5, 2.0])
def test_the_starved_bacterium_now_leaves_more_malate_at_both_doses(dose):
    """The repair UN-MASKS the physical signal, and this is the one case nothing else runs.

    Every shipped MLF and Brett test isolates the uptake competitor in both arms, which is why
    the artefact could ride so long: nothing exercised a co-inoculated MLF with uptake LIVE.

    Yeast uptake starves the bacterium, so the starved arm must leave MORE malate unconverted.
    Before D-250 that held at 0.5 g/L (4.00x) and **reversed** at 2 g/L (0.853x) — the pH
    excursion lifted MLF's own pH logistic `1/(1 + 10^(pH_half_mlf - pH))` enough to
    over-compensate for the catalyst it had lost. The reading depended on the DOSE, which is the
    defect, not the ratio at either point. It now reads 4.66x and 1.18x: same sign, both doses.
    """
    must = _dosed_must(dose)
    must["mlf_pitch_gpl"] = 0.05
    _, on = _run(dict(must), uptake=True)
    _, off = _run(dict(must), uptake=False)

    malic_on = float(np.interp(72.0, on.t, on.series("malic")))
    malic_off = float(np.interp(72.0, off.t, off.series("malic")))
    assert malic_on > malic_off, (
        f"at {dose} g/L the starved arm left {malic_on:.4f} g/L malate against the control's "
        f"{malic_off:.4f} — LESS, so losing its nitrogen made the bacterium convert more. Some "
        "second channel is over-compensating in this observable; before D-250 it was a pH "
        "excursion feeding MLF's own pH gate"
    )


# -- 2. The defect was CHARGE, never mass ---------------------------------------------------


@pytest.mark.parametrize("uptake", [True, False])
def test_nitrogen_is_conserved_in_both_arms_so_the_slot_rise_was_never_creation(uptake):
    """``N`` exceeding its own start read like nitrogen appearing. It never was.

    The store is on ``total_nitrogen`` at weight 1.0, exactly like ``N``, so the ledger closes
    through the uptake transfer as well as through growth's draw back out of it.
    """
    compiled, traj = _run(_dosed_must(2.0), uptake=uptake)
    # The store inherited no nonnegativity check from anywhere: it is a new slot. The draw out
    # of it is proportional to its own share, so it is self-limiting, but that is an argument
    # rather than a check.
    #
    # THE atol IS RAISED FROM THE 1e-9 DEFAULT AND THE REASON IS MEASURED (decision D-253).
    # D-251 recorded this guard as within a factor of two of firing on solver noise; at D-253's
    # shipped capacity it fired here, at -2.36e-9 on the 2 g/L dosed must. It is noise, and the
    # discriminator is not the size but the SCALING: tightening the solver two orders takes the
    # same dip to -1.4e-12 and four orders takes it to -3.4e-14, while the store's own peak is
    # unchanged to five decimals. An over-draw does not do that -- it is a property of the
    # derivative and survives any tolerance. So the loosened bound below is paired with the
    # tightened run that follows, and it is that pair, not either number, that forbids an
    # over-draw. Do not raise this further without re-measuring the scaling.
    assert_nonnegative(traj, (STORED_NITROGEN_KEY,), atol=1e-8)

    tight = _run(_dosed_must(2.0), uptake=uptake, rtol=1e-8, atol=1e-11)[1]
    assert_nonnegative(tight, (STORED_NITROGEN_KEY,))  # the DEFAULT bound, not loosened
    worst_default = float(np.min(traj.series(STORED_NITROGEN_KEY)))
    worst_tight = float(np.min(tight.series(STORED_NITROGEN_KEY)))
    assert abs(worst_tight) < abs(worst_default) / 100.0 or worst_default >= 0.0, (
        f"the store's dip went {worst_default:.3e} -> {worst_tight:.3e} across two orders of "
        "solver tolerance. A dip that does NOT collapse is an over-draw in the proportional "
        "split -- a physics defect, and the loosened atol above would then be hiding it"
    )
    assert float(np.max(tight.series(STORED_NITROGEN_KEY))) == pytest.approx(
        float(np.max(traj.series(STORED_NITROGEN_KEY))), rel=1e-4
    ), (
        "the store's PEAK moved with solver tolerance, so the two runs are not the same "
        "trajectory read at two accuracies and the comparison above is not about noise"
    )
    assert_conserved(
        traj,
        total_nitrogen(
            compiled.schema,
            biomass_nitrogen_fraction=compiled.param_values["biomass_N_fraction"],
        ),
        rtol=1e-8,
        label="total_nitrogen",
    )


def test_the_store_is_in_no_charge_balance_at_all(wine_params):
    """Structural, not a magnitude: loading the store must not move the solved pH by one bit.

    This is the property the slot exists for, asserted directly rather than inferred from a run.
    The same mutation applied to ``N`` MUST move pH — otherwise this test would pass in a build
    where the charge balance had stopped reading nitrogen at all, and would be pinning nothing.
    """
    schema = wine_schema()
    base = schema.pack(
        {"X": 1.0, "S": [100.0], "E": 40.0, "N": 0.15, "T": 293.15, "CO2": 5.0, "tartaric": 3.0}
    )
    ph_base = acidbase.ph_of_state(base, schema, wine_params)

    stored = base.copy()
    stored[schema.slice(STORED_NITROGEN_KEY)] = 0.30  # twice the ammonium, held intracellularly
    assert acidbase.ph_of_state(stored, schema, wine_params) == ph_base, (
        "loading the intracellular store moved the solved pH — it is in a charge balance, which "
        "is precisely the D-248 defect this slot was minted to remove"
    )

    # The non-vacuity arm: the same nitrogen in `N` DOES move pH.
    extracellular = base.copy()
    extracellular[schema.slice("N")] = 0.45
    assert acidbase.ph_of_state(extracellular, schema, wine_params) != ph_base, (
        "adding nitrogen to `N` did not move pH either, so the test above pins nothing: this "
        "build's charge balance does not read nitrogen and D-209 has been undone"
    )


# -- 3. The split: growth's draw and the swap's refund must agree ---------------------------


def test_the_swap_cannot_refund_to_a_pool_growth_did_not_debit(wine_params):
    """The second door the charge artefact could come back through.

    Growth draws ``f_N·base_dx`` split across ``N`` and the store. The D-32 swap refunds a
    fraction ``ψ·gate`` of that same draw. If the refund were booked wholly to ``N`` while the
    draw was split, net ``dN/dt`` could go POSITIVE and the must's ammonium would rise again —
    with the uptake Process innocent. Driven at a state where the store holds most of the
    nitrogen, which is where an unsplit refund is worst.
    """
    schema = wine_schema()
    y = schema.pack(
        {
            "X": 2.0,
            "S": [150.0],
            "E": 40.0,
            "N": 0.02,  # nearly all the reachable nitrogen is INSIDE the cells
            STORED_NITROGEN_KEY: 0.40,
            "T": 293.15,
            "CO2": 5.0,
            "amino_acids": 0.4,
            "amino_acids_generic": 0.4,
        }
    )
    assert biomass_growth_rate(y, schema, wine_params) > 0.0, (
        "growth is stopped at this state, so neither Process contributes and the test is vacuous"
    )

    pset = ProcessSet(schema, [GrowthNitrogenLimited(), AminoAcidAssimilation()], strict=True)
    d = pset.total_derivatives(0.0, y, wine_params)
    assert float(d[schema.slice("N")][0]) <= 0.0, (
        f"net dN/dt is {float(d[schema.slice('N')][0]):.6e} > 0 with only growth and the swap "
        "active: the swap is refunding to a pool growth did not debit in the same proportion, "
        "and the must's ammonium can rise again"
    )


def test_the_split_is_proportional_and_totals_the_undivided_flux():
    """The helper's two properties: it conserves the flux, and it splits it by the holdings.

    Proportional rather than store-first is deliberate — a preferential draw puts a C0 kink
    exactly where the store empties, and the BDF Jacobian probe straddles that kind of gate.
    """
    schema = wine_schema()
    y = schema.zeros()
    y[schema.slice("N")] = 0.03
    y[schema.slice(STORED_NITROGEN_KEY)] = 0.09  # 25 % / 75 %

    d = schema.zeros()
    add_assimilable_nitrogen(d, y, schema, -0.4)
    to_n = float(d[schema.slice("N")][0])
    to_store = float(d[schema.slice(STORED_NITROGEN_KEY)][0])

    assert to_n + to_store == pytest.approx(-0.4, rel=1e-12), "the split lost or created nitrogen"
    assert to_n == pytest.approx(-0.1, rel=1e-12)
    assert to_store == pytest.approx(-0.3, rel=1e-12)

    # An empty store books the whole flux to `N` — the pre-D-250 form, term for term.
    y[schema.slice(STORED_NITROGEN_KEY)] = 0.0
    d = schema.zeros()
    add_assimilable_nitrogen(d, y, schema, -0.4)
    assert float(d[schema.slice("N")][0]) == pytest.approx(-0.4, rel=1e-12)
    assert float(d[schema.slice(STORED_NITROGEN_KEY)][0]) == 0.0


def test_growth_is_limited_by_both_pools_together(wine_params):
    """A cell with a full store and an empty must must still grow.

    If growth's Monod had been left reading ``N`` alone, D-248's whole repair would unwind the
    moment the surplus moved out of ``N``: uptake would fill a store nothing could spend.
    """
    schema = wine_schema()
    y = schema.pack({"X": 1.0, "S": [150.0], "E": 20.0, "N": 0.0, "T": 293.15, "CO2": 5.0})
    assert biomass_growth_rate(y, schema, wine_params) == 0.0, "no nitrogen anywhere ⇒ no growth"

    y[schema.slice(STORED_NITROGEN_KEY)] = 0.20
    assert biomass_growth_rate(y, schema, wine_params) > 0.0, (
        "growth is still zero with a full intracellular store and an empty must — the Monod is "
        "reading `N` alone, so D-248's uptake now fills a pool nothing can spend"
    )
    assert assimilable_nitrogen_pools(y, schema) == (0.0, 0.20)


# -- 4. Scope: what stays blind, and what the beer schema never gains -----------------------


def test_the_bacteria_are_blind_to_the_store_on_purpose():
    """MLF and Brett must NOT read a pool that is inside a yeast cell.

    D-248 recorded the bacterial starvation as an open residue (``Flags: D-100``) on the reading
    that ``N`` was extracellular ammonium the bacteria could not see. Under D-250 that reading is
    withdrawn: the pool is intracellular, so the blindness is correct and the residue closes.
    What the model actually lacks is a bacterial nitrogen source yeast do not compete for
    (peptides). Pinned so a later beat cannot "fix" it back.
    """
    for process in (MalolacticGrowth(), BrettGrowth()):
        assert STORED_NITROGEN_KEY not in process.touches, (
            f"{process.name} now draws the yeast's intracellular nitrogen store. Bacteria cannot "
            "reach inside a yeast cell; the missing substrate is peptides, not this store"
        )
        assert STORED_NITROGEN_KEY not in process.touches_where_present


def test_beer_never_gains_the_slot_and_the_exemption_is_scoped():
    """The store mirrors the wine-only amino-acid ledger, and ``touches_where_present`` is why.

    Growth is the one primary-fermentation Process wired into both media. Declaring the store in
    ``touches`` outright would make every beer ProcessSet raise on an unknown variable; declaring
    it nowhere would make every wine ProcessSet raise on a leak. The medium-conditional
    declaration is the third option, and it must not have quietly become a blanket exemption.
    """
    assert STORED_NITROGEN_KEY not in beer_schema()
    assert STORED_NITROGEN_KEY in wine_schema()
    assert GrowthNitrogenLimited.touches == ("X", "S", "N")
    assert GrowthNitrogenLimited.touches_where_present == (STORED_NITROGEN_KEY,)

    class _Leaky(Process):
        name = "leaky"
        tier = Tier.SPECULATIVE
        touches = ("X",)

        def derivatives(self, t, y, schema, params):
            d = schema.zeros()
            d[schema.slice("X")] = 1.0
            d[schema.slice(STORED_NITROGEN_KEY)] = 1.0
            return d

    schema = wine_schema()
    pset = ProcessSet(schema, [_Leaky()], strict=True)
    with pytest.raises(ValueError, match="undeclared variables"):
        pset.total_derivatives(0.0, schema.zeros(), {})


# -- 5. The reading D-250 makes possible for the first time ---------------------------------


def test_the_extracellular_reading_is_now_separable_and_d249s_verdict_survives_it():
    """Crépin sampled the MEDIUM. Before D-250 the model could not report that separately.

    ``tests/test_defined_media._assimilable_n_mgl`` keeps counting the store, because that is
    what D-248's "40.8 % → 0.62 % residual" meant and dropping it would redefine the number
    rather than migrate it. But the narrower quantity — nitrogen still OUTSIDE the cells, which
    is what her Data Set S1 measures — is newly readable, and it exhausts sooner. D-249's
    conclusion is that the nitrogen channel is slower than the fermentation containing it; that
    survives the narrower reading, which is the point of pinning it here.

    **At the sourced pitch this is the reading that moved most (decision D-253).** The gap was
    1.71× when the fixture carried the house 0.25 g/L inoculum; on the sourced 0.04 it is
    **0.98×** — the medium's nitrogen exhausts at 28.6 h against Crépin's measured N_T of 28 h.
    So on the one frame she sampled, at the one inoculum her lab published, the exhaustion clock
    agrees to 2 %. D-249's verdict is unaffected: 0.98× is still comfortably under the run's own
    1.54×, and the comparison below is what carries that claim, not the number above it.
    """
    compiled = compile_scenario(commensurate_scenario("crepin"))
    for name in _OTHER_PRECURSOR_CONSUMERS:
        if name in compiled.process_set:
            compiled.process_set.disable(name)
    traj = simulate_scheduled(
        compiled.process_set,
        compiled.param_values,
        compiled.y0,
        compiled.t_span_h,
        events=compiled.events,
        param_tiers=compiled.parameters.tier_map(),
    )
    assert traj.success, traj.message
    schema = compiled.schema

    total = np.array([_assimilable_n_mgl(traj, schema, i) for i in range(traj.y.shape[1])])
    store = np.maximum(traj.y[schema.slice(STORED_NITROGEN_KEY), :][0], 0.0) * 1000.0
    outside = total - store

    def _hours_to(series: np.ndarray, fraction: float) -> float:
        """First crossing, INTERPOLATED between the bracketing samples.

        The same idiom ``test_nitrogen_timing_attribution._cross`` uses, and not a detail: taking
        the sample index outright reads this crossing late by however coarse the output grid is
        there, which is what makes the number incomparable with the one the record it is scored
        against derived. Error scales with the local slope, and the nitrogen curve is steep
        exactly here. At the sourced pitch (D-253) the two idioms differ by 0.10 h — 28.61 h
        interpolated against 28.70 h off the index, a ratio of 0.979 against 0.975 — where at the
        house pitch the same choice was worth 0.5 h and moved the ratio 1.71 → 1.66. The
        sensitivity shrank because the run is slower there, not because the idiom stopped
        mattering; keep the interpolation.
        """
        consumed = 1.0 - series / series[0]
        hit = np.nonzero(consumed >= fraction)[0]
        assert hit.size, "the run never consumes that much nitrogen"
        i = int(hit[0])
        if i == 0:
            return float(traj.t[0])
        return float(
            np.interp(
                fraction,
                [float(consumed[i - 1]), float(consumed[i])],
                [float(traj.t[i - 1]), float(traj.t[i])],
            )
        )

    sugar = traj.y[schema.slice("S"), :].sum(axis=0)
    dry = np.nonzero(sugar <= 2.0)[0]
    assert dry.size, "the run never reaches dryness, so there is no clock to divide out"
    j = int(dry[0])
    # Sugar falls, so the bracketing pair is reversed to give np.interp an increasing x.
    run_h = (
        float(traj.t[0])
        if j == 0
        else float(
            np.interp(
                2.0,
                [float(sugar[j]), float(sugar[j - 1])],
                [float(traj.t[j]), float(traj.t[j - 1])],
            )
        )
    )

    outside_h, total_h = _hours_to(outside, 0.9), _hours_to(total, 0.9)
    assert outside_h < total_h, (
        "the medium's nitrogen does not exhaust before the model's total — the store is holding "
        "nothing, so this reading is not separable after all and the test is vacuous"
    )

    # D-249's verdict, re-derived on the narrower quantity: the nitrogen channel is still SLOWER
    # than the run containing it (Crépin's N_T = 28 h against her EF = 150 h).
    nitrogen_gap = 28.0 / outside_h
    clock_gap = 150.0 / run_h
    # D-250 measured 1.71x here at the house pitch of 0.25 g/L. D-253 moved the fixture onto the
    # sourced inoculum and this is the reading that moved most: **0.98x**, i.e. the medium's
    # nitrogen now exhausts within 2 % of Crepin's own N_T, on the frame she actually sampled.
    # The pin is the number, not the verdict; the verdict is the comparison below, and it holds
    # at both pitches (tests/test_nitrogen_timing_attribution.py measures both).
    assert nitrogen_gap == pytest.approx(0.98, abs=0.05), (
        f"the extracellular nitrogen gap reads {nitrogen_gap:.2f}x, not the 0.98 D-253 measured "
        "at the sourced pitch. A jump back toward 1.7x means the fixture is pitched at the "
        "house 0.25 again, which reverses D-253 rather than loosening a tolerance"
    )
    assert nitrogen_gap < clock_gap, (
        f"read outside the cells the nitrogen gap is {nitrogen_gap:.2f}x against the run's own "
        f"{clock_gap:.2f}x. D-249's attribution rests on the channel being slower than its own "
        "fermentation, and on this reading it no longer is"
    )
