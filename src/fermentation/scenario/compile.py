"""The scenario → core compile seam.

A :class:`~fermentation.scenario.schema.Scenario` is declarative and expressed in
industry units (degrees Brix, mg/L of nitrogen, degrees C, days).
:func:`compile_scenario` turns it into everything the runtime needs to integrate:

    * ``y0``         — the initial state vector in canonical units (g/L, K),
    * ``process_set`` — the medium's Processes, assembled against its schema,
    * ``parameters``  — the provenance-backed parameter set for that medium/strain.

This is the *only* place industry units cross into the canonical internal
representation (decision D-3); the core never sees a degree Brix. Physics does not
live here — it stays in the core's Processes — so this module is pure plumbing:
look up the medium, convert the initial composition, load the parameters, and
assemble the Process set.

The accepted ``Scenario.initial`` keys are validated here (the schema deliberately
leaves them as a free ``dict`` so the vocabulary can live at this boundary):

    wine: brix, yan_mgl, pitch_gpl, [ethanol_gpl]
    beer: glucose_gpl, maltose_gpl, maltotriose_gpl, yan_mgl, pitch_gpl, [ethanol_gpl]

Beer's three sugars are given explicitly rather than split from a single original
gravity: that wort spectrum is a provenance-backed parameter (Milestone 1's
sourcing task), not a magic constant to bury in the compile step.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from fermentation.core import acidbase
from fermentation.core.chemistry import sugar_species
from fermentation.core.kinetics import (
    AcetaldehydeBridgedCondensation,
    AceticAcidOverflow,
    AminoAcidAssimilation,
    AnthocyaninFading,
    AntioxidantBurstOxidation,
    AutolyticHydrogenSulfide,
    AutolyticMercaptan,
    BiomassCarryingCapacity,
    BoundHydrogenSulfideRelease,
    BoundMethanethiolRelease,
    BrettDeath,
    BrettDecarboxylation,
    BrettEthanolToxicity,
    BrettGrowth,
    BrettVinylphenolReduction,
    Caramelization,
    ClosureOxygenIngress,
    EllagitanninOxidation,
    EsterHydrolysis,
    EthylAcetateEsterification,
    EthylHexanoateHydrolysis,
    FuselAminoAcidReroute,
    IsoAlphaAcidLoss,
    MaillardBrowning,
    MaillardStrecker,
    MalolacticCitrateMetabolism,
    MalolacticConversion,
    MalolacticDeath,
    MalolacticGrowth,
    MalolacticSenescence,
    OakExtraction,
    OenococcusDiacetylReduction,
    OrganicAcidExcretion,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    PrecursorNonEhrlichFates,
    SMMHydrolysis,
    StreckerDegradation,
    SulfiteOxidation,
    TanninAnthocyaninCondensation,
    TanninEthylTanninCondensation,
    TanninSelfPolymerization,
    ThermalAnthocyaninFade,
    WortAcidRemoval,
    WortOxygenUptake,
    YeastAutolysis,
    YeastPOFDecarboxylation,
)
from fermentation.core.kinetics.amino_acid_pools import (
    AMINO_ACID_SPECS,
    GENERIC_POOL,
    AminoAcidSpec,
)
from fermentation.core.kinetics.carbon_routing import ESTER_SPECS, FUSEL_SPECS
from fermentation.core.kinetics.hops import iso_alpha_fraction
from fermentation.core.kinetics.oxidative_cascade import (
    OxygenActivation,
    PeroxideEthanolOxidation,
    PeroxideSulfiteOxidation,
    QuinoneAnthocyaninFading,
    QuinoneAscorbateReduction,
    QuinoneEllagitanninOxidation,
    QuinoneHydrogenSulfideCapture,
    QuinonePolymerization,
    QuinoneStreckerDegradation,
    QuinoneSulfonation,
)
from fermentation.core.kinetics.temperature import RAMP_RATE
from fermentation.core.media import Medium, get_medium
from fermentation.core.process import ProcessSet
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier, combine
from fermentation.parameters.schema import Parameter, Provenance, Uncertainty
from fermentation.parameters.store import ParameterSet, default_data_dir, load_parameters
from fermentation.runtime.ensemble import Ensemble, simulate_ensemble
from fermentation.runtime.schedule import ScheduledEvent, ScheduledTrajectory, simulate_scheduled
from fermentation.scenario.schema import Intervention, Scenario
from fermentation.units.convert import (
    brix_to_sugar_gpl,
    celsius_to_kelvin,
    days_to_hours,
    mgl_to_gpl,
    ugl_to_gpl,
)

#: Coleman Y_X/N regression coefficients (decision D-14). Present iff a medium
#: ships the nitrogen-dependent biomass yield; gates the compile-time override.
_N_YIELD_COEFFS = ("biomass_N_yield_log_intercept", "biomass_N_yield_log_slope")

#: The malolactic Processes gated on an *Oenococcus oeni* pitch (decisions D-23, D-31, D-39):
#: malate→lactate conversion, the citrate co-metabolism feeding the diacetyl reservoir, the
#: bacterial diacetyl reduction, and bacterial death/decay. They are wired into the wine medium
#: but contribute nothing until bacteria are present, so the compile step DISABLES them when
#: unpitched and the ``pitch_mlf`` intervention (decision D-36) re-enables *exactly* this set at
#: its breakpoint — a single source of truth so the compile-time gate and the mid-run pitch cannot
#: drift apart. :class:`MalolacticDeath` (D-39) and :class:`MalolacticSenescence` (MLF v2, D-41) are
#: pitch-gated too (bacteria die/age whether or not amino acids were dosed), unlike
#: :class:`MalolacticGrowth`, which is amino-acid-gated below.
_MLF_GATED_PROCESSES = (
    MalolacticConversion,
    MalolacticCitrateMetabolism,
    OenococcusDiacetylReduction,
    MalolacticDeath,
    MalolacticSenescence,
)

#: The *Brettanomyces* Processes gated on a Brett pitch (decision D-40): hydroxycinnamate
#: decarboxylation, vinylphenol reduction, the SO₂-driven :class:`BrettDeath` (D-40 pt3), and the
#: ethanol-toxicity :class:`BrettEthanolToxicity` (D-58). Wired into the wine medium but
#: contributing nothing until Brett is present, so compile DISABLES them (unpitched) and the
#: ``pitch_brett`` intervention re-enables *exactly* this set at its breakpoint — one source of
#: truth so the compile-time gate and the mid-run pitch cannot drift apart (the
#: ``_MLF_GATED_PROCESSES`` pattern). :class:`BrettDeath`/:class:`BrettEthanolToxicity` are
#: pitch-gated too (Brett dies whether or not amino acids were dosed), unlike :class:`BrettGrowth`,
#: which is amino-acid-gated below (the exact :class:`MalolacticDeath` vs
#: :class:`MalolacticGrowth` split).
_BRETT_GATED_PROCESSES = (
    BrettDecarboxylation,
    BrettVinylphenolReduction,
    BrettDeath,
    BrettEthanolToxicity,
)

#: The aging Processes ``begin_aging`` enables (decisions D-70/D-71/D-72/D-74/D-75/D-77/D-78/D-79):
#: :class:`EsterHydrolysis` (the banana-acetate fade), :class:`EthylHexanoateHydrolysis` (the D-126
#: apple-ethyl-ester fade — the sibling hydrolysis, Makhotkina & Kilmartin 2012),
#: :class:`OxidativeAcetaldehyde` (the O₂-driven ethanol
#: oxidation), :class:`PhenolicBrowning` (the O₂-driven browning, D-74, accumulating ``A420``) and
#: :class:`SulfiteOxidation` (the O₂-driven SO₂ scavenging, D-72) and :class:`StreckerDegradation`
#: (the O₂/amino-acid-driven Strecker aldehydes, D-75). The first three are medium-agnostic
#: (wired into both media); :class:`SulfiteOxidation`, :class:`StreckerDegradation`,
#: :class:`OakExtraction` (the NON-oxidative barrel/chip aroma extraction, D-77 — a separate axis
#: drawing no O2) and :class:`EllagitanninOxidation` (the D-78 oak-tannin O₂ scavenging — oak
#: PROTECTION, the bridge from the oak axis to the O₂ sub-axis) are wine-only (they read wine-only
#: ``so2_total``/pH, ``amino_acids``/``N``, the oak ceiling/extractive slots and the
#: ``ellagitannin`` pool respectively), as are :class:`TanninAnthocyaninCondensation` (the D-79
#: red-wine colour-stabilization + astringency softening — grape ``anthocyanin`` + grape ``tannin``
#: condense to stable polymeric pigment; a NON-oxidative grape axis drawing neither O₂ nor oak) and
#: :class:`AcetaldehydeBridgedCondensation` (the D-80 SPLIT-LEDGER beat — dissolved-O₂ acetaldehyde
#: bridges grape tannin to anthocyanin, the first link from the oxidative sub-axis to red-wine
#: colour; its ``ethyl_bridge`` slot captures the acetaldehyde carbon on the ledger) and
#: :class:`AnthocyaninFading` (the D-81 O₂-coupled bleaching loss — dissolved O₂ fades free grape
#: anthocyanin to the colourless ``faded_anthocyanin`` slot, so colour genuinely declines and SO₂
#: protection is emergent via the shared o2 pool) and :class:`ThermalAnthocyaninFade` (the D-83
#: O₂-INDEPENDENT thermal fade — heat degrades free anthocyanin to the same ``faded_anthocyanin``
#: slot with NO oxygen, so a sealed/sulfited/anaerobic red still fades and SO₂ gives no protection)
#: and :class:`TanninSelfPolymerization` (the D-84 tannin–tannin axis — grape ``tannin`` condenses
#: with itself, ``[tannin]²``, into a soft polymer, softening astringency WITHOUT anthocyanin) and
#: :class:`TanninEthylTanninCondensation` (the D-85 acetaldehyde-bridged tannin–ethyl–tannin route —
#: dissolved-O₂ acetaldehyde bridges two flavanols, an O₂-driven softener that captures its carbon
#: in the shared ``ethyl_bridge`` slot and deposits no pigment) and :class:`MaillardStrecker` (the
#: D-87 NON-oxidative THERMAL Strecker route — residual sugar + heat, no O₂, degrade amino acids to
#: the sweet-wine/Madeira aldehyde suite; wine-only, reads ``amino_acids`` + deaminates to ``N``)
#: and
#: :class:`Caramelization` (the D-88 sugar-only THERMAL browning — residual sugar browns to the
#: on-ledger ``melanoidin`` carbon-park by heat with no O₂, raising the shared ``A420``;
#: MEDIUM-AGNOSTIC per D-90 — beer's residual dextrins caramelize too, the vectorized draw
#: apportions across beer's 3-slot ``S``; the first aging Process to consume core ``S``) and
#: :class:`MaillardBrowning` (the D-89
#: amino-acid-incorporating THERMAL browning — residual sugar + amino acids brown to the on-ledger
#: N-bearing ``maillard_melanoidin`` carbon+nitrogen-park by heat with no O₂, raising the same
#: ``A420``; wine-only, the first aging Process on the nitrogen ledger), so
#: on beer they are simply absent from the ProcessSet — both the compile-disable and the
#: begin_aging-enable loops guard with ``name in process_set``, so listing them here is beer-safe.
#: All are DISABLED unconditionally at compile (aging is inherently post-ferment); the
#: ``begin_aging`` verb re-enables exactly this tuple at its breakpoint and the compile seam
#: disables exactly this tuple — one list,
#: so the enable/disable stay symmetric as the aging axis grows. Their shared
#: aging.yaml/oak.yaml/polymerization.yaml parameters are guarded together at the verb boundary.
_AGING_GATED_PROCESSES = (
    EsterHydrolysis,
    EthylHexanoateHydrolysis,
    EthylAcetateEsterification,
    OxidativeAcetaldehyde,
    PhenolicBrowning,
    SulfiteOxidation,
    AntioxidantBurstOxidation,
    StreckerDegradation,
    OakExtraction,
    EllagitanninOxidation,
    TanninAnthocyaninCondensation,
    AcetaldehydeBridgedCondensation,
    AnthocyaninFading,
    ThermalAnthocyaninFade,
    TanninSelfPolymerization,
    TanninEthylTanninCondensation,
    MaillardStrecker,
    Caramelization,
    MaillardBrowning,
    SMMHydrolysis,
    BoundHydrogenSulfideRelease,
    BoundMethanethiolRelease,
    ClosureOxygenIngress,
    # The D-141 oxidative cascade. Both alternatives are listed: the enable loop guards on
    # ``name in process_set``, so whichever set a given build wired is switched on and the other
    # is skipped silently — which is exactly what makes the two isolable at this seam too.
    # (That same guard is why AntioxidantBurstOxidation above has never run: it is listed here
    # and wired into no medium, so the loop skips it every time. Pinned at D-140, still open.)
    OxygenActivation,
    PeroxideEthanolOxidation,
    PeroxideSulfiteOxidation,
    QuinoneSulfonation,
    QuinoneStreckerDegradation,
    QuinoneAnthocyaninFading,
    QuinoneEllagitanninOxidation,
    QuinonePolymerization,
    # D-201's sulfide sink. Registering it here is NOT bookkeeping: every Process in the cascade
    # is disabled at the compile seam and re-enabled by ``begin_aging`` off THIS tuple, so a
    # cascade Process omitted from it is never disabled and therefore runs from t = 0 — through
    # fermentation, against a quinone pool that is zero then, but with no gate saying so. D-200's
    # throwaway probe hit exactly this and had to gate itself on ``quinone_sulfonation`` instead.
    QuinoneHydrogenSulfideCapture,
    # D-202's ascorbate route. It needs the aging gate for the same reason, and the reason is NOT
    # weakened by its pool defaulting to 0: a scenario that doses ascorbate BEFORE ``begin_aging``
    # (the natural way to model an addition at crush or at bottling) would otherwise have it
    # scavenging quinone during fermentation.
    QuinoneAscorbateReduction,
)

#: A name → value(s) mapping ready for :meth:`StateSchema.pack`.
_Initial = dict[str, float | list[float]]


@dataclass(frozen=True, eq=False)
class CompiledScenario:
    """Everything the runtime needs to integrate one scenario.

    Realises the documented compile seam ``(y0, ProcessSet, params)`` as a named
    record, plus the schema and time span that travel with them. ``param_values``
    is the plain ``{name: float}`` mapping ``simulate`` and ``Process.derivatives``
    consume; ``parameters`` retains the full provenance and tier information for
    honest downstream reporting.
    """

    scenario: Scenario
    schema: StateSchema
    y0: FloatArray
    process_set: ProcessSet
    parameters: ParameterSet
    t_span_h: tuple[float, float]
    #: Timed interventions compiled from the scenario, in canonical hours, ready to hand to
    #: :func:`fermentation.runtime.simulate_scheduled`: the temperature-schedule slope-change
    #: events (decision D-35) merged with the discrete winemaking verbs — ``add_dap`` / ``add_so2``
    #: / ``rack`` / ``pitch_mlf`` (decision D-36) and the ``begin_aging`` aging-phase switch
    #: (decision D-70). Empty ⇒ an un-scheduled run (plain :func:`simulate` suffices).
    events: tuple[ScheduledEvent, ...] = field(default_factory=tuple)

    @property
    def param_values(self) -> dict[str, float]:
        """Resolved ``{name: value}`` mapping for the integration hot loop."""
        return self.parameters.resolve()

    def run(self, **kwargs: object) -> ScheduledTrajectory:
        """Integrate this scenario, **honouring its timed events** (decision D-35).

        The single "run a compiled scenario" entry point. It always dispatches through
        :func:`~fermentation.runtime.simulate_scheduled` with the compiled ``events``, so a
        temperature ramp (and, from D-36, a dosing/pitching schedule) is applied correctly —
        a multi-knot ramp changes slope at its breakpoints, a hold holds. With no events this
        is byte-for-byte a plain :func:`~fermentation.runtime.simulate` (an empty schedule is a
        single segment), so it is always the right call.

        This exists because a hand-wired ``simulate(cs.process_set, cs.param_values, cs.y0,
        cs.t_span_h)`` **silently ignores** ``events``: the injected ``temperature_ramp_rate``
        would then apply the *first* segment's slope for the whole run — correct only for a
        single-slope ramp. ``param_tiers`` defaults to the scenario's own tier map for honest
        D-1 reporting; ``t_eval``/solver kwargs pass straight through.

        **Calling this twice on one compiled scenario does not repeat the run** — the events'
        reconfigure persists, so the second call starts with the first's enables live from
        ``t = 0``. See :func:`~fermentation.runtime.simulate_scheduled` for why that is the
        contract and what to bracket a reused set with (decision D-206).

        The stochastic sibling is :meth:`run_ensemble`, which threads the same ``events`` into
        :func:`~fermentation.runtime.simulate_ensemble` (decision D-37).
        """
        kwargs.setdefault("param_tiers", self.parameters.tier_map())
        return simulate_scheduled(
            self.process_set,
            self.param_values,
            self.y0,
            self.t_span_h,
            events=self.events,
            **kwargs,  # type: ignore[arg-type]
        )

    def run_ensemble(self, **kwargs: object) -> Ensemble:
        """Run a stochastic ensemble of this scenario, **honouring its timed events** (D-37).

        The uncertainty-band counterpart to :meth:`run`: it hands the compiled ``events`` to
        :func:`~fermentation.runtime.simulate_ensemble` so every sampled member is integrated
        through the same schedule (temperature ramp + dosing/pitching), then reports the nominal
        run plus the median + spread over the parameters' provenance bands. Passes the full
        :class:`~fermentation.parameters.store.ParameterSet` (the ensemble needs the *bands*, not
        just resolved floats) and defaults ``param_tiers`` to the scenario's own tier map. Sampling
        scope, per-member Process-set isolation, and the per-member external-flow ledger are all
        handled by :func:`simulate_ensemble`; ``n_members``/``seed``/``sampler``/``t_eval`` and the
        solver kwargs pass straight through.
        """
        kwargs.setdefault("param_tiers", self.parameters.tier_map())
        kwargs.setdefault("y0_for_member", self.y0_for_member())
        return simulate_ensemble(
            self.process_set,
            self.parameters,
            self.y0,
            self.t_span_h,
            events=self.events,
            **kwargs,  # type: ignore[arg-type]
        )

    def y0_for_member(
        self,
    ) -> Callable[[Mapping[str, float]], FloatArray] | None:
        """Per-member ``y0`` builder for the parts of the seed a **parameter** derives.

        Was ``reanchor_for_member`` and did only the pH anchor (D-233); D-236 added the copper
        seed and the name followed the scope, and D-238 the peptide buffer capacity. It is still
        **not** a re-run of the initial
        builder: it rebuilds exactly the slots whose compile-time value is a sampled parameter,
        one measured rule at a time, and leaves everything else at the compiled array. A full
        rebuild would move seeds no beat has measured, which is the trap D-233 declined and this
        one does not re-open.

        ``None`` when no rule applies, so a caller stays byte-identical to the fixed-``y0`` path.

        **Rule 1 — the pH anchor (D-233).** Applies when the schema has ``cation_charge`` and the
        scenario gave an ``initial_ph``. ``initial_ph`` is an *input* and its contract is that the
        model reproduces it at t=0 (D-18 for wine, D-179/D-180 for beer). The engine honours it by
        back-solving ``cation_charge`` at COMPILE — and that back-solve reads the **pKa map**,
        which the ensemble samples (every ``pKa_*`` reaches the sampled set through the ``reads``
        of any active pH-reading Process, D-160). So a member drew its own pKas and then started
        from a cation fitted to somebody else's: measured, beer members began at pH 5.5062-5.7778
        against an anchor of 5.65 (worst miss **0.1438**), and wine at 3.4208-3.5780 against 3.50
        (worst **0.0792**), where the correct spread is **exactly zero**. Re-solving restores it
        to 2.3e-11 pH. :func:`~fermentation.core.acidbase.cation_charge_for_ph` is the exact
        inverse of ``ph_of_state`` (D-186), and at t=0 — ``Byp`` 0, no evolved CO2 — it reduces
        term for term to the compile seam's own ``solve_cation_charge``, so at the nominal draw it
        reproduces the compiled slot rather than competing with it (pinned as a test).

        *Scope, stated because the number is smaller than it first looks.* The reported *band* is
        essentially unchanged — 1.008x at day 14 across all 83 sampled parameters, not the 1.287x
        a single-parameter sweep suggests (that figure is ``pKa_peptide_buffer``'s own
        contribution with everything else nominal and must never be quoted as the band's). What it
        repairs is each member's own trajectory: worst per-member day-14 shift 0.0346 pH. The case
        for it is the t=0 contract, never the spread. The **second** pH anchor, ``set_ph``'s, is
        not here at all — it is a mutation, repaired at D-235 by handing mutations the running map.

        **Rule 2 — the copper seed (decision D-236).** Applies to a wine that did **not** name
        ``copper_gpl``, where the seam seeds the ``copper`` slot from ``copper_typical`` itself.
        D-134 made those two numerically identical on purpose so that ``PhenolicBrowning``'s
        mean-centred ``f(Cu) = 1 + k·(copper − copper_typical)`` is **exactly** 1 for an
        un-overridden wine. A drawn member moved the reference and not the seed, so the design
        invariant held at the nominal draw and nowhere else: with ``copper_typical`` the only name
        sampled, aged ``A420`` moved **16.65 %** across 12 members — and D-234 measured the
        coherent channel of the same parameter at **exactly zero**, so every bit of that was the
        broken cancellation, carrying the sign the parameter's name argues against (a wine whose
        *typical* copper is drawn higher browned *less*). Re-seeding restores ``f(Cu) == 1`` per
        member, and with it bit-identical members.

        *Why the condition is a branch and not a blanket.* A scenario that names ``copper_gpl`` is
        stating this wine's copper, which is a genuinely independent quantity from the reference
        the multiplier is centred on — there, drawing the reference alone is **correct**, and
        re-seeding would overwrite a scenario input and breach D-24's surviving exclusion. Both
        halves of the branch are guarded. The condition is checked twice over, on the scenario key
        *and* on the compiled slot still holding the nominal ``copper_typical``, so a future seam
        that stops deriving this slot silently stops the rule instead of silently overwriting.

        **Rule 3 — the peptide buffer capacity (decision D-238), and it runs FIRST.**
        ``peptide_buffer_capacity_beer`` is back-solved OFFLINE against Peyer's published wort
        BC = 1.18 at the *nominal* ``pKa_peptide_buffer``, then shipped as a compile-time seed —
        while that pKa is read at runtime and is drawn. So a member carried a wort whose buffering
        capacity was 1.1161-1.180 rather than the 1.18 the constant exists to reproduce (D-214,
        re-measured at D-233: 0.0100 pH at day 14 on the low-pKa arm, 21 % of that arm's defect).
        :func:`~fermentation.core.acidbase.peptide_capacity_for_wort_bc` re-roots it on the
        member's own map, and at the nominal draw it returns the shipped literal **bit-for-bit**.

        *Why D-233's reason for declining this no longer applies.* It declined the repair because
        moving the root-find into ``src`` would make
        ``test_the_peptide_capacity_still_reproduces_peyers_published_wort_bc`` compare the
        root-finder against itself. That is true of *deriving the shipped constant*, and this does
        not: the YAML literal is untouched and still scored against Peyer's published 1.18 by a
        titration neither side produces. The teeth that forced the D-180 and D-181 re-anchors are
        in those two literals, not in owning a second copy of the arithmetic.

        *It runs before rule 1, and that ordering is load-bearing.* ``peptide_buffer`` is an acid
        slot the t=0 cation back-solve reads, so an anchor solved against the *old* capacity would
        put the member at a pH the scenario never asked for — the very defect rule 1 exists to
        close. Rule 1 therefore reads the array under construction rather than the compiled one.
        A wrong order needs no new guard: it turns
        ``test_every_sampled_member_starts_at_the_ph_the_scenario_anchored[beer]`` red, with wine
        (no peptide slot) staying green as the free control.

        *The condition, and it is the same shape as rule 2's.* A scenario that names
        ``peptide_buffer_gpl`` is stating this wort's buffering protein — a scenario INPUT, which
        D-24 excludes from sampling — so the rule does not fire, and the compiled slot is checked
        for the nominal capacity too so a future seam that stops deriving it stops the rule.

        *Scope: the whole running map feeds the root, deliberately.* The eight ``*_typical_wort``
        acid LEVELS are not in the sampled set at all (measured: 0 of 8), but the other acids'
        **pKas** are, and they move the computed BC. They are fed in, because Peyer's 1.18 is a
        measurement of a real wort and a back-solve against a measurement should see every constant
        that enters it. Priced rather than assumed: worst member 0.40 % of the capacity against a
        peptide-pKa-only root. Filtering the map would also have shipped a deliberately half-pinned
        read *inside* the fix for half-pinning, which is this census's own defect class.
        """
        base = self.y0
        rules: list[Callable[[FloatArray, Mapping[str, float]], None]] = []

        target = self.scenario.initial.get("initial_ph")
        if (
            "peptide_buffer" in self.schema
            and "peptide_buffer_gpl" not in self.scenario.initial
            and "peptide_buffer_capacity_beer" in self.parameters
            and "wort_buffering_capacity_peyer" in self.parameters
            and "pKa_peptide_buffer" in self.parameters
        ):
            peptide_slot = self.schema.slice("peptide_buffer")
            nominal_pka = self.parameters["pKa_peptide_buffer"].value
            target_bc = self.parameters["wort_buffering_capacity_peyer"].value
            if (
                float(base[peptide_slot][0])
                == self.parameters["peptide_buffer_capacity_beer"].value
            ):

                def recapacitate(out: FloatArray, values: Mapping[str, float]) -> None:
                    # The exact-nominal skip is structural, not numerical. The root-find DOES
                    # return the shipped literal bit-for-bit here (asserted as a test, and D-233
                    # §1's one-ULP disagreement was a looser root, not a floor) — but resting
                    # D-24's byte-for-byte nominal claim on a root-finder's tolerance surviving a
                    # scipy upgrade is not the same as guaranteeing it.
                    if values["pKa_peptide_buffer"] == nominal_pka:
                        return
                    # Reads the acid slots off `out` and overwrites only its own; the titration
                    # inverse-solves its own sample cation, so this never reads `cation_charge`
                    # and can safely run before the anchor below.
                    out[peptide_slot] = acidbase.peptide_capacity_for_wort_bc(
                        out, self.schema, values, target_bc
                    )

                rules.append(recapacitate)

        if "cation_charge" in self.schema and target is not None:
            cation_slot = self.schema.slice("cation_charge")
            target_ph = float(target)

            def reanchor(out: FloatArray, values: Mapping[str, float]) -> None:
                # Reads the acid slots off `out` — the array RULE 3 has already re-capacitated —
                # and never the cation slot it is solving for, so it is not circular. Reading
                # `base` here (as this did before D-238) would anchor the member against a
                # peptide pool it no longer carries.
                out[cation_slot] = acidbase.cation_charge_for_ph(
                    out, self.schema, values, target_ph
                )

            rules.append(reanchor)

        if (
            "copper" in self.schema
            and "copper_gpl" not in self.scenario.initial
            and "copper_typical" in self.parameters
        ):
            copper_slot = self.schema.slice("copper")
            seeded = float(base[copper_slot][0])
            if seeded == self.parameters["copper_typical"].value:

                def reseed_copper(out: FloatArray, values: Mapping[str, float]) -> None:
                    out[copper_slot] = values["copper_typical"]

                rules.append(reseed_copper)

        if not rules:
            return None

        def build(values: Mapping[str, float]) -> FloatArray:
            out = base.copy()
            for rule in rules:
                rule(out, values)
            return out

        return build


# -- initial-composition vocabulary (the industry-unit boundary) --------------


def _amino_acid_override_key(spec: AminoAcidSpec) -> str:
    """The scenario key overriding one speciated amino-acid pool's dose (decision D-100).

    ``<species>_gpl`` — named for the **molecule**, not the slot, so a scenario says
    ``arginine_gpl`` rather than leaking the historical ``amino_acids``-slot-is-arginine detail
    (D-100 kept that slot name to avoid touching every consumer twice; the dose API need not
    inherit the compromise). The generic bucket is the one pool with no single molecule, so it
    keeps its slot name: ``amino_acids_generic_gpl``.
    """
    return f"{spec.pool if spec.pool == GENERIC_POOL else spec.species}_gpl"


#: Every per-species amino-acid override key the wine seam accepts (decision D-100).
_AMINO_ACID_OVERRIDE_KEYS: tuple[str, ...] = tuple(
    _amino_acid_override_key(spec) for spec in AMINO_ACID_SPECS
)

#: Keys accepted in ``Scenario.initial`` per medium. Validated at compile time so
#: a typo ("brixx") fails loudly instead of being silently ignored.
_ALLOWED_KEYS: dict[str, frozenset[str]] = {
    "wine": frozenset(
        # tartaric_gpl/malic_gpl/initial_ph are the optional pH-solver inputs (D-18);
        # lactic is produced-only (MLF product) so it is not an input, and the
        # strong cation is back-solved from initial_ph, not given. so2_total_mgl is the
        # optional total-SO₂ dose for the free/bound + molecular-SO₂ readout (D-22/D-28);
        # mlf_pitch_gpl is the optional Oenococcus oeni dose driving malolactic conversion (D-23).
        # carrying_capacity_gpl is the optional opt-in biomass cap K that enables the
        # residual-nitrogen floor (D-30); absent ⇒ the cap modifier is disabled (core untouched).
        # citrate_gpl is the optional citric-acid must input O. oeni co-metabolises into
        # MLF-derived diacetyl (D-31); absent ⇒ 0 (no citrate, the diacetyl branch is silent).
        # amino_acids_gpl is the optional assimilable amino-acid dose the AminoAcidAssimilation
        # swap funds biomass from (D-32); absent ⇒ 0 (the swap Process is disabled, core untouched).
        # autolysis_rate_per_h is the optional opt-in first-order autolysis rate (1/h) that both
        # enables YeastAutolysis and overrides its k_autolysis reference (D-34); absent/0 ⇒ the
        # autolysis Process is disabled (an undosed wine run is byte-for-byte the validated core).
        {
            "brix",
            "yan_mgl",
            "pitch_gpl",
            "ethanol_gpl",
            "tartaric_gpl",
            "malic_gpl",
            "initial_ph",
            "so2_total_mgl",
            # 5-oxofructose botrytis carbonyl SO₂-binder must input (decision D-130); mg/L,
            # absent ⇒ 0 (a non-botrytis must — the binding equilibrium is byte-for-byte the
            # D-51 3-carbonyl form). Set for a botrytized/sweet wine (~80-150 mg/L typical).
            "oxofructose_mgl",
            # Aggregate medium-chain fatty acids (octanoic+decanoic, octanoic-equivalent) that
            # inhibit MLF (decision D-131); mg/L, absent ⇒ 0 (byte-for-byte the pre-D-131 MLF).
            # ~3 mg/L is a normal-AF wine (barely bites); ~15-22 mg/L a stressed/N-deficient AF
            # (halves→stalls MLF — the classic yeast-strain-dependent stuck MLF, Lonvaud-Funel).
            "mcfa_mgl",
            "mlf_pitch_gpl",
            "carrying_capacity_gpl",
            "citrate_gpl",
            "amino_acids_gpl",
            # Per-species amino-acid overrides (decision D-100): each speciated pool's dose can be
            # set directly, overriding its must-spectrum share of amino_acids_gpl. Spread from the
            # canonical registry so a ninth amino acid needs no edit here and cannot be silently
            # unreachable from a scenario.
            *_AMINO_ACID_OVERRIDE_KEYS,
            "autolysis_rate_per_h",
            "hydroxycinnamic_gpl",
            # Ferulic-acid must precursor (decision D-55) — the second, genuinely distinct
            # volatile-phenol branch split out from hydroxycinnamic_gpl (which is booked as
            # p-coumaric acid specifically). Absent ⇒ 0 (that branch is inert), exactly mirroring
            # hydroxycinnamic_gpl's own isolability.
            "ferulic_acid_gpl",
            "brett_pitch_gpl",
            # anthocyanin_gpl / tannin_gpl are the optional GRAPE must inputs driving red-wine
            # tannin–anthocyanin condensation (decision D-79): free monomeric anthocyanin
            # (bleachable
            # red pigment) and condensed grape tannin (harsh young astringency), which
            # TanninAnthocyaninCondensation combines into stable polymeric pigment during aging.
            # Both absent/0 ⇒ a white / no-tannin wine — the Process is byte-for-byte inert (doubly
            # substrate-gated), so an undosed run is unchanged. Off every ledger (grape-derived).
            "anthocyanin_gpl",
            "tannin_gpl",
            # pof_positive is a binary POF+ strain opt-in (decision D-40 pt4): present/>0 enables
            # YeastPOFDecarboxylation (the yeast cinnamate decarboxylase filling vinylphenols during
            # AF), WHOLLY INDEPENDENT of brett_pitch_gpl. Absent/0 ⇒ a POF-negative wine — the
            # Process stays disabled and the run is byte-for-byte the validated core. Not a state
            # slot: it is a compile-time gate only, never packed into y0.
            "pof_positive",
            # The grape's DMS potential in DMS-EQUIVALENT µg/L (decision D-102) — the precursor
            # SMMHydrolysis converts to aged-wine DMS. UNLIKE every other optional key here, absent
            # does NOT mean 0: it falls back to the sourced must level `dms_potential_initial`,
            # because DMSp is a property of the GRAPE that every must carries rather than a
            # winemaking dose (a 0 default would assert aged wine makes no DMS — the D-45 hard-zero
            # defect). Scenarios SHOULD set it: DMSp is strongly variety-dependent and the sourced
            # default is Syrah's, which over-predicts a low-DMSp variety (see dms.yaml). Explicit 0
            # is still honoured, and makes the Process byte-for-byte inert.
            "dms_potential_ugl",
            # The grape's initial-burst antioxidant charge in g/L (decision D-133) — the finite,
            # unidentified, non-SO2 pool AntioxidantBurstOxidation scavenges. Unlike its
            # dms_potential_ugl/copper_gpl siblings, the D-45 "absent does not mean 0" fallback is
            # CONDITIONAL ON THE CONSUMER BEING WIRED (decision D-147): it falls back to the sourced
            # burst_antioxidant_initial under ``oxidative="direct_burst"``, and to 0.0 otherwise.
            # The D-45 argument INVERTS without the Process — with nothing able to draw the pool, a
            # non-zero seed does not assert "this wine has a day-1 burst", it asserts an antioxidant
            # that exists and is never spent, which is strictly worse than 0 and visible in output.
            # Dosing it into a build that cannot consume it is an ERROR, not a silent no-op.
            "burst_antioxidant_gpl",
            # Dissolved copper in g/L-equivalent mg/L input (decision D-134) — the mean-centered
            # multiplier PhenolicBrowning reads. Like dms_potential_ugl/burst_antioxidant_gpl,
            # absent does NOT mean 0: it falls back to the sourced copper_typical, because copper
            # is a property every must/wine carries rather than a winemaking dose (a 0 default
            # would understate a typical wine's browning rate, not merely leave the multiplier
            # isolable at zero). Scenarios SHOULD override it to explore atypical copper wines;
            # explicit 0 is still honoured (an explicit, documented deviation, not silent).
            "copper_gpl",
            # The wine's metal-complexed ("bonded") sulfide reservoirs at bottling, in µg/L
            # (decision D-135) — what BoundHydrogenSulfideRelease / BoundMethanethiolRelease empty
            # into the free h2s / methanethiol pools during anaerobic bottle aging. Like
            # dms_potential_ugl / burst_antioxidant_gpl / copper_gpl, absent does NOT mean 0: they
            # fall back to the sourced bound_h2s_initial / bound_methanethiol_initial, because 94 %
            # and 62 % of a red wine's H2S and MeSH are bonded at bottling (Franco-Luesma & Ferreira
            # 2016) and defaulting to 0 would assert that no bottled wine ever turns reductive —
            # the D-45 hard-zero defect. Scenarios SHOULD override for a WHITE: whites carry MORE
            # bonded H2S (23.7 µg/L) and about HALF the bonded MeSH (0.8 µg/L) of the red defaults.
            # Explicit 0 is still honoured and makes both Processes byte-for-byte inert.
            "bound_h2s_ugl",
            "bound_methanethiol_ugl",
        }
    ),
    "beer": frozenset(
        # initial_ph is the OPT-IN GATE for beer's whole pH system (decision D-179), exactly as
        # it is for wine (D-18): absent ⇒ every acid slot and the cation stay 0, the charge
        # balance is empty and the run is byte-for-byte the pre-D-179 beer. Present ⇒ the acid
        # slots are seeded from their sourced levels and the strong cation is back-solved to
        # reproduce the given pH (inverse anchoring).
        #
        # The eight acid keys override those sourced levels individually — the "hold the must,
        # spike the leucine" shape of the D-100 amino-acid overrides. Each is in g/L.
        # peptide_buffer_gpl overrides the lumped protein buffer; set it to 0 to run the
        # organic-acids-only arm (which measures ~1/6th of real wort's buffering — the honest
        # under-prediction the peptide term exists to close).
        #
        # pyruvic/formic/oxalic (decision D-181) are the three the ferment REMOVES. Setting one
        # to 0 is how the "what does the missing base cost" arm is run without touching code —
        # and it is the only way to recover the pre-D-181 one-sided acidification, which is why
        # they are overridable rather than fixed at their wort levels.
        {
            "glucose_gpl",
            "maltose_gpl",
            "maltotriose_gpl",
            "yan_mgl",
            "pitch_gpl",
            "ethanol_gpl",
            "initial_ph",
            "lactic_gpl",
            "acetic_gpl",
            "citric_gpl",
            "malic_gpl",
            "succinic_gpl",
            "pyruvic_gpl",
            "formic_gpl",
            "oxalic_gpl",
            "peptide_buffer_gpl",
            # Wort aeration (decision D-213), in mg/L to match brewery practice and the source's
            # own units — the ONE key here that is not an acid and not gated on `initial_ph`.
            # Absent falls back to the sourced `o2_wort_aeration_beer`, because a cast wort is
            # always aerated: this is D-147's "the quantity every wort carries" case, the same
            # call as dms_potential/copper, NOT the acids' opt-in. Set it to 0 to recover the
            # pre-D-213 anaerobic-from-t0 beer, which is how the isolability arm is run.
            "o2_mgl",
        }
    ),
}


def _nonneg(value: float, key: str) -> float:
    if value < 0.0:
        raise ValueError(f"scenario.initial[{key!r}] must be >= 0, got {value}")
    return value


def _require(values: Mapping[str, float], key: str, medium: str) -> float:
    if key not in values:
        raise ValueError(f"{medium} scenario.initial is missing required key {key!r}")
    return _nonneg(float(values[key]), key)


def _optional(values: Mapping[str, float], key: str, default: float) -> float:
    return _nonneg(float(values[key]), key) if key in values else default


def _wine_amino_acids(values: Mapping[str, float], parameters: ParameterSet) -> dict[str, float]:
    """Split the assimilable amino-acid dose across the eight speciated pools (decision D-100).

    **Fixed spectrum + per-species overrides (the owner's D-100 dose API).** ``amino_acids_gpl``
    stays the one knob a scenario normally turns: it is apportioned by the sourced
    ``must_aa_fraction_*`` spectrum (recorded from published *Vitis vinifera* must profiles before
    any wiring — the D-99 anti-tuning discipline), so the default composition is a real must's,
    not a modelling convenience. Any pool can then be overridden individually with
    ``<species>_gpl`` — which is what a study of a *specific* precursor needs (spike the leucine,
    hold the rest) and what keeps the fixed spectrum from becoming a straitjacket.

    The fractions are **normalized** rather than asserted to sum to 1: they are eight independent
    provenance entries, each free to be re-sourced on its own, and an ensemble sampling their
    uncertainty bands would otherwise break a sum-to-one assertion on nearly every draw. So
    ``amino_acids_gpl`` is exactly conserved into the pools whatever the fractions are, and the
    dose means "this much assimilable amino acid" — never "this much times whatever the spectrum
    happens to add up to". An override **replaces** that pool's share outright (absolute, not
    additive) and does not re-normalize the others, so it raises or lowers the total assimilable
    dose accordingly — which is what "hold the must, spike the leucine" requires.

    An absent/zero dose leaves every pool at 0 — the isolability guarantee (D-32): every gate
    reads exactly 0 and the run is byte-for-byte the validated core.
    """
    dose = _optional(values, "amino_acids_gpl", 0.0)
    fractions = {spec.pool: parameters[spec.fraction_param].value for spec in AMINO_ACID_SPECS}
    total = sum(fractions.values())
    pools = {pool: dose * fraction / total for pool, fraction in fractions.items()}
    for spec in AMINO_ACID_SPECS:
        key = _amino_acid_override_key(spec)
        if key in values:
            pools[spec.pool] = _nonneg(float(values[key]), key)
    return pools


def _wine_initial(
    values: Mapping[str, float], temperature_k: float, parameters: ParameterSet
) -> _Initial:
    # Brix measures *total* dissolved solids; only ~90-95% of ripe-must solids are
    # fermentable hexose (the rest is acids/minerals/phenolics). The sourced
    # must_fermentable_fraction corrects brix_to_sugar_gpl so a 24 Brix must loads
    # realistic fermentable sugar (~245 g/L, not 264) and the wine ABV is realistic
    # (decision D-16). Absent ⇒ 1.0 (no correction), so older parameter sets still
    # compile. Produced-only pools (X_dead, Gly, Byp, esters, fusels) default to 0
    # (see VarSpec) and so start empty at pitch.
    fermentable_fraction = (
        parameters["must_fermentable_fraction"].value
        if "must_fermentable_fraction" in parameters
        else 1.0
    )
    sugar_gpl = brix_to_sugar_gpl(_require(values, "brix", "wine")) * fermentable_fraction
    # pH-solver acid inputs (decision D-18), all optional so acid-free scenarios still
    # compile (slots default to 0, inert). tartaric/malic are must inputs in g/L;
    # lactic is produced-only (MLF product), 0 at pitch. The net strong cation is
    # back-solved from the measured initial_ph so the modelled pH reproduces it at t=0
    # (inverse anchoring): D-18 predicts pH *changes*, not absolute initial pH.
    tartaric = _optional(values, "tartaric_gpl", 0.0)
    malic = _optional(values, "malic_gpl", 0.0)
    initial: _Initial = {
        "X": _require(values, "pitch_gpl", "wine"),
        "S": [sugar_gpl],
        "E": _optional(values, "ethanol_gpl", 0.0),
        "N": mgl_to_gpl(_require(values, "yan_mgl", "wine")),
        "T": temperature_k,
        "CO2": 0.0,
        "X_dead": 0.0,  # no inactivated biomass at pitch
        "Gly": 0.0,  # no byproducts at pitch (decision D-16)
        "Byp": 0.0,
        # Produced-only aroma pools, empty at pitch (decision D-19). The three ester pools and
        # their headspace twins (D-96) are spread from the canonical registry, so a fourth ester
        # needs no edit here — and cannot be silently omitted from a pitch (pack() would raise).
        **{spec.pool: 0.0 for spec in ESTER_SPECS},
        # Volatilized-ester bookkeeping pools, empty at pitch (decisions D-20/D-96)
        **{spec.gas_pool: 0.0 for spec in ESTER_SPECS},
        # The five Ehrlich higher-alcohol pools, empty at pitch (decision D-99) — spread from
        # the canonical registry for the same reason as the esters. No gas twins: higher
        # alcohols are not stripped.
        **{spec.pool: 0.0 for spec in FUSEL_SPECS},
        "tartaric": tartaric,
        "malic": malic,
        "lactic": 0.0,
        "cation_charge": 0.0,  # back-solved below iff initial_ph is given
        # Total-SO₂ dose for the free/bound + molecular-SO₂ readout (D-22/D-28); mg/L→g/L,
        # default 0 (no dose). Inert/conserved state (readout-only, not in the charge
        # balance), so it does NOT enter the cation back-solve below — SO₂'s minor bisulfite
        # charge is a scoped omission the inverse anchoring would absorb at t=0 anyway (D-22).
        # Free/bound are derived from this total + acetaldehyde at the solved pH (D-28).
        "so2_total": mgl_to_gpl(_optional(values, "so2_total_mgl", 0.0)),
        # 5-Oxofructose — the dominant BOTRYTIS-specific SO2-binding carbonyl (decision D-130);
        # mg/L→g/L, default 0 (a non-botrytis must). Botrytis oxidises must fructose ON THE BERRY
        # (pre-crush), so it is a must-composition INPUT like the SO2 dose above, not a
        # fermentation product — carried as an inert slot no Process touches. Yeast-inert (Handbook
        # Vol 1 via Barbe 2000: "not altered by alcoholic fermentation"), so it persists into the
        # bottle and drives a botrytized wine's high SO2-combining power; typical botrytis load
        # ~80-150 mg/L. Off the charge balance (like so2_total), so it does NOT enter the cation
        # back-solve. At 0 the binding equilibrium is byte-for-byte the D-51 3-carbonyl form.
        "oxofructose": mgl_to_gpl(_optional(values, "oxofructose_mgl", 0.0)),
        # Medium-chain fatty acids (octanoic+decanoic, octanoic-equivalent) — the aggregate
        # yeast-secreted MLF inhibitor (decision D-131, Lonvaud-Funel 1988); mg/L→g/L, default 0.
        # A wine-composition-at-MLF INPUT, not a produced quantity (v1 defers the yeast-synthesis
        # production layer), so like so2_total/oxofructose it is a dosed inert slot no Process
        # touches. Read only by malolactic_environmental_gate (g_FA); at 0 the MLF is byte-for-byte
        # the pre-D-131 form. Off the charge balance (a weak acid, mostly protonated at wine pH;
        # its trace dissociation is a scoped omission like SO2's bisulfite charge, D-22).
        "mcfa": mgl_to_gpl(_optional(values, "mcfa_mgl", 0.0)),
        # DMS potential — the grape-borne precursor of aged-wine DMS (decision D-102), in
        # DMS-EQUIVALENT µg/L. Unlike every other optional above, this does NOT default to 0: it
        # defaults to the SOURCED must level (dms_potential_initial). The distinction is real —
        # so2_total/oak/anthocyanin are winemaking DOSES, and 0 is a true statement about a
        # scenario that made no addition, whereas DMSp is a property of the GRAPE that every must
        # carries. Defaulting it to 0 would assert that aged wine develops no DMS, which is the
        # D-45 hard-zero defect (a Process that silently never fires). Scenarios override via
        # `dms_potential_ugl` — and should, since DMSp is strongly variety-dependent and the
        # sourced default is Syrah's (see dms.yaml's notes on the Amarone over-prediction).
        # Absent from the ParameterSet ⇒ 0.0, so older parameter sets still compile inertly.
        "dms_potential": ugl_to_gpl(
            _optional(
                values,
                "dms_potential_ugl",
                (
                    parameters["dms_potential_initial"].value
                    if "dms_potential_initial" in parameters
                    else 0.0
                ),
            )
        ),
        "dms": 0.0,  # produced-only: no DMS at pitch, it accumulates over bottle aging
        # The metal-complexed sulfide reservoirs (decision D-135), µg/L. Like dms_potential above
        # these do NOT default to 0 but to their sourced levels: Franco-Luesma & Ferreira 2016
        # measured free AND total H2S/MeSH in 24 wines and found 94 %/62 % of a red's is already
        # bound at bottling, so 0 would assert that a sealed wine never turns reductive (the D-45
        # hard-zero defect — the very gap D-101 recorded as unmodellable). Seeded at PITCH rather
        # than at bottling as a v1 simplification: the reservoir is a property of the *finished*
        # wine, but no PROCESS reads or writes either slot until `begin_aging` enables the release
        # Processes, so carrying it inertly through fermentation is observationally identical to
        # seeding it at the aging boundary — and it is the oxofructose/mcfa "dosed inert slot"
        # idiom. D-193 NARROWS what follows: a `add_copper` fining now DOES build the reservoir
        # mid-run, from the stoichiometry that verb already computes, so "the model cannot build
        # bonded sulfide" is true only of the SPONTANEOUS route — bonded sulfide forming as
        # fermentative H2S meets must copper with no intervention, which still needs a binding
        # equilibrium and a binding constant nobody has published (see bound_sulfides.yaml).
        # Absent from the ParameterSet ⇒ 0.0, so older parameter sets still compile inertly.
        "bound_h2s": ugl_to_gpl(
            _optional(
                values,
                "bound_h2s_ugl",
                (
                    parameters["bound_h2s_initial"].value
                    if "bound_h2s_initial" in parameters
                    else 0.0
                ),
            )
        ),
        "bound_methanethiol": ugl_to_gpl(
            _optional(
                values,
                "bound_methanethiol_ugl",
                (
                    parameters["bound_methanethiol_initial"].value
                    if "bound_methanethiol_initial" in parameters
                    else 0.0
                ),
            )
        ),
        # Oenococcus oeni dose driving malolactic conversion (D-23); g/L, default 0 (no
        # MLF). Inert catalyst in v1 (no Process grows/kills it) and carbon-free, so an
        # undosed run is byte-for-byte the validated core; the compile step below disables
        # the MLF Processes entirely when this is 0 (tier + perf isolability).
        "X_mlf": _optional(values, "mlf_pitch_gpl", 0.0),
        # Citric acid must input (decision D-31); g/L, default 0 (no citrate ⇒ no MLF-derived
        # diacetyl). O. oeni co-metabolises it into α-acetolactate feeding the shared VDK
        # reservoir. Carbon-active (weighted in total_carbon) but not charge-active (kept out
        # of the D-18 pH balance in v1); inert at 0, so an un-dosed run is unchanged.
        "citrate": _optional(values, "citrate_gpl", 0.0),
        # Assimilable amino-acid dose (decisions D-32, SPECIATED at D-100); g/L, default 0 (no
        # amino-acid ledger). One ``amino_acids_gpl`` dose is split across the eight speciated
        # pools by the sourced must spectrum, with optional per-species override keys — see
        # :func:`_wine_amino_acids`. Carbon- AND nitrogen-bearing (every pool is weighted in both
        # ledgers); inert at 0 and the compile step below disables the amino-acid Processes
        # entirely when the dose is 0, so an undosed run is byte-for-byte the validated core
        # (tier + perf isolability, the MLF/carrying pattern).
        **_wine_amino_acids(values, parameters),
        # Non-assimilable cell-wall debris pool (decision D-34); produced-only, empty at pitch.
        # YeastAutolysis routes the carbon-rich remainder of autolysed dead biomass here after
        # releasing the nitrogen-rich amino acids; inert (weight 0) until autolysis is opted in.
        "debris": 0.0,
        # Lumped hydroxycinnamic-acid must precursors (decision D-40); g/L, default 0. Real must
        # carries ~10-200 mg/L (p-coumaric + ferulic); defaulted 0 for isolability (an undosed
        # run is byte-for-byte the validated core), so a Brett scenario doses it. Carbon-active
        # (weighted in total_carbon as p-coumaric); Brettanomyces decarboxylates it to vinylphenols.
        "hydroxycinnamics": _optional(values, "hydroxycinnamic_gpl", 0.0),
        "vinylphenols": 0.0,  # shared decarboxylase→reductase intermediate, empty at pitch (D-40)
        "ethylphenols": 0.0,  # terminal Brett volatile-phenol readout, empty at pitch (D-40)
        # Ferulic-acid branch (decision D-55): a genuinely distinct second precursor pool, split
        # out from hydroxycinnamics because ferulic acid (10 C) is a different molecule from
        # p-coumaric (9 C) whose decarboxylation cannot be a fixed-ratio split of the p-coumaric
        # flow without breaking carbon closure. Same isolability shape as the p-coumaric branch.
        "ferulic_acid": _optional(values, "ferulic_acid_gpl", 0.0),
        "vinylguaiacols": 0.0,  # ferulic-branch decarboxylase→reductase intermediate (D-55)
        "ethylguaiacols": 0.0,  # ferulic-branch terminal Brett volatile-phenol readout (D-55)
        # Brettanomyces dose driving the volatile-phenol spoilage (decision D-40); g/L, default 0
        # (no Brett). Constant inert catalyst in pt1 and carbon-free, so an undosed run is
        # byte-for-byte the validated core; the compile step below disables the Brett Processes
        # entirely when this is 0 (tier + perf isolability, the mlf_pitch_gpl pattern).
        "X_brett": _optional(values, "brett_pitch_gpl", 0.0),
        "X_brett_dead": 0.0,  # non-viable Brett lees, empty until BrettDeath (D-40 pt3)
        # Grape anthocyanin + condensed tannin must inputs (decision D-79); g/L, default 0 (a white
        # /
        # no-tannin wine). TanninAnthocyaninCondensation combines them into stable polymeric pigment
        # during aging (colour stabilization + astringency softening). Off every ledger (grape-
        # derived, the iso_alpha/ellagitannin precedent), so an undosed run is byte-for-byte the
        # validated core; the Process is doubly substrate-gated AND disabled at compile until
        # begin_aging (aging is post-ferment), so present-but-un-aged wine carries them inertly too.
        "anthocyanin": _optional(values, "anthocyanin_gpl", 0.0),
        "tannin": _optional(values, "tannin_gpl", 0.0),
        # Initial-burst antioxidant pool (decision D-133), g/L. Like dms_potential, this is NOT a
        # winemaking dose but a GRAPE-composition property every must carries — so, unlike
        # so2_total/tannin/anthocyanin above, it does NOT default to 0: it falls back to the
        # sourced burst_antioxidant_initial. Defaulting to 0 would assert every wine's Ferreira-
        # measured day-1 O2-burst is zero (the D-45 hard-zero defect). Scenarios override via
        # `burst_antioxidant_gpl` — and should, since Ferreira found >15x between-wine spread
        # (0.54-8.2 mg/L/day initial rate, Cu-driven, untracked here). Absent from the
        # ParameterSet ⇒ 0.0, so older parameter sets still compile inertly.
        #
        # D-147: the fallback above is CONDITIONAL and this builder cannot see the condition — it
        # does not know which oxidative set was selected. `_resolve_burst_antioxidant_seed` in
        # `compile_scenario` applies that condition after the pack, the `iso_alpha`/hops precedent
        # below it. Read the two together: the D-45 fallback is right only where the consumer is.
        "burst_antioxidant": _optional(
            values,
            "burst_antioxidant_gpl",
            (
                parameters["burst_antioxidant_initial"].value
                if "burst_antioxidant_initial" in parameters
                else 0.0
            ),
        ),
        # Dissolved copper (decision D-134), g/L. Like burst_antioxidant/dms_potential, this is NOT
        # a winemaking dose but a must-composition property every wine carries — so it does NOT
        # default to 0: it falls back to the sourced copper_typical, the SAME level D-132's
        # k_browning_phenolic is implicitly calibrated at (so an un-overridden wine reproduces the
        # D-132/D-133 rate byte-for-byte via PhenolicBrowning's mean-centered f(Cu) = 1). Scenarios
        # override via `copper_gpl` to explore atypical-copper wines. Absent from the ParameterSet
        # ⇒ 0.0 (an un-seeded/older ParameterSet compiles, but is NOT byte-for-byte D-132 — see
        # the `copper` VarSpec docstring).
        "copper": _optional(
            values,
            "copper_gpl",
            parameters["copper_typical"].value if "copper_typical" in parameters else 0.0,
        ),
    }
    if "initial_ph" in values:
        # Byp = 0 at pitch, so the anchoring cation reproduces initial_ph from the named
        # acids alone; as Byp accumulates during the ferment, pH drifts emergently.
        # `phosphate` (D-210) is deliberately absent: nothing seeds it, so it is 0 at every
        # anchor, and an absent key contributes nothing to `charge_residual`'s sum — which keeps
        # this seam bitwise what it was. It is also the whole reason the dosed phosphate MOVES
        # anything: a species present here is absorbed by the fitted cation (D-178's result), and
        # one dosed later is not. If a native must phosphate is ever seeded, it belongs here and
        # will correctly become a near no-op.
        acid_gpl = {"tartaric": tartaric, "malic": malic, "lactic": 0.0}
        totals_molar = {n: g / acidbase.ACID_STATE[n].molar_mass for n, g in acid_gpl.items()}
        try:
            # The assimilable-nitrogen pool carries net POSITIVE charge and is present at the
            # anchor, so the cation the target pH needs is partly supplied by ``N`` (decision
            # D-209). Subtract that share, or the slot would double-count it and the must would
            # not read back at ``initial_ph``. Must nitrogen is arginine-dominated, hence a
            # wine-specific z-bar rather than beer's.
            total_cation = acidbase.solve_cation_charge(
                totals_molar,
                byp_succinic_molar=0.0,
                # Dissolved CO2 is 0 at the anchor and that is STRUCTURAL, not a shortcut
                # (decision D-182): a must has not fermented, so the ``CO2`` slot this term
                # is derived from is 0. It is also exactly why the term moves the finished pH
                # instead of vanishing — a species present at the anchor gets absorbed into
                # the fitted cation and becomes a near no-op (D-178's phosphate result).
                carbonic_molar=0.0,
                pka_map=acidbase.build_pka_map(parameters.resolve()),
                target_ph=float(values["initial_ph"]),
            )
            initial["cation_charge"] = acidbase.cation_slot_after_nitrogen(
                total_cation,
                acidbase.nitrogen_charge_from_gpl(
                    mgl_to_gpl(_require(values, "yan_mgl", "wine")), "wine", parameters.resolve()
                ),
                float(values["initial_ph"]),
                "wine",
            )
        except ValueError as exc:  # initial_ph below the acid load's intrinsic pH
            raise ValueError(f"wine scenario.initial['initial_ph'] is unphysical: {exc}") from exc
        except KeyError as exc:  # acidbase.yaml pKa parameters not loaded
            raise ValueError(
                "wine scenario gives 'initial_ph' but the pKa parameters are missing "
                f"({exc}); include acidbase.yaml in parameter_paths (the default lookup "
                "merges it automatically)."
            ) from exc
    return initial


#: Beer's charge-active acid slots and the sourced level each is seeded from when the pH
#: system is opted into (decision D-179). Slot → (scenario override key, parameter name).
#: Ordered so the back-solve below is deterministic. ``peptide_buffer`` rides in the same
#: table because it enters the charge balance identically — it is simply not an organic acid.
#:
#: **These are WORT levels since D-180, and that change was forced, not chosen.** Until then
#: they were the ``*_typical_beer`` levels — a FINISHED BEER's acid composition dosed into a
#: wort at pitch, which was harmless only for as long as nothing produced acid. D-180 gives
#: four of these slots a producer, and a producer on top of finished-beer seeds finishes at
#: roughly twice the measured beer (verified: pH 4.26, below any real beer). So "give beer an
#: acid producer" and "start beer from a wort" are one decision.
#:
#: ``initial_ph``'s own contract is UNCHANGED — it still means "reproduce this pH at t=0". What
#: changed is which pH a caller should hand it: a beer scenario compiles at PITCH, so the
#: physical input is now the wort pH (~5.4-5.7), not the finished beer's, and the finished
#: beer's pH becomes an OUTPUT. A scenario that still passes 4.4 is not an error and is not
#: rejected — it describes a wort that is already at beer pH, and it will finish lower still.
_BEER_ACID_SEEDS: tuple[tuple[str, str, str], ...] = (
    ("lactic", "lactic_gpl", "lactic_typical_wort"),
    ("acetic", "acetic_gpl", "acetic_typical_wort"),
    ("citrate", "citric_gpl", "citric_typical_wort"),
    ("malic", "malic_gpl", "malic_typical_wort"),
    ("succinic", "succinic_gpl", "succinic_typical_wort"),
    # The three that FALL (decision D-181). Seeded from the SAME wort as the five above —
    # Tyrell's Figs 6/7/11 against their 9/10/12/13/14 — and drained by ``WortAcidRemoval``.
    # They enter the cation back-solve below like any other acid, which is what makes the
    # anchored start pH unchanged and the FINISHED pH higher: more acid at t=0 means a larger
    # back-solved cation to still hit ``initial_ph``, and that cation stays when the acid goes.
    ("pyruvic", "pyruvic_gpl", "pyruvic_typical_wort"),
    ("formic", "formic_gpl", "formic_typical_wort"),
    ("oxalic", "oxalic_gpl", "oxalic_typical_wort"),
    ("peptide_buffer", "peptide_buffer_gpl", "peptide_buffer_capacity_beer"),
    # `phosphate` (D-210) is NOT here, and its absence is the reason it does anything. Malt
    # phosphate really is in a wort — D-178 measured it and refused it as beer's buffer, because
    # a near-constant charge PRESENT AT THE ANCHOR is absorbed by the back-solved cation. Only
    # the DOSED phosphate `add_dap` writes is post-anchor, so only that one acidifies. Seeding a
    # wort phosphate here would reproduce D-178's near-no-op, which is why it is not a gap.
)


def _beer_acids(values: Mapping[str, float], parameters: ParameterSet) -> dict[str, float]:
    """Beer's acid doses (decision D-179) — sourced levels, individually overridable.

    **``initial_ph`` is the gate, and absent means zero.** This is the OPPOSITE call to
    ``dms_potential`` / ``copper`` / ``burst_antioxidant``, where absent falls back to a
    sourced level because the quantity is a property every must carries (the D-45 hard-zero
    argument). The distinction is D-147's: a fallback is only right where the CONSUMER is
    wired, and here the consumer is a charge balance that is meaningless without its
    counter-cation. Seeding acids into a beer with no ``initial_ph`` would give an acid load
    with no strong cation — the very configuration D-178 measured at pH 4.47 and called
    "plausible-looking, produced by nothing real". So the whole system opts in together.

    With ``initial_ph`` given, each slot takes its sourced level (Tyrell 2013 for the organic
    acids, the back-solve for the peptide buffer) unless the scenario overrides it. An
    override REPLACES that slot outright — including an explicit ``0.0``, which is how the
    organic-acids-only arm is run.
    """
    if "initial_ph" not in values:
        return {slot: 0.0 for slot, _, _ in _BEER_ACID_SEEDS}
    out: dict[str, float] = {}
    for slot, key, param in _BEER_ACID_SEEDS:
        default = parameters[param].value if param in parameters else 0.0
        out[slot] = _optional(values, key, default)
    return out


def _beer_cation(
    values: Mapping[str, float], acids: Mapping[str, float], parameters: ParameterSet
) -> float:
    """Back-solve beer's net strong cation from its measured ``initial_ph`` (decision D-179).

    The same inverse anchoring wine has used since D-18, now that beer has an acid load to
    anchor against: the claim is that the model predicts pH *changes*, not absolute initial
    pH, and the absolute is an input. ``Byp`` is 0 at pitch (beer has no producer for it), so
    the anchor is set by the named acids alone.

    Returns 0.0 when no ``initial_ph`` is given — with every acid slot also 0 (see
    :func:`_beer_acids`) that is an empty charge balance, byte-for-byte the pre-D-179 beer.
    """
    if "initial_ph" not in values:
        return 0.0
    totals_molar = {slot: gpl / acidbase.BEER_ACIDS[slot].molar_mass for slot, gpl in acids.items()}
    try:
        # Wort nitrogen is on the cation side and present at the anchor, so its charge is
        # subtracted off the fitted slot (decision D-209) — the same correction wine's anchor
        # makes, with wort's own z-bar. Wort's amino-acid halves nearly cancel, so this is
        # essentially the wort ammonium; the term is what lets beer's pH FALL as YAN is taken up.
        total = acidbase.solve_cation_charge(
            totals_molar,
            byp_succinic_molar=0.0,
            # 0 for the same structural reason as wine's anchor above (decision D-182): a
            # wort has not fermented, so its evolved-CO2 slot is 0 and the carbonic term is
            # absent from the anchor — which is what leaves it free to move the finished pH.
            carbonic_molar=0.0,
            pka_map=acidbase.build_pka_map(parameters.resolve()),
            target_ph=float(values["initial_ph"]),
        )
        return acidbase.cation_slot_after_nitrogen(
            total,
            acidbase.nitrogen_charge_from_gpl(
                mgl_to_gpl(_require(values, "yan_mgl", "beer")), "beer", parameters.resolve()
            ),
            float(values["initial_ph"]),
            "beer",
        )
    except ValueError as exc:  # initial_ph below the acid load's intrinsic pH
        raise ValueError(f"beer scenario.initial['initial_ph'] is unphysical: {exc}") from exc
    except KeyError as exc:  # acidbase.yaml / beer_acids.yaml pKa parameters not loaded
        raise ValueError(
            "beer scenario gives 'initial_ph' but the pKa parameters are missing "
            f"({exc}); include acidbase.yaml in parameter_paths (the default lookup "
            "merges it automatically)."
        ) from exc


def _beer_initial(
    values: Mapping[str, float], temperature_k: float, parameters: ParameterSet
) -> _Initial:
    acids = _beer_acids(values, parameters)
    return {
        "X": _require(values, "pitch_gpl", "beer"),
        # Beer's charge-active acids + the inverse-anchored strong cation (decision D-179).
        #
        # THIS COMMENT WAS STALE AND IS CORRECTED AT D-213. It used to say every slot here is
        # inert "because beer still has no organic-acid producer (D-16, open)", and that beer's
        # pH "does not FALL during fermentation the way a real beer's does". Both statements
        # were overtaken and neither is true now: `OrganicAcidExcretion` and `WortAcidRemoval`
        # (D-181) and `AceticAcidOverflow` (D-183) all move these slots, and beer's falling pH
        # is a PREDICTION that agrees with Tyrell's measured course (D-207, re-framed to the
        # decarbonated frame at D-208, and closed on day-1 timing at D-211). Found while adding
        # the wort-oxygen seed below — the same class as the two prose numbers D-212 repaired.
        **acids,
        "cation_charge": _beer_cation(values, acids, parameters),
        "S": [
            _require(values, "glucose_gpl", "beer"),
            _require(values, "maltose_gpl", "beer"),
            _require(values, "maltotriose_gpl", "beer"),
        ],
        "E": _optional(values, "ethanol_gpl", 0.0),
        "N": mgl_to_gpl(_require(values, "yan_mgl", "beer")),
        "T": temperature_k,
        "CO2": 0.0,
        "X_dead": 0.0,  # no inactivated biomass at pitch
        # Wort aeration (decision D-213): a cast wort is aerated on its way to the fermenter, so
        # beer's `o2` starts at the sourced level rather than the slot default of 0. Stripped by
        # `WortOxygenUptake` during the lag phase. `o2_mgl` overrides it; 0 recovers the
        # pre-D-213 anaerobic-from-t0 beer. NOTE it feeds nothing in the default set — the three
        # O2 consumers are aging-gated — so this moves the `o2` column and no other.
        #
        # `if in parameters else 0.0` follows `_beer_acids`, and the fallback is deliberate: a
        # parameter set that omits `beer_generic.yaml` is a deliberately reduced core (several
        # tests compile against exactly that), and 0.0 is the pre-D-213 behaviour — i.e. the
        # isolable default this Process is required to reduce to. Raising instead would make a
        # reduced core uncompilable for the sake of a term that moves nothing.
        "o2": mgl_to_gpl(
            _optional(
                values,
                "o2_mgl",
                parameters["o2_wort_aeration_beer"].value
                if "o2_wort_aeration_beer" in parameters
                else 0.0,
            )
        ),
        "Gly": 0.0,  # beer carries zero byproduct diversion in M1 (decision D-16)
        "Byp": 0.0,
        # Produced-only aroma pools, empty at pitch (decision D-19). The three ester pools and
        # their headspace twins (D-96) are spread from the canonical registry, so a fourth ester
        # needs no edit here — and cannot be silently omitted from a pitch (pack() would raise).
        **{spec.pool: 0.0 for spec in ESTER_SPECS},
        # Volatilized-ester bookkeeping pools, empty at pitch (decisions D-20/D-96)
        **{spec.gas_pool: 0.0 for spec in ESTER_SPECS},
        # The five Ehrlich higher-alcohol pools, empty at pitch (decision D-99) — spread from
        # the canonical registry for the same reason as the esters. No gas twins: higher
        # alcohols are not stripped.
        **{spec.pool: 0.0 for spec in FUSEL_SPECS},
    }


_INITIAL_BUILDERS: dict[str, Callable[[Mapping[str, float], float, ParameterSet], _Initial]] = {
    "wine": _wine_initial,
    "beer": _beer_initial,
}


#: The closure menu (decision D-136), in ascending order of steady oxygen transmission. What Lopes
#: et al. 2007 actually establishes is the THREE-TIER grouping — "low in screw-caps and 'technical'
#: corks, intermediate in conventional natural cork stoppers, and high in synthetic closures" — and
#: that is the defensible claim; the two WITHIN-TIER orderings this tuple additionally asserts
#: (technical cork vs screwcap, Nomacorc vs SupremeCorq) come off Table I's nominals, not off that
#: sentence, and they are not equally strong. D-162 measured each against the shipped bands; the
#: scoping lives in ``closure.yaml``'s header. Each name maps to ``otr_<name>`` in ``closure.yaml``,
#: so the value carries full provenance instead of being an inlined constant (prime directive #2).
#:
#: NOTE the order: **technical cork sits BELOW screwcap**, which contradicts the widely repeated
#: "screwcaps are the least permeable closure" (Godden et al. 2005, quoted in Oliveira et al. 2013's
#: own introduction). Both are true on their own terms — the screwcap's famous number is `<500
#: µL/day` *at the moment of bottling*, headspace air trapped at sealing rather than transmission
#: through the liner, which dominates any comparison made on a total-including-burst basis. At
#: steady state, which is what this axis models, technical cork wins. Do not "fix" this ordering.
_CLOSURES: tuple[str, ...] = (
    "hermetic",
    "technical_cork",
    "screwcap",
    "natural_cork",
    "synthetic_nomacorc",
    "synthetic_supremecorq",
)


def _closure_otr(closure: str, parameters: ParameterSet) -> float:
    """The named closure's steady oxygen transmission rate [g/L/h] (decision D-136).

    Resolves ``scenario.closure`` against :data:`_CLOSURES` and returns the sourced ``otr_<name>``
    value from ``closure.yaml``. Both failure modes are loud and name the menu, in the
    ``_ALLOWED_KEYS`` spirit: an unknown closure name is a scenario error (not a silent fallback to
    no ingress, which would quietly age the wine wrongly), and a missing parameter file is a
    configuration error surfaced HERE at compile rather than as a bare ``KeyError`` mid-run.
    """
    if closure not in _CLOSURES:
        raise ValueError(
            f"unknown scenario.closure {closure!r}; expected one of {', '.join(_CLOSURES)} "
            "(decision D-136)"
        )
    name = f"otr_{closure}"
    if name not in parameters:
        raise ValueError(
            f"scenario.closure {closure!r} needs parameter {name!r}, which is not loaded; "
            "closure.yaml is missing from the parameter set (decision D-136)"
        )
    return float(parameters[name].value)


def _bottling_burst(closure: str, parameters: ParameterSet) -> float:
    """The named closure's one-off bottling charge [g/L] (decision D-187).

    :func:`_closure_otr`'s sibling, resolving ``bottling_burst_<name>`` instead of ``otr_<name>``,
    with the same two loud failure modes for the same reasons. The two quantities are deliberately
    separate parameters rather than one entry with a transient: the OTR is *permeation through* the
    closure and this is the closure's *own trapped air* decompressing out of it — different
    mechanisms, measured in different columns of Lopes et al. 2007's Table I, and only the second
    is a finite stock.
    """
    if closure not in _CLOSURES:
        raise ValueError(
            f"unknown scenario.closure {closure!r}; expected one of {', '.join(_CLOSURES)} "
            "(decision D-136)"
        )
    name = f"bottling_burst_{closure}"
    if name not in parameters:
        raise ValueError(
            f"scenario.closure {closure!r} needs parameter {name!r}, which is not loaded; "
            "closure.yaml is missing from the parameter set (decision D-187)"
        )
    return float(parameters[name].value)


def _iso_alpha_at_pitch(scenario: Scenario, parameters: ParameterSet) -> float:
    """Iso-alpha-acids [g/L] delivered to the fermenter from the boil (decision D-64).

    For each hop addition, runs the Malowicki closed-form isomerization
    (:func:`~fermentation.core.kinetics.hops.iso_alpha_fraction`) at the scenario's boil
    temperature, weights it by that addition's *dissolved* alpha concentration (hop mass /
    ``batch_volume_liters``; full dissolution assumed — extraction incompleteness is folded into
    ``hop_utilization_efficiency``), sums the additions, and applies the kettle->fermenter
    utilization efficiency. Evaluated once here at the compile boundary because the boil is a
    wort-side input (373 K, no yeast), not a fermentation phase — running it through the
    integrator would drive the yeast-free wort at boiling temperature. The result seeds the
    ``iso_alpha`` state; :class:`~fermentation.core.kinetics.hops.IsoAlphaAcidLoss` then reduces
    it during fermentation. ``batch_volume_liters`` is guaranteed present by the scenario
    validator whenever ``hops`` is non-empty.
    """
    volume = scenario.batch_volume_liters
    if volume is None:  # defensive; the schema validator already enforces this
        raise ValueError("hop bittering needs 'batch_volume_liters' (decision D-64)")
    boil_temp_k = celsius_to_kelvin(scenario.boil_celsius)
    resolved = parameters.resolve()
    total_iso_gpl = 0.0
    for hop in scenario.hops:
        alpha0_gpl = (hop.alpha_acid_percent / 100.0) * hop.grams / volume
        total_iso_gpl += alpha0_gpl * iso_alpha_fraction(hop.boil_minutes, boil_temp_k, resolved)
    return total_iso_gpl * resolved["hop_utilization_efficiency"]


def _resolve_burst_antioxidant_seed(
    scenario: Scenario, medium: Medium, process_set: ProcessSet, y0: FloatArray
) -> None:
    """Make the ``burst_antioxidant`` seed follow its consumer (decision D-147).

    ``_wine_initial`` seeds the slot from the sourced ``burst_antioxidant_initial`` on the D-45
    "absent does not mean 0" reasoning D-133 borrowed from ``dms_potential``: *a 0 default would
    silently assert that every wine's Ferreira-measured day-1 O2-burst is absent*. That reasoning
    is sound **only where something can draw the pool**, and D-133 shipped
    :class:`~fermentation.core.kinetics.aging.AntioxidantBurstOxidation` wired into no medium at
    all. D-140 found the discrepancy; D-147 measured what it costs.

    **Without the consumer the D-45 argument inverts.** With nothing able to draw it, the model
    already asserts the burst never happens — the missing Process says that, not the seed — and the
    non-zero seed *additionally* asserts an antioxidant that is present and is never spent, for the
    whole life of the wine. That is strictly worse than 0, and unlike a quiet modelling choice it
    is **visible in output**: a 2 y run of the default build emits ``burst_antioxidant`` constant at
    3.3e-3 g/L with ``ptp == 0.0`` exactly.

    So the seed is conditional on the wiring, and dosing the pool into a build that cannot consume
    it **raises** rather than silently seeding an inert number — the ``hops``-without-a-bitterness-
    model precedent one block up in :func:`compile_scenario`. A build that *does* wire the Process
    (``oxidative="direct_burst"``) keeps D-133's fallback unchanged, because there the argument for
    it is the one D-133 actually made.
    """
    if "burst_antioxidant" not in medium.schema:
        return  # beer: no slot, nothing to resolve (Ferreira's dataset is red wine only)
    if AntioxidantBurstOxidation.name in process_set:
        return  # the consumer is wired ⇒ D-133's sourced fallback stands, as written
    if "burst_antioxidant_gpl" in scenario.initial:
        raise ValueError(
            "scenario dosed 'burst_antioxidant_gpl' but the compiled oxidative set wires no "
            f"{AntioxidantBurstOxidation.name!r} Process, so nothing can draw the pool and the "
            "dose would sit in the output unspent (decision D-147). Compile with "
            "oxidative='direct_burst' to wire the consumer, or drop the dose."
        )
    y0[medium.schema.slice("burst_antioxidant")] = 0.0


def _validate_initial_keys(scenario: Scenario) -> None:
    allowed = _ALLOWED_KEYS.get(scenario.medium)
    if allowed is None:
        raise ValueError(
            f"medium {scenario.medium!r} has no initial-composition vocabulary defined"
        )
    unknown = set(scenario.initial) - allowed
    if unknown:
        raise ValueError(
            f"scenario.initial has unknown key(s) {sorted(unknown)} for medium "
            f"{scenario.medium!r}; allowed: {sorted(allowed)}"
        )


def _initial_temperature_kelvin(scenario: Scenario) -> float:
    schedule = scenario.temperature_schedule
    if not schedule:
        raise ValueError(
            f"scenario {scenario.name!r}: temperature_schedule needs at least one point "
            "to seed the initial temperature"
        )
    earliest = min(schedule, key=lambda point: point.day)
    return celsius_to_kelvin(earliest.celsius)


def _load_parameters(
    scenario: Scenario,
    parameter_paths: Sequence[str | Path] | None,
    data_dir: str | Path | None,
) -> ParameterSet:
    if parameter_paths is not None:
        # Caller-controlled override: a caller wanting the pH solver must include
        # acidbase.yaml in their paths (the pKa set the charge balance reads, D-18).
        return load_parameters(*parameter_paths)
    base = Path(data_dir) if data_dir is not None else default_data_dir()
    path = base / f"{scenario.medium}_{scenario.strain}.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"no parameter file for medium={scenario.medium!r} strain={scenario.strain!r}: "
            f"expected {path}. Pass parameter_paths=... or add the YAML "
            "(see the Milestone 1 parameter-sourcing task)."
        )
    # Merge the shared, medium-agnostic parameter files alongside the medium file so every
    # default-lookup scenario can compute pH (acidbase.yaml, decision D-18), run the diacetyl
    # pathway (vicinal_diketones.yaml, decision D-26 — the load-bearing decarb step is
    # non-enzymatic, so its constants are medium-agnostic), the acetaldehyde buffer
    # (acetaldehyde.yaml, decision D-27 — main-pathway yeast metabolism, likewise generic) and
    # H₂S production (hydrogen_sulfide.yaml, decision D-29 — the sulfate-reduction sequence,
    # generic yeast metabolism). The names are collision-free with the per-medium kinetic
    # parameters; load_parameters merges left-to-right.
    # additions.yaml carries the industry-unit → canonical conversion constants the
    # discrete-intervention verb registry reads at this boundary (decision D-36); like the
    # others it is medium-agnostic and collision-free.
    shared_files = [
        base / "acidbase.yaml",
        base / "vicinal_diketones.yaml",
        base / "acetaldehyde.yaml",
        base / "keto_acids.yaml",
        base / "hydrogen_sulfide.yaml",
        base / "additions.yaml",
        # Hop bittering kinetics (decision D-64): the Malowicki boil isomerization constants and
        # the iso-alpha loss/utilization parameters. Beer-only in effect (only beer carries an
        # iso_alpha slot and the hop Process/boil calc), but loaded universally like the other
        # shared files — collision-free names, inert for wine.
        base / "hops.yaml",
        # Aging chemistry (decision D-70): the ester-hydrolysis constants (k_ester_hydrolysis,
        # E_a_ester_hydrolysis, isoamyl_acetate_eq) the post-fermentation EsterHydrolysis Process
        # reads — and, since D-126, the ethyl-hexanoate-hydrolysis constants
        # (k_ethyl_hexanoate_hydrolysis, E_a_ethyl_hexanoate_hydrolysis, ethyl_hexanoate_eq) its
        # sibling EthylHexanoateHydrolysis reads.
        # Medium-agnostic (acid-catalysed hydrolysis is a molecule/pH property, the
        # vicinal_diketones.yaml pattern) and collision-free, so loaded universally like the other
        # shared files; INERT until a begin_aging intervention enables the Process (which is
        # disabled at compile), so an un-aged scenario carries the params but never reads them.
        base / "aging.yaml",
        # Oak extraction (decision D-77): the barrel/chip aroma-extractive constants —
        # k_oak_extraction, the weak diffusion E_a_oak_extraction, and the 12 toast-specific
        # yields the add_oak verb reads to set each ceiling. Wine-only in effect (only wine carries
        # the oak slots + wires OakExtraction), but loaded universally like the other shared files —
        # collision-free names, inert for beer; INERT until an add_oak dose + begin_aging enable.
        base / "oak.yaml",
        # Tannin–anthocyanin condensation (decision D-79): the red-wine colour-stabilization +
        # astringency-softening rate/E_a/yield the TanninAnthocyaninCondensation Process reads.
        # Wine-only in effect (only wine carries the anthocyanin/tannin slots + wires the Process),
        # but loaded universally like the other shared files — collision-free names, inert for beer;
        # INERT until an anthocyanin+tannin must dose + begin_aging enable.
        base / "polymerization.yaml",
        # Non-oxidative THERMAL aging axis (decisions D-87/D-88/D-89): the sugar+heat-driven,
        # O2-independent Strecker aldehydes (MaillardStrecker), sugar-only caramelization browning
        # (Caramelization), and amino-acid-incorporating Maillard browning (MaillardBrowning — the
        # N-bearing melanoidin branch). Caramelization is MEDIUM-AGNOSTIC (D-90: both media wire it
        # + carry the melanoidin slot — beer's residual dextrins caramelize); MaillardStrecker /
        # MaillardBrowning stay wine-only in effect (only wine wires them + carries amino_acids /
        # the thermal aroma slots). Loaded universally like the other shared files — collision-free
        # names; every rate is INERT until a begin_aging enable.
        base / "thermal.yaml",
        # DMS via SMM hydrolysis (decision D-102): the k/E_a of the grape-borne precursor's
        # bottle-aging hydrolysis + the sourced must DMS-potential level that seeds the pool.
        # Wine-only in effect (only wine carries the dms_potential/dms slots + wires
        # SMMHydrolysis), but loaded universally like the other shared files — collision-free
        # names, inert for beer. Beer's DMS is real but arrives by OTHER routes this file's
        # wine-anchored constants must not be used for (see dms.yaml / the Process docstring);
        # INERT until a begin_aging enable.
        base / "dms.yaml",
        # Bottle-reduction sulfides (decision D-135): the release rates of the metal-complexed
        # H2S/MeSH reservoirs Franco-Luesma & Ferreira 2016 measured, the sourced at-bottling
        # reservoir levels that seed them, and the mass-balance shares recording how much of real
        # sulfide aging the release-only model accounts for. Wine-only in effect (only wine carries
        # the bound_h2s/bound_methanethiol slots + wires the two Processes), but loaded universally
        # like the other shared files — collision-free names, inert for beer. INERT until a
        # begin_aging enable, and doubly so on an unseeded reservoir.
        base / "bound_sulfides.yaml",
        # Closure oxygen ingress (decision D-136): the steady per-closure OTRs Lopes et al.
        # 2007 measured, converted to g/L/h for a 750 mL bottle. Loaded universally like the
        # other shared files (collision-free names, inert for beer) but WINE-ONLY in effect —
        # only wine carries the closure_otr slot and wires ClosureOxygenIngress. UNLIKE every
        # other file here these are read at COMPILE time (to seed the state slot) rather than
        # per-RHS-step, so a scenario naming a `closure` without this file fails loudly in
        # _closure_otr. INERT until a begin_aging enable, and at closure=hermetic/absent.
        base / "closure.yaml",
        # Beer's acid COMPOSITION (decision D-179): the sourced Tyrell 2013 levels the beer
        # charge balance is dosed with, plus the back-solved lumped peptide-buffer capacity.
        # Loaded universally like the other shared files (collision-free names) but BEER-ONLY
        # in effect — only beer's schema carries these slots and only beer's acid registry
        # reads them, so every value is inert for wine. Like closure.yaml these are read at
        # COMPILE time (to seed the state slots), not per-RHS-step.
        base / "beer_acids.yaml",
    ]
    return load_parameters(path, *(f for f in shared_files if f.exists()))


def _apply_nitrogen_dependent_yield(scenario: Scenario, parameters: ParameterSet) -> ParameterSet:
    """Override ``biomass_N_fraction`` from Coleman's ``Y_X/N(N_init)`` regression.

    Coleman, Fish & Block (2007) found the cell-mass-per-nitrogen yield to depend
    on the *initial* nitrogen (Fig 4 / Table A2): ``ln(Y_X/N) = a0 + a1·YAN``
    (YAN in mg N/L). This is the one parameter that cannot be pre-evaluated into
    the YAML the way the temperature regressions are — the evaluation point is the
    scenario's nitrogen, not a fixed reference — so it is computed here at the
    compile boundary and nowhere else (decision D-14). Because every assimilated
    gram of nitrogen enters biomass in our model, ``Y_X/N = 1/f_N`` identically;
    setting ``biomass_N_fraction = 1/Y_X/N`` leaves the nitrogen balance exact (the
    ``total_nitrogen`` check reads this same per-run constant).

    Gated on the regression coefficients being present, so a medium without them
    (beer) keeps the static elemental ``biomass_N_fraction`` untouched.
    """
    if not all(name in parameters for name in _N_YIELD_COEFFS):
        return parameters
    yan_mgl = scenario.initial.get("yan_mgl")
    if yan_mgl is None:
        return parameters

    a0, a1 = (parameters[name] for name in _N_YIELD_COEFFS)
    y_xn = math.exp(a0.value + a1.value * float(yan_mgl))  # g cell / g N
    f_n = 1.0 / y_xn
    override = Parameter(
        name="biomass_N_fraction",
        value=f_n,
        unit="g/g",
        tier=combine((a0.tier, a1.tier)),
        uncertainty=Uncertainty(
            # Bracketing metadata, not a tuned value: f_N = 1/Y_X/N ranges
            # ~0.039-0.107 across Coleman's 70-350 mg N/L treatment span
            # (Y_X/N ~25.7 down to ~9.4); [0.03, 0.15] brackets that with margin.
            low=0.03,
            high=0.15,
            note="nitrogen-status-dependent; brackets f_N across Coleman's 70-350 mg N/L range",
        ),
        provenance=Provenance(
            source=a0.provenance.source,
            doi=a0.provenance.doi,
            conditions=(
                f"computed at compile from Coleman Y_X/N regression at YAN={float(yan_mgl):g} mg/L"
            ),
            notes=(
                f"Y_X/N = exp({a0.value} + {a1.value}*{float(yan_mgl):g}) = {y_xn:.2f} g cell/g N; "
                f"f_N = 1/Y_X/N = {f_n:.4f} g N/g cell. Overrides the static elemental "
                "biomass_N_fraction so a nitrogen-limited must builds realistically little "
                "biomass (decision D-14)."
            ),
        ),
    )
    return parameters.merge(ParameterSet([override]), override=True)


def _override_in_band(base: Parameter, value: float, knob: str) -> float:
    """Forbid a scenario override that would land outside the reference parameter's band.

    Decision D-164. The bound is real and load-bearing, but it was previously enforced
    only as a SIDE EFFECT: an override is built as a new :class:`Parameter` carrying the
    base's ``uncertainty``, so ``Parameter._value_in_range`` rejected it — with a message
    naming the *parameter's epistemic band*, never the *scenario knob the user actually
    set*. This states the same bound up front, named for what it forbids.

    Why the band is the admissible range **today**: the ensemble sampler draws this
    parameter as ``triangular(low, value, high)`` (``runtime.ensemble.sample_parameters``,
    and ``_inverse_cdf`` via ``scipy.stats.triang`` with ``c = (val-lo)/(hi-lo)``), which
    *requires* ``low <= value <= high``. Both knobs' parameters are confirmed present in
    the sampled-name set of the very run that enables them, so an out-of-band override
    would break a live call site, not merely an abstract invariant.

    This deliberately does NOT decouple the epistemic band from the admissible range —
    that needs its own provenance-bearing field, and D-164 flags it. Note the asymmetry it
    leaves: to sweep further you must widen the band in the YAML *with provenance*, which
    is the intended cost.
    """
    u = base.uncertainty
    if not (u.low <= value <= u.high):
        raise ValueError(
            f"scenario knob {knob!r}={value:g} is outside the admissible override range "
            f"[{u.low:g}, {u.high:g}] {base.unit} for parameter {base.name!r} "
            f"(decision D-164). A scenario override must stay inside the reference "
            f"parameter's stated band: the ensemble sampler draws it as "
            f"triangular(low, value, high), which requires low <= value <= high. "
            f"To sweep beyond this, widen {base.name!r}'s uncertainty band in the "
            f"parameter YAML with provenance — do not widen it here."
        )
    return value


def _override_carrying_capacity(parameters: ParameterSet, cap_gpl: float) -> ParameterSet:
    """Override the reference ``biomass_carrying_capacity`` with a scenario opt-in value.

    Only reached when a wine scenario passes ``carrying_capacity_gpl > 0`` (decision D-30),
    so the biomass cap modifier is enabled and its cap ``K`` is the scenario's value rather
    than the YAML reference — letting a demonstration sweep the cap **within the reference
    band** (:func:`_override_in_band`, decision D-164 — the sweep is bounded, and the YAML
    band is what bounds it). Keeps the reference parameter's units/tier/uncertainty (the
    form and confidence are unchanged; only the operating point moves) and records the
    override in provenance for the audit trail.
    """
    base = parameters["biomass_carrying_capacity"]
    _override_in_band(base, cap_gpl, "carrying_capacity_gpl")
    override = Parameter(
        name="biomass_carrying_capacity",
        value=cap_gpl,
        unit=base.unit,
        tier=base.tier,
        uncertainty=base.uncertainty,
        provenance=Provenance(
            source=base.provenance.source,
            doi=base.provenance.doi,
            conditions=(
                f"scenario opt-in override (decision D-30): carrying_capacity_gpl={cap_gpl:g} g/L, "
                f"replacing the {base.value:g} g/L YAML reference"
            ),
            notes=base.provenance.notes,
        ),
    )
    return parameters.merge(ParameterSet([override]), override=True)


def _override_autolysis_rate(parameters: ParameterSet, rate_per_h: float) -> ParameterSet:
    """Override the reference ``k_autolysis`` with a scenario opt-in value (decision D-34).

    Reached only when a wine scenario passes ``autolysis_rate_per_h > 0`` (which enables the
    autolysis Process), so the rate is the scenario's value rather than the YAML reference —
    letting a demonstration sweep the *sur lie* timescale **within the reference band**
    (:func:`_override_in_band`, decision D-164). Keeps the reference parameter's
    units/tier/uncertainty (only the operating point moves) and records the override in
    provenance, mirroring :func:`_override_carrying_capacity`.
    """
    base = parameters["k_autolysis"]
    _override_in_band(base, rate_per_h, "autolysis_rate_per_h")
    override = Parameter(
        name="k_autolysis",
        value=rate_per_h,
        unit=base.unit,
        tier=base.tier,
        uncertainty=base.uncertainty,
        provenance=Provenance(
            source=base.provenance.source,
            doi=base.provenance.doi,
            conditions=(
                f"scenario opt-in override (decision D-34): autolysis_rate_per_h={rate_per_h:g} "
                f"1/h, replacing the {base.value:g} 1/h YAML reference"
            ),
            notes=base.provenance.notes,
        ),
    )
    return parameters.merge(ParameterSet([override]), override=True)


def _temperature_ramp_schedule(
    scenario: Scenario, t_end_h: float
) -> tuple[float, tuple[ScheduledEvent, ...]]:
    """Compile ``temperature_schedule`` into an initial slope + slope-change events.

    A temperature schedule is a piecewise-*linear* ramp between its ``(day, celsius)``
    knots (decision D-35): between two knots ``dT/dt`` is a single constant, and the
    :class:`~fermentation.core.kinetics.temperature.TemperatureRamp` Process writes that
    constant into ``dT/dt``. This converts the knots (industry units) into canonical hours
    and Kelvin and returns ``(initial_slope, events)`` where ``events`` restart the
    integrator only at the interior knots where the slope **changes** — so collinear knots
    (a straight ramp given by three points) produce a single segment, and a flat or
    single-knot schedule produces no events and a zero initial slope (isothermal). ``T`` is
    held at the nearest knot's value outside the schedule's span (slope 0 before the first
    knot and after the last).
    """
    knots = sorted(scenario.temperature_schedule, key=lambda p: p.day)
    times = [days_to_hours(p.day) for p in knots]
    temps = [celsius_to_kelvin(p.celsius) for p in knots]

    def slope_after(t: float) -> float:
        # Slope of the segment starting at t. Held flat outside the schedule's span.
        if t < times[0] or t >= times[-1]:
            return 0.0
        for i in range(len(times) - 1):
            if times[i] <= t < times[i + 1]:
                dt = times[i + 1] - times[i]
                return (temps[i + 1] - temps[i]) / dt if dt > 0.0 else 0.0
        return 0.0

    initial_slope = slope_after(0.0)
    events: list[ScheduledEvent] = []
    prev = initial_slope
    for bt in sorted({t for t in times if 0.0 < t < t_end_h}):
        s = slope_after(bt)
        # Only a genuine slope change opens a new segment (collinear knots do not); isclose
        # absorbs float noise so a straight multi-point ramp stays one segment.
        if not math.isclose(s, prev, rel_tol=1e-12, abs_tol=1e-15):
            events.append(
                ScheduledEvent(
                    time_h=bt,
                    label=f"temperature_ramp@{bt / 24.0:g}d",
                    param_update={RAMP_RATE: s},
                )
            )
            prev = s
    return initial_slope, tuple(events)


def _inject_temperature_ramp_rate(parameters: ParameterSet, slope_k_per_h: float) -> ParameterSet:
    """Register the scenario's initial temperature-ramp slope as a provenance-backed parameter.

    Reached only when the schedule actually ramps (decision D-35). The rate is a
    scenario-*exact* set-point forcing, not an empirical kinetic constant, so it is
    VALIDATED with a zero-width band (never swept by the ensemble) — but prime directive #2
    still requires it to travel as a :class:`Parameter` with provenance, so it is minted
    here at the boundary (the D-14/D-30/D-34 pattern) rather than inlined. This value seeds
    the first segment; later segments' slopes are supplied by the ``simulate_scheduled``
    events. ``TemperatureRamp`` reads it with a ``0.0`` default, so an un-ramped scenario
    needs no such parameter at all.
    """
    param = Parameter(
        name=RAMP_RATE,
        value=slope_k_per_h,
        unit="K/h",
        tier=Tier.VALIDATED,
        uncertainty=Uncertainty(
            low=slope_k_per_h,
            high=slope_k_per_h,
            note="scenario-exact temperature set-point schedule; not an uncertain parameter",
        ),
        provenance=Provenance(
            source="scenario temperature_schedule",
            conditions=(
                f"initial piecewise-linear ramp slope {slope_k_per_h:g} K/h (decision D-35)"
            ),
            notes="later intervals' slopes are supplied by simulate_scheduled events",
        ),
    )
    return parameters.merge(ParameterSet([param]), override=True)


# -- discrete-intervention verb registry (the winemaking vocabulary boundary) -----------------
#
# ``scenario.interventions`` is a declarative timeline of winemaking verbs in industry units
# (decision D-36). This registry is where each verb's *meaning* lives — which canonical state
# slot a dose lands on, which unit conversion applies, which Processes a pitch enables — exactly
# the layer that owns the initial-composition vocabulary and the temperature-schedule compile
# above (decision D-3). The runtime driver (``simulate_scheduled``) stays verb-agnostic: a verb
# compiles to an opaque :class:`ScheduledEvent` and the driver just segments-and-restarts around
# it, booking each state jump as an :class:`~fermentation.runtime.schedule.ExternalFlow` for the
# conservation ledger. New verbs are added here and nowhere else.


def _iv_check_keys(iv: Intervention, allowed: frozenset[str], verb: str) -> None:
    unknown = set(iv.params) - allowed
    if unknown:
        raise ValueError(
            f"intervention {verb!r} at day {iv.day:g} has unknown param(s) {sorted(unknown)}; "
            f"allowed: {sorted(allowed)}"
        )


def _iv_float(iv: Intervention, key: str, verb: str) -> float:
    """Read a required numeric intervention param, non-negative (the ``_nonneg`` discipline)."""
    if key not in iv.params:
        raise ValueError(
            f"intervention {verb!r} at day {iv.day:g} is missing required param {key!r}"
        )
    try:
        value = float(iv.params[key])
    except (TypeError, ValueError):
        raise ValueError(
            f"intervention {verb!r} param {key!r} must be a number, got {iv.params[key]!r}"
        ) from None
    return _nonneg(value, key)


def _iv_str(iv: Intervention, key: str, verb: str) -> str:
    """Read a required string intervention param (the categorical sibling of :func:`_iv_float`)."""
    if key not in iv.params:
        raise ValueError(
            f"intervention {verb!r} at day {iv.day:g} is missing required param {key!r}"
        )
    value = iv.params[key]
    if not isinstance(value, str):
        raise ValueError(f"intervention {verb!r} param {key!r} must be a string, got {value!r}")
    return value


def _verb_add_dap(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_dap`` — dose diammonium phosphate, BOTH of its ions (decisions D-36, D-210).

    Doses DAP by mass (``dap_gpl``, faithful to the commercial additive) and converts to the
    assimilable-nitrogen jump on the lumped ``N`` slot via the sourced ``dap_nitrogen_fraction``
    (exact (NH₄)₂HPO₄ stoichiometry, VALIDATED). The D-36 headline consequence is a *timing*
    effect the static D-29 lever could not produce: restoring N mid-ferment momentarily closes
    the inverse H₂S gate ``K_h2s_n/(K_h2s_n+N)`` while sugar (hence the flux the gate multiplies)
    is still present.

    **D-210 makes the dose a SALT rather than a nitrogen injection.** Until then this verb wrote
    one slot, and the charge consequences of the other two thirds of the compound were both
    missing:

    * the dosed nitrogen was charged at its medium's *composition average* (D-209's ``z̄``),
      when DAP's nitrogen is pure ammonium and carries **+1 per mole N** — 3× wine's average.
      ``nitrogen_charge_excess`` now carries the difference, re-mixed against whatever the pool
      already held by :func:`~fermentation.core.acidbase.remix_nitrogen_charge_excess`;
    * the **phosphate was dropped entirely**, on the correct premise that no phosphorus pool
      exists and the incorrect inference that therefore nothing was owed. The charge balance
      needs a total, not a pool, and ~0.95 equivalents per mole of anion charge were going
      unbooked — so the verb was adding a cation with no counter-anion, i.e. booking a salt as
      a base. It now writes the ``phosphate`` slot too.

    The two are ONE change, not two independent ones, and the ordering matters for what may be
    claimed: the ammonium half moves a *dry* wine's final pH by exactly zero (the dosed nitrogen
    leaves either way) while tripling the dose-time excursion, and the phosphate half is what
    permanently acidifies (−0.159 pH at 1.1 g/L). Shipping the nitrogen half alone would have
    credited it with +0.242 pH of dose-time rise that the phosphate cancels 72 % of.

    Both new writes are guarded on slot presence rather than assumed, because a caller-supplied
    schema predating D-210 must keep working exactly as it did.

    **The two charge writes ride D-179's opt-in gate, and that was found in review rather than
    designed.** Without it the halves were gated *differently*: ``nitrogen_charge_molar`` checks
    ``charge_balance_is_populated``, while ``phosphate`` is an ordinary acid slot that
    ``_totals_molar`` reads unconditionally. Two things went wrong on a wine with no ``initial_ph``
    and no dosed acids — an empty balance, which is the shape of the Palma 2012 benchmark:

    * the dose booked the salt's ANION with its CATION suppressed, the mirror of the strong-base
      artefact D-179's gate exists to prevent, and more acidic than the chemistry allows;
    * worse, writing an acid slot **opened the gate**, because that is what the gate tests. So a
      *nutrient addition* switched on the whole D-209 nitrogen term for a beverage whose pH no
      scenario had supplied: pH 3.10 → **4.53** at the dose instant and **−0.647 pH** at the end,
      reaching every pH-reading Process. A dose is not pH information (D-179's own principle).

    Gating both writes together fixes both: an unanchored wine's balance gains nothing, the gate
    cannot be opened by a dose, and the ``N`` jump D-36 shipped for the H₂S gate is untouched.
    Evaluated on ``y`` (pre-mutation) and inside the event, so it is not a discontinuity in the
    RHS — an event already is one, which is why the same gate would be forbidden in a Process
    ([[feedback-a-gate-is-a-discontinuity-the-solver-probes]]).
    """
    _iv_check_keys(iv, frozenset({"dap_gpl"}), "add_dap")
    dap_gpl = _iv_float(iv, "dap_gpl", "add_dap")
    try:
        n_fraction = parameters["dap_nitrogen_fraction"].value
        phosphate_fraction = parameters["dap_phosphate_fraction"].value
        dosed_charge = parameters["dap_nitrogen_charge"].value
    except KeyError as exc:  # additions.yaml not loaded (caller-supplied parameter_paths)
        raise ValueError(
            "intervention 'add_dap' needs 'dap_nitrogen_fraction', 'dap_phosphate_fraction' and "
            f"'dap_nitrogen_charge' but one is missing ({exc}); include additions.yaml in "
            "parameter_paths (the default lookup merges it automatically)."
        ) from None
    added_n_gpl = dap_gpl * n_fraction
    added_phosphate_gpl = dap_gpl * phosphate_fraction
    n_slice = schema.slice("N")
    phosphate_slice = schema.slice("phosphate") if "phosphate" in schema else None
    excess_key = acidbase.NITROGEN_CHARGE_EXCESS_KEY
    excess_slice = schema.slice(excess_key) if excess_key in schema else None

    def mutate(
        mutate_schema: StateSchema, y: FloatArray, params: Mapping[str, float]
    ) -> FloatArray:
        out = y.copy()
        out[n_slice] += added_n_gpl
        # The two CHARGE writes are atomic and ride D-179's gate together (decision D-210). Both,
        # or neither: a state must never be able to say "this pool is ammonium-rich" while the
        # phosphate that arrived with the ammonium is missing.
        if acidbase.charge_balance_is_populated(y, mutate_schema):
            if excess_slice is not None:
                # Re-mixed against the pool BEFORE the nitrogen lands — the weighting is
                # pool-before against dose, so `y` and not `out`.
                # The RUNNING map (decision D-235), not the compile-time one. This call reads
                # `nitrogen_uptake_charge_<medium>` — the z̄ an excess is measured RELATIVE TO —
                # and the charge balance reads the member's own z̄ back at every later step. A
                # nominal z̄ here would store the pool's excess against one reference and have it
                # read against another: 0.0108 pH at the dose on the worst of 24 members, ~7×
                # smaller than `set_ph`'s and the same defect. `dosed_charge` stays compile-time
                # on purpose — `dap_nitrogen_charge` is a property of the COMPOUND, is not in the
                # sampled set (it reaches no Process's `reads`), and so has no member value.
                out[excess_slice] = acidbase.remix_nitrogen_charge_excess(
                    float(y[n_slice][0]),
                    float(y[excess_slice][0]),
                    added_n_gpl,
                    float(dosed_charge),
                    mutate_schema.medium,
                    params,
                )
            if phosphate_slice is not None:
                out[phosphate_slice] += added_phosphate_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_dap@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_add_so2(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_so2`` — dose total SO₂ onto the conserved ``so2_total`` slot (decision D-36).

    Doses total sulfur dioxide by the industry unit (``so2_mgl``, mg/L) and converts to the
    canonical g/L jump on ``so2_total`` — the same slot the initial ``so2_total_mgl`` addition
    lands on (D-22/D-28). Free/bound/molecular SO₂ are then re-derived at the solved pH from that
    total (D-28), so a mid-ferment addition raises the antimicrobial molecular fraction from that
    time forward. SO₂ carries neither carbon nor nitrogen, so this flow perturbs neither elemental
    ledger — the single-run carbon and nitrogen balances still close with no correction term.
    """
    _iv_check_keys(iv, frozenset({"so2_mgl"}), "add_so2")
    so2_mgl = _iv_float(iv, "so2_mgl", "add_so2")
    if "so2_total" not in schema:
        raise ValueError(
            f"intervention 'add_so2' at day {iv.day:g} needs a 'so2_total' slot, but medium "
            f"{schema!r} has none (SO₂ is a wine-only pool, decision D-22)"
        )
    added_gpl = mgl_to_gpl(so2_mgl)
    so2_slice = schema.slice("so2_total")

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[so2_slice] += added_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_so2@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_add_copper(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_copper`` — copper-fine reductive sulfur (H₂S + mercaptans) out of the wine (D-44/D-45).

    The remediation half of the reductive-fault beat. Copper (Cu²⁺, dosed as copper sulfate) binds
    the dissolved reductive-sulfur compounds into copper-sulfhydryl complexes that are no longer
    volatile — so the *odour* goes, while the complexes themselves stay dispersed in the wine
    (D-191/D-193; the settle-out-with-the-lees account is the retracted one) — the standard fix for
    the sur-lie "reduction"
    :class:`~fermentation.core.kinetics.hydrogen_sulfide.AutolyticHydrogenSulfide` (D-44) and
    :class:`~fermentation.core.kinetics.mercaptans.AutolyticMercaptan` (D-45) build up un-stripped
    post-dryness. Doses copper by the industry unit (``copper_mgl``, mg/L Cu) and binds, in **order
    of affinity**:

    1. **H₂S first** — copper sulfide (Cu²⁺ + H₂S → CuS↓ + 2 H⁺, **1:1 mol**), CuS being far more
       insoluble (Ksp ~10⁻³⁶) than the mercaptide, so sulfide is bound preferentially. Capacity
       ``copper·copper_h2s_binding``; removes ``min(h2s, capacity)``.
    2. **Mercaptans with the leftover copper** — copper mercaptide (Cu²⁺ + 2 RSH → Cu(SR)₂↓ + 2 H⁺,
       **1:2 mol**, so a gram of Cu binds ~2.8× the thiol mass it does sulfide), capacity
       ``copper_left·copper_mercaptan_binding``; removes ``min(mercaptans, capacity)``.

    Copper in excess simply clears all dissolved reductive sulfur (the real outcome). Copper is
    also imperfect on mercaptans and useless on the disulfides they oxidise to (see
    ``copper_mercaptan_binding``). The ``methanethiol`` slot is wine-only, so on a medium without it
    copper binds H₂S alone.

    **Ledger (D-193 changed this).** Both bound species now go to the wine's complexed reservoirs
    rather than out of the system, so **the verb moves no carbon at all**: ``total_carbon`` closes
    across the fining jump to machine precision (measured −1.19e-23 g/L, against −1.09e-06 g/L
    before) and the external flow it books is elementally empty. H₂S was always carbon-free (D-29,
    the ``add_so2`` precedent); the thiol's carbon used to be booked as a **negative external
    flow** leaving the wine, which is the half of D-45 that D-191's ``Flags:`` marker identified as
    resting on a retracted mechanism.

    **3. The copper stays behind (D-191), which closes D-149's "the two coppers never meet".**
    Until D-191 this verb wrote ``h2s``/``methanethiol`` and **never** the ``copper`` slot D-134
    gave wine, so fining could not raise a wine's oxidation rate — though residual copper is the
    commonest way real wine's copper level rises. The stated reason for leaving it was that
    "nothing sources a residual-copper fraction"; **that was false, and the source settles the
    mechanism as well as the number.** Understanding Wine Chemistry 2nd ed. §26.2.4.1 reports white
    wines dosed to 1 mg/L Cu with equimolar H₂S retaining >95 % of their copper through *filtering*
    or *5 days settling plus racking*, and concludes copper additions "are not necessarily a
    'fining' operation, as the copper remains in solution, albeit in a different form" — Ch. 24
    gives that form as dispersed Cu(I)-sulfhydryl nanoparticles, and names the settles-with-the-lees
    account this docstring used to give as an incorrect assumption of older textbooks. So the verb
    now credits ``copper_fining_residual_fraction · dose`` to the ``copper`` slot.

    The sulfide arithmetic above is **untouched** — the odour still goes. What changes is that a
    fined wine now carries the oxidative cost UWC Table 26.2.1 prints for this treatment: at the
    0.5 mg/L dose, ``f_copper`` 1.000 → 1.285, i.e. ~29 % faster phenolic browning. **Ledger:**
    ``copper`` is off every ledger (a trace metal), so the credit weighs exactly zero in carbon and
    nitrogen — the flow is booked and balances, the ``adjust_cations`` shape.

    **4. The sulfur stays behind too (D-193), which spends D-191's ``Flags:`` on D-45.**
    Until D-193 the bound sulfur was **annihilated**: ``h2s``/``methanethiol`` were decremented and
    nothing received them, so a fining permanently destroyed sulfur that the source says is still
    in the bottle. The same paragraph that licenses the copper credit says the complexes
    **regenerate H₂S during storage**, and the model has owned that reservoir since D-135. So the
    bound mass is now **transferred, not deleted**::

        removed_h2s  -> bound_h2s          (1:1, the same molecule in a different binding state)
        removed_merc -> bound_methanethiol (likewise)

    **The whole removed mass moves, with no retention fraction applied.** The 0.95 above is a
    *printed lower bound* on copper retention measured after filtering or racking — operations this
    verb does not model and cannot know happened (the repo has a separate ``rack`` verb, and this
    parameter's own band note says those cellar sinks are "separate operations on longer
    timescales, not part of this event"). Scaling the sulfur by it was measured and rejected: it
    would convert a bound into a loss fraction and leave an invented 5 % as the only carbon
    outflow, re-asserting the retracted mechanism at 1/20th size instead of discharging it. The
    two shares therefore differ on purpose — copper 0.95 because that share is *sourced*, sulfur
    1.0 because no source measures any loss at this event.

    **The free pools are untouched at the event**: the removal arithmetic in 1–2 is byte-for-byte
    what it was, so the odour fix is unchanged and only the destination differs.

    **This is NOT the copper-coupled release D-135 refused.** That refusal is about the release
    *rate constant* — a PLS regression coefficient "is not a binding stoichiometry" — and it stands
    untouched: ``k_bound_h2s_release`` still reads no copper. What moves here is the *reservoir*,
    seeded from a stoichiometry this verb already computes (``copper_h2s_binding``, 1:1 CuS).

    **Named extrapolation on the thiol half.** Franco-Luesma & Ferreira 2016 give bonded MeSH a
    *negative* copper coefficient and conclude copper's role as a MeSH trapping agent "is not
    really important" — a statement about *natural* wine copper, not about a deliberate dose, which
    forms Cu(SR)₂ by the very stoichiometry step 2 uses. So the destination is sound (UWC Ch. 24's
    dispersed sulfhydryl nanoparticles cover the mercaptide — D-191's ``Flags:`` says so in those
    words), but the *rate* at which fined thiol comes back out is unmeasured, and this reuses the
    natural pool's 8.1 %/yr. Deleting it instead is the strictly worse claim: an instant, permanent,
    total removal that the source contradicts outright.

    Both reservoirs are wine-only, so the routing is guarded on the slots exactly as the copper
    credit is: a medium without them keeps the removal-only behaviour.
    """
    _iv_check_keys(iv, frozenset({"copper_mgl"}), "add_copper")
    copper_mgl = _iv_float(iv, "copper_mgl", "add_copper")
    if "h2s" not in schema:  # both current media carry h2s; guard for symmetry with add_so2
        raise ValueError(
            f"intervention 'add_copper' at day {iv.day:g} needs an 'h2s' slot, but medium "
            f"{schema!r} has none"
        )
    try:
        binding_h2s = parameters["copper_h2s_binding"].value
    except KeyError as exc:  # additions.yaml not loaded (caller-supplied parameter_paths)
        raise ValueError(
            "intervention 'add_copper' needs 'copper_h2s_binding' but it is missing "
            f"({exc}); include additions.yaml in parameter_paths (the default lookup merges "
            "it automatically)."
        ) from None
    copper_gpl = mgl_to_gpl(copper_mgl)  # g/L Cu dosed
    h2s_slice = schema.slice("h2s")
    # The thiol pool is wine-only; bind it with leftover copper iff the slot exists (D-45). The
    # slot is `methanethiol` since D-110 (it was the false plural `mercaptans` through D-109);
    # `copper_mercaptan_binding` keeps its name deliberately — copper mercaptide precipitation
    # Cu2+ + 2 RSH -> Cu(SR)2 is real class chemistry, general over thiols, and asserts no lump.
    has_methanethiol = "methanethiol" in schema
    binding_merc = parameters["copper_mercaptan_binding"].value if has_methanethiol else 0.0
    merc_slice = schema.slice("methanethiol") if has_methanethiol else None
    # D-191: the retained share of the dose, credited to the `copper` slot below. Guarded (not
    # gated) on the slot the same way the thiol half is: `copper` is wine-only, and a medium
    # without it simply keeps the removal behaviour. Missing parameter => hard error rather than a
    # silent 0.0 default, because a 0.0 here would restore exactly the pre-D-191 behaviour while
    # looking like a configured choice (the D-45 hard-zero defect).
    has_copper = "copper" in schema
    copper_slice = schema.slice("copper") if has_copper else None
    # D-193: where the bound sulfur goes. Both reservoirs are wine-only, so each is guarded on its
    # own slot the way the copper credit and the thiol half are — a medium without them keeps the
    # removal-only behaviour rather than erroring. No parameter is read: the transfer is 1:1 by
    # construction (one molecule changing binding state), which is exactly why no retention
    # fraction is applied to it (see the docstring's part 4).
    bound_h2s_slice = schema.slice("bound_h2s") if "bound_h2s" in schema else None
    bound_merc_slice = (
        schema.slice("bound_methanethiol") if "bound_methanethiol" in schema else None
    )
    if has_copper:
        try:
            residual_fraction = parameters["copper_fining_residual_fraction"].value
        except KeyError as exc:
            raise ValueError(
                "intervention 'add_copper' needs 'copper_fining_residual_fraction' but it is "
                f"missing ({exc}); include additions.yaml in parameter_paths (the default lookup "
                "merges it automatically)."
            ) from None
    else:
        residual_fraction = 0.0

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        # 1. H₂S first (higher affinity). Clamp present ≥ 0 so a solver undershoot is not "removed".
        h2s_present = max(float(out[h2s_slice][0]), 0.0)
        removed_h2s = min(h2s_present, copper_gpl * binding_h2s)
        out[h2s_slice] = float(out[h2s_slice][0]) - removed_h2s
        # ...and it is COMPLEXED, not destroyed (D-193): the bound mass lands in the reservoir the
        # aging release Process empties, so a fined wine can turn reductive again in the bottle the
        # way the source says it does. Whole mass, no fraction — see the docstring's part 4.
        if bound_h2s_slice is not None:
            out[bound_h2s_slice] = float(out[bound_h2s_slice][0]) + removed_h2s
        # 2. Mercaptans with the copper left after binding H₂S (its stoichiometric share).
        if merc_slice is not None:
            copper_left = max(copper_gpl - removed_h2s / binding_h2s, 0.0)
            merc_present = max(float(out[merc_slice][0]), 0.0)
            removed_merc = min(merc_present, copper_left * binding_merc)
            out[merc_slice] = float(out[merc_slice][0]) - removed_merc
            # Same transfer for the thiol (D-193), and this one also fixes the CARBON: the two
            # slots carry the identical carbon weight, so the fining stops booking a carbon
            # outflow that rested on the retracted precipitation mechanism.
            if bound_merc_slice is not None:
                out[bound_merc_slice] = float(out[bound_merc_slice][0]) + removed_merc
        # 3. The copper STAYS (D-191). The retained fraction applies to the WHOLE dose, not to
        # whatever is left after binding: the source's wines were dosed alongside equimolar H₂S, so
        # >95 % retention is already a post-reaction figure — the copper-sulfhydryl complexes are
        # themselves what stays dispersed. Accumulates onto the slot's current value rather than
        # setting it, so repeated finings add up and the grape/must background (copper_typical) is
        # preserved.
        if copper_slice is not None:
            out[copper_slice] = float(out[copper_slice][0]) + residual_fraction * copper_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_copper@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_add_oxygen(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_oxygen`` — dose dissolved oxygen onto the ``o2`` aging substrate (decision D-71).

    The oxidative-aging substrate lever: doses dissolved O₂ by the industry unit (``o2_mgl``, mg/L)
    and converts to the canonical g/L jump on the ``o2`` slot — the ingress a finished wine/beer
    takes up in bottle, under micro-oxygenation, or across a barrel. One dose models a single
    exposure (a bottle's total ingress); repeated doses model continuous micro-ox / barrel aging.
    The dosed O₂ is then drawn down by the oxidative aging Processes (once ``begin_aging`` has
    enabled them), each taking its own share of the shared ``o2`` pool (D-73/D-74):
    :class:`~fermentation.core.kinetics.aging.PhenolicBrowning` (medium-agnostic, the DOMINANT sink)
    oxidises phenolics to brown pigment, accumulating the ``A420`` browning index;
    :class:`~fermentation.core.kinetics.aging.OxidativeAcetaldehyde` (medium-agnostic) oxidises
    ethanol → acetaldehyde at its per-O₂ yield; and
    :class:`~fermentation.core.kinetics.aging.SulfiteOxidation` (wine) diverts O₂ to spend free SO₂
    —
    so a dose **browns** the finished wine/beer and raises its acetaldehyde ('sherry'/oxidised),
    with
    browning taking most of the O₂ (suppressing the acetaldehyde) and SO₂ intercepting O₂ while it
    lasts; any acetaldehyde formed is further mopped up by dosed SO₂ via the D-47 binding
    equilibrium
    for free.

    **The add_so2 pattern exactly** (a carbon-free dosed pool): O₂ carries neither carbon nor
    nitrogen and the ``o2`` slot is off every conservation ledger (``total_carbon``/``total_mass``/
    ``total_nitrogen`` weight only their named pools), so this flow perturbs no elemental balance —
    the single-run carbon and nitrogen ledgers still close with **no** external-flow correction
    term (unlike the carbon-bearing ``add_acid``/``add_sugar`` doses). Concentration model: no
    volume change on the addition (the shared verb caveat).

    Medium-agnostic (``o2`` is in ``_common_specs``, so both media carry it). Ordering note: dosing
    O₂ *without* a ``begin_aging`` leaves it inert in the slot — the oxidation Process stays
    disabled until the aging phase begins — so the natural usage is ``begin_aging`` at the
    ferment/aging boundary plus ``add_oxygen`` for each exposure over the aging tail.
    """
    _iv_check_keys(iv, frozenset({"o2_mgl"}), "add_oxygen")
    o2_mgl = _iv_float(iv, "o2_mgl", "add_oxygen")
    if "o2" not in schema:  # both current media carry o2; guard for symmetry with add_so2
        raise ValueError(
            f"intervention 'add_oxygen' at day {iv.day:g} needs an 'o2' slot, but medium "
            f"{schema!r} has none (the dissolved-oxygen aging substrate, decision D-71)"
        )
    added_gpl = mgl_to_gpl(o2_mgl)
    o2_slice = schema.slice("o2")

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[o2_slice] += added_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_oxygen@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_add_ascorbate(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_ascorbate`` — dose ascorbic acid (vitamin C) as an antioxidant (decision D-202).

    **The verb exists because the pool's default is 0, and the default is 0 because the source
    says so.** *Understanding Wine Chemistry* 2nd ed. §24.4.3.2: *"There is a small amount in
    grapes that is quickly depleted during fermentation, such that new wine has a negligible
    ascorbic acid content"*, and ascorbic acid is *"a permitted winemaking additive in most
    wine-producing countries"* (Ch. 27). So ascorbate is a **dose**, never a must-composition
    property — the opposite call to D-134's ``copper`` (where 0 was unphysical and the slot had to
    be seeded from the typical level) and the same call as ``add_so2``/``add_oxygen``.

    **The add_so2 / add_oxygen pattern exactly** (a dosed pool that is off every ledger): the
    ``ascorbate`` slot carries no conservation weight — the carbon is exogenous, the oxidation
    product is untracked, and the o-diphenol the reaction regenerates is off-ledger by fork D2 (see
    ``M_ASCORBIC``) — so this flow perturbs no elemental balance and needs no external-flow
    correction term, unlike the carbon-bearing ``add_acid``/``add_sugar`` doses. Concentration
    model: no volume change on the addition (the shared verb caveat).

    Wine-only, because the ``ascorbate`` slot is (the consumer,
    :class:`~fermentation.core.kinetics.oxidative_cascade.QuinoneAscorbateReduction`, is a wine
    cascade node). **Ordering note, and it differs from ``add_oxygen``'s:** a dose is *not* inert
    before ``begin_aging`` in the way dosed O2 is. The consumer is aging-gated, so it will not run
    early — but the ascorbate sits in the slot and is fully available the moment aging starts.
    That is the faithful reading of an addition at crush or at bottling, which is when a winemaker
    actually makes it. UWC's reference dose, and the one D-200/D-202 measure against, is
    **60 mg/L**; the EU permits up to 250 mg/L. Neither is enforced here — the verb takes the
    winemaker's number, as every other dosing verb does.
    """
    _iv_check_keys(iv, frozenset({"ascorbate_mgl"}), "add_ascorbate")
    ascorbate_mgl = _iv_float(iv, "ascorbate_mgl", "add_ascorbate")
    if "ascorbate" not in schema:
        raise ValueError(
            f"intervention 'add_ascorbate' at day {iv.day:g} needs an 'ascorbate' slot, but "
            f"medium {schema!r} has none (the wine-only dosed antioxidant, decision D-202)"
        )
    added_gpl = mgl_to_gpl(ascorbate_mgl)
    ascorbate_slice = schema.slice("ascorbate")

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[ascorbate_slice] += added_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_ascorbate@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_seal_bottle(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet, scenario: Scenario
) -> ScheduledEvent:
    """``seal_bottle`` — dose the closure's OWN trapped-air charge at bottling (decision D-187).

    **The gap this closes.** D-136 turned closure oxygen from a dosed stock into a continuous flow,
    but deliberately shipped only the STEADY column of Lopes et al. 2007's Table I. Its first-month
    column — 10–150× higher, and a different mechanism — was left to the author as a bare
    ``add_oxygen`` number, which made it the one oxygen input on this axis with no provenance: the
    quantity that decides whether a freshly bottled wine is modelled at all had to be invented at
    the call site. D-136's own Next said the data was already in hand ("P1/P2 both quantify it per
    closure … left out deliberately to keep the D-133 line clean, **not** for lack of data"). This
    verb spends it.

    **It is ``add_oxygen``'s dose with the number taken out of the author's hands.** The mutation is
    identical — one bolus onto the ``o2`` slot, carbon- and nitrogen-free, no volume change — and
    the ONLY difference is where the magnitude comes from: ``bottling_burst_<closure>`` in
    ``closure.yaml``, resolved from ``scenario.closure``. That is the whole point, so the verb takes
    **no params at all**: a scenario cannot say "seal it with a natural cork but dose 4 mg/L", which
    would be an ``add_oxygen`` wearing a sourced verb's name. An author who wants their own number
    still has ``add_oxygen``, and an author whose bottling line adds headspace oxygen on top of the
    closure's charge should dose BOTH — the line's contribution is not a closure property and is
    deliberately not in these parameters (see ``closure.yaml``'s bottling-burst header).

    **Why it needs the scenario, and is therefore in :data:`_SCENARIO_INTERVENTION_VERBS`.** The
    dose is a function of a *scenario-level* field, not of the intervention. Every other verb is
    handed only its own ``Intervention``, the schema and the parameters — the ``set_ph`` gate
    (D-186) already had to reach for ``scenario.initial`` and did it from
    :func:`_compile_interventions` because a verb cannot see the scenario at all. This one needs the
    scenario for its *value*, not merely for a gate, so it takes it as a fourth argument rather than
    duplicating ``closure`` into ``params`` where the two could disagree.

    **Two gates, both in :func:`_compile_interventions` and both loud** (see there): the scenario
    must name a ``closure``, and the seal must not precede ``begin_aging``. The second is not
    fastidiousness — the shipped charge is P1's first-month total *less the steady permeation over
    the same 30 days*, an anti-double-count that assumes
    :class:`~fermentation.core.kinetics.aging.ClosureOxygenIngress` is switched on across that
    month. Sealing before the aging boundary would subtract a flux nobody paid. A third gate
    (wine-only) was written and deleted because it provably cannot fire — D-136 rejects a
    ``closure`` on beer before any intervention compiles, so the first gate is already the whole
    medium check.

    **``hermetic`` doses exactly 0.0**, so the verb is isolable in the prime-directive-3 sense: a
    hermetic seal is byte-for-byte the run without it, and that zero is *named* in ``closure.yaml``
    (the D-45 lesson) rather than expressed by leaving the intervention out.
    """
    _iv_check_keys(iv, frozenset(), "seal_bottle")
    if "o2" not in schema:  # unreachable via the medium gate below; the add_oxygen symmetry guard
        raise ValueError(
            f"intervention 'seal_bottle' at day {iv.day:g} needs an 'o2' slot, but medium "
            f"{schema!r} has none (decision D-187)"
        )
    assert scenario.closure is not None  # guaranteed by the _compile_interventions gate
    added_gpl = _bottling_burst(scenario.closure, parameters)
    o2_slice = schema.slice("o2")

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[o2_slice] += added_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"seal_bottle@{iv.day:g}d",
        mutate=mutate,
    )


#: The oak toast levels :func:`_verb_add_oak` accepts (decision D-77) and the extractives it doses:
#: the FIVE aroma extractives (four D-77 + ``furaneol`` caramel D-94) plus ``ellagitannin`` (the
#: D-78 taste/O₂-scavenging bridge).
#: The categorical ``toast`` selects the per-gram yield set (``oak_yield_<compound>_<toast>`` in
#: oak.yaml); the compound → ceiling-slot pairing mirrors ``aging._OAK_COMPOUND_CEILINGS``. So one
#: ``add_oak`` dose sets all six saturation ceilings (aroma + tannin) from a single ``oak_gpl``/
#: ``toast`` choice.
_OAK_TOASTS = ("light", "medium", "heavy")
_OAK_COMPOUNDS = (
    "whiskey_lactone",
    "vanillin",
    "guaiacol",
    "eugenol",
    "furaneol",  # caramel/toffee — the caramel furanone (decision D-94)
    "ellagitannin",
)
#: The ex-spirit barrel types :func:`_verb_add_oak` accepts for the D-92 residual-spirit soak-back:
#: a first-fill ex-``spirit`` barrel donates residual ethanol (``spirit_soak_ethanol_<spirit>`` in
#: oak.yaml) into the beverage, raising ABV. ``bourbon`` this beat (whiskey/rum extensible); the
#: categorical is how the caller asserts a soaked barrel (soak-back is a barrel, not chips).
_OAK_SPIRITS = ("bourbon",)
#: The oak-aroma extractives an ex-spirit barrel's residual spirit BUMPS the ceiling of (D-93/D-94):
#: the bourbon-barrel aroma soak-back. A DELIBERATE subset of :data:`_OAK_COMPOUNDS` — vanilla +
#: coconut + char are bourbon's signature (``vanillin``/``whiskey_lactone``/``guaiacol``, D-93) and
#: ``furaneol`` its prominent CARAMEL/toffee note (D-94, bourbon matures in charred new oak);
#: ``eugenol`` (clove) is not a bourbon note and ``ellagitannin`` is a wood taste tannin, so both
#: are excluded.
#: Each is bumped by ``spirit_soak_<compound>_<spirit>`` (toast- and ``oak_gpl``-independent, ×
#: ``spirit_scale``); ``OakExtraction`` then leaches it in gradually — a CEILING bump is the only
#: form additive with the wood pool (a bolus into the pool is erased by the extraction gate, D-93).
_OAK_SPIRIT_AROMAS = ("vanillin", "whiskey_lactone", "guaiacol", "furaneol")


def _verb_add_oak(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_oak`` — put the beverage in oak, setting each extractive's ceiling (D-77/D-78/D-86).

    The oak-extraction substrate lever, the aging-axis sibling of ``add_oxygen``: ``params`` names
    the oak-contact dose ``oak_gpl`` (the generalized chips-g/L / barrel surface-to-volume dose) and
    the categorical ``toast`` (``light``/``medium``/``heavy`` — the ``add_acid`` string-param move).
    For each of the five extractives — the four aroma compounds (D-77) plus ``ellagitannin`` (the
    D-78 taste/O₂-scavenging tannin) — it computes the **saturation ceiling** ``oak_gpl ·
    oak_yield_<compound>_<toast>`` (the provenance-backed toast-specific per-gram yields in
    ``oak.yaml``) and writes it to that compound's **set-and-hold** ceiling state slot. The
    :class:`~fermentation.core.kinetics.aging.OakExtraction` Process (enabled by ``begin_aging``)
    then rises the extracted pools toward those ceilings — so the toast selects the aroma *profile*
    (light → coconut-dominant, medium → vanilla, heavy → smoky/clove) and ``oak_gpl`` scales the
    ceilings linearly. The ellagitannin ceiling is set the same way; the D-78
    :class:`~fermentation.core.kinetics.aging.EllagitanninOxidation` sink then draws that tannin
    down as it scavenges O₂ (oak protection), so oaking a beverage both flavours it and buffers its
    redox.

    **The add_oxygen pattern (a dosed off-ledger substrate), NOT begin_aging.** Like ``add_oxygen``,
    this verb only doses — it does **not** enable the Process (``begin_aging`` does, alongside the
    other aging Processes). So the natural usage is ``begin_aging`` at the ferment/aging boundary,
    plus ``add_oak`` for the oak charge; a second ``add_oak`` **raises** the ceilings (a fresh
    charge / more chips — the ``+=`` dose idiom). Note ``+=`` and ``fill_number`` (D-91, below) are
    **orthogonal** levers, not coarse-vs-fine versions of one: ``+=`` adds *more oak contact* and
    **raises** the ceiling, whereas ``fill_number`` models the *same oak, more depleted* by prior
    fills and **lowers** it. The
    ceiling slots are **off every ledger** (wood-derived, the ``iso_alpha`` precedent), so the
    jump perturbs no elemental balance — the run-wide carbon/nitrogen ledgers close with **no**
    correction term (like ``add_oxygen``; unlike carbon-bearing ``add_acid``/``add_sugar``).
    Concentration model: no volume change on the addition (the shared verb caveat).

    **Wine + barrel-beer** (decision D-86): the oak slots are carried by both ``wine_schema`` and
    ``beer_schema`` (via ``core.media._oak_specs``), so ``add_oak`` works on either medium. Only a
    bare/other medium with no ``whiskey_lactone`` slot raises. Guards that ``oak.yaml`` is loaded
    (the ``add_dap`` discipline) so a caller-supplied ``parameter_paths`` without it fails loudly
    HERE at compile, not as a bare ``KeyError`` when the verb reads a yield.

    **Barrel fill-number depletion (decision D-91).** An OPTIONAL ``fill_number`` (int ≥ 1, default
    1) counts the barrel's use: a reused barrel has a depleted accessible extractable pool, so it
    sets LOWER ceilings than a fresh first-fill one at the same ``oak_gpl``/``toast``. Every ceiling
    is scaled by ``oak_fill_retention ** (fill_number − 1)``, so ``fill_number = 1`` is UNSCALED
    (``r**0 == 1.0`` exactly ⇒ byte-for-byte the pre-D-91 dose) and each prior fill geometrically
    discounts the extractables — the signature lever of barrel-aged BEER programs (a first-fill
    bourbon barrel for the imperial stout, then the neutralised barrel for a sour). This is an
    ACROSS-FILL dose input (barrel history known at charge time), NOT a within-fill dynamic
    reservoir; the mechanistic finite-reservoir model and per-compound retention are documented
    refinements deferred here. ``fill_number`` is validated int-valued ≥ 1 (a "zeroth
    fill" is meaningless), and ``oak_fill_retention`` is read only when it bites (``fill_number ≠
    1``), so a fresh fill stays inert even against a partial ``oak.yaml``.

    **Bourbon-barrel spirit soak-back (decision D-92).** An OPTIONAL categorical ``spirit`` (v1:
    ``"bourbon"``) marks the barrel as an ex-spirit cask: its staves are soaked with residual
    high-ABV spirit that leaches back into the beverage, DONATING ethanol and RAISING ABV (the
    "a bourbon-barrel imperial stout gains ~1% ABV from the barrel" effect). A DISCRETE ethanol
    dose to the core ``E`` slot (the ``add_oxygen`` precedent), ``spirit_soak_ethanol_<spirit> ×
    spirit_soak_retention ** (fill_number − 1)`` g/L — a SEPARATE contribution from the wood
    extractives (the ethanol is from the SPIRIT, not the wood), so it does NOT touch the ceilings.
    Ethanol is ON the carbon+mass ledger, but the scheduler books this dose as a POSITIVE external
    flow (the ``add_sugar`` precedent), so the run-wide ``final == initial + Σ flows`` still closes.
    DECOUPLED from ``oak_gpl`` (soak-back is a barrel, not a chips/S:V, phenomenon) and anchored
    straight to the observed ABV gain. Residual spirit depletes with reuse via its OWN steep
    ``spirit_soak_retention`` (spirit ~gone by fill 2–3, far faster than the wood's
    ``oak_fill_retention``), read only when it bites. ``spirit`` DEFAULTS ABSENT ⇒ no ethanol dose
    ⇒ **byte-for-byte** the pre-D-92 charge.

    **Bourbon-barrel aroma soak-back (decision D-93).** The same ``spirit`` also carries the
    residual spirit's own aroma **congeners** — bourbon matures in **charred new oak**, so its
    residual spirit reads vanilla/coconut/char-forward. So a ``spirit`` dose ALSO **bumps the
    ceilings** of the bourbon-signature aroma extractives (:data:`_OAK_SPIRIT_AROMAS`:
    ``vanillin``/``whiskey_lactone``/``guaiacol`` — not clove ``eugenol`` or the taste tannin
    ``ellagitannin``) by ``spirit_soak_<compound>_<spirit> × spirit_scale`` g/L, and
    :class:`~fermentation.core.kinetics.aging.OakExtraction` then leaches them in **gradually** on
    top of the wood diffusion. A CEILING bump — NOT a bolus into the extracted pool, which the
    extraction gate (``gap = ceiling − conc``) would ERASE, giving ``max(wood, spirit)`` not the
    sum; bumping the ceiling is the **only** wood + spirit **additive** form. Legal because the
    aroma ceilings are **off the carbon/mass ledger** (wood-derived, ``iso_alpha`` precedent) — so
    unlike the on-ledger ethanol (FORCED to a discrete dose lest a gradual leach create carbon
    within-segment), the aroma leach is gradual for free, the **more faithful** form. Toast- and
    ``oak_gpl``-INDEPENDENT (the congener profile is set by the bourbon's char, not the cooper's
    toast, and residual spirit is a barrel not a chips/S:V property), depleting with reuse by the
    SAME ``spirit_scale`` as the ethanol. NOT double-counting: one shared pool bumped, not a
    parallel pool (the D-77 yields stay generic new-oak wood; the ex-bourbon barrel's *depleted
    wood* is the orthogonal ``fill_number`` effect, D-91). ``spirit`` absent ⇒ no bump ⇒
    byte-for-byte the pre-D-92 charge on the aroma ceilings too.

    **Bourbon-barrel CARAMEL soak-back (decision D-94).** The caramel/toffee note D-93 deferred is
    now modelled as ``furaneol`` (HDMF), a fifth oak aroma extractive (:data:`_OAK_COMPOUNDS`) with
    its own toast yields (``oak_yield_furaneol_<toast>``, RISING with toast — a thermal
    sugar-degradation furanone of toasted/charred oak) and, being in :data:`_OAK_SPIRIT_AROMAS`, a
    ``spirit_soak_furaneol_<spirit>`` ceiling bump exactly like the three D-93 congeners. The D-93
    collision worry with the D-88 caramelization/``A420`` axis is DISSOLVED, not relocated:
    ``furaneol`` is on the OAK axis — off every ledger (wood/spirit-derived, ``iso_alpha``), never
    touching core ``S`` or the on-ledger ``melanoidin`` — so it cannot perturb D-88's sugar→
    melanoidin carbon closure. ``melanoidin`` is caramelization's *colour body* (on-ledger,
    ``A420``); ``furaneol`` the *volatile aroma* of the same browning chemistry (off-ledger, OAV).
    The genuinely deferred beat is caramel aroma from the *beverage's own* thermal caramelization
    (on-ledger — it would divert a sliver of sugar carbon out of the melanoidin park); this D-94
    pool is oak/spirit-derived only.
    """
    _iv_check_keys(iv, frozenset({"oak_gpl", "toast", "fill_number", "spirit"}), "add_oak")
    oak_gpl = _iv_float(iv, "oak_gpl", "add_oak")
    toast = _iv_str(iv, "toast", "add_oak")
    if toast not in _OAK_TOASTS:
        raise ValueError(
            f"intervention 'add_oak' at day {iv.day:g}: unknown toast {toast!r}; the oak toast "
            f"levels are {sorted(_OAK_TOASTS)} (decision D-77)"
        )
    if "whiskey_lactone" not in schema:  # wine + beer carry the oak slots (D-86); a bare one won't
        raise ValueError(
            f"intervention 'add_oak' at day {iv.day:g} needs a 'whiskey_lactone' slot, but "
            f"medium {schema!r} has none (oak needs the oak-axis slots; wine and beer carry "
            f"them, decisions D-77/D-86)"
        )
    # Barrel fill-number depletion (D-91): a reused barrel extracts LESS. fill_number is an OPTIONAL
    # int >= 1 counting the barrel's use (1 = a fresh first-fill barrel, the default). Each prior
    # fill geometrically discounts every ceiling by oak_fill_retention; fill_number = 1 is UNSCALED
    # (retention ** 0 == 1.0 exactly ⇒ byte-for-byte the pre-D-91 behaviour). Validate int-valued
    # >= 1 (a "zeroth fill" is meaningless — brewers count first/second/third), the toast-string
    # rejection pattern. Read oak_fill_retention only when it BITES (fill_number != 1), so a fresh
    # fill stays inert even against a partial oak.yaml.
    fill_number = 1
    if "fill_number" in iv.params:
        raw = iv.params["fill_number"]
        try:
            fill_f = float(raw)
        except (TypeError, ValueError):
            raise ValueError(
                f"intervention 'add_oak' param 'fill_number' must be a number, got {raw!r}"
            ) from None
        if fill_f < 1.0 or fill_f != int(fill_f):
            raise ValueError(
                f"intervention 'add_oak' at day {iv.day:g}: fill_number must be an integer >= 1 "
                f"(1 = a fresh first-fill barrel), got {raw!r} (decision D-91)"
            )
        fill_number = int(fill_f)
    fill_scale = 1.0
    if fill_number != 1:
        try:
            retention = parameters["oak_fill_retention"].value
        except KeyError:  # oak.yaml not loaded (caller-supplied parameter_paths)
            raise ValueError(
                "intervention 'add_oak' with fill_number > 1 needs 'oak_fill_retention' but it "
                "is missing; include oak.yaml in parameter_paths (decision D-91)."
            ) from None
        fill_scale = retention ** (fill_number - 1)
    # Bourbon-barrel spirit soak-back (D-92): an ex-spirit barrel donates residual ETHANOL, raising
    # ABV — the "barrel-aged stout gains ~1% ABV" effect. An OPTIONAL categorical `spirit` (default
    # absent ⇒ no soak-back ⇒ byte-for-byte the pre-D-92 dose): when given, add a DISCRETE ethanol
    # bolus (the add_oxygen dose precedent) to the core E slot. UNLIKE the off-ledger wood ceilings,
    # ethanol is ON the carbon+mass ledger, but the scheduler books this dose's delta as a POSITIVE
    # external flow (add_sugar precedent), so the run-wide identity final == initial + Σ flows still
    # closes. Anchored straight to the g/L ABV gain and DECOUPLED from oak_gpl (soak-back is a
    # barrel, not a chips/S:V, effect). Depletes with fill_number via its OWN steep
    # spirit_soak_retention (spirit ~gone by fill 2-3, far faster than the wood's
    # oak_fill_retention), read only when it
    # BITES (fill_number != 1). ETHANOL (ABV) only — bourbon AROMA congeners (overlapping the D-77
    # oak aroma pools) and a gradual reservoir leach are deferred refinements.
    ethanol_soak_delta = 0.0
    spirit_aroma_bumps: dict[str, float] = {}  # D-93: compound -> ceiling bump from residual spirit
    if "spirit" in iv.params:
        spirit = _iv_str(iv, "spirit", "add_oak")
        if spirit not in _OAK_SPIRITS:
            raise ValueError(
                f"intervention 'add_oak' at day {iv.day:g}: unknown spirit {spirit!r} "
                f"(known ex-spirit barrels: {', '.join(_OAK_SPIRITS)}, decision D-92)"
            )
        if "E" not in schema:
            raise ValueError(
                f"intervention 'add_oak' at day {iv.day:g} with spirit={spirit!r} needs an 'E' "
                f"(ethanol) slot for the soak-back, but medium {schema!r} has none (decision D-92)"
            )
        soak_name = f"spirit_soak_ethanol_{spirit}"
        try:
            soak_gpl = parameters[soak_name].value
        except KeyError:
            raise ValueError(
                f"intervention 'add_oak' with spirit={spirit!r} needs {soak_name!r} but it is "
                "missing; include oak.yaml in parameter_paths (decision D-92)."
            ) from None
        spirit_scale = 1.0
        if fill_number != 1:  # residual spirit depletes with reuse, its OWN steep retention (D-92)
            try:
                spirit_retention = parameters["spirit_soak_retention"].value
            except KeyError:
                raise ValueError(
                    "intervention 'add_oak' with spirit and fill_number > 1 needs "
                    "'spirit_soak_retention' but it is missing; include oak.yaml (decision D-92)."
                ) from None
            spirit_scale = spirit_retention ** (fill_number - 1)
        ethanol_soak_delta = soak_gpl * spirit_scale
        # Bourbon AROMA soak-back (D-93): the residual spirit also BUMPS the ceilings of the
        # bourbon-signature aroma extractives (vanilla/coconut/char), which OakExtraction then
        # leaches in gradually — a CEILING bump, NOT a bolus into the pool (the extraction gate
        # would erase it). Toast- and oak_gpl-INDEPENDENT flat g/L bumps, × the SAME spirit_scale as
        # ethanol (one residual spirit, one depletion). Off-ledger like the ceilings they raise.
        for compound in _OAK_SPIRIT_AROMAS:
            bump_name = f"spirit_soak_{compound}_{spirit}"
            try:
                bump_val = parameters[bump_name].value
            except KeyError:
                raise ValueError(
                    f"intervention 'add_oak' with spirit={spirit!r} needs {bump_name!r} but it is "
                    "missing; include oak.yaml in parameter_paths (decision D-93)."
                ) from None
            spirit_aroma_bumps[compound] = bump_val * spirit_scale
    ceiling_deltas: dict[str, float] = {}
    for compound in _OAK_COMPOUNDS:
        yield_name = f"oak_yield_{compound}_{toast}"
        try:
            yield_val = parameters[yield_name].value
        except KeyError:  # oak.yaml not loaded (caller-supplied parameter_paths)
            raise ValueError(
                f"intervention 'add_oak' needs {yield_name!r} but it is missing; include oak.yaml "
                "in parameter_paths (the default lookup merges it automatically, decision D-77)."
            ) from None
        # fill_scale discounts the fresh-barrel ceiling by the barrel's use history (D-91); the D-93
        # spirit-aroma bump (0.0 for compounds not in _OAK_SPIRIT_AROMAS, or when no spirit) adds
        # the ex-bourbon barrel's residual-spirit congeners on top of the wood diffusion into the
        # SAME ceiling — additive, so OakExtraction rises the pool to wood + spirit, not the max.
        ceiling_deltas[f"{compound}_ceiling"] = (
            oak_gpl * yield_val * fill_scale + spirit_aroma_bumps.get(compound, 0.0)
        )
    slices = {name: schema.slice(name) for name in ceiling_deltas}
    ethanol_slice = schema.slice("E") if ethanol_soak_delta else None

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        for name, delta in ceiling_deltas.items():
            out[slices[name]] += delta  # += so a second oak charge raises the ceiling (refill)
        if ethanol_slice is not None:  # D-92 spirit soak-back: residual ethanol raises ABV
            out[ethanol_slice] += ethanol_soak_delta
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_oak@{iv.day:g}d",
        mutate=mutate,
    )


#: The lees-associated pools racking removes: inactivated yeast biomass ``X_dead`` and (if autolysis
#: is opted in) the non-assimilable cell-wall ``debris`` (decision D-36); plus **both** *O. oeni*
#: pools — settled dead ``X_mlf_dead`` **and viable ``X_mlf``** (decision D-39). Racking viable
#: bacteria is the deliberate **asymmetry with yeast**: a rack leaves viable yeast ``X`` untouched
#: because it works in *suspension* during AF (racking gross lees leaves it fermenting), but
#: *O. oeni* carries out MLF *on the lees* and is drawn off with them — so racking removes the
#: bacteria that clear diacetyl, the physical twin of the SO₂ kill (the other half of the D-31
#: "rack/SO₂ locks in diacetyl" lever, D-39). Both bacterial pools carry biomass carbon and nitrogen
#: (weighted since D-38), so their removal books a negative C/N external flow like ``X_dead``.
#: Dissolved species (sugar, ethanol, acids, glycerol, byproducts, SO₂, YAN) stay with the
#: racked-off liquid — a concentration model has no volume change on racking, so touching them
#: would be physically wrong (decision D-36).
#:
#: **Both *Brettanomyces* pools** — viable ``X_brett`` and settled dead ``X_brett_dead`` — are
#: racked off too (decision D-40), the same lees-organism asymmetry as *O. oeni*: Brett colonises
#: the lees, so drawing the wine off them removes the spoilage catalyst and halts volatile-phenol
#: production — the physical twin of the SO₂ kill (:class:`~fermentation.core.kinetics.brett.\
#: BrettDeath`, pt3). ``X_brett`` is carbon-free in pt1 (constant catalyst, unweighted), so its
#: removal books no C/N flow yet; :class:`~fermentation.core.kinetics.brett.BrettGrowth` (pt2)
#: promotes both pools to weighted biomass, and the same ExternalFlow machinery then books their
#: removal like ``X_mlf`` (D-38/D-39).
_LEES_SLOTS = ("X_dead", "debris", "X_mlf", "X_mlf_dead", "X_brett", "X_brett_dead")


def _verb_rack(iv: Intervention, schema: StateSchema, parameters: ParameterSet) -> ScheduledEvent:
    """``rack`` — draw the wine off a fraction of its lees (decisions D-36, D-39).

    Removes ``fraction`` ∈ [0, 1] of each lees-associated pool (:data:`_LEES_SLOTS`: inactivated
    yeast ``X_dead`` and, when autolysis is opted in, the cell-wall ``debris``; plus both *O. oeni*
    pools — viable ``X_mlf`` and settled dead ``X_mlf_dead``), booking the negative jump as an
    :class:`~fermentation.runtime.schedule.ExternalFlow` (the ledger's removal side). Viable
    **yeast** ``X`` and every dissolved species are left untouched — a normal post-AF rack settles
    dead yeast, and a concentration model has no volume change on racking.

    **Racking removes viable *O. oeni* (the D-39 asymmetry, the D-31 lever's physical half).**
    Unlike viable yeast — which ferments in suspension, so a rack leaves it working — *O. oeni*
    carries out MLF on the lees and is drawn off with them. So racking removes the bacteria that
    clear diacetyl, the physical twin of the SO₂ kill (:class:`~fermentation.core.kinetics.\
    malolactic.MalolacticDeath`): the deferred D-31 "rack early ⇒ diacetyl locked in" case. Both
    bacterial pools carry biomass carbon and nitrogen (weighted since D-38), so — like ``X_dead`` —
    their removal shows up as a negative term in the run-wide carbon and nitrogen ledgers
    (``X_dead``/``X_mlf``/``X_mlf_dead`` all carry N; every racked pool carries C).
    """
    _iv_check_keys(iv, frozenset({"fraction"}), "rack")
    fraction = _iv_float(iv, "fraction", "rack")
    if fraction > 1.0:
        raise ValueError(
            f"intervention 'rack' at day {iv.day:g}: fraction must be in [0, 1], got {fraction:g}"
        )
    slices = [schema.slice(name) for name in _LEES_SLOTS if name in schema]
    retained = 1.0 - fraction

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        for sl in slices:
            out[sl] *= retained
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"rack@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_pitch_mlf(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``pitch_mlf`` — inoculate *Oenococcus oeni* mid-run, enabling malolactic fermentation.

    The verb that exercises the driver's third effect (in-place reconfiguration): it both
    **mutates** ``X_mlf`` (adds the bacterial catalyst dose, ``pitch_gpl`` g/L) and
    **reconfigures** the Process set to enable :data:`_MLF_GATED_PROCESSES` — *exactly* the set an
    unpitched compile disables, so a sequential mid-run pitch is symmetric with an initial
    co-inoculation. Since D-38 ``X_mlf`` is real bacterial biomass carrying carbon and nitrogen (no
    longer the inert catalyst of D-23), so the pitch's state jump adds biomass C/N — booked as an
    :class:`~fermentation.runtime.schedule.ExternalFlow`, exactly like the ``add_dap`` dose, so the
    run-wide ledgers still close.

    Because the Processes are enabled only from the breakpoint onward, ``simulate_scheduled``
    min-combines the per-segment tier maps (D-35): the malate/lactate/citrate slots the enabled
    speculative Processes touch report speculative for the *whole* run, not just the back half.

    Honest scope (decisions D-23, D-31): a *post-AF* pitch lands past the Luong ethanol wall
    (~110 g/L; a 24-Brix wine finishes near ~135), so the environmental gate keeps conversion
    near zero — malolactic must be co-inoculated or pitched early to complete. This verb makes
    that timing a *scenario* choice; it does not change the kinetics.
    """
    _iv_check_keys(iv, frozenset({"pitch_gpl"}), "pitch_mlf")
    pitch_gpl = _iv_float(iv, "pitch_gpl", "pitch_mlf")
    if "X_mlf" not in schema:
        raise ValueError(
            f"intervention 'pitch_mlf' at day {iv.day:g} needs an 'X_mlf' slot, but medium "
            f"{schema!r} has none (malolactic fermentation is wine-only, decision D-23)"
        )
    x_mlf_slice = schema.slice("X_mlf")
    gated_names = tuple(p.name for p in _MLF_GATED_PROCESSES)

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[x_mlf_slice] += pitch_gpl
        return out

    def reconfigure(ps: ProcessSet) -> None:
        for name in gated_names:
            if name in ps:
                ps.enable(name)

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"pitch_mlf@{iv.day:g}d",
        mutate=mutate,
        reconfigure=reconfigure,
    )


def _verb_pitch_brett(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``pitch_brett`` — inoculate *Brettanomyces* mid-run, enabling volatile-phenol spoilage.

    The Brett twin of :func:`_verb_pitch_mlf`: it both **mutates** ``X_brett`` (adds the spoilage
    dose, ``pitch_gpl`` g/L) and **reconfigures** the Process set to enable
    :data:`_BRETT_GATED_PROCESSES` — *exactly* the set an unpitched compile disables, so a mid-run
    contamination is symmetric with an initial co-inoculation. The realistic Brett scenario is a
    *post-AF* contamination in the cellar/barrel, which this verb expresses as a scenario choice.

    Because the Processes are enabled only from the breakpoint onward, ``simulate_scheduled``
    min-combines the per-segment tier maps (D-35): the ``vinylphenols``/``ethylphenols`` slots the
    enabled speculative Processes touch report speculative for the whole run. ``X_brett`` is a
    carbon-free catalyst in pt1, so the pitch's state jump adds no biomass carbon/nitrogen (no
    ExternalFlow needed); :class:`~fermentation.core.kinetics.brett.BrettGrowth` (pt2) makes it
    weighted biomass, at which point the pitch books an :class:`~fermentation.runtime.schedule.\
    ExternalFlow` like ``pitch_mlf`` (D-38).
    """
    _iv_check_keys(iv, frozenset({"pitch_gpl"}), "pitch_brett")
    pitch_gpl = _iv_float(iv, "pitch_gpl", "pitch_brett")
    if "X_brett" not in schema:
        raise ValueError(
            f"intervention 'pitch_brett' at day {iv.day:g} needs an 'X_brett' slot, but medium "
            f"{schema!r} has none (Brettanomyces spoilage is wine-only, decision D-40)"
        )
    x_brett_slice = schema.slice("X_brett")
    gated_names = tuple(p.name for p in _BRETT_GATED_PROCESSES)

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[x_brett_slice] += pitch_gpl
        return out

    def reconfigure(ps: ProcessSet) -> None:
        for name in gated_names:
            if name in ps:
                ps.enable(name)

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"pitch_brett@{iv.day:g}d",
        mutate=mutate,
        reconfigure=reconfigure,
    )


def _verb_add_acid(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_acid`` — dose a charge-active organic acid onto its slot (decision D-65, §3.3).

    The general acidulation verb over the D-18 charge-active acids
    (:data:`~fermentation.core.acidbase.ACID_STATE` — tartaric/malic/lactic): ``params`` names the
    ``acid`` and its dose ``gpl``, and the whole mass lands on that acid's state slot. Because
    those slots are wine-only (D-18), this is wine-only by slot presence — a beer scenario has no
    ``tartaric``/``malic``/``lactic`` slot and raises. The dose is the pure acid (it brings its own
    protons, no counter-cation), so it is NOT added to ``cation_charge``; the pH charge balance
    then re-solves the SAME back-anchored strong cation against MORE diprotic/monoprotic anion, so
    pH **drops** and titratable acidity **rises** — the standard acidulation outcome, *emergent*
    from the D-18 keystone rather than scripted (potassium bitartrate, which adds a counter-cation,
    would be a different verb). Each acid carries carbon (tartaric/malic C4, lactic C3, all weighted
    in ``total_carbon``), so the dose is a POSITIVE carbon external flow (the ``add_dap`` +N
    precedent, opposite sign to the copper mercaptan −C removal) and nitrogen-free; the run-wide
    ledger ``final == initial + Σ flows`` still closes to machine precision. Concentration model:
    no volume change on the addition (the shared verb caveat). The acid slot is inert (no Process),
    so no tier moves — pH's tier is already the PLAUSIBLE-floored pKa tier (D-18).
    """
    _iv_check_keys(iv, frozenset({"acid", "gpl"}), "add_acid")
    acid = _iv_str(iv, "acid", "add_acid")
    gpl = _iv_float(iv, "gpl", "add_acid")
    if acid not in acidbase.ACID_STATE:
        raise ValueError(
            f"intervention 'add_acid' at day {iv.day:g}: unknown acid {acid!r}; the charge-active "
            f"acids are {sorted(acidbase.ACID_STATE)} (decision D-18)"
        )
    if acid not in schema:
        raise ValueError(
            f"intervention 'add_acid' at day {iv.day:g} needs a {acid!r} slot, but medium "
            f"{schema!r} has none (the organic-acid slots are wine-only, decision D-18)"
        )
    acid_slice = schema.slice(acid)

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[acid_slice] += gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_acid@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_set_ph(iv: Intervention, schema: StateSchema, parameters: ParameterSet) -> ScheduledEvent:
    """``set_ph`` — re-anchor the strong cation so the beverage sits at ``ph`` from this day on.

    **The gap this closes (decision D-186, closing D-150's open item).** ``initial_ph`` anchors
    t=0 and nothing else: the back-solved cation is computed with ``Byp`` = 0 at pitch, and the
    ferment then drags pH somewhere the scenario never chose — D-150 measured ``initial_ph``
    3.26/3.61 arriving at the oxygen dose as 3.2084/3.5135, a span of 0.3052 where 0.35 was
    asked for. Since every oxidative and SO₂ rate in the aging phase reads pH, an aging study
    at a *stated* pH was simply not writable. This verb makes it writable: put it at the
    ferment/aging boundary (beside ``begin_aging``) and the wine ages at the pH you name.

    **It is a cation-moving verb, which is what makes it physical — this is NOT a pH dial.**
    D-65 anticipated the CATEGORY and deferred it as v1 scope ("potassium bitartrate / K-tartrate
    **additions** — deacidification via a counter-cation, a different, cation-moving verb"). Note
    what that names: the **dose** form. This ships the **target** form, which is the one an aging
    study needs — you know the pH you want to compare at, not the gram-per-litre that reaches it.
    So the axis was anticipated and the API was not. Both directions are real cellar operations
    acting on the same quantity:

    - **raising pH** = deacidification with potassium/calcium carbonate — base in, strong
      cation up. The classic move on an over-acid wine.
    - **lowering pH** = cation-exchange resin, which strips K⁺ and acidifies the wine without
      adding anything. (Lowering pH by *acid addition* is a different operation and already
      has a verb: ``add_acid`` doses tartaric/malic/lactic onto their own slots.)

    So the model change — move ``cation_charge``, touch nothing else — corresponds to a real
    treatment in either direction, and the verb is stated as the adjustment rather than as a
    setter. **Two scope limits, stated rather than left implicit:**

    1. **Cold stabilisation is NOT this verb.** KHT precipitation removes potassium *and*
       tartrate; only the cation is booked here, so the tartaric slot is untouched. A scenario
       wanting the tartrate loss too must also schedule the acid change itself.
    2. **The ``CO2`` slot is deliberately untouched**, even though carbonate deacidification
       evolves gas. That slot is the cumulative *evolved* integral, and the charge balance reads
       ``min(evolved, C_sat(T))`` (decision D-182) — already saturated many times over by the
       time anyone ages a wine, so adding the carbonate's CO₂ moves the dissolved term by
       exactly nothing. Writing it anyway would be bookkeeping theatre.

    **No external flow is booked, unlike ``add_acid`` and ``add_copper``.** A strong cation is
    K⁺/Ca²⁺ charge, not a carbon or nitrogen species: ``cation_charge`` carries weight 0 in
    ``total_carbon`` and is absent from ``total_nitrogen``, so the run-wide identity
    ``final == initial + Σ flows`` closes across this jump with a zero contribution. (The
    scheduler still records the state difference in its ledger, as it does for every mutation;
    it simply weighs nothing.) No tier moves either — pH's tier is already the PLAUSIBLE-floored
    pKa tier (D-18), and this writes a state slot rather than any parameter.

    **Where the validation happens, which is the one place this verb differs from its
    siblings.** Every other verb validates entirely at compile. The reachable pH range depends
    on the *state* — the acid load and accumulated ``Byp`` at the moment of the adjustment — so
    the "is this target reachable" check cannot be made until the event fires. What *is* checked
    at compile: the param keys, the pH bracket, the ``cation_charge`` slot, and that the pKa and
    CO₂-solubility parameters are loaded (the ``begin_aging`` discipline). The state-dependent
    check runs inside the mutation and raises a ``ValueError`` naming the achievable floor.
    That surfaces cleanly rather than as a mid-integration traceback, because
    ``simulate_scheduled`` applies mutations *between* segments — the error comes out of the
    event application with the verb's label, and no partial state is committed.

    **It requires the scenario to have opted into the pH system** (``initial_ph`` present);
    :func:`_compile_interventions` enforces that, and the reason differs per medium, so both
    are stated there rather than one covering for the other.

    **Under an ensemble the anchor is PER-MEMBER, and so is ``initial_ph``'s (decision D-235).**
    This verb's mutation re-solves the cation from the **running** parameter map that
    ``simulate_scheduled`` hands it at the breakpoint, which under an ensemble is the member's own
    draw — so a member sits
    at ``ph`` after this event, not near it. That was not always so, and the history is the point:
    D-186 measured the nominal-only anchor here (24 members spreading **0.1292** pH after the event
    against ``initial_ph``'s **0.1273** at t=0, ratio **1.015**) and declined to repair it *alone*,
    because a lone repair would leave the two anchors disagreeing about what a member's pH means.
    D-233 then repaired the t=0 one and D-234 measured the asymmetry that produced (worst member
    miss **0.07896**, spread **0.13202**, against **2.03e-11** at t=0). D-235 repairs this half, so
    the pair agrees again — which is what D-186's instruction asked for, in the direction it did not
    take. **Both anchors move together: a change to either is a change to the pair.**
    """
    _iv_check_keys(iv, frozenset({"ph"}), "set_ph")
    target_ph = _iv_float(iv, "ph", "set_ph")
    if not 0.0 < target_ph < 14.0:
        raise ValueError(
            f"intervention 'set_ph' at day {iv.day:g}: ph must lie strictly inside (0, 14), got "
            f"{target_ph:g} — solve_ph reports the bracket ends as saturating answers, not roots "
            "(decision D-46), so anchoring to one would not be invertible"
        )
    if "cation_charge" not in schema:
        raise ValueError(
            f"intervention 'set_ph' at day {iv.day:g} needs a 'cation_charge' slot, but medium "
            f"{schema!r} has none — there is nothing to anchor with (decision D-18)"
        )
    resolved = parameters.resolve()
    try:
        # Evaluated for its KeyError, not its value: building the whole pKa lookup here makes a
        # missing acidbase.yaml a scenario error NOW rather than when the event fires (the
        # add_dap/begin_aging discipline). The map the mutation actually uses is the RUNNING one
        # handed to it at the breakpoint, which since D-235 CAN differ from this one — it is the
        # member's draw. What this guard still establishes is that the NAMES are present, which
        # is a property of the ParameterSet the sampler draws from and so holds for every member.
        acidbase.build_pka_map(resolved)
    except KeyError as exc:
        raise ValueError(
            "intervention 'set_ph' needs the pKa parameters but they are missing "
            f"({exc}); include acidbase.yaml in parameter_paths (the default lookup merges "
            "it automatically)."
        ) from exc
    for name in (
        # Read only at mutation time, inside dissolved_co2_molar — so unlike the pKa map above
        # they cannot be guarded by evaluating them, and are named explicitly (decision D-182).
        "H_co2_beverage",
        "T_ref_co2_solubility",
        "vant_hoff_co2_solubility",
    ):
        if name not in parameters:
            raise ValueError(
                f"intervention 'set_ph' at day {iv.day:g} needs {name!r} but it is missing; "
                "include acidbase.yaml in parameter_paths (the default lookup merges it "
                "automatically, decision D-182)."
            )
    cation_slice = schema.slice("cation_charge")
    label = f"set_ph@{iv.day:g}d"

    def mutate(_schema: StateSchema, y: FloatArray, params: Mapping[str, float]) -> FloatArray:
        # `params` — the RUNNING map — and never the compile-time `resolved` (decision D-235).
        # Under an ensemble this is the member's own drawn pKa set, so the member re-anchors to
        # its own numbers exactly as `y0_for_member`'s anchor rule does at t=0. On a nominal run
        # maps are the same object's resolve(), so the jump is bit-identical to what shipped.
        out = y.copy()
        try:
            out[cation_slice] = acidbase.cation_charge_for_ph(y, schema, params, target_ph)
        except acidbase.NitrogenExceedsCationDemandError:
            # A DIFFERENT failure with a different remedy, so it passes through with its own
            # message (decision D-210). Rewriting it as the floor case below blamed the acid load
            # for a state whose nitrogen pool was the cause — reachable once `add_dap` puts a
            # supplement's +1/mol N into the pool, and misleading exactly there.
            raise
        except ValueError:
            # Below the achievable floor: no cation addition reaches it, and removing cation
            # bottoms out at zero. Name that floor — it is what the caller has to work with, and
            # solve_ph is total (D-46) so computing it cannot itself raise.
            floor_state = y.copy()
            floor_state[cation_slice] = 0.0
            # The member's own map again: a floor computed at the nominal pKas would name a
            # number this member's state does not have, which is the reported-value half of the
            # very defect the line above repairs.
            floor = acidbase.ph_of_state(floor_state, schema, params)
            raise ValueError(
                f"intervention {label!r}: target pH {target_ph:g} is below this state's "
                f"intrinsic pH {floor:.4f} — the acid load, plus whatever cation charge the "
                "assimilable-nitrogen pool carries (D-209, and NOT the caller's to remove), "
                "holds it there with zero strong cation, so no deacidification or cation "
                "exchange reaches the target. Raise the target, or lower the acid load."
            ) from None
        return out

    return ScheduledEvent(time_h=days_to_hours(iv.day), label=label, mutate=mutate)


def _verb_add_sugar(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``add_sugar`` — chaptalize: dose sucrose, inverted to fermentable hexose (decision D-65).

    Chaptalization (and beer priming/adjunct) doses SUCROSE by mass (``sugar_gpl``, the commercial
    additive). Sucrose is not fermented as such: yeast invertase hydrolyses it near-instantly into
    glucose + fructose, so the verb inverts it AT THE DOSE (a state mutation, not a kinetic pool —
    invertase is fast vs the ferment) via the exact ``sucrose_inversion_mass_ratio`` (~1.0526; the
    +5.26 % over the sucrose mass is hydrolysis water, the same di-/tri-saccharide mass gain beer's
    wort sugars carry, D-8). The hexose-equivalent lands on the fermentable sugar slot: wine's
    single lumped hexose ``S``, or beer's **glucose** component specifically (found by name, not
    broadcast across the maltose/maltotriose slots) — fructose is lumped as glucose-equivalent,
    exact on carbon and mass since the two are isomers. More fermentable sugar ⇒ a higher finished
    ethanol/ABV once it ferments out (emergent, not imposed). Carbon is conserved through inversion
    (water is carbon-free), so the flow books exactly the sucrose carbon (a POSITIVE carbon external
    flow, the ``add_dap`` precedent) and nitrogen-free; the run-wide ledger still closes to machine
    precision. Concentration model: no volume change (the shared verb caveat).
    """
    _iv_check_keys(iv, frozenset({"sugar_gpl"}), "add_sugar")
    sugar_gpl = _iv_float(iv, "sugar_gpl", "add_sugar")
    if "S" not in schema:
        raise ValueError(
            f"intervention 'add_sugar' at day {iv.day:g} needs an 'S' slot, but medium "
            f"{schema!r} has none"
        )
    try:
        ratio = parameters["sucrose_inversion_mass_ratio"].value
    except KeyError as exc:  # additions.yaml not loaded (caller-supplied parameter_paths)
        raise ValueError(
            "intervention 'add_sugar' needs 'sucrose_inversion_mass_ratio' but it is missing "
            f"({exc}); include additions.yaml in parameter_paths (the default lookup merges "
            "it automatically)."
        ) from None
    hexose_gpl = sugar_gpl * ratio
    # Target the glucose/hexose slot by name (wine's lumped S is treated as glucose), never a
    # broadcast across beer's 3-wide S — fructose from the inversion lumps as glucose-equivalent.
    species = sugar_species(schema)
    glucose_offset = species.index("glucose")
    glucose_index = schema.slice("S").start + glucose_offset

    def mutate(_schema: StateSchema, y: FloatArray, _params: Mapping[str, float]) -> FloatArray:
        out = y.copy()
        out[glucose_index] += hexose_gpl
        return out

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"add_sugar@{iv.day:g}d",
        mutate=mutate,
    )


def _verb_begin_aging(
    iv: Intervention, schema: StateSchema, parameters: ParameterSet
) -> ScheduledEvent:
    """``begin_aging`` — start the post-fermentation aging phase (decisions D-70/D-71, §4.1).

    The aging-axis wiring: it **reconfigures** the Process set to enable the aging Processes
    (:data:`_AGING_GATED_PROCESSES` — :class:`~fermentation.core.kinetics.aging.EsterHydrolysis`,
    :class:`~fermentation.core.kinetics.aging.OxidativeAcetaldehyde`,
    :class:`~fermentation.core.kinetics.aging.PhenolicBrowning` and the wine-only
    :class:`~fermentation.core.kinetics.aging.SulfiteOxidation`) from its ``day`` onward —
    the ``pitch_mlf`` reconfigure pattern MINUS the state mutation (aging inoculates nothing; it
    just switches on the spontaneous chemistry the compile seam left off). All are wired into their
    media but DISABLED at compile (aging is inherently post-ferment — there is no aging at t0), so
    this verb is the *only* way to turn them on; before the breakpoint the run is byte-for-byte the
    pre-aging model and after it the young fruity acetate esters hydrolyse back toward equilibrium
    (fading the ester OAV, raising the fusel OAV, drifting VA/pH up) and — if oxygen has been dosed
    (``add_oxygen``) — dissolved O₂ **browns** the wine/beer (raising the ``A420`` index, D-74) and
    oxidises ethanol to acetaldehyde (the 'sherry'/oxidised note). With no oxygen dosed the
    O₂-driven
    Processes are inert (``o2 = 0``), so ``begin_aging`` alone is purely *reductive* aging —
    byte-for-byte the ester-hydrolysis-only case (D-71/D-74).

    **The aging span is expressed by ``duration_days``** (this is a pure reconfigure with no
    "how long" of its own): put ``begin_aging`` at the ferment/aging boundary day and set
    ``duration_days`` to cover the aging tail. The §7 slow-phase concern (do not integrate years
    at ferment resolution) is answered for free by ``simulate_scheduled``'s segment restart — the
    BDF solver re-initialises its order at the breakpoint and, with the fermentative flux gone
    (``S ≈ 0``), takes large steps across the quiescent aging segment (default ``max_step=∞``); no
    new integration machinery. Every other producer of ``esters``/``fusels``/``Byp`` is
    flux-gated and silent at dryness, so the aging signal is unconfounded (Stance A, D-70).

    Because the Process is enabled only from the breakpoint, ``simulate_scheduled`` min-combines
    the per-segment tier maps (D-35): the speculative ``EsterHydrolysis`` drags ``esters`` /
    ``fusels`` / ``Byp`` to speculative for the WHOLE run, not just the aging back half — a run is
    only as trustworthy as its least-trustworthy segment.

    Takes no params (a pure phase switch). Guards that the aging parameters are loaded (the
    ``add_dap`` discipline) so a caller-supplied ``parameter_paths`` without ``aging.yaml`` fails
    loudly HERE at compile, not as a bare ``KeyError`` mid-integration when the Process reads
    ``k_ester_hydrolysis``.
    """
    _iv_check_keys(iv, frozenset(), "begin_aging")
    # No schema-slot requirement — the aging Processes are medium-agnostic (esters/fusels/Byp/
    # acetaldehyde/o2/A420 exist in both media). Guard the aging params are present (the add_dap/
    # additions.yaml pattern): the reconfigure takes effect at runtime, so an absent aging.yaml
    # would otherwise surface as a KeyError deep in an aging Process's derivatives rather than a
    # clear compile-time scenario error. Guards ALL aging Processes' params (D-70 hydrolysis +
    # D-71 ethanol oxidation + D-72 SO₂ oxidation + D-74 phenolic browning), since begin_aging
    # enables all of them, incl. the D-77 non-oxidative oak axis (k_oak_extraction/E_a). The D-72/
    # D-74 aging.yaml + D-77 oak.yaml params ride in every medium's shared files, so guarding is
    # beer-safe (present in every medium) even though SulfiteOxidation/StreckerDegradation/
    # OakExtraction are wine-only.
    for name in (
        "k_ester_hydrolysis",
        "E_a_ester_hydrolysis",
        "isoamyl_acetate_eq",
        "k_o2_depletion_total",
        "f_ethanol_o2_share",
        "E_a_ethanol_oxidation",
        "y_acetaldehyde_per_o2",
        "k_so2_oxidation",
        "E_a_so2_oxidation",
        "k_browning_phenolic",
        "E_a_browning",
        "y_a420_per_o2",
        # Copper multiplier on the browning rate (D-134): a mean-centered f(Cu) boost, guarded
        # before reading in-Process (like the tannin/anthocyanin phenolic-boost reads), so these
        # can never be missing-when-needed in practice; guarded here for parity.
        "copper_typical",
        "k_copper_multiplier",
        # Initial-burst antioxidant pool (D-133): the finite, fast-reacting non-SO2 sink that
        # produces Ferreira's day-1 O2-consumption spike. Substrate-gated on burst_antioxidant
        # (guarded before reading params in-Process, like Strecker/EllagitanninOxidation), so these
        # can never be missing-when-needed in practice; guarded here for parity.
        "k_burst_oxidation",
        "E_a_burst_oxidation",
        "y_burst_per_o2",
        # Oak extraction (D-77): the non-oxidative barrel/chip axis begin_aging also enables. Only
        # the shared rate + activation energy are guarded here (the 15 toast-specific yields — 4
        # aroma + ellagitannin — are guarded at the add_oak verb, which is the only reader that
        # needs them); k_oak_extraction/E_a_oak_extraction are read by OakExtraction on every
        # enabled aging segment.
        "k_oak_extraction",
        "E_a_oak_extraction",
        # Ellagitannin O₂ scavenging (D-78): the oak-tannin protection sink. Substrate-gated on the
        # ellagitannin pool (guarded before reading params in-Process, like Strecker), so these can
        # never be missing-when-needed (add_oak, the only way to get ellagitannin, already requires
        # oak.yaml); guarded here for parity with the oak-extraction params it ships alongside.
        "k_ellagitannin_oxidation",
        "E_a_ellagitannin_oxidation",
        "y_ellag_per_o2",
    ):
        if name not in parameters:
            raise ValueError(
                f"intervention 'begin_aging' at day {iv.day:g} needs {name!r} but it is missing; "
                "include aging.yaml in parameter_paths (the default lookup merges it "
                "automatically, decisions D-70/D-71)."
            )

    def reconfigure(ps: ProcessSet) -> None:
        for aging_process in _AGING_GATED_PROCESSES:
            if aging_process.name in ps:
                ps.enable(aging_process.name)
        # THE FIRST THING begin_aging HAS EVER DISABLED (decision D-213), and it is the other
        # half of the same switch rather than an exception to it. `WortOxygenUptake` is
        # *pitched* yeast building membrane sterol in the lag phase; past the ferment/aging
        # breakpoint the yeast is settled and dormant, and the O2 chemistry that matters is the
        # aging set enabled just above. Left on, it competes with those sinks for a dosed
        # `add_oxygen` and eats roughly half of it — which is exactly what
        # `test_beer_depletes_its_packaging_oxygen` caught: packaging oxygen must be consumed by
        # the oxidative sinks whose rates are calibrated against measured O2 depletion, not by a
        # fermentation-phase term with an author-estimate constant.
        if WortOxygenUptake.name in ps:
            ps.disable(WortOxygenUptake.name)

    return ScheduledEvent(
        time_h=days_to_hours(iv.day),
        label=f"begin_aging@{iv.day:g}d",
        reconfigure=reconfigure,
    )


#: action verb → compiler turning one :class:`Intervention` into a :class:`ScheduledEvent`.
_INTERVENTION_VERBS: dict[
    str, Callable[[Intervention, StateSchema, ParameterSet], ScheduledEvent]
] = {
    "add_dap": _verb_add_dap,
    "add_so2": _verb_add_so2,
    "add_copper": _verb_add_copper,
    "add_acid": _verb_add_acid,
    "set_ph": _verb_set_ph,
    "add_sugar": _verb_add_sugar,
    "add_oxygen": _verb_add_oxygen,
    "add_ascorbate": _verb_add_ascorbate,
    "add_oak": _verb_add_oak,
    "rack": _verb_rack,
    "pitch_mlf": _verb_pitch_mlf,
    "pitch_brett": _verb_pitch_brett,
    "begin_aging": _verb_begin_aging,
}

#: action verb → compiler that additionally needs the whole :class:`Scenario` (decision D-187).
#:
#: A SECOND table rather than a fourth parameter on every verb above, because the two kinds are
#: genuinely different and the split says which is which: a verb in :data:`_INTERVENTION_VERBS`
#: is a function of its own ``Intervention`` alone — it can be read, tested and reasoned about
#: without knowing anything else about the run — while a verb here draws its magnitude from a
#: scenario-level choice made elsewhere. ``seal_bottle`` is the first: its dose is whatever
#: ``scenario.closure`` names. Widening all twelve signatures to carry a ``scenario`` almost none
#: of them may look at would have made that distinction invisible instead of explicit.
#:
#: :func:`_compile_interventions` dispatches on this table first and both are searched before an
#: action is called unknown, so the two are one namespace with one error message.
_SCENARIO_INTERVENTION_VERBS: dict[
    str, Callable[[Intervention, StateSchema, ParameterSet, Scenario], ScheduledEvent]
] = {
    "seal_bottle": _verb_seal_bottle,
}

#: Verbs whose dose is resolved from ``scenario.closure`` and which therefore require one
#: (decision D-187). Separate from :data:`_SCENARIO_INTERVENTION_VERBS` because *needing the
#: scenario* and *needing a closure* are different claims — a later scenario-level verb reading,
#: say, ``batch_volume_liters`` would join the first set and not this one.
_CLOSURE_SOURCED_VERBS: frozenset[str] = frozenset({"seal_bottle"})

#: Verbs that may only be scheduled when the scenario opted into the pH system by giving
#: ``initial_ph`` (decision D-186, riding D-179's gate). The check lives in
#: :func:`_compile_interventions` rather than in the verb because a verb is handed only its own
#: ``Intervention``, the schema and the parameters — never ``scenario.initial``.
#:
#: **The reason it is needed differs by medium, and neither reason covers for the other.**
#: For BEER the failure is structural: without ``initial_ph`` every acid slot is 0 (D-179's
#: opt-in), so anchoring would write a strong cation into an *empty* acid load — a charge
#: balance with a counter-cation and nothing to counter, which solves to the top of the bracket
#: and is not a beverage. For WINE the balance is not empty at all — ``tartaric``/``malic`` are
#: seeded from their own scenario keys regardless, and D-182 measured an un-anchored wine at pH
#: 2.92 off ``Byp`` alone — so the objection is epistemic rather than structural: re-anchoring a
#: wine whose pH the scenario never supplied would let one intervention manufacture the pH
#: information the whole D-18 inverse-anchoring design says must be an input.
_PH_SYSTEM_VERBS: frozenset[str] = frozenset({"set_ph"})


def _compile_interventions(
    scenario: Scenario, schema: StateSchema, parameters: ParameterSet, t_end_h: float
) -> tuple[ScheduledEvent, ...]:
    """Compile ``scenario.interventions`` into timed :class:`ScheduledEvent`\\ s (decision D-36).

    Each verb is looked up in :data:`_INTERVENTION_VERBS`, then in
    :data:`_SCENARIO_INTERVENTION_VERBS` for the few whose *magnitude* comes from a scenario-level
    field (D-187's ``seal_bottle``, whose dose is whatever ``scenario.closure`` names); an action in
    neither raises loudly, naming both (the ``_ALLOWED_KEYS`` discipline). An intervention at or
    after the run duration is rejected here with a scenario-level message rather than deferred to
    ``simulate_scheduled``'s window check, so the error names the scenario and the verb.

    **Cross-cutting gates live here, not in the verbs**, because a verb is handed only its own
    ``Intervention``, the schema and the parameters: D-186's pH-system opt-in reads
    ``scenario.initial``, and D-187's three closure gates read ``scenario.closure``,
    ``scenario.medium`` and the *other* interventions (a ``seal_bottle`` may not precede
    ``begin_aging``).
    """
    events: list[ScheduledEvent] = []
    aging_days = [iv.day for iv in scenario.interventions if iv.action == "begin_aging"]
    for iv in scenario.interventions:
        verb = _INTERVENTION_VERBS.get(iv.action)
        scenario_verb = _SCENARIO_INTERVENTION_VERBS.get(iv.action)
        if verb is None and scenario_verb is None:
            raise ValueError(
                f"scenario {scenario.name!r}: unknown intervention action {iv.action!r}; "
                f"known verbs: "
                f"{sorted({*_INTERVENTION_VERBS, *_SCENARIO_INTERVENTION_VERBS})}"
            )
        if iv.action in _CLOSURE_SOURCED_VERBS:
            # Gate 1 — the dose IS the closure, so there is nothing to dose without one. Silently
            # skipping (or dosing zero) would let a scenario ask for a sourced bottling charge and
            # receive none, which is the failure `_closure_otr` refuses for the same reason.
            if scenario.closure is None:
                raise ValueError(
                    f"scenario {scenario.name!r}: intervention {iv.action!r} at day {iv.day:g} "
                    "takes its dose from 'closure', which this scenario does not name; give one "
                    f"of {', '.join(_CLOSURES)} (decision D-187). To dose an oxygen charge you "
                    "supply yourself, use 'add_oxygen'"
                )
            # THERE IS DELIBERATELY NO WINE-ONLY GATE HERE, and its absence is measured rather
            # than assumed. One was written, and it could not fire: this verb requires a
            # `closure` (above), and `compile_scenario` rejects a `closure` on a medium without
            # the `closure_otr` slot BEFORE it compiles any intervention (D-136). So a beer
            # scenario reaching this line is already impossible — with a closure it died at the
            # compile seam, without one it dies on the gate above. Adding a third message would
            # have documented a path no scenario can take (the D-186 lesson: falsify a guard
            # against the case it names, and if it cannot fire, do not ship it).
            #
            # Gate 2 — the anti-double-count in `bottling_burst_*` subtracts 30 days of steady
            # permeation, which only runs once `begin_aging` has enabled ClosureOxygenIngress.
            # Sealing earlier would net off a flux the run never paid. Bottling IS the start of
            # bottle aging, so this rejects nothing an author would want to write.
            if not aging_days or min(aging_days) > iv.day:
                when = (
                    f"the earliest 'begin_aging' is day {min(aging_days):g}"
                    if aging_days
                    else "this scenario never calls 'begin_aging'"
                )
                raise ValueError(
                    f"scenario {scenario.name!r}: intervention {iv.action!r} at day {iv.day:g} "
                    f"must be at or after 'begin_aging' ({when}); the sourced charge is net of "
                    "the steady ingress over the same first month, which is not running until "
                    "aging begins (decision D-187)"
                )
        if iv.action in _PH_SYSTEM_VERBS and "initial_ph" not in scenario.initial:
            raise ValueError(
                f"scenario {scenario.name!r}: intervention {iv.action!r} at day {iv.day:g} "
                "needs the pH system, which is opted into by giving 'initial_ph' in "
                "scenario.initial (decisions D-179/D-186); without it there is no anchored "
                "pH for this verb to re-anchor"
            )
        if days_to_hours(iv.day) >= t_end_h:
            raise ValueError(
                f"scenario {scenario.name!r}: intervention {iv.action!r} at day {iv.day:g} is "
                f"at or beyond the run duration ({scenario.duration_days:g} d); interventions "
                "must fall within the run"
            )
        if scenario_verb is not None:
            events.append(scenario_verb(iv, schema, parameters, scenario))
        else:
            assert verb is not None  # one of the two tables matched, checked above
            events.append(verb(iv, schema, parameters))
    return tuple(events)


def compile_scenario(
    scenario: Scenario,
    *,
    parameter_paths: Sequence[str | Path] | None = None,
    data_dir: str | Path | None = None,
    strict: bool = False,
    oxidative: str = "direct",
) -> CompiledScenario:
    """Compile a declarative scenario into an integrable :class:`CompiledScenario`.

    Industry units in ``scenario.initial`` are converted to canonical units here
    and nowhere else. ``parameter_paths`` overrides the default lookup of
    ``<medium>_<strain>.yaml`` under ``data_dir`` (or the packaged data dir);
    ``strict=True`` enables the Process ``touches`` contract on the returned set.

    ``oxidative`` selects which oxidative alternative is wired (decisions D-141/D-147):
    ``"direct"`` — the default — is the six calibrated pre-cascade sinks that each draw straight on
    ``o2``; ``"cascade"`` routes them all behind one Fe(II)+O2 activation node; ``"direct_burst"``
    is the direct six plus D-133's ``AntioxidantBurstOxidation``. It is passed through to
    :func:`~fermentation.core.media.get_medium`; see there for why neither non-default alternative
    is default. Under ``"direct_burst"`` the ``burst_antioxidant`` slot is seeded from the sourced
    ``burst_antioxidant_initial``; under the other two it is seeded 0.0 and dosing
    ``burst_antioxidant_gpl`` raises (see :func:`_resolve_burst_antioxidant_seed`).

    Raises ``KeyError`` for an unknown medium, ``ValueError`` for an invalid
    initial composition or missing temperature, and ``FileNotFoundError`` when the
    medium/strain has no parameter file yet.
    """
    medium = get_medium(scenario.medium, oxidative=oxidative)
    _validate_initial_keys(scenario)

    builder = _INITIAL_BUILDERS.get(scenario.medium)
    if builder is None:
        raise ValueError(f"no initial-composition builder for medium {scenario.medium!r}")

    temperature_k = _initial_temperature_kelvin(scenario)

    # Parameters are loaded *before* y0 because the wine initial sugar applies a
    # sourced must_fermentable_fraction (decision D-16), mirroring how the
    # nitrogen-dependent yield (D-14) is also resolved at this boundary.
    parameters = _load_parameters(scenario, parameter_paths, data_dir)
    parameters = _apply_nitrogen_dependent_yield(scenario, parameters)

    y0 = medium.schema.pack(builder(scenario.initial, temperature_k, parameters))
    process_set = medium.build_process_set(strict=strict)

    # Hop bittering (decision D-64): the boil isomerization is a wort-side calc, run once here and
    # wired into ``iso_alpha`` at t=0 (like the measured ``initial_ph`` back-solve, D-18). When
    # hops are scheduled, seed the state; when they are NOT, DISABLE the fermentation loss so the
    # empty ``iso_alpha`` slot keeps its VALIDATED tier (an enabled speculative Process touching it
    # would drag ``tier_of`` even with a zero contribution) and no flux is paid — the MLF/Brett
    # isolability pattern. Guard that hops are only given for a medium that HAS a bitterness model
    # (beer): a wine scenario with hops is a user error, not a silently-ignored field.
    if scenario.hops:
        if "iso_alpha" not in medium.schema:
            raise ValueError(
                f"scenario has 'hops' but medium {scenario.medium!r} has no bitterness model "
                "(no 'iso_alpha' state); hop bittering is beer-only (decision D-64)"
            )
        y0[medium.schema.slice("iso_alpha")] = _iso_alpha_at_pitch(scenario, parameters)
    elif IsoAlphaAcidLoss.name in process_set:
        process_set.disable(IsoAlphaAcidLoss.name)

    _resolve_burst_antioxidant_seed(scenario, medium, process_set, y0)

    # Closure oxygen ingress (decision D-136): seed the `closure_otr` state slot from the named
    # closure's sourced OTR. The `iso_alpha` pattern above exactly — a scenario-level choice
    # resolved once here at the compile boundary rather than read per-RHS-step, because a closure
    # does not change during a run. It rides in STATE rather than as a parameter the Process reads
    # because the scenario layer has no parameter-override seam, so a per-run choice has nowhere
    # else to live (the `copper`/`bound_h2s` precedent, D-134/D-135).
    #
    # Absent `closure` ⇒ the slot keeps its 0.0 VarSpec default ⇒ ClosureOxygenIngress contributes
    # byte-for-byte zero and the whole pre-D-136 aging axis is bit-identical. UNLIKE the D-134
    # copper case, 0 is the correct neutral here: it is an additive source, and a zero-ingress
    # bottle is a real measured case (Lopes et al. 2007 found only their flame-sealed control fully
    # air-tight), not an unphysical setting. Nothing reads the slot until `begin_aging` enables the
    # Process, so a scenario that names a closure but never ages is inert rather than wrong.
    if scenario.closure is not None:
        if "closure_otr" not in medium.schema:
            raise ValueError(
                f"scenario has 'closure' but medium {scenario.medium!r} has no closure-ingress "
                "model (no 'closure_otr' state); closure oxygen ingress is wine-only "
                "(decision D-136)"
            )
        y0[medium.schema.slice("closure_otr")] = _closure_otr(scenario.closure, parameters)

    # MLF isolability (decisions D-23, D-31): the malolactic Processes are wired into the wine
    # medium but contribute nothing until Oenococcus oeni is pitched. When it is not, DISABLE
    # them all so (a) the inert ``malic``/``lactic``/``citrate`` slots keep their VALIDATED
    # tier — an *enabled* Process that touches them drops them to speculative even with a zero
    # contribution, since ``tier_of`` counts enabled, not nonzero, Processes — and (b) no
    # per-RHS pH ``brentq`` solve is paid on an undosed run. When pitched, MalolacticConversion
    # is the first RHS consumer of the D-18 pH solver / D-22 molecular-SO₂ readout, and the two
    # D-31 Processes co-metabolise citrate into diacetyl and reduce it on the lees.
    # An initial ``mlf_pitch_gpl`` co-inoculates at t0; a mid-run ``pitch_mlf`` intervention
    # (decision D-36) instead leaves this 0 and re-enables the same _MLF_GATED_PROCESSES at its
    # breakpoint. Either way, an unpitched compile disables them here.
    mlf_pitch_gpl = float(scenario.initial.get("mlf_pitch_gpl", 0.0) or 0.0)
    if mlf_pitch_gpl <= 0.0:
        for mlf_process in _MLF_GATED_PROCESSES:
            if mlf_process.name in process_set:
                process_set.disable(mlf_process.name)

    # Brett isolability (decision D-40, D-55): the volatile-phenol Processes are wired into the
    # wine medium but contribute nothing until Brettanomyces is pitched. When it is not, DISABLE
    # them so (a) the inert ``hydroxycinnamics``/``vinylphenols``/``ethylphenols`` slots AND their
    # D-55 ferulic-branch counterparts (``ferulic_acid``/``vinylguaiacols``/``ethylguaiacols``)
    # keep their VALIDATED tier (an *enabled* zero-contribution Process still drags them to
    # speculative via ``tier_of``) and (b) no per-RHS pH ``brentq`` is paid on an unpitched run. An
    # initial ``brett_pitch_gpl`` co-inoculates at t0; a mid-run ``pitch_brett`` intervention
    # instead leaves this 0 and re-enables the same _BRETT_GATED_PROCESSES at its breakpoint (the
    # MLF pattern).
    brett_pitch_gpl = float(scenario.initial.get("brett_pitch_gpl", 0.0) or 0.0)
    if brett_pitch_gpl <= 0.0:
        for brett_process in _BRETT_GATED_PROCESSES:
            if brett_process.name in process_set:
                process_set.disable(brett_process.name)

    # POF+ yeast decarboxylase isolability (decision D-40 pt4): YeastPOFDecarboxylation is the yeast
    # cinnamate decarboxylase filling ``vinylphenols`` from must ``hydroxycinnamics`` during AF.
    # POF+ is a BINARY STRAIN TRAIT, gated on its own opt-in ``pof_positive`` and WHOLLY INDEPENDENT
    # of the Brett pitch (a POF+ ferment need not have Brett; a POF-negative wine must make no
    # vinylphenol). Absent/<=0 => DISABLE it so (a) the empty ``vinylphenols`` slot keeps its
    # VALIDATED tier - an *enabled* zero-contribution Process would drag it to speculative via
    # ``tier_of`` - and (b) no wasted flux/Monod recompute is paid on a POF- run, which is then
    # byte-for-byte the validated core (the Brett-unpitched pattern). Opted in => the Process runs;
    # ``vinylphenols`` honestly reports speculative, while ``ethylphenols`` stays VALIDATED at 0
    # unless Brett (the only reductase) is also present - the emergent stranding.
    pof_positive = float(scenario.initial.get("pof_positive", 0.0) or 0.0)
    if pof_positive <= 0.0 and YeastPOFDecarboxylation.name in process_set:
        process_set.disable(YeastPOFDecarboxylation.name)

    # Beer's organic-acid producer (decision D-180) rides the SAME opt-in gate as the acid
    # slots it fills — ``initial_ph``, D-179's gate for beer's whole pH system. Absent ⇒
    # DISABLE, and this one is a correctness gate rather than the usual tier/cost argument.
    #
    # Without ``initial_ph`` every acid slot AND the counter-cation are 0, which is an EMPTY
    # charge balance — the state ``acidbase.charge_balance_is_populated`` reports as "this
    # beverage does not claim a pH", so that the aging rate laws hold their pH factor at 1
    # instead of aging against pure water's 7.0 (the defect D-179 shipped an amendment for).
    # An ENABLED producer would fill those slots from sugar carbon as the ferment ran, so a
    # beer that supplied no pH would ACQUIRE a populated charge balance mid-run and start
    # aging against a pH nobody gave it. That is the same defect wearing a producer's hat, and
    # it would arrive silently — the run stays green, the number just becomes fiction.
    # Disabled, an un-anchored beer is byte-for-byte the pre-D-179 beer, and the empty acid
    # slots keep their tier (``tier_of`` counts enabled, not nonzero, Processes).
    #
    # ``WortAcidRemoval`` (D-181) rides the identical gate, and its argument is the mirror
    # image: absent ``initial_ph`` its three slots are 0, already below their floors, so the
    # rate law is a no-op anyway — but an ENABLED Process holds those empty slots' tier below
    # VALIDATED (``tier_of`` counts enabled Processes, not nonzero ones), so leaving it on
    # would change what an un-anchored beer REPORTS about acids it does not carry.
    #
    # ``AceticAcidOverflow`` (D-183) rides it too, and for the FIRST of those two reasons rather
    # than the second: it is a genuine producer, so left enabled it would fill the ``acetic``
    # slot from sugar on a beer that supplied no pH — the correctness case, not the tier one.
    # It is here rather than beside the growth Processes it reads because what gates it is the
    # acid it makes, not the growth it tracks.
    for _acid_process in (OrganicAcidExcretion, AceticAcidOverflow, WortAcidRemoval):
        if _acid_process.name in process_set and "initial_ph" not in scenario.initial:
            process_set.disable(_acid_process.name)

    # Residual-nitrogen floor (decision D-30): the biomass carrying-capacity cap is a
    # deliberate DEPARTURE from the validated Coleman anchor (which caps nothing and strips
    # YAN to zero at every dose), so it ships OPT-IN. Absent ``carrying_capacity_gpl`` ⇒
    # DISABLE the modifier so (a) growth's whole contribution is unscaled (factor 1) and the
    # run is byte-for-byte the validated core, and (b) the enabled-but-inert modifier does not
    # drag growth's X/S/N outputs from PLAUSIBLE to speculative (``tier_of`` counts enabled,
    # not nonzero, modifiers — the exact MLF *tier* isolability argument above). Opted in ⇒
    # enable it and override the reference cap with the scenario's value so demonstrations can
    # sweep K; growth's outputs then honestly report speculative.
    if BiomassCarryingCapacity.name in process_set:
        raw_cap = scenario.initial.get("carrying_capacity_gpl")
        # A negative cap is a typo, not an intent — raise loudly like every other initial key
        # (the _nonneg gate), rather than silently disabling. Absent or 0 ⇒ opt out (disable).
        cap_gpl = _nonneg(float(raw_cap), "carrying_capacity_gpl") if raw_cap is not None else 0.0
        if cap_gpl <= 0.0:
            process_set.disable(BiomassCarryingCapacity.name)
        else:
            parameters = _override_carrying_capacity(parameters, cap_gpl)

    # Amino-acid ledger isolability (decisions D-32, D-33): the AminoAcidAssimilation swap and
    # the FuselAminoAcidReroute are wired into the wine medium but contribute nothing until amino
    # acids are dosed. When they are not, DISABLE them so (a) the empty ``amino_acids`` slot keeps
    # its VALIDATED tier — an *enabled* speculative Process touching ``S``/``N`` would drag those
    # outputs down even with a zero contribution (``tier_of`` counts enabled, not nonzero,
    # Processes) — and (b) no rate recompute is paid on an undosed run. Dosed, the swap funds a
    # fraction of biomass from amino acids (refunding sugar carbon and ammonium N, scaled alongside
    # growth by the wine Arrhenius/carrying modifiers so it never creates sugar, D-32), and the
    # re-route sources a fraction of Ehrlich fusel carbon from amino acids, deaminating the nitrogen
    # to ammonium (D-33). The re-route is paired with FuselAlcoholsEhrlich (it refunds sugar that
    # producer drew), which is always enabled in the wine set, so it is safe to enable here.
    amino_acids_gpl = float(scenario.initial.get("amino_acids_gpl", 0.0) or 0.0)
    if amino_acids_gpl <= 0.0:
        for aa_process in (
            AminoAcidAssimilation,
            FuselAminoAcidReroute,
            PrecursorNonEhrlichFates,
        ):
            if aa_process.name in process_set:
                process_set.disable(aa_process.name)

    # MLF-growth isolability (decision D-38, the deferred growth beat). MalolacticGrowth builds
    # bacterial biomass X_mlf from the amino-acid pool, so the FEATURE it represents (amino-acid-fed
    # bacterial growth) is keyed on amino acids being dosed — the SAME gate as the swap/re-route
    # above. Disable it when amino_acids_gpl ≤ 0 so (a) the empty amino_acids slot / X_mlf keep
    # their tier — an enabled speculative Process touching them would drag tier_of even with a zero
    # contribution — and (b) no rate recompute is paid undosed. This alone prevents the D-23/D-31
    # tier regression: those tests pitch O. oeni but dose NO amino acids, so growth stays disabled.
    # NOT additionally gated on the pitch: "bacteria present" is runtime state the Process's own
    # ``X_mlf ≤ 0`` guard handles (zero until a co-inoculation dose or a mid-run pitch_mlf mutation
    # adds X_mlf), and — mirroring how MalolacticConversion trusts its ethanol gate rather than a
    # compile rule — whether post-pitch bacteria then GROW is left to the emergent environmental
    # gate (g_EtOH·γ(T)·…). So co-inoculation-dominance is emergent, not hard-coded: a high-ABV must
    # arrests growth via the ethanol wall, a normal-ABV sequential MLF can still grow (D-38).
    if MalolacticGrowth.name in process_set and amino_acids_gpl <= 0.0:
        process_set.disable(MalolacticGrowth.name)

    # Brett-growth isolability (decision D-40 pt2). BrettGrowth builds X_brett from the amino-acid
    # pool (drawing the carbon shortfall from ethanol, so Brett grows in a dry wine), so — exactly
    # like MalolacticGrowth — it is keyed on amino acids being dosed, NOT on the Brett pitch (the
    # Process's own ``X_brett ≤ 0`` guard handles "Brett present"; whether it then grows is left to
    # the emergent SO₂/temperature gate). Disable it when amino_acids_gpl ≤ 0 so the empty
    # amino_acids / X_brett slots keep their tier and no rate recompute is paid undosed. This keeps
    # every pitched-but-not-aa-dosed Brett run (e.g. the pt1 headline) with growth disabled.
    if BrettGrowth.name in process_set and amino_acids_gpl <= 0.0:
        process_set.disable(BrettGrowth.name)

    # Autolytic source (decisions D-34, D-44, D-45): YeastAutolysis refills the amino-acid pool from
    # dead biomass (X_dead) post-AF — the second MLF-with-growth prerequisite — AutolyticHydrogen-
    # Sulfide (D-44) feeds the shared h2s pool the sulfide those self-digesting cells release, and
    # AutolyticMercaptan (D-45) feeds the mercaptans pool their thiols (drawing carbon from
    # amino_acids, deaminating N) — all three yields on the SAME autolysis flux. They *consume/read*
    # core state gated on autolysis, so they ship OPT-IN TOGETHER: absent ``autolysis_rate_per_h`` ⇒
    # DISABLE all three so (a) the X_dead/amino_acids/debris/h2s/mercaptans/N columns are untouched
    # and the run is byte-for-byte the validated core, and (b) the inert
    # Processes do not drag those outputs to speculative (``tier_of`` counts enabled, not nonzero,
    # Processes — the MLF/carrying *tier* isolability argument). Opted in ⇒ enable them and override
    # k_autolysis with the scenario's rate so demonstrations can sweep the sur-lie timescale — the
    # override drives BOTH the peptide refill and the sulfide yield (they read one k_autolysis).
    if YeastAutolysis.name in process_set:
        raw_rate = scenario.initial.get("autolysis_rate_per_h")
        rate_per_h = (
            _nonneg(float(raw_rate), "autolysis_rate_per_h") if raw_rate is not None else 0.0
        )
        if rate_per_h <= 0.0:
            for autolysis_process in (
                YeastAutolysis,
                AutolyticHydrogenSulfide,
                AutolyticMercaptan,
            ):
                if autolysis_process.name in process_set:
                    process_set.disable(autolysis_process.name)
        else:
            parameters = _override_autolysis_rate(parameters, rate_per_h)

    # Aging isolability (decisions D-70/D-71): the aging Processes (EsterHydrolysis +
    # OxidativeAcetaldehyde) are wired into both media but aging is INHERENTLY post-ferment — there
    # is no aging at t0 — so unlike the pitch-gated MLF/Brett tuples (which can co-inoculate at t0)
    # they are DISABLED unconditionally here. The ONLY way to turn them on is a ``begin_aging``
    # intervention, which re-enables exactly this tuple at its breakpoint (the pitch_mlf reconfigure
    # pattern). Disabled ⇒ skipped by ``active``/``tier_of``/strict, so an un-aged scenario is
    # byte-for-byte the pre-aging core and the esters/fusels/Byp/acetaldehyde/o2 pools keep their
    # pre-aging tier (prime directive #3). aging.yaml's params ride in every ParameterSet
    # (shared_files) but are read by nothing until a begin_aging event fires.
    for aging_process in _AGING_GATED_PROCESSES:
        if aging_process.name in process_set:
            process_set.disable(aging_process.name)

    t_span_h = (0.0, days_to_hours(scenario.duration_days))

    # Temperature schedule (decision D-35): compile the piecewise-linear ramp into the
    # TemperatureRamp's initial slope + slope-change events. Only when it actually ramps do
    # we mint the provenance-backed rate parameter and emit events — a flat/single-knot
    # schedule leaves ``temperature_ramp_rate`` absent, so the always-enabled TemperatureRamp
    # reads its 0.0 default and an isothermal run is byte-for-byte the pre-ramp core.
    initial_slope, ramp_events = _temperature_ramp_schedule(scenario, t_span_h[1])
    if initial_slope != 0.0 or ramp_events:
        parameters = _inject_temperature_ramp_rate(parameters, initial_slope)

    # Discrete winemaking interventions (decision D-36): compile the declarative timeline of
    # verbs into timed events and merge with the ramp's slope-change events into the single
    # ``events`` tuple ``simulate_scheduled`` sorts by time. Empty ⇒ the temp-only path is
    # unchanged (byte-for-byte core when there is no ramp either).
    intervention_events = _compile_interventions(scenario, medium.schema, parameters, t_span_h[1])
    events = (*ramp_events, *intervention_events)

    return CompiledScenario(
        scenario=scenario,
        schema=medium.schema,
        y0=y0,
        process_set=process_set,
        parameters=parameters,
        t_span_h=t_span_h,
        events=events,
    )
