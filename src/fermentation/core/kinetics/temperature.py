"""Temperature as a driven state variable — the schedule-following ramp (decision D-35).

Every kinetic rate in the model reads temperature from the state vector ``T`` (Kelvin),
and the Arrhenius modifier is already written to consume a *time-varying* ``T`` (see
:mod:`fermentation.core.kinetics.arrhenius`, "isothermal in M1 … already correct for the
non-isothermal temperature dynamics of a later tier"). Until now nothing *drove* ``T``:
with no Process touching it, ``dT/dt = 0`` and a run sat at its pitch temperature. This
Process closes that gap.

**A constant-slope forcing, segmented at slope changes.** A cellar temperature schedule
is a piecewise-*linear* ramp between set-points, not a staircase (real temperature moves
over hours; it does not teleport). Between two schedule knots the slope ``dT/dt`` is a
single constant, so this Process simply writes that constant into the ``T`` derivative:

    dT/dt = temperature_ramp_rate   (K/h)

The per-segment slope is supplied by the runtime event loop
(:func:`fermentation.runtime.simulate_scheduled`), which restarts the integrator at each
knot where the slope *changes* (compiled from the scenario's ``temperature_schedule`` at
the D-3 unit boundary). Because ``dT/dt`` is then constant *within* every segment, the
implicit BDF solver integrates ``T`` **exactly** — a first-degree polynomial is integrated
to round-off by BDF of any order — so ``T`` as an integrated state carries no meaningful
numerical error versus an analytic line, and every other Process sees the true
instantaneous ramped temperature at each solver substep.

**Isothermal default, byte-for-byte.** The rate is read with a ``0.0`` default, so with
no schedule ramp (or when a caller has not injected the parameter) this Process
contributes exactly ``0.0`` to ``dT/dt`` — adding nothing to an isothermal run. It is
wired into both media and always enabled; an isothermal run is therefore byte-for-byte
the pre-ramp core (``0.0 + 0.0 == 0.0`` on the ``T`` slot) and ``T`` stays VALIDATED.

**Tier / provenance.** ``temperature_ramp_rate`` is a *scenario-exact* forcing — a
controlled set-point schedule, known without uncertainty — not an empirical kinetic
parameter. So the Process is VALIDATED and, deliberately, declares **no** ``reads``: the
``reads`` mechanism exists for D-1 credibility propagation (a speculative parameter
dragging an output's tier down), and a value that is exact by construction borrows no
credibility. Declaring it would only force ``temperature_ramp_rate`` into every
``param_tiers`` map with no honesty benefit (and it would be pointlessly swept by the
stochastic ensemble). The scenario compile boundary still records it as a provenance-
backed :class:`~fermentation.parameters.schema.Parameter` (source: the temperature
schedule) when a ramp is present, per prime directive #2.
"""

from __future__ import annotations

from collections.abc import Mapping

from fermentation.core.process import Process
from fermentation.core.state import FloatArray, StateSchema
from fermentation.core.tiers import Tier

#: Parameter name carrying the current segment's temperature slope, in K/h. Supplied by
#: the scenario compile boundary / event loop; absent ⇒ isothermal (the ``0.0`` default).
RAMP_RATE = "temperature_ramp_rate"


class TemperatureRamp(Process):
    """Drive ``T`` along the scenario's piecewise-linear temperature schedule.

    ``dT/dt = params[temperature_ramp_rate]`` (K/h), constant within each event-loop
    segment. Touches only ``T``; reads the rate with a ``0.0`` isothermal default so an
    un-ramped run is a no-op. VALIDATED (a set-point schedule is an exact input); see the
    module docstring for why it declares no ``reads``.
    """

    name = "temperature_ramp"
    tier = Tier.VALIDATED
    touches = ("T",)

    def derivatives(
        self, t: float, y: FloatArray, schema: StateSchema, params: Mapping[str, float]
    ) -> FloatArray:
        d = schema.zeros()
        # Absent ⇒ 0.0: the isothermal no-op that keeps an un-ramped run byte-for-byte the
        # pre-ramp core (and shields bare-built Process sets whose hand-made param maps do
        # not carry the scenario-derived rate).
        d[schema.slice("T")] = params.get(RAMP_RATE, 0.0)
        return d
