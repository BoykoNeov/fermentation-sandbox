"""Three different things get called "fidelity". They are kept apart here on purpose.

Rolling them into one quality slider would be the kind of quiet mixing the project's own
rules forbid, because only the third of them has anything to do with being right:

1. **How carefully the maths is done** — the solver and its tolerances. Turn this up and the
   answer settles onto a fixed value. It does not become more like a real fermentation; it
   becomes more like this model. The only question it answers is whether the calculation has
   stopped moving, and :func:`tightened` plus one re-run is how you find out.
2. **How many points are kept** — how much of the run is stored for charting. This changes
   nothing about the answer. It is *not* cosmetic for the uncertainty band, though: the
   re-runs can only be combined if they all land on the same points, so the same count has
   to be used for both.
3. **What is being modelled** — which chemistry is switched on, and whether the model's
   internal contract is enforced. This is the only setting here that can move the result
   towards or away from reality.

``RK45`` is deliberately absent from every preset. Fermentation always mixes very fast
chemistry with very slow chemistry, and the solvers built for that are the implicit ones;
an explicit method will crawl or fail outright. It is reachable under "custom" only, and
labelled as the diagnostic it is. Capping the step size is custom-only for the same kind of
reason: the solver already chooses its own steps, so a cap almost never changes the answer
and always costs time.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

#: Precision preset → (rtol, atol). ``standard`` is the engine's own default, so a
#: default-preset run in the app is byte-for-byte a default-argument run in the library.
PRECISION_PRESETS: dict[str, tuple[float, float]] = {
    "draft": (1e-4, 1e-7),
    "standard": (1e-6, 1e-9),
    "high": (1e-8, 1e-11),
}

#: Loosest → tightest. :func:`tightened` walks this, and the convergence check compares a
#: run against its next entry.
PRECISION_ORDER: tuple[str, ...] = ("draft", "standard", "high")

#: One-line explanation per preset, shown next to the control.
PRECISION_BLURB: dict[str, str] = {
    "draft": "Quick and rough. Good for shaping a batch; do not quote a number off it.",
    "standard": "What the engine uses by default and what every test in the project runs at.",
    "high": (
        "A hundred times stricter. Use it to confirm the standard answer rather than to get "
        "a better one — if the two disagree, the standard one was not settled."
    ),
    "custom": (
        "Your own settings. RK45 is included only as a diagnostic: it is not built for a "
        "system where fast and slow chemistry run side by side, which this always is."
    ),
}

#: Implicit methods appropriate to a stiff system, plus the non-stiff diagnostic last.
METHODS: tuple[str, ...] = ("BDF", "Radau", "LSODA", "RK45")

#: Stored-point counts offered for the output grid.
POINT_CHOICES: tuple[int, ...] = (100, 200, 500, 1000, 2000)

#: ``oxidative`` set → what selecting it means. Keys must match
#: ``fermentation.core.media.get_medium``'s accepted values.
OXIDATIVE_BLURB: dict[str, str] = {
    "direct": (
        "Six ways oxygen gets used up, each drawing on the dissolved oxygen directly. This "
        "is the default and the only version the project's benchmarks were fitted against."
    ),
    "cascade": (
        "The same six, but with oxygen having to be activated by iron first before any of "
        "them can run. Fully built and tested, and closer to the textbook chemistry — but "
        "the numbers on the new step have no published source yet, so it is not the default."
    ),
    "direct_burst": (
        "The default six plus an extra one: a pool of antioxidant that mops up oxygen early "
        "and then runs out. It reproduces one of the two things it was fitted to and not "
        "the other, which is why it is off unless you ask for it."
    ),
}


@dataclass(frozen=True)
class Fidelity:
    """A complete, hashable description of *how* a scenario is to be computed.

    Every field is a knob about computation, never about the scenario itself — which is why
    a run is cached on ``(scenario, fidelity)`` and why swapping one changes nothing about
    what is being simulated, only how hard the machine works at it.
    """

    #: One of :data:`PRECISION_PRESETS` or ``"custom"``.
    precision: str = "standard"
    method: str = "BDF"
    rtol: float = 1e-6
    atol: float = 1e-9
    #: Hours. ``inf`` (the engine default) lets the adaptive solver choose freely.
    max_step: float = math.inf
    #: Stored output points. Threaded identically into the nominal run and the ensemble.
    points: int = 200
    #: Which oxidative mechanism set is wired — a model-scope choice, not a solver one.
    oxidative: str = "direct"
    #: Enforce the Process ``touches`` contract during the run. Catches a Process writing a
    #: slot it never declared, at a real cost per step.
    strict: bool = False

    @classmethod
    def preset(cls, name: str, **overrides: object) -> Fidelity:
        """Build a Fidelity at one of the named precision presets."""
        rtol, atol = PRECISION_PRESETS[name]
        return cls(precision=name, rtol=rtol, atol=atol, **overrides)  # type: ignore[arg-type]

    @property
    def label(self) -> str:
        if self.precision == "custom":
            return f"custom {self.method} rtol={self.rtol:g}"
        return f"{self.precision} ({self.method})"

    def solver_kwargs(self, t_span_h: tuple[float, float]) -> dict[str, object]:
        """Keyword arguments for ``run`` / ``run_ensemble``, including the shared grid."""
        return {
            "method": self.method,
            "rtol": self.rtol,
            "atol": self.atol,
            "max_step": self.max_step,
            "t_eval": np.linspace(t_span_h[0], t_span_h[1], self.points),
        }

    def cache_key(self) -> tuple[object, ...]:
        return (
            self.precision,
            self.method,
            self.rtol,
            self.atol,
            self.max_step,
            self.points,
            self.oxidative,
            self.strict,
        )


def tightened(fidelity: Fidelity) -> Fidelity | None:
    """The next preset up, for the convergence check. ``None`` if there is none.

    A custom fidelity is tightened by two decades on both tolerances rather than by moving
    to a preset, so the comparison stays a comparison of *the same* configuration with the
    solver leaned on harder — the method and every model-scope choice are carried over.
    """
    if fidelity.precision == "custom":
        return replace(fidelity, rtol=fidelity.rtol / 100.0, atol=fidelity.atol / 100.0)
    i = PRECISION_ORDER.index(fidelity.precision)
    if i + 1 >= len(PRECISION_ORDER):
        return None
    rtol, atol = PRECISION_PRESETS[PRECISION_ORDER[i + 1]]
    return replace(fidelity, precision=PRECISION_ORDER[i + 1], rtol=rtol, atol=atol)
