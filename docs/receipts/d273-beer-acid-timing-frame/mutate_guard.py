"""D-273 mutation arm -- does the tightened inter-acid guard actually catch a per-acid shape?

The guard's docstring claims a future Process giving these acids a shape of their own would
break it. That claim is only worth what the tolerance is worth. Measured inter-acid spread is
1.1e-13 points; the tolerance is now 1e-9. This gives succinic ALONE an ethanol-dependent
factor -- a shape, not a level -- and checks the guard goes red at a size the old 1e-6 would
have waved through.
"""

import sys

import numpy as np

sys.path.insert(0, "M:/claud_projects/Fermentation")

import fermentation.core.kinetics.organic_acids as oa  # noqa: E402

_original = oa.organic_acid_rates
SCALE = 1.0e-10  # per g/L of ethanol; ~1.2e-8 relative at 12 % ABV -- INVISIBLE to the old 1e-6 tolerance


def _patched(y, schema, params):
    ethanol = float(np.asarray(y)[schema.slice("E")][0])
    return [
        (spec, rate * (1.0 + SCALE * ethanol) if spec.slot == "succinic" else rate)
        for spec, rate in _original(y, schema, params)
    ]


def rise_fractions(patched: bool) -> dict[str, dict[int, float]]:
    oa.organic_acid_rates = _patched if patched else _original
    from tests.test_organic_acids import TYRELL_ACID_COURSE_PPM, TYRELL_SCENARIO, _daily_ppm, _run

    _, res = _run(dict(TYRELL_SCENARIO), days=7.0)
    got = {slot: _daily_ppm(res, slot) for slot in TYRELL_ACID_COURSE_PPM}
    return {
        slot: {
            day: 100.0 * (c[day] - c[0]) / (c[7] - c[0]) for day in range(1, 8)
        }
        for slot, c in got.items()
    }


for label, patched in (("BASELINE", False), ("MUTANT ", True)):
    fr = rise_fractions(patched)
    worst = max(
        abs(fr[slot][day] - fr["lactic"][day]) for slot in fr for day in range(1, 8)
    )
    verdict_new = "RED" if worst > 1e-9 else "green"
    verdict_old = "RED" if worst > 1e-6 else "green"
    print(
        f"{label}: worst |succinic-lactic| = {worst:.3e} points"
        f"   -> tolerance 1e-9: {verdict_new}   old tolerance 1e-6: {verdict_old}"
    )
