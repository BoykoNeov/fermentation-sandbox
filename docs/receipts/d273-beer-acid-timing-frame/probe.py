"""D-273 probe: beer's three flux-linked acid courses, scored in three frames.

Every number in the record comes from here. The point of the probe is the FRAME, not a new
measurement: `TYRELL_ACID_COURSE_PPM` and `TYRELL_FLUX_FRACTION` were both transcribed at
D-215 off the same four ferments of the same wort, and D-215 §3 only ever compared the first
against the MODEL's flux. Compared against the source's own flux they say something different.

Three frames:

* **A, calendar** -- measured acid vs modelled acid by day. What D-215 §3 used. Carries the
  model's own ~3x day-2 speed deficit (D-215 §4, parked at D-223) inside every number.
* **B, source-internal** -- measured acid vs Tyrell's measured flux. No model quantity at all.
* **C, model-internal** -- modelled acid vs the model's OWN fermented fraction. Immune to the
  speed defect, because the model's clock never enters.

Run: ``uv run python docs/receipts/d273-beer-acid-timing-frame/probe.py``
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tests.test_organic_acids import (  # noqa: E402
    TYRELL_ACID_COURSE_PPM,
    TYRELL_ACID_COURSE_READ_TOL,
    TYRELL_FLUX_FRACTION,
    TYRELL_SCENARIO,
    _run,
)

OUT = Path(__file__).resolve().parent
DAYS = 12.0


def main() -> dict[str, object]:
    compiled, res = _run(dict(TYRELL_SCENARIO), days=DAYS)
    grid = np.linspace(0.0, DAYS * 24.0, int(DAYS * 24) + 1)
    sugar = np.asarray(res.y, dtype=float)[compiled.schema.slice("S"), :]
    total = np.vstack([np.interp(grid, res.t, row) for row in sugar]).sum(axis=0)
    ferm = (float(total[0]) - total) / float(total[0])
    acid = {
        slot: np.interp(grid, res.t, np.asarray(res.series(slot), dtype=float)) * 1000.0
        for slot in TYRELL_ACID_COURSE_PPM
    }
    out: dict[str, object] = {}

    # -- frame C: the model against its own fermented fraction ------------------------------
    print("FRAME C -- model acid-rise fraction vs the model's OWN fermented fraction")
    print("(both normalised on day 7; the acids are one curve, so one column)")
    print(f"{'day':>4} {'acid %':>9} {'ferm %':>9} {'residual':>10}")
    frame_c = {}
    for day in range(1, 8):
        hour, hour7 = day * 24, 7 * 24
        series = acid["lactic"]
        acid_pct = 100.0 * (series[hour] - series[0]) / (series[hour7] - series[0])
        ferm_pct = 100.0 * ferm[hour] / ferm[hour7]
        frame_c[day] = {"acid": acid_pct, "ferm": ferm_pct, "residual": acid_pct - ferm_pct}
        print(f"{day:>4} {acid_pct:>9.2f} {ferm_pct:>9.2f} {acid_pct - ferm_pct:>+10.3f}")
    out["frame_c_model_own"] = frame_c

    spread = max(
        abs(
            100.0 * (s[d * 24] - s[0]) / (s[7 * 24] - s[0])
            - 100.0 * (acid["lactic"][d * 24] - acid["lactic"][0])
            / (acid["lactic"][7 * 24] - acid["lactic"][0])
        )
        for s in acid.values()
        for d in range(1, 8)
    )
    print(f"\nlargest disagreement BETWEEN the three modelled acids: {spread:.2e} points")
    out["model_inter_acid_spread_points"] = spread

    # -- frame B: the source against itself -------------------------------------------------
    print("\nFRAME B -- measured acid vs Tyrell's OWN flux (lag, points of the day-0->7 rise)")
    frame_b = {}
    for day in (1, 2, 5):
        frame_b[day] = {}
        for slot, course in TYRELL_ACID_COURSE_PPM.items():
            done = 100.0 * (course[day] - course[0]) / (course[7] - course[0])
            frame_b[day][slot] = 100.0 * TYRELL_FLUX_FRACTION[day] - done
        line = "  ".join(f"{k}={v:+.1f}" for k, v in frame_b[day].items())
        print(f"day {day} (flux {100 * TYRELL_FLUX_FRACTION[day]:5.1f} %): {line}")
    out["frame_b_lag_points"] = frame_b

    print("\n  read floor per acid (+-2 ppm as a share of that acid's whole rise):")
    floors = {}
    for slot, course in TYRELL_ACID_COURSE_PPM.items():
        floors[slot] = 100.0 * TYRELL_ACID_COURSE_READ_TOL / (course[7] - course[0])
        print(f"    {slot:>9}: {floors[slot]:.1f} points")
    out["read_floor_points"] = floors

    # -- why the calendar frame changes sign ------------------------------------------------
    print("\nWHY FRAME A CHANGES SIGN -- the model's flux deficit against the span of the lags")
    deficits = {}
    for day in (1, 2, 5):
        deficit = 100.0 * TYRELL_FLUX_FRACTION[day] - 100.0 * float(ferm[day * 24])
        low, high = min(frame_b[day].values()), max(frame_b[day].values())
        deficits[day] = {"deficit": deficit, "low": low, "high": high, "inside": low < deficit < high}
        print(
            f"day {day}: deficit {deficit:+6.1f}   lags span [{low:+.1f}, {high:+.1f}]"
            f"   INSIDE = {low < deficit < high}"
        )
    out["deficit_vs_lag_span"] = deficits

    # -- the post-attenuation rise ----------------------------------------------------------
    print("\nPOST-ATTENUATION -- share of the day-0->7 rise arriving after the wort is out")
    print(f"  Tyrell (his flux reads {TYRELL_FLUX_FRACTION[5]:.3f} at day 5):")
    post = {"measured": {}}
    for slot, course in TYRELL_ACID_COURSE_PPM.items():
        after5 = 100.0 * (course[7] - course[5]) / (course[7] - course[0])
        after4 = 100.0 * (course[7] - course[4]) / (course[7] - course[0])
        post["measured"][slot] = {"after_day5": after5, "after_day4": after4}
        print(f"    {slot:>9}: {after5:5.1f} % after day 5   {after4:5.1f} % after day 4")

    attenuated = int(np.argmax(ferm >= 0.99))
    print(
        f"\n  model (99 % attenuated at hour {attenuated}, day {attenuated / 24:.2f}; "
        f"run ends at {100 * float(ferm[-1]):.2f} %):"
    )
    post["model"] = {}
    for slot, series in acid.items():
        share = 100.0 * float(series[-1] - series[attenuated]) / float(series[-1] - series[0])
        post["model"][slot] = share
        print(f"    {slot:>9}: {share:5.3f} % of its full rise after its own 99 % attenuation")
    post["model_attenuated_hour"] = attenuated
    out["post_attenuation_share_pct"] = post

    (OUT / "findings.json").write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\nwrote {OUT / 'findings.json'}")
    return out


if __name__ == "__main__":
    main()
