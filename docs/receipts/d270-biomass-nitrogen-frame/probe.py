"""D-270 probe: what nitrogen fraction does each medium's compiled run actually use?

Reads the realised ``biomass_N_fraction`` out of a COMPILED scenario (not off the YAML),
across every declared assimilable-nitrogen level the suite's wine musts use, plus beer.
No model is integrated: the override is a compile-time pure function of the declared YAN.
"""

from __future__ import annotations

import json
from pathlib import Path

from fermentation.scenario import (
    Scenario,
    TemperaturePoint,
    amino_acid_dose_nitrogen_mgl,
    compile_scenario,
)

OUT = Path(r"M:\claud_projects\temp\ferment\d270-nitrogen-frame")

# The two sourced wine-yeast statements, as total nitrogen mass fraction of dry weight (D-267 §6).
SOURCED_LO = 40.0 / 6.25 / 100.0   # Concise Encyclopedia: 40-45 % crude protein (N x 6.25)
SOURCED_MID = 45.0 / 6.25 / 100.0
SOURCED_HI = 50.0 / 6.0 / 100.0    # Understanding Wine Chemistry: protein ~50 % DW, 1/6 N


def wine_f_n(yan_mgl: float, *, amino_acids_gpl: float | None = None) -> dict[str, float]:
    initial: dict[str, float] = {"brix": 24.0, "yan_mgl": yan_mgl, "pitch_gpl": 0.25}
    if amino_acids_gpl is not None:
        initial["amino_acids_gpl"] = amino_acids_gpl
        initial["yan_mgl"] += amino_acid_dose_nitrogen_mgl(initial)
    sc = Scenario(
        name=f"probe-wine-{yan_mgl:g}",
        medium="wine",
        initial=initial,
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=20.0)],
        duration_days=10.0,
    )
    compiled = compile_scenario(sc, strict=True)
    p = compiled.parameters["biomass_N_fraction"]
    return {
        "declared_yan_mgl": initial["yan_mgl"],
        "f_N": p.value,
        "tier": p.tier.value,
        "conditions": p.provenance.conditions,
    }


def beer_f_n(yan_mgl: float = 200.0) -> dict[str, float]:
    sc = Scenario(
        name="probe-beer",
        medium="beer",
        initial={
            "glucose_gpl": 15.0,
            "maltose_gpl": 60.0,
            "maltotriose_gpl": 15.0,
            "yan_mgl": yan_mgl,
            "pitch_gpl": 0.3984,
        },
        temperature_schedule=[TemperaturePoint(day=0.0, celsius=15.0)],
        duration_days=10.0,
    )
    compiled = compile_scenario(sc, strict=True)
    p = compiled.parameters["biomass_N_fraction"]
    return {
        "declared_yan_mgl": yan_mgl,
        "f_N": p.value,
        "tier": p.tier.value,
        "conditions": p.provenance.conditions,
    }


def main() -> None:
    rows = []
    # Every bare declared YAN the suite's wine scenarios use, low to high.
    for yan in (50.0, 80.0, 100.0, 150.0, 250.0, 300.0, 330.0, 350.0, 400.0, 500.0):
        rows.append({"arm": f"wine yan={yan:g}", **wine_f_n(yan)})
    # The migrated form: a dosed fixture declares the sum, which is what the fit reads.
    for yan, dose in ((250.0, 2.0), (250.0, 4.0), (80.0, 2.0)):
        rows.append({"arm": f"wine yan={yan:g}+aa{dose:g}", **wine_f_n(yan, amino_acids_gpl=dose)})
    rows.append({"arm": "beer yan=200 (static, no override)", **beer_f_n()})

    lo = min(r["f_N"] for r in rows if r["arm"].startswith("wine"))
    hi = max(r["f_N"] for r in rows if r["arm"].startswith("wine"))
    summary = {
        "sourced_range_g_N_per_g_DW": [SOURCED_LO, SOURCED_HI],
        "sourced_points": {"crude_protein_40": SOURCED_LO, "crude_protein_45": SOURCED_MID,
                           "uwc_arithmetic": SOURCED_HI},
        "wine_compiled_f_N_span": [lo, hi],
        "wine_arms_inside_sourced_range": sum(
            1 for r in rows if r["arm"].startswith("wine") and SOURCED_LO <= r["f_N"] <= SOURCED_HI
        ),
        "wine_arm_count": sum(1 for r in rows if r["arm"].startswith("wine")),
        "static_yaml_value": 0.114,
        "static_band": [0.08, 0.14],
        "compile_bracket": [0.03, 0.15],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "findings.json").write_text(
        json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8"
    )
    w = max(len(r["arm"]) for r in rows)
    print(f"{'arm':<{w}}  {'declared YAN':>12}  {'f_N':>8}  {'vs sourced':>12}  tier")
    for r in rows:
        if r["f_N"] < SOURCED_LO:
            verdict = "BELOW"
        elif r["f_N"] > SOURCED_HI:
            verdict = "ABOVE"
        else:
            verdict = "inside"
        print(
            f"{r['arm']:<{w}}  {r['declared_yan_mgl']:>12.1f}  {r['f_N']:>8.4f}  "
            f"{verdict:>12}  {r['tier']}"
        )
    print()
    print(f"sourced range        : {SOURCED_LO:.4f} - {SOURCED_HI:.4f} g N / g DW")
    print(f"wine compiled span   : {lo:.4f} - {hi:.4f}")
    print(f"static YAML value    : 0.1140   band [0.08, 0.14]   compile bracket [0.03, 0.15]")


if __name__ == "__main__":
    main()


def repricing() -> None:
    """What the sourced composition does to convert.py's own independent per-cell check."""
    import math

    a0, a1 = 3.50, -3.61e-3
    y_xn = math.exp(a0 + a1 * 330.0)          # g cell / g N, Coleman Fig. 4 at 330 mg N/L
    cells_per_g_n = y_xn / 4.0e-11
    print(f"\nY_X/N(330)        = {y_xn:.4f} g cell/g N")
    print(f"cells per g N     = {cells_per_g_n:.4e}")
    for label, f_n in (
        ("Roels 0.114 (shipped static)", 0.114),
        ("sourced hi 0.0833 (UWC)", 50.0 / 6.0 / 100.0),
        ("sourced lo 0.0640 (N x 6.25, 40 %)", 40.0 / 6.25 / 100.0),
        ("sourced 0.0720 (N x 6.25, 45 %)", 45.0 / 6.25 / 100.0),
        ("compiled at 330 mg N/L", 1.0 / y_xn),
        ("static band low 0.08", 0.08),
        ("static band high 0.14", 0.14),
    ):
        pg = (1.0 / f_n) / cells_per_g_n * 1e12
        print(f"  f_N = {f_n:.4f}  ({label:<34}) -> {pg:7.2f} pg/cell")

    print("\nBeer: cell nitrogen Tyrell's counted crop would demand (D-230): 0.202-0.262 g N/g DW")
    for label, ref in (
        ("static band top 0.14", 0.14),
        ("shipped static 0.114", 0.114),
        ("sourced hi 0.0833", 50.0 / 6.0 / 100.0),
        ("sourced lo 0.0640", 40.0 / 6.25 / 100.0),
    ):
        print(f"  vs {label:<22}: {0.202 / ref:.2f}-{0.262 / ref:.2f}x outside")

    print("\nBeer biomass ceiling dX = YAN/f_N, relative to the shipped static 0.114:")
    for label, f_n in (
        ("sourced hi 0.0833", 50.0 / 6.0 / 100.0),
        ("sourced mid 0.0720", 45.0 / 6.25 / 100.0),
        ("sourced lo 0.0640", 40.0 / 6.25 / 100.0),
    ):
        print(f"  {label:<20}: x{0.114 / f_n:.3f} biomass  "
              f"(extent 5.378x -> {5.378 * 0.114 / f_n:.2f}x)")
