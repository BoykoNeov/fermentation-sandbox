"""D-271 probe: price the count-plus-weight pairing hunt's arithmetic.

No model is integrated here. Every number is an identity or a unit conversion on
figures transcribed from named sources.
"""

from __future__ import annotations

import json
import math
import pathlib

# --- sourced constants -----------------------------------------------------
# Klis, de Koster & Brul 2014, Eukaryot. Cell (PMC3910951), Table 1 footnote a:
# "The biomass (dry weight) was calculated by multiplying the volume with the
# density (1.11) to obtain the biomass (wet weight) of the cell and multiplying
# the obtained value with the dry weight fraction (0.34) of the wet weight."
KLIS_DENSITY_G_PER_ML = 1.11
KLIS_DRY_FRACTION = 0.34
KLIS_HAPLOID_VOLUME_FL = 44.0
KLIS_DIPLOID_VOLUME_FL = 83.0
KLIS_HAPLOID_PG = 16.5
KLIS_DIPLOID_PG = 31.2

# D-219's geometric cross-check, its own assumed cell-volume range for a
# wine/ale cell, and the band it printed off an UNSOURCED "~30 % dry matter".
D219_ASSUMED_VOLUME_FL = (100.0, 150.0)
D219_PRINTED_CROSSCHECK_PG = (30.0, 57.0)
D219_SETTLED_PG = 40.0
D219_SETTLED_BAND_PG = (28.0, 50.0)

# D-230 branch 1: Tyrell's counted crop read as a per-cell dry mass. UNCHANGED
# by D-270 -- D-270 sec 7 re-priced the ENGINE-side elemental estimate, not this.
BRANCH1_DEMAND_PG = (70.9, 91.9)
D270_ELEMENTAL_REPRICED_PG = (47.71, 62.12)

# Rule of thumb, BNID 101795, from Physical Biology of the Cell Table 1.1.
RULE_OF_THUMB_PG = 60.0

# Okada et al. 2023, Sci. Rep. (PMC9883461), Table 2: APPARENT diameters [um] from
# flow-cytometry forward scatter calibrated on size beads -- a DIFFERENT measurement
# frame from Klis's microscopy volumes, which is the point of comparing them.
OKADA_APPARENT_DIAMETER_UM = {
    "haploid_BY4741L": 7.3,
    "haploid_BY4742L": 7.1,
    "diploid_BY4743L": 9.4,
    "triploid": 12.0,
    "tetraploid": 14.5,
    "brewing_diploid_K7A": 12.6,
}


def volume_fl_to_dry_pg(volume_fl: float, dry_fraction: float = KLIS_DRY_FRACTION) -> float:
    """Cell volume [fL] -> dry mass [pg], by Klis's stated route.

    1 fL = 1e-12 mL, so wet mass [g] = volume_fl * 1e-12 * rho, and 1 g = 1e12 pg.
    """
    wet_g = volume_fl * 1e-12 * KLIS_DENSITY_G_PER_ML
    return wet_g * dry_fraction * 1e12


def dry_pg_to_volume_fl(dry_pg: float, dry_fraction: float = KLIS_DRY_FRACTION) -> float:
    """Inverse of :func:`volume_fl_to_dry_pg`."""
    wet_g = dry_pg * 1e-12 / dry_fraction
    return wet_g / KLIS_DENSITY_G_PER_ML * 1e12


def sphere_diameter_um(volume_fl: float) -> float:
    """Equivalent-sphere diameter [um] for a volume in fL (1 fL == 1 um^3)."""
    return (6.0 * volume_fl / math.pi) ** (1.0 / 3.0)


def main() -> dict[str, object]:
    # 1. Klis's own two figures reproduce from his stated route -- the transcription check.
    reproduced = {
        "haploid_pg": volume_fl_to_dry_pg(KLIS_HAPLOID_VOLUME_FL),
        "diploid_pg": volume_fl_to_dry_pg(KLIS_DIPLOID_VOLUME_FL),
        "printed_haploid_pg": KLIS_HAPLOID_PG,
        "printed_diploid_pg": KLIS_DIPLOID_PG,
    }

    # 2. D-219's cross-check, re-run on the SOURCED dry fraction instead of "~30 %".
    sourced_crosscheck = tuple(volume_fl_to_dry_pg(v) for v in D219_ASSUMED_VOLUME_FL)
    at_thirty = tuple(volume_fl_to_dry_pg(v, 0.30) for v in D219_ASSUMED_VOLUME_FL)
    # what dry fraction each printed edge implies, at the matching volume edge
    implied_fraction = (
        D219_PRINTED_CROSSCHECK_PG[0] * 1e-12 / (D219_ASSUMED_VOLUME_FL[0] * 1e-12 * KLIS_DENSITY_G_PER_ML),
        D219_PRINTED_CROSSCHECK_PG[1] * 1e-12 / (D219_ASSUMED_VOLUME_FL[1] * 1e-12 * KLIS_DENSITY_G_PER_ML),
    )

    # 3. Branch 1's demand, expressed as a cell SIZE through the sourced constants.
    demand_volume_fl = tuple(dry_pg_to_volume_fl(m) for m in BRANCH1_DEMAND_PG)
    demand_diameter_um = tuple(sphere_diameter_um(v) for v in demand_volume_fl)
    demand_vs_klis_diploid = tuple(v / KLIS_DIPLOID_VOLUME_FL for v in demand_volume_fl)
    demand_vs_d219_assumed = (
        demand_volume_fl[0] / D219_ASSUMED_VOLUME_FL[1],
        demand_volume_fl[1] / D219_ASSUMED_VOLUME_FL[0],
    )

    # 4. Where each non-pairing literature figure lands against branch 1's demand.
    def lands(value: float) -> str:
        lo, hi = BRANCH1_DEMAND_PG
        return "below" if value < lo else ("above" if value > hi else "inside")

    landings = {
        "klis_haploid": lands(KLIS_HAPLOID_PG),
        "klis_diploid": lands(KLIS_DIPLOID_PG),
        "rule_of_thumb_60pg": lands(RULE_OF_THUMB_PG),
        "d219_settled_40pg": lands(D219_SETTLED_PG),
        "d270_elemental_low": lands(D270_ELEMENTAL_REPRICED_PG[0]),
        "d270_elemental_high": lands(D270_ELEMENTAL_REPRICED_PG[1]),
        "sourced_crosscheck_low": lands(sourced_crosscheck[0]),
        "sourced_crosscheck_high": lands(sourced_crosscheck[1]),
    }

    # 5. THE SIZE ROUTE IS FRAME-BROKEN: the two size sources disagree on the SAME
    #    lab strains, so neither can be carried onto the other's scale.
    def sphere_volume_fl(diameter_um: float) -> float:
        return math.pi / 6.0 * diameter_um**3

    okada_haploid_volume_fl = sphere_volume_fl(OKADA_APPARENT_DIAMETER_UM["haploid_BY4741L"])
    okada_diploid_volume_fl = sphere_volume_fl(OKADA_APPARENT_DIAMETER_UM["diploid_BY4743L"])
    frame_gap = {
        "klis_haploid_volume_fl": KLIS_HAPLOID_VOLUME_FL,
        "klis_haploid_equivalent_diameter_um": sphere_diameter_um(KLIS_HAPLOID_VOLUME_FL),
        "okada_haploid_apparent_diameter_um": OKADA_APPARENT_DIAMETER_UM["haploid_BY4741L"],
        "okada_haploid_implied_volume_fl": okada_haploid_volume_fl,
        "haploid_volume_ratio_okada_over_klis": okada_haploid_volume_fl / KLIS_HAPLOID_VOLUME_FL,
        "klis_diploid_volume_fl": KLIS_DIPLOID_VOLUME_FL,
        "okada_diploid_implied_volume_fl": okada_diploid_volume_fl,
        "diploid_volume_ratio_okada_over_klis": okada_diploid_volume_fl / KLIS_DIPLOID_VOLUME_FL,
        # branch 1's demanded volume is INSIDE Okada's frame and ABOVE Klis's --
        # which is exactly why the size route cannot adjudicate.
        "branch1_demand_vs_okada_haploid_volume": (
            demand_volume_fl[0] / okada_haploid_volume_fl,
            demand_volume_fl[1] / okada_haploid_volume_fl,
        ),
    }

    out: dict[str, object] = {
        "size_route_is_frame_broken": frame_gap,
        "klis_route_reproduces_its_own_printed_values": reproduced,
        "d219_crosscheck_at_sourced_0_34_pg": sourced_crosscheck,
        "d219_crosscheck_at_0_30_pg": at_thirty,
        "d219_printed_crosscheck_pg": D219_PRINTED_CROSSCHECK_PG,
        "dry_fraction_d219_printed_edges_imply": implied_fraction,
        "branch1_demand_pg": BRANCH1_DEMAND_PG,
        "branch1_demanded_volume_fl": demand_volume_fl,
        "branch1_demanded_sphere_diameter_um": demand_diameter_um,
        "branch1_demand_over_klis_diploid_volume": demand_vs_klis_diploid,
        "branch1_demand_over_d219_assumed_volume": demand_vs_d219_assumed,
        "where_each_non_pairing_figure_lands_vs_branch1": landings,
    }
    return out


if __name__ == "__main__":
    result = main()
    text = json.dumps(result, indent=2)
    pathlib.Path(__file__).with_name("findings.json").write_text(text, encoding="utf-8")
    print(text)
