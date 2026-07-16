"""Tests for the speciated amino-acid pool registry — the D-100 shared idioms.

The lump `amino_acids` became eight single-molecule pools at D-100, and the arithmetic every
consumer shares lives in :mod:`fermentation.core.kinetics.amino_acid_pools`. This suite pins the
properties that make the split safe rather than merely plausible:

* **the reduction property** — at must-spectrum composition every per-species relative-depletion
  gate is ALGEBRAICALLY the pre-split lumped gate, so the split does not silently move the dosed
  baseline (the advisor's must-pin, and the reason every closed-form suite could keep asserting
  its old numbers);
* **the decoupling** — no Ehrlich/Strecker precursor is drawn by an identity-agnostic consumer and
  vice versa, which is the structural claim D-100 rests on;
* the C:N orderings the D-32 no-sugar-creation, D-34 positive-debris and D-89 positive-denominator
  guarantees depend on — now for BLENDS, not one species.
"""

import pytest

from fermentation.core.chemistry import carbon_mass_fraction, nitrogen_mass_fraction
from fermentation.core.kinetics.amino_acid_pools import (
    AMINO_ACID_SPECS,
    ASSIMILABLE_SPECS,
    SPEC_BY_SPECIES,
    amino_acid_pool,
    assimilable_carbon_per_nitrogen,
    depletion_gate,
    draw_assimilable_nitrogen,
    draw_precursor_carbon,
    release_spectrum_nitrogen,
    spectrum_carbon_per_nitrogen,
)
from fermentation.core.media import wine_schema
from fermentation.parameters.store import default_data_dir, load_parameters
from tests.conftest import seed_amino_acids


@pytest.fixture
def params():
    return load_parameters(default_data_dir() / "wine_generic.yaml").resolve()


@pytest.fixture
def schema():
    return wine_schema()


# -- THE reduction property (the advisor's must-pin) --------------------------


@pytest.mark.parametrize("total", [0.05, 0.3, 0.8, 2.0, 10.0])
def test_every_per_species_gate_reduces_to_the_lumped_gate_at_spectrum_composition(
    schema, params, total
):
    """At must-spectrum composition, gate_i == the pre-split lumped gate — for EVERY species.

    THE property that makes the D-100 gate rule honest rather than a guess. ``K_i = K·f_i`` is
    derived from constants already sourced (``K_amino_acids`` + the must spectrum), and it means
    ``aa_i/(K·f_i + aa_i) = f_i·aa/(f_i·K + f_i·aa) = aa/(K + aa)`` exactly. So speciation is a
    provable no-op on the rate at the composition a dose creates, and only the pool DRIFTING away
    from that composition — the actual D-100 pathology — changes anything.

    Without this, every closed-form assertion the eight consumer suites carry would have been
    silently re-baselined by the split rather than preserved.
    """
    y = schema.zeros()
    seed_amino_acids(y, schema, params, total)
    lumped = total / (params["K_amino_acids"] + total)
    for spec in AMINO_ACID_SPECS:
        assert depletion_gate(y, schema, params, (spec,)) == pytest.approx(lumped, rel=1e-12)
    # ...and for any SUBSET, since the rule scales K by the substrate's summed share. This is what
    # lets the identity-agnostic consumers gate on {arginine, generic} without shrinking to 0.81x.
    assert depletion_gate(y, schema, params, ASSIMILABLE_SPECS) == pytest.approx(lumped, rel=1e-12)
    assert depletion_gate(y, schema, params, AMINO_ACID_SPECS) == pytest.approx(lumped, rel=1e-12)


def test_the_gate_bites_per_species_once_the_pool_leaves_spectrum_composition(schema, params):
    """Draining ONE species throttles only ITS gate — the D-100 pathology, now expressible.

    The inverse of the reduction property, and the reason the split is not cosmetic. Under the
    lump, a pool 38% arginine reported "amino acids available" while leucine sat at zero, so the
    re-route kept drawing leucine carbon that did not exist. Each gate now reads its own molecule.
    """
    y = schema.zeros()
    seed_amino_acids(y, schema, params, 0.8)
    before = depletion_gate(y, schema, params, (SPEC_BY_SPECIES["leucine"],))
    others = {
        spec.species: depletion_gate(y, schema, params, (spec,))
        for spec in AMINO_ACID_SPECS
        if spec.species != "leucine"
    }
    y[schema.slice("leucine")] = 0.0  # the Ehrlich pathway ate it
    assert depletion_gate(y, schema, params, (SPEC_BY_SPECIES["leucine"],)) == 0.0
    assert before > 0.0
    for spec in AMINO_ACID_SPECS:
        if spec.species == "leucine":
            continue
        assert depletion_gate(y, schema, params, (spec,)) == pytest.approx(others[spec.species])


def test_empty_and_untracked_pools_gate_to_exactly_zero(schema, params):
    # The isolability guarantee every consumer leans on (prime directive #3): an undosed run has
    # every gate at EXACTLY 0, so the validated core is byte-for-byte untouched.
    y = schema.zeros()
    for spec in AMINO_ACID_SPECS:
        assert depletion_gate(y, schema, params, (spec,)) == 0.0
    assert depletion_gate(y, schema, params, AMINO_ACID_SPECS) == 0.0
    # ...and a solver undershoot below zero cannot resurrect it (or flip a draw's sign).
    y[schema.slice("leucine")] = -1e-12
    assert amino_acid_pool(y, schema, "leucine") == 0.0
    assert depletion_gate(y, schema, params, (SPEC_BY_SPECIES["leucine"],)) == 0.0


# -- the decoupling: the structural claim D-100 rests on ----------------------


def test_no_precursor_is_also_an_identity_agnostic_pool():
    """The two draw idioms address DISJOINT pools — the D-100 decoupling, at the registry level.

    Fusel production drains precursors; bacterial/yeast growth drains {arginine, generic}. Because
    those sets do not intersect, no amount of Ehrlich flux can starve MLF, Brett, the yeast swap or
    Maillard browning. That was the D-99 finding: three unrelated subsystems broke through one
    lumped substrate. If a future beat ever makes a precursor identity-agnostic (or vice versa),
    this test fails and the decoupling argument must be re-made rather than silently lost.
    """
    assimilable = {spec.pool for spec in ASSIMILABLE_SPECS}
    precursors = {spec.pool for spec in AMINO_ACID_SPECS} - assimilable
    assert assimilable == {"amino_acids", "amino_acids_generic"}
    assert precursors == {
        "leucine",
        "isoleucine",
        "valine",
        "threonine",
        "phenylalanine",
        "methionine",
    }
    assert not (assimilable & precursors)


def test_every_pool_is_weighted_on_both_conservation_ledgers():
    # Prime directive: a pool carrying carbon and nitrogen must be weighted in both ledgers, or the
    # split would leak mass. Every species must resolve in the chemistry source of truth.
    for spec in AMINO_ACID_SPECS:
        assert carbon_mass_fraction(spec.species) > 0.0
        assert nitrogen_mass_fraction(spec.species) > 0.0


# -- the draw idioms close their own ledgers ---------------------------------


@pytest.mark.parametrize("species", ["leucine", "threonine", "methionine", "phenylalanine"])
def test_precursor_draw_is_carbon_exact_and_returns_its_own_nitrogen(schema, params, species):
    # The per-precursor idiom: carbon out of the pool == carbon the caller asked for, and the
    # nitrogen returned is what THAT molecule carries (not arginine's — the retired D-33 lump
    # over-released ~4x because it deaminated a 4-nitrogen molecule for a 1-nitrogen job).
    d = schema.zeros()
    carbon = 1.0e-4
    nitrogen = draw_precursor_carbon(d, schema, species, carbon)
    mass = -float(d[schema.slice(species)][0])
    assert mass * carbon_mass_fraction(species) == pytest.approx(carbon, rel=1e-12)
    assert nitrogen == pytest.approx(mass * nitrogen_mass_fraction(species), rel=1e-12)


def test_precursor_draws_accumulate_so_shared_precursors_are_not_dropped(schema, params):
    # `+=`, not `=`: threonine feeds BOTH propanol (Ehrlich) and sotolon (D-87), and the
    # MaillardStrecker route
    # draws five precursors in one pass. A plain assignment would silently drop all but the last
    # draw — carbon would still "close" against a smaller draw, so no conservation test would catch
    # it. This is the guard.
    d = schema.zeros()
    draw_precursor_carbon(d, schema, "threonine", 1.0e-4)
    draw_precursor_carbon(d, schema, "threonine", 3.0e-4)
    mass = -float(d[schema.slice("threonine")][0])
    assert mass * carbon_mass_fraction("threonine") == pytest.approx(4.0e-4, rel=1e-12)


def test_assimilable_draw_splits_by_nitrogen_and_matches_its_advertised_ratio(schema, params):
    # The identity-agnostic idiom: the demand is met from {arginine, generic} in proportion to the
    # nitrogen each holds, and the carbon returned is EXACTLY `nitrogen * the advertised ratio`.
    # MaillardBrowning sizes its melanoidin from that ratio BEFORE drawing, so any drift between
    # the two would break its carbon ledger silently.
    y = schema.zeros()
    seed_amino_acids(y, schema, params, 0.8)
    d = schema.zeros()
    nitrogen = 1.0e-4
    carbon = draw_assimilable_nitrogen(d, y, schema, nitrogen)
    ratio = assimilable_carbon_per_nitrogen(y, schema)
    assert carbon == pytest.approx(nitrogen * ratio, rel=1e-12)
    # nitrogen closes: the two pools give up exactly the demand
    drawn_n = sum(
        -float(d[schema.slice(spec.pool)][0]) * nitrogen_mass_fraction(spec.species)
        for spec in ASSIMILABLE_SPECS
    )
    assert drawn_n == pytest.approx(nitrogen, rel=1e-12)
    # and no precursor was touched (the decoupling, at the arithmetic level)
    for spec in AMINO_ACID_SPECS:
        if spec in ASSIMILABLE_SPECS:
            continue
        assert float(d[schema.slice(spec.pool)][0]) == 0.0


def test_the_blend_ratio_stays_between_its_two_members_for_any_composition(schema, params):
    """The C:N ordering the D-32/D-89 guarantees rest on — now for BLENDS, not one species.

    D-32's no-sugar-creation and D-89's positive-denominator both turned on arginine's C:N (1.29)
    sitting below biomass's (4.3) and melanoidin's (8). D-100 replaces that fixed number with a
    state-dependent blend, so the guarantees now need the blend to be BOUNDED — which it is, by its
    two members, for any pool composition. That is what keeps both proofs structural (no clamp, no
    C-zero kink for the BDF solver) rather than empirical.
    """
    lo = carbon_mass_fraction("arginine") / nitrogen_mass_fraction("arginine")
    hi = carbon_mass_fraction("glutamine") / nitrogen_mass_fraction("glutamine")
    assert lo < hi < 4.0  # both below biomass's C:N ~4.3 — the D-32 guarantee's load-bearing fact
    for arg, gen in [(1.0, 0.0), (0.0, 1.0), (0.5, 0.5), (0.9, 0.1), (0.05, 0.95)]:
        y = schema.zeros()
        y[schema.slice("amino_acids")] = arg
        y[schema.slice("amino_acids_generic")] = gen
        assert lo - 1e-12 <= assimilable_carbon_per_nitrogen(y, schema) <= hi + 1e-12


# -- the release idiom (autolysis, the model's only amino-acid source) --------


def test_spectrum_release_is_nitrogen_exact_and_holds_the_must_composition(schema, params):
    # The inverse of the draws: autolysis (D-34) deposits the dead-cell nitrogen across all eight
    # pools AT MUST-SPECTRUM COMPOSITION. Under the lump this released pure arginine — refilling a
    # pool without restoring a single aroma precursor. Releasing the spectrum is what makes
    # sur-lie thermal aroma work (test_thermal_aroma_from_drained_precursors_requires_autolysis).
    d = schema.zeros()
    nitrogen = 1.0e-4
    carbon = release_spectrum_nitrogen(d, schema, params, nitrogen)
    released_n = sum(
        float(d[schema.slice(spec.pool)][0]) * nitrogen_mass_fraction(spec.species)
        for spec in AMINO_ACID_SPECS
    )
    assert released_n == pytest.approx(nitrogen, rel=1e-12)
    assert carbon == pytest.approx(nitrogen * spectrum_carbon_per_nitrogen(params), rel=1e-12)
    # every pool gains, in must-spectrum mass ratios
    total_mass = sum(float(d[schema.slice(spec.pool)][0]) for spec in AMINO_ACID_SPECS)
    denom = sum(params[spec.fraction_param] for spec in AMINO_ACID_SPECS)
    for spec in AMINO_ACID_SPECS:
        share = float(d[schema.slice(spec.pool)][0]) / total_mass
        assert share > 0.0
        assert share == pytest.approx(params[spec.fraction_param] / denom, rel=1e-12)


def test_the_spectrum_ratio_keeps_autolysis_debris_structurally_positive(params):
    """D-34's no-clamp guarantee, re-verified for the spectrum (its margin narrowed at D-100).

    Autolysis routes the dead-cell carbon its released amino acids CANNOT carry to debris. That
    remainder is `r·(f_C − f_N·R)`, so it stays positive only while biomass C:N (4-11 across
    Coleman's whole nitrogen range) exceeds the released blend's R. D-100 raised R from arginine's
    1.29 to the spectrum's ~1.9 — a real narrowing, and worth pinning rather than assuming, since a
    negative remainder would put a C-zero kink under the BDF solver and silently create carbon.
    """
    ratio = spectrum_carbon_per_nitrogen(params)
    assert (
        1.5 < ratio < 2.5
    )  # the must spectrum's mass C:N — between arginine's and the C-rich AAs'
    assert ratio < 4.0  # below the FLOOR of biomass C:N (~4 at Coleman's N-richest) — the guarantee


# -- the registry is the single source of truth ------------------------------


def test_spectrum_fractions_are_all_sourced_and_positive(params):
    # Prime directive #2: every fraction is a provenance-backed parameter, not an inline number.
    # They are deliberately NOT asserted to sum to 1 — they are eight independent entries, each
    # re-sourceable on its own, and an ensemble sampling their bands would break such an assertion
    # on nearly every draw. The compile seam and the release idiom both NORMALIZE instead.
    for spec in AMINO_ACID_SPECS:
        assert params[spec.fraction_param] > 0.0
    assert sum(params[spec.fraction_param] for spec in AMINO_ACID_SPECS) == pytest.approx(
        1.0, abs=0.05
    )


def test_spec_by_species_covers_every_pool():
    assert set(SPEC_BY_SPECIES) == {spec.species for spec in AMINO_ACID_SPECS}
    assert len(AMINO_ACID_SPECS) == 8
    # proline is EXCLUDED and must stay so: it is not assimilated anaerobically (excluded from YAN
    # by definition), so `amino_acids_gpl` is honestly an *assimilable* dose and must's ~48% proline
    # never enters the model.
    assert "proline" not in SPEC_BY_SPECIES
