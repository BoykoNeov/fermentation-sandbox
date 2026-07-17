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

from collections.abc import Mapping

import pytest

from fermentation.core.chemistry import (
    MOLAR_MASS,
    carbon_mass_fraction,
    nitrogen_mass_fraction,
)
from fermentation.core.kinetics.aging import (
    _CO2_PER_STRECKER_ALDEHYDE,
    _MAILLARD_PRODUCTS,
    _STRECKER_ROUTES,
    MaillardStrecker,
    StreckerDegradation,
)
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
from fermentation.core.kinetics.byproducts import FuselAminoAcidReroute
from fermentation.core.kinetics.carbon_routing import FUSEL_SPECS
from fermentation.core.media import wine_schema
from fermentation.core.state import FloatArray, StateSchema
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
    ``aa_i/(K·f_i + aa_i) = f_i·aa/(f_i·K + f_i·aa) = aa/(K + aa)``. So speciation is a no-op on the
    rate at the composition a dose creates, and only the pool DRIFTING away from that composition —
    the actual D-100 pathology — changes anything. Without it, every closed-form assertion the
    eight consumer suites carry would have been silently re-baselined by the split.

    TWO CLAIMS OF DIFFERENT STANDING, asserted separately below. That all species and subsets agree
    is STRUCTURAL (any fractions). That the shared value equals the PRE-SPLIT LUMPED gate is
    CONTINGENT on the spectrum summing to 1: a dose is apportioned `f_i·D/Σf`, so a spectrum summing
    to F seeds `D/F` and gives `(D/F)/(K + D/F)`. The sourced fractions sum to exactly 1.000 today,
    so it holds — but an ENSEMBLE sampling their uncertainty bands has F != 1 on nearly every draw
    and shifts the baseline slightly (acceptable: the fractions are speculative). This test pins the
    CURRENT values and will fail if a re-source breaks the sum — which is the point.
    """
    y = schema.zeros()
    seed_amino_acids(y, schema, params, total)

    # (a) STRUCTURAL: every species and every subset agree, because the rule scales K by the
    # substrate's summed share. This is what lets the identity-agnostic consumers gate on
    # {arginine, generic} without silently shrinking to 0.81x.
    gates = [depletion_gate(y, schema, params, (spec,)) for spec in AMINO_ACID_SPECS]
    for gate in gates:
        assert gate == pytest.approx(gates[0], rel=1e-12)
    assert depletion_gate(y, schema, params, ASSIMILABLE_SPECS) == pytest.approx(
        gates[0], rel=1e-12
    )
    assert depletion_gate(y, schema, params, AMINO_ACID_SPECS) == pytest.approx(gates[0], rel=1e-12)

    # (b) CONTINGENT on the sourced fractions summing to 1 (they do, exactly): the shared value is
    # the PRE-SPLIT lumped gate, so the split does not move the dosed baseline.
    assert sum(params[spec.fraction_param] for spec in AMINO_ACID_SPECS) == pytest.approx(
        1.0, abs=1e-12
    )
    lumped = total / (params["K_amino_acids"] + total)
    assert gates[0] == pytest.approx(lumped, rel=1e-12)


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


# -- the stoichiometric signature: a closed ledger is NOT a correct draw ------


#: Every precursor→product route in the tree, with the **carbon each one charges to its precursor's
#: draw beyond the named product**: the ``CO₂`` it books, plus any **tracked co-product** it books
#: (decision D-105, generalised at D-107). Together those decide whether a carbon-sized draw is also
#: the real stoichiometry. Read off the route tables: :data:`_STRECKER_ROUTES` and
#: :data:`_MAILLARD_PRODUCTS` each charge 1 CO₂ and no co-product; the Ehrlich re-route charges its
#: decarboxylation CO₂ since D-106; and the mercaptan charges **no CO₂ but a C4 co-product** —
#: α-ketobutyrate — since D-107.
#:
#: **The co-product column is what D-105 predicted and could not yet express.** Its wording was that
#: a failing route "is either sourcing carbon it does not name (sotolon's acetaldehyde) or
#: **discarding carbon it should charge for** (the mercaptan's 2-oxobutyrate, the Ehrlich CO₂)".
#: Charging is only possible into a pool that exists — CO₂ already did (so D-106 could fix the
#: Ehrlich routes immediately), and α-ketobutyrate did not (so the mercaptan waited for D-107 to
#: build it). The signature did not change; the model finally has the slots to satisfy it.
_ROUTES: tuple[tuple[str, str, float, tuple[str, ...], str], ...] = (
    *((pool, prec, 1.0, (), "D-75 oxidative Strecker") for pool, _f, prec in _STRECKER_ROUTES),
    *((pool, prec, 1.0, (), "D-87 thermal Strecker") for pool, _m, _w, prec in _MAILLARD_PRODUCTS),
    *(
        (spec.species, spec.precursor_amino_acid, 1.0, (), "D-33 Ehrlich re-route")
        for spec in FUSEL_SPECS
    ),
    # Demethiolation: 1 methionine -> 1 methanethiol + 1 alpha-ketobutyrate + NH3. No CO2 (nothing
    # decarboxylates), but the C4 keto-acid is now charged — which is the whole of the D-107 fix.
    ("methanethiol", "methionine", 0.0, ("alpha_ketobutyrate",), "D-45 demethiolation"),
)

#: The routes whose carbon-sized draw is **knowingly not** their molar stoichiometry, each with the
#: carbon that makes it so (decision D-105). This is an allow-list, not a waiver: a route may only
#: sit here with a reason, and anything NEW that drifts fails the test below rather than joining
#: them silently.
#:
#: **IT IS EMPTY (decision D-107), and that is the beat's win condition.** D-105 wrote that this
#: list
#: "is exactly the keto-acid node's work-list, and nothing else" — both remaining entries were
#: blocked on one missing molecule, 2-ketobutyrate, which the mercaptan **produced and discarded**
#: while sotolon **consumed it and invented it from sugar**. Building that pool
#: (:mod:`~fermentation.core.kinetics.keto_acids`) closed both from opposite sides:
#:
#: * ``("methanethiol", "methionine")`` — the 5× under-draw: fixed. The draw is 1 mol methionine per
#:   mol thiol and the C4 goes to ``alpha_ketobutyrate``, so ``5 == 1 + 4`` and the route charges
#:   every carbon it consumes.
#: * ``("sotolon", "threonine")`` — the 1.5× over-draw: **dissolved rather than fixed**, and the
#:   distinction matters. Sotolon does not appear in :data:`_ROUTES` any more because it is not a
#:   carbon-sized draw off an amino acid at all: it is an aldol of two tracked pools
#:   (:class:`~fermentation.core.kinetics.aging.SotolonAldolCondensation`), drawing 1 mol of each
#:   substrate because that is what is written. **A route with no carbon-sized draw has no D-105
#:   blind spot to check** — the ratio it was failing was an artifact of asking a Strecker question
#:   about a molecule that was never a Strecker product. Its replacement is pinned by
#:   ``test_the_sotolon_aldol_draws_one_mole_of_each_substrate_when_driven``, which is a *stronger*
#:   test: it reads the real stoichiometry off ``dy/dt`` rather than off a declared table.
#:
#: **Keep it empty.** The stale-waiver arm below has nothing to guard now, so the *only* thing
#: standing between a new carbon-sized draw and a silent wrong mole count is the empty dict and the
#: driven tests. An entry added here must carry the carbon that explains it and a reason to believe
#: it cannot be charged — which, as D-107 demonstrates, usually means "the pool does not exist yet",
#: i.e. a work-list item rather than a permanent exemption.
_KNOWN_NON_STOICHIOMETRIC: dict[tuple[str, str], str] = {}


def test_a_carbon_sized_draw_equals_real_stoichiometry_only_where_it_charges_the_co2():
    """The signature that separates a true degradation from a carbon-sized stand-in (D-105).

    Every draw in this tree closes the carbon ledger **by construction** — the mass is sized so C
    out of the pool equals C into the products. That makes conservation blind here: a draw can
    conserve carbon perfectly while consuming the wrong number of moles of precursor, and no
    conservation test can tell. This is the test that can.

    A product that is a genuine degradation of its precursor satisfies ``C(precursor) ==
    C(product) + C(CO2 charged)`` — every precursor carbon is accounted for, so sizing the draw by
    carbon lands on **exactly 1 mol precursor per mol product**, the real stoichiometry, with no
    freedom left over. Where that identity fails the draw is a stand-in: the route is either
    sourcing carbon it does not name (sotolon's acetaldehyde) or discarding carbon it should charge
    for (the mercaptan's 2-oxobutyrate, the Ehrlich CO2).

    **This is the D-104 error class, made mechanical.** D-104 found sotolon rooted in a must amino
    acid where reality uses a mostly-de-novo keto acid, and the audit of the D-75 oxidative route
    it prescribed (D-105) found the same signature already sitting on the carbon ledger — visible
    with no literature at all. D-75 and D-87's five true Strecker aldehydes pass **exactly**; they
    are the only routes in the tree that charge their decarboxylation CO2, and that term is why.
    """
    m_c = MOLAR_MASS["CO2"] * carbon_mass_fraction("CO2")  # g C per mol carbon

    def carbons(species: str) -> float:
        return MOLAR_MASS[species] * carbon_mass_fraction(species) / m_c

    for product, precursor, n_co2, co_products, route in _ROUTES:
        charged = carbons(product) + n_co2 + sum(carbons(c) for c in co_products)
        implied = charged / carbons(precursor)  # mol precursor / mol product
        known = _KNOWN_NON_STOICHIOMETRIC.get((product, precursor))
        if known is None:
            assert implied == pytest.approx(1.0, abs=1e-12), (
                f"{route}: {precursor} -> {product} charges {n_co2:g} CO2 and co-products "
                f"{co_products or '()'}, implying {implied:.4f} mol {precursor} per mol {product} "
                f"where a true degradation demands 1.0. Either this route is not a degradation of "
                f"{precursor} (the D-104 error -- it needs a de-novo/keto-acid source, not a hard "
                f"gate on the pool), or it is failing to charge carbon it really releases: its "
                f"CO2, or a co-product pool. Fix it, or add it to _KNOWN_NON_STOICHIOMETRIC with "
                f"the carbon that explains it."
            )
        else:
            assert implied != pytest.approx(1.0, abs=1e-12), (
                f"{route}: {precursor} -> {product} is listed in _KNOWN_NON_STOICHIOMETRIC but "
                f"now draws stoichiometrically ({implied:.4f}). If it was fixed, delete the "
                f"entry -- a stale waiver hides the next regression.\nRecorded reason: {known}"
            )


def test_the_known_non_stoichiometric_allow_list_is_empty():
    """The keto-acid node's work-list, closed (decision D-107).

    D-105 found the signature above and left two routes on the allow-list, writing that it "is now
    exactly the keto-acid node's work-list, and nothing else". Both were blocked on **one missing
    molecule**: the mercaptan *produced* 2-ketobutyrate and discarded it (under-drawing methionine
    5×), while sotolon *consumed* it and invented it from sugar (over-drawing threonine 1.5×) —
    producer and consumer of the same untracked pool, in the same wine, on the same aging phase.
    D-107 built the pool and both closed, from opposite sides.

    This test exists so the list cannot quietly refill. It is not the same claim as the loop above:
    that one permits an entry *with a reason*, and this one says the current correct number of
    reasons is zero. If a future beat needs to add one, delete this test **deliberately** and say
    why in DECISIONS — do not let it happen as a side effect.
    """
    assert _KNOWN_NON_STOICHIOMETRIC == {}, (
        "The keto-acid node's work-list was emptied at D-107. A new entry means a route is "
        "knowingly drawing the wrong number of moles while conservation stays green -- which is "
        "exactly the D-104/D-105 error class. Justify it in DECISIONS, or fix the route."
    )


@pytest.fixture
def aging_params():
    # The Strecker routes read their rate constants from three files; the pool spectrum lives in
    # wine_generic. Driving a Process needs all of them resolved together.
    return {
        **load_parameters(default_data_dir() / "wine_generic.yaml").resolve(),
        **load_parameters(default_data_dir() / "aging.yaml").resolve(),
        **load_parameters(default_data_dir() / "thermal.yaml").resolve(),
    }


def _driveable_state(schema: StateSchema, params: Mapping[str, float]) -> FloatArray:
    # A warm, oxygenated, mid-ferment-composition state: enough for BOTH the oxidative (needs o2)
    # and the thermal (needs sugar) routes to fire on one seeded pool.
    y = schema.zeros()
    y[schema.slice("T")] = 298.15
    y[schema.slice("S")] = 100.0
    y[schema.slice("X")] = 2.0
    y[schema.slice("E")] = 50.0
    y[schema.slice("o2")] = 8.0e-3
    y[schema.slice("N")] = 0.1
    seed_amino_acids(y, schema, params, 1.0)
    return y


#: The seven routes whose draw ratio is **gate-independent** — the Process gates its product and its
#: precursor draw together, so the gate cancels and the measured ratio is a structural constant that
#: one seeded state pins exactly (decision D-105). The Ehrlich re-route is not here for a different
#: reason than de-novo sotolon: it never touches ``fusels`` (production stays in the producer), so
#: there is no product rate on its own ``dy/dt`` to divide by. **D-106 gives it one anyway** — its
#: CO₂, one mole per alcohol re-sourced, which is what
#: :func:`test_the_ehrlich_reroute_charges_one_co2_per_precursor_mole_when_driven` reads.
_EXACT_DRIVEN_ROUTES: tuple[tuple[str, str, str], ...] = (
    ("strecker", "methional", "methionine"),
    ("strecker", "phenylacetaldehyde", "phenylalanine"),
    ("maillard", "methional", "methionine"),
    ("maillard", "phenylacetaldehyde", "phenylalanine"),
    ("maillard", "2_methylbutanal", "isoleucine"),
    ("maillard", "3_methylbutanal", "leucine"),
    ("maillard", "2_methylpropanal", "valine"),
)


@pytest.mark.parametrize(("which", "product", "precursor"), _EXACT_DRIVEN_ROUTES)
def test_every_true_strecker_route_draws_1_to_1_when_the_process_is_actually_driven(
    schema, aging_params, which, product, precursor
):
    """The same claim, but **against the code** rather than the route table (decision D-105).

    The arithmetic test above reads ``n_co2`` out of ``_ROUTES`` — a literal *this file* declares.
    That pins the **declaration** layer and nothing else: delete the CO₂ term from
    :meth:`StreckerDegradation.derivatives` and methional starts drawing 0.8 mol methionine per mol
    methional — **the D-104 error, reintroduced** — while carbon still closes (4 C out, 4 C in) and
    the table still *says* ``n_co2=1.0``. Conservation would stay green and so would that test.
    **Both blind to it.** This test drives the Process and reads the debit off ``dy/dt``, so the
    deletion fails here — which is the only place it fails.

    The distinction is the beat's own lesson turned on its own tripwire: **a calculation that can
    only produce the answer you expect is not a check**, and `4 + 1 == 5` is that calculation.
    """
    proc = StreckerDegradation() if which == "strecker" else MaillardStrecker()
    y = _driveable_state(schema, aging_params)
    d = proc.derivatives(0.0, y, schema, aging_params)
    product_rate = float(d[schema.slice(product)][0])
    precursor_rate = -float(d[schema.slice(precursor)][0])
    assert product_rate > 0.0, "the route must actually fire, or this asserts nothing"
    assert precursor_rate > 0.0
    mol_product = product_rate / MOLAR_MASS[product]
    mol_precursor = precursor_rate / MOLAR_MASS[precursor]
    # A true Strecker degradation consumes exactly one precursor per aldehyde: the aldehyde IS the
    # amino acid minus its carboxyl (charged as CO2) and its amino group (deaminated to N).
    assert mol_precursor / mol_product == pytest.approx(1.0, abs=1e-9)


def test_the_ehrlich_reroute_charges_one_co2_per_precursor_mole_when_driven(schema, aging_params):
    """The re-route's draw is 1:1 **against the code**, not the table (decision D-106).

    D-105 measured this route at ``(n-1)/n``: it charged the precursor for the alcohol's carbon but
    not for the CO2 the same decarboxylation releases, so leucine went out at 5/6 mol per mol
    isoamyl alcohol. D-106 charges it. The claim is now that **every mole of precursor drawn emits
    exactly one mole of CO2** — the Ehrlich decarboxylation, one per alcohol.

    This route cannot be checked the way the seven Strecker routes are: it never touches ``fusels``
    (production stays in :class:`FuselAlcoholsEhrlich`, which is the whole D-33 swap design), so
    there is no product rate on its ``dy/dt``. Its CO2 **is** the product rate — one mole per mole
    of alcohol re-sourced — and it is gate-independent for the same reason the others are: the
    draw and the CO2 both scale with ``gate x fusel_carbon``, so the gate cancels and any seeded
    state pins the ratio exactly.

    Deleting the CO2 term fails **here**: the precursor debit drops by 1/n, the moles stop matching,
    and — exactly as at D-105 — carbon still closes and no conservation test notices.
    """
    y = _driveable_state(schema, aging_params)
    d = FuselAminoAcidReroute().derivatives(0.0, y, schema, aging_params)
    co2_mol = float(d[schema.slice("CO2")][0]) / MOLAR_MASS["CO2"]
    precursor_mol = sum(
        -float(d[schema.slice(spec.precursor_amino_acid)][0])
        / MOLAR_MASS[spec.precursor_amino_acid]
        for spec in FUSEL_SPECS
    )
    assert precursor_mol > 0.0, "the re-route must actually fire, or this asserts nothing"
    assert co2_mol > 0.0, "no CO2 emitted -- the D-106 decarboxylation term is gone"
    # One CO2 per precursor consumed: the Ehrlich pathway decarboxylates exactly once per alcohol.
    assert co2_mol / precursor_mol == pytest.approx(1.0, abs=1e-9)


def test_the_d75_oxidative_strecker_routes_draw_at_exact_stoichiometry():
    """The D-104-prescribed audit of D-75, pinned as the acquittal it produced (D-105).

    D-104 left the oxidative Strecker route un-audited for its own error while noting that D-75
    shares methional and phenylacetaldehyde with the thermal route. It does not have it: methional
    IS methionine minus its carboxyl and amino groups (C5 -> C4 + CO2), phenylacetaldehyde IS
    phenylalanine's (C9 -> C8 + CO2). Both charge that CO2, so both land on 1.0000 mol precursor
    per mol aldehyde and no de-novo share is even expressible -- the mass balance forbids it. The
    hard gate on the amino-acid pool is therefore the CORRECT topology here, not the sotolon
    mistake, and this pins that rather than leaving it to be re-litigated.
    """
    m_c = MOLAR_MASS["CO2"] * carbon_mass_fraction("CO2")
    for pool, _fparam, precursor in _STRECKER_ROUTES:
        c_product = MOLAR_MASS[pool] * carbon_mass_fraction(pool) / m_c
        c_precursor = MOLAR_MASS[precursor] * carbon_mass_fraction(precursor) / m_c
        # the whole claim in one line: aldehyde + its carboxyl CO2 == the amino acid, exactly
        assert c_product + _CO2_PER_STRECKER_ALDEHYDE == pytest.approx(c_precursor, abs=1e-9)


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
