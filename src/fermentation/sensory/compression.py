"""Stevens compression — per-compound odor-intensity curves (decision D-98, beat 1b slice 2).

Slice 2 of the handoff §4.2 descriptor beat, deferred at D-95 and built here **narrower than
D-95 scoped it**. D-95 bundled four items under "slice 2": weights, a compression exponent,
masking/suppression, and matrix effects. Three of those do not survive contact with the
literature, and saying so is most of this module's content:

* **Matrix effects — ALREADY DONE, one layer down.** Beat 1a's thresholds are already
  matrix-specific (``threshold_<pool>_wine`` vs ``_beer``), which is where a matrix effect
  belongs: the reported water→wine threshold shifts span factors of ~2 to ~23000. Nothing to
  build.
* **Weights — NOT A THING.** D-95 cited ``thermal.yaml``'s relative weights as the precedent,
  but those are *production-flux* weights and that file says outright that "perceptibility
  lives in threshold_phenylacetaldehyde_wine". A perceptual weight is just the coefficient
  ``k`` in ``I = k * C**n``, and calibrating on the threshold fixes ``k`` (see below). There is
  no third quantity to source.
* **Masking / suppression — BLOCKED, and deliberately not built.** The vector model needs a
  per-pair interaction coefficient ``cosα``, obtained by putting panels in front of pairs; the
  one published shortcut (a constant ``cosα ≈ -0.129``) is explicitly scoped to compounds "with
  similar functional groups", and this medium's pools span esters, aldehydes, ketones, thiols,
  phenols, furanones and lactones. For 21 pools that is 210 pair coefficients nobody has
  measured. Inventing them is not on the table.

What is left — and what this module is — is the **per-compound compression curve**.

**The additivity through-line survives intact, and this is the load-bearing claim.** D-67
refused to sum OAVs; D-95's max rule refused again; the plan requires that any weighted or
compressed rule "answer that objection head-on before it ships". Here is the answer:
compression is applied **per compound, before the combination rule**, and the combination rule
is **still max**. ``I_i = OAV_i ** n_i`` is a monotone transform of one pool's own OAV — it
says nothing whatever about how two compounds combine. :class:`StevensProjector` therefore
asserts no additivity, exactly as :class:`~fermentation.sensory.descriptors.MaxRuleProjector`
does not. *We never assume additivity, at any layer* — still true after slice 2.

That max stays is not merely convenient. Mixture perception is generally **hypoadditive**: the
true combined intensity sits *between* max and sum (with suppression sometimes pushing it below
even the stronger component). So max is a rough **lower bound** — an honest under-claim, which
is what a speculative layer owes. The truthful middle is precisely the region that needs the
``cosα`` coefficients that do not exist, so the honest rule is the bound, not a guess at the
middle.

**Why the threshold cannot calibrate these curves.** The natural idea — calibrate each curve on
the compound's measured detection threshold, which beat 1a already holds for all 33 keys — is
provably empty. ``I = k * C**n`` has two unknowns; a threshold ``T`` is one point, and by
definition the point where every compound is equally just-detectable (``I(T) == I_thr``). So
``k = I_thr / T**n`` and ``I(C) = I_thr * (C/T)**n = I_thr * OAV**n``: the threshold pins ``k``
**given** ``n`` and never pins ``n``. It reproduces OAV — which beat 1a already computes — and
leaves the curve's *shape* untouched. A threshold says **where** a smell becomes detectable and
is silent on **how fast** it grows after that. With ``I_thr`` normalized to 1 this layer's
intensity is therefore exactly ``OAV ** n``, and ``n`` is the part no measurement here supplies.

**Why the exponents are author estimates, and why a global one is not offered.** The wine-aroma
literature states plainly that the exponent is unknown for most compounds and that the common
workaround is to assign a global value of one — which *is* OAV, i.e. what beat 1a already does.
A **global** exponent is a provable **no-op** on every observable this layer reports: a monotone
transform preserves argmax, so ``dominant`` is unchanged, and ``I > 1`` iff ``OAV > 1``, so
``above_threshold`` is unchanged. It would mint a parameter and change nothing. **Per-compound**
exponents can flip ``dominant``, and that flip is this module's entire payload — so the numbers
are author estimates (``psychophysics.yaml``, ordered by aqueous solubility per Cain 1969, whose
citation justifies the *direction* and rough *scale* of the spread and **not** any single value).

**Which makes the sensitivity pass the real deliverable, not the projector.** A guess that
changes the answer is worth nothing on its own; a guess whose *consequences you have mapped* is
worth something. :func:`dominant_flip_sensitivity` samples the exponents across their honest
uncertainty bands and reports, per descriptor axis, whether ``dominant`` is stable or whether
the answer is an artefact of the guess. A knife-edge flip is not a finding about wine — it is a
finding about ``psychophysics.yaml``, and it is reported as one. See D-98.

**Isolation.** ``psychophysics.yaml`` loads standalone here (never through the compile seam —
no RHS reads an exponent), and slice 1 neither imports nor knows about this module. Delete
this file and beat 1a + slice 1 are byte-for-byte unaffected (prime directive #3).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from fermentation.parameters.store import ParameterSet, default_data_dir, load_parameters
from fermentation.sensory.descriptors import (
    DescriptorProfile,
    DescriptorReading,
    axes_for_medium,
    descriptor_tier,
)
from fermentation.sensory.oav import SensoryProfile


def load_exponents() -> ParameterSet:
    """Load ``psychophysics.yaml`` standalone — the per-compound Stevens exponents.

    Deliberately **not** routed through the compile seam's ``shared_files`` (no RHS reads an
    exponent), mirroring :func:`~fermentation.sensory.oav.load_thresholds`. Kept out of
    ``sensory.yaml`` for two reasons: a threshold is a *measured* quantity with a real
    citation and these are 21 author estimates that must not borrow its credibility by
    adjacency (§4.3); and slice 2 must stay togglable off (prime directive #3).

    NB (D-24): like the thresholds, these load standalone and therefore sit **outside** the
    ensemble parameter sweep — ``simulate_ensemble`` will not propagate these bands. That is
    what makes :func:`dominant_flip_sensitivity` a *manual* Monte Carlo rather than a free
    consequence of the existing machinery.
    """
    return load_parameters(default_data_dir() / "psychophysics.yaml")


def _exponent_key(pool: str) -> str:
    """The exponent parameter name for ``pool``, e.g. ``stevens_n_diacetyl``.

    Matrix-**agnostic**, unlike :func:`~fermentation.sensory.oav._threshold_key`. Exponents
    are not measured for these compounds in *any* matrix, so minting ``_wine``/``_beer``
    variants would double the guesses for exactly zero information. The matrix-dependence of
    compression is unmodelled and flagged in every entry's ``conditions``.
    """
    return f"stevens_n_{pool}"


def compressed_intensity(oav: float, exponent: float) -> float:
    """Perceived intensity in threshold units: ``OAV ** n``.

    The whole of Stevens compression, once the threshold has pinned ``k`` (see the module
    docstring): ``I = I_thr * OAV**n`` with ``I_thr`` normalized to 1. Properties that matter
    and are pinned by tests:

    * **1 at threshold** for every ``n`` — so ``I > 1`` iff ``OAV > 1``, and ``above_threshold``
      carries the identical meaning under compression. Compression can never invent or silence
      a detectable smell; it only re-scales how loud a detectable one is.
    * **Monotone increasing** in ``oav`` for ``n > 0`` — more compound never smells like less.
    * **0 at 0** — a clean run raises no false aroma.
    * **Compressive** for ``n < 1``: above threshold it pulls values *down* toward 1, below
      threshold it pulls them *up* toward 1. Both are the same statement — the nose squashes
      ratios — and both are why an OAV of 761 was never a claim about loudness.
    """
    return float(oav**exponent)


class StevensProjector:
    """Descriptor projector applying a **per-compound** compression curve, then **max**.

    The slice-2 counterpart to :class:`~fermentation.sensory.descriptors.MaxRuleProjector`,
    and a drop-in at the same :class:`~fermentation.sensory.descriptors.DescriptorProjector`
    seam — which is what that seam was built for. It differs in exactly one respect: each
    contributor's OAV is first mapped to a perceived intensity ``OAV ** n`` with ``n`` read
    per-compound from ``psychophysics.yaml``, *then* the axis takes the max as before. The
    combination rule is untouched, so no additivity is assumed (module docstring).

    **What it can change, and what it cannot.** It cannot change ``above_threshold`` — ``I > 1``
    iff ``OAV > 1`` for any ``n > 0``. The *only* new observable is ``dominant``: because
    exponents differ per compound, the loudest-by-OAV compound need not be the
    loudest-by-intensity one. That is the documented critique of OAV made executable (a
    compound can be detectable at a lower concentration than another and still be the weaker
    smell at realistic levels), and it is this projector's entire reason to exist.

    **Reading its output honestly.** :attr:`DescriptorReading.oav` holds a compressed intensity
    here, not an OAV, and :attr:`DescriptorReading.rule` says so (``"stevens"``). Every number
    it reports rides on 21 author estimates; consult :func:`dominant_flip_sensitivity` before
    believing any ``dominant`` it names. ``MaxRuleProjector`` remains the default precisely
    because its claim is weaker.
    """

    #: Names the rule in every :class:`DescriptorReading` this projector emits, so a reading
    #: is self-identifying rather than depending on the caller to recall which seam it used.
    rule: str = "stevens"

    def __init__(self, exponents: ParameterSet) -> None:
        self._exponents = exponents

    def exponent_of(self, pool: str) -> float:
        """The Stevens exponent for ``pool``, or a clear error naming the missing key."""
        key = _exponent_key(pool)
        try:
            return self._exponents.value(key)
        except KeyError as exc:  # pragma: no cover - defensive, pinned by a test
            raise ValueError(
                f"no Stevens exponent for aroma pool {pool!r}: psychophysics.yaml has no "
                f"{key!r}. Every aroma pool of every medium needs one — see D-98."
            ) from exc

    def project(self, profile: SensoryProfile) -> DescriptorProfile:
        readings: dict[str, DescriptorReading] = {}
        for axis in axes_for_medium(profile.medium):
            contributors = {
                pool: compressed_intensity(profile.readings[pool].oav, self.exponent_of(pool))
                for pool in axis.pools
            }
            dominant = max(contributors, key=lambda p: contributors[p])
            intensity = contributors[dominant]
            readings[axis.name] = DescriptorReading(
                descriptor=axis.name,
                contributors=contributors,
                dominant=dominant,
                oav=intensity,
                # identical to `profile.readings[dominant].oav > 1.0` — compression is
                # threshold-preserving, and a test pins that equivalence rather than trusting it.
                above_threshold=intensity > 1.0,
                lumped=profile.readings[dominant].lumped,
                tier=descriptor_tier(profile.readings[p].tier for p in axis.pools),
                rule=self.rule,
            )
        return DescriptorProfile(
            medium=profile.medium, time_index=profile.time_index, readings=readings
        )


@dataclass(frozen=True)
class FlipVerdict:
    """Whether one descriptor axis' ``dominant`` survives the exponents' uncertainty bands.

    ``nominal`` is the compound named at the file's face values; ``share`` maps each compound
    to the fraction of Monte-Carlo draws in which it won. ``robust`` is True when ``nominal``
    won *every* draw **and the axis actually smells of anything** — i.e. the attribution does
    not depend on the guess anywhere in its honest band. Anything else is knife-edge, and its
    ``dominant`` is a statement about ``psychophysics.yaml``, not about the drink.

    **``silent`` is not a detail — it is the guard against a lie.** An axis whose contributors
    are all zero (an un-oaked wine's ``vanilla_oak``) has no dominant compound at all: since
    ``0 ** n == 0`` for every ``n``, every draw ties, the tie breaks to the first-listed pool,
    and a naive reading would report "vanillin, wins every draw" — maximal confidence in an
    aroma that does not exist. That is the "clean run raises no false descriptor" sin slice 1
    tests against, arriving one layer up wearing a statistic. A silent axis is reported silent
    and is never called robust.
    """

    descriptor: str
    nominal: str
    share: Mapping[str, float]
    silent: bool
    _unanimous: bool

    @property
    def robust(self) -> bool:
        """The attribution holds across the whole band **and** there is an aroma to attribute."""
        return self._unanimous and not self.silent

    @property
    def contested(self) -> bool:
        """True when more than one compound wins somewhere in the band."""
        return len([c for c, s in self.share.items() if s > 0.0]) > 1

    def summary(self) -> str:
        """A one-line human verdict, for the DECISIONS record and reports."""
        if self.silent:
            return f"{self.descriptor}: silent (no contributor present — nothing to attribute)"
        if not self.contested:
            return f"{self.descriptor}: {self.nominal} (robust — wins every draw)"
        ranked = sorted(self.share.items(), key=lambda kv: -kv[1])
        parts = ", ".join(f"{c} {s:.0%}" for c, s in ranked if s > 0.0)
        return f"{self.descriptor}: CONTESTED — {parts}"


def dominant_flip_sensitivity(
    profile: SensoryProfile,
    exponents: ParameterSet,
    *,
    draws: int = 4000,
    seed: int = 0,
) -> dict[str, FlipVerdict]:
    """Does each axis' ``dominant`` survive the exponents' uncertainty bands? (D-98)

    **The actual deliverable of slice 2.** The projector reports a ``dominant`` computed from 21
    author estimates; on its own that number is worth nothing, because a guess that changes the
    answer is indistinguishable from a guess that fabricates it. This maps the consequence:
    sample every exponent independently and uniformly across its ``psychophysics.yaml``
    uncertainty band, re-project, and count how often each compound wins its axis.

    Sampling is **independent** per compound because the estimates *are* independent guesses —
    there is no covariance to model, and pretending otherwise would invent structure. Uniform
    rather than normal for the same reason: the bands are honest ignorance, not panel spreads,
    and a uniform draw declines to claim the centre is likelier than the edge.

    Read the verdicts as follows. ``robust`` means the attribution holds across the whole band,
    so it is real conditional knowledge — the OAV gap is wide enough that no plausible exponent
    ratio closes it. ``contested`` means the answer is an artefact of the guess and must be
    reported as "cannot say", never as a sensory claim. ``silent`` means the axis has no
    contributor present at all, so there is nothing to attribute and *neither* verdict applies
    (see :class:`FlipVerdict` — this case otherwise masquerades as perfect robustness).
    Single-contributor axes are trivially robust (nothing to flip) and are included for
    completeness.

    Deterministic in ``seed``. Only axes of ``profile``'s medium are covered.
    """
    rng = np.random.default_rng(seed)
    axes = axes_for_medium(profile.medium)
    pools = sorted({p for a in axes for p in a.pools})

    bands = {p: exponents[_exponent_key(p)].uncertainty for p in pools}
    samples = {p: rng.uniform(bands[p].low, bands[p].high, size=draws) for p in pools}

    nominal = StevensProjector(exponents).project(profile)
    wins: dict[str, dict[str, int]] = {a.name: dict.fromkeys(a.pools, 0) for a in axes}
    for d in range(draws):
        for axis in axes:
            best, best_i = axis.pools[0], -np.inf
            for pool in axis.pools:
                i = compressed_intensity(profile.readings[pool].oav, float(samples[pool][d]))
                if i > best_i:
                    best, best_i = pool, i
            wins[axis.name][best] += 1

    verdicts: dict[str, FlipVerdict] = {}
    for axis in axes:
        share = {p: wins[axis.name][p] / draws for p in axis.pools}
        nom = nominal.readings[axis.name].dominant
        verdicts[axis.name] = FlipVerdict(
            descriptor=axis.name,
            nominal=nom,
            share=share,
            # every contributor absent => no aroma to attribute; see FlipVerdict.silent
            silent=all(profile.readings[p].oav == 0.0 for p in axis.pools),
            _unanimous=share[nom] == 1.0,
        )
    return verdicts
