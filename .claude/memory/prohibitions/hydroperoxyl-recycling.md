# Hydroperoxyl recycling limb — REFUSED on measurement (D-198)

Detail for the ledger row in [[project-fermentation-sandbox]]. Reached BY PATH; no `MEMORY.md` row.

## What it is

`CH₃·C(OH)H + O₂ → CH₃CHO + HO₂·` — the ethanol limb of the Fenton branch takes an O₂ (built at
D-196) and leaves a hydroperoxyl radical behind. If that radical hands an oxidising equivalent
back, the model is charging too much oxygen for the step.

## Why D-197's expected refusal would have been a MISATTRIBUTION

D-197 §8 predicted the refusal would rest on *Understanding Wine Chemistry*'s "an attempt to
detect this species was not successful". **That sentence is about a different reaction.** The
undetected HO₂· is the one-electron intermediate of the **initiation** node
`Fe(II) + O₂ → Fe(III) + HO₂·` (§24.4.1, Fig. 24.9, index p. 327-328). The ethanol-limb
co-product is four pages later (index p. 332). Same chemical species, different role.

The model already declines the initiation route: `_H2O2_PER_O2` implements the two-Fe(II)
alternative UWC itself calls "more probable". So the failed detection was never evidence about
this limb at all.

## The four grounds for refusing

1. **No source for this limb's co-product.** 101 hits across 17 of 24 texts, skips counted; zero
   state the ethanol-limb fate. UWC Fig. 24.13 draws the step with **no co-product** — silent,
   not negative.

2. **Three sourced fates disagree**, all in other contexts: → H₂O₂ (Carrascón 2018, Marrufo 2018,
   Danilewicz 2007 scheme 8); dismutation `2 HO₂· → H₂O₂ + O₂` (beer texts); Haber-Weiss.

3. **A pole with no published terminator.** H₂O₂ is quasi-steady-state, so a returning limb feeds
   its own production: `F = A/(1−s)`. `s` = **exactly 1.0** unsulfited (measured) — diverges, and
   the arm is **not integrable**. Terminator = the hydroxycinnamate quench; UWC cites Gislason
   et al., on disk twice, both citations. D-196 §2's "absorbed into the calibrated rate" **does
   not extend**: the quench removes the returned H₂O₂ too.

4. **Nothing here can separate the candidates.** Routes 1 and 2 give the *identical* SO₂:O₂ ratio
   (`2 − 1.5s`); the hermetic O₂ budget spans **0.078 %** across all of them, that pool being
   supply-limited (D-196 §3). Even a sourced fate could not be told apart here.

Route 3 is additionally **not expressible**: Haber-Weiss consumes an H₂O₂ and regenerates the ·OH,
making it a *third consumer* of the H₂O₂ node. `h2o2_branch_fraction` partitions over exactly two
weights; a third needs `k_HW·[HO₂·]`, with no slot and no source.

## What D-198 corrected in D-196

`_O2_PER_ACETALDEHYDE = 1.0` was glossed as "an **UPPER BOUND** on the net draw". For the fate
three wine texts describe it is a **lower** bound: `A/(1−s)` vs the shipped `A(1+s)`, excess
`s²/(1−s)` at every `s > 0`. D-196's *other* §5 claim — that the limb pushes SO₂:O₂ the other way
— is **CONFIRMED**; D-198's own pre-registration predicted it false and was wrong, because it set
the recycled radical's quinone yield to zero. Derived from the module's `2 Fe(III) → o-quinone`
rule, `q = 0.5` and the crossover sits at `s = 1.0`, outside the interior.

## The guard, and why one was owed

`test_the_fenton_limb_returns_no_h2o2_to_the_node_it_drew_from` — asserts the limb's acetaldehyde
flux equals `share × activation`, the identity that holds only while the limb is a pure *consumer*
of the H₂O₂ node.

D-196's 1:1 test asserts a **ratio**, and every recycling term scales O₂ and acetaldehyde by the
same factor — **measured GREEN** under a capped route-1 arm. The uncapped arm's red is the
*fixture crashing on the pole*, a red that names nothing.
See [[feedback-a-ratio-guard-cannot-see-a-common-factor]].

## Re-opening it

Needs the **quench ratio**, not another sweep of the same 24 texts. And even with it, this
scenario cannot discriminate the fates — a discriminating measurement would have to be
non-hermetic (continuous O₂ ingress), where the pool is not supply-limited.

Receipts: `M:\claud_projects\temp\ferment\d198-hydroperoxyl\`.
