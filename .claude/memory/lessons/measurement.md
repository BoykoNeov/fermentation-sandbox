---
name: lessons-measurement
description: "Reading a number off the model and deciding what caused it: baselines, controls, censuses, sampling grids, summary statistics and null results"
metadata:
  node_type: memory
  type: feedback
---

**Lessons — Measuring a run & attributing the number.** Split out of `.claude/memory/MEMORY.md` on 2026-08-26; that file's
index points here by path. These rows carry **no `MEMORY.md` row of their own**, so they cost
nothing until this file is read — the same arrangement `prohibitions/` has had since D-185.

**Read this file before you write the code, not after the review.** Each row is *the trap* +
*what to do instead*; the measurement that earned it is in the linked file. If a row looks
like pedantry, open its file — every one of them cost a beat.

- [Compute the clean fix before adopting it](../feedback-compute-the-clean-fix-before-adopting-it.md) — the structurally-cleanest candidate came out 3.71× off where the shipped one was 0.67×; and a consequential number needs its own baseline before you attribute it to your mechanism
- [Measure which side before building](../feedback-measure-which-side-before-building.md) — a one-directional corrective only helps in one direction; check the sign, the reach, and whether it's a rate knob on a supply-limited quantity
- [Count and print your skips](../feedback-count-and-print-your-skips.md) — a harness that silently drops what it can't parse reports "5 of 5 clean" on a denominator it never measured; D-157's live defect was in the dropped two
- [A majority is not a direction](../feedback-a-majority-is-not-a-direction.md) — 6-of-8 seeds "widened" flipped to 10-of-24 when the seeds tripled; fix the stopping rule before the run lands, and label rows with the run that produced them
- [Pre-register the cheap prediction](../feedback-pre-register-the-cheap-prediction.md) — write the free static pass to disk before the expensive campaign; it's what makes a miss diagnosable and a surprise claimable
- [Name the field your predicate read](../feedback-name-the-field-your-predicate-read.md) — a census counts whatever FIELD the predicate touched; D-163 tested `provenance.source` and reported it about band EDGES, turning 110 citation-backed edges into 4
- [A screen isn't idempotent under its own repair](../feedback-a-screen-is-not-idempotent-under-its-own-repair.md) — repairing a note made it score WORSE, so 15→15 hid four fixes; run the denominator on the pre-repair tree and diff MEMBERSHIP, not the count
- [A tautology can smuggle an attribution](../feedback-a-tautology-can-smuggle-an-attribution.md) — "outside the band either way" conceded a caveat then declared it non-binding; the control landed 6 figures from the mutation, so the SHAPE did nothing and the finding got STRONGER (D-195)
- [A magnitude defence hides in its smallest member](../feedback-a-magnitude-defence-can-be-scoped-to-the-smallest-member.md) — "only ~2 %" was TRUE and about the smallest of FIVE terms; the lump was 62-76%. Enumerate from the REGISTRY, and price removal against data first (D-197)
- [A hit can be two errors cancelling](../feedback-a-hit-can-be-two-errors-cancelling.md) — a prediction that LANDS is evidence only if its inputs were; re-derive a hit
- [A summary statistic is not the curve](../feedback-a-summary-statistic-is-not-the-curve.md) — an endpoint metric passed at 87 % while day 1 sat 8.1× outside the band; an OUT call needs a MEASURED tolerance
- [A scope note can size the wrong half](../feedback-a-scope-note-can-size-the-wrong-half.md) — D-209 named two omissions and sized the one moving a dry endpoint by **exactly 0.0**; the unsized half owned all of it, and they OPPOSE at the dose. Ask which one the observable depends on at the HORIZON you report (D-210)
- [A margin can be borrowed from a defect](../feedback-a-margin-can-be-borrowed-from-a-defect.md) — D-183 was scored while growth ran 2.88× too fast; re-pin the number, never relax the threshold
- [A uniform rescale can't test composition](../feedback-a-uniform-rescale-cannot-test-a-composition-claim.md) — D-246 scaled ONE constant and read the gain as a composition defect; that probe is bit-identical to scaling every share, so nothing moved relative. The real one: −6.9 % of the gap, wrong way (D-247)
- [A pool size is not a flux](../feedback-a-pool-size-is-not-a-flux.md) — "pool 32.2 vs 47.7 needed, so it dies on mass" was UNSOUND for an intermediate: throughput was ~3038 and the gap needed 1.57 %. A refusal's value is its enumeration of causes, and a wrong cause is what a later beat overturns
- [Read a fast curve on a fixed grid](../feedback-read-a-fast-curve-on-a-fixed-grid.md) — `argmin` over adaptive output landed 22 min early on the steepest limb: 0.034 vs 0.0274, a 26 % error, stable enough to look measured. Error scales with LOCAL SLOPE, so it hides everywhere but the number you want (D-214)
- [A null needs the predicate the CLAIM makes](../feedback-a-null-needs-the-predicate-the-claim-makes.md) — 75 hits over 10 patterns looked exhaustive, but the claim also said "no two-temperature rate pair", a shape NONE of them match; a second census found 28 more hits. Count what the sentence asserts
- [Two measured quantities don't locate a model defect](../feedback-two-measured-quantities-do-not-locate-a-model-defect.md) — I diagnosed "the model front-loads" from two SOURCE curves without running the model; it tracks its own, 2.8× slower, driver. If a diagnosis needs no model run, it's about the literature
- [A baseline log goes stale between edits](../feedback-a-baseline-log-goes-stale-between-edits.md) — diffing 16→15 against a log from a different tree credited another edit to my fix, which repaired nothing; diff NAMES, re-run the baseline, and treat an invisible fix as one owing a guard
- [Reproduce a published number first](../feedback-reproduce-a-published-number-before-trusting-the-new-column.md) — my A/B harness ran ONE tree's code twice (`sys.path` can't un-import) and printed a plausible WRONG table; the baseline column failing to reproduce D-226's 0.9549 caught it. One process per tree
- [A limitation can belong to its frame](../feedback-a-limitation-can-belong-to-its-frame.md) — "growth stops dead at day 1" rode 3 surfaces as the next thing to BUILD; it was the test frame's unsourced pitch, and the frame with data already agreed
- [A discriminator can flip, not just go inert](../feedback-a-discriminator-can-flip-across-a-fitted-band.md) — the day1-vs-day3 settling test REVERSED its verdict between `mu_max`'s band edges (0.348 span at day 1, 0.003 at day 3). Inert = no information; flipping = information about the FIT. Run both edges
- [A one-parameter sweep is not the band](../feedback-a-one-parameter-sweep-is-not-the-band.md) — swept ONE parameter, got 1.287x, wrote it as "the band is 29 % wider"; over all 83 it is 1.008x — wrong 36x, flattering the beat. Vary what the SAMPLER varies, or scope the claim
- [Size a set with an instrument that can see it](../feedback-size-a-set-with-the-instrument-that-can-see-it.md) — sized an unenumerated set with the grep whose blind spot CREATED it: predicted 12-20, got 32, and the 21 invisible names held the whole prior defect. Different KIND of instrument, or a floor (D-234)
- [A blast radius counted at definitions misses call sites](../feedback-a-blast-radius-counted-at-definitions-misses-call-sites.md) — priced a signature widening at 13 defs + 13 handovers; the break landed on 28 CALL sites in tests no grep matched. Let the type checker enumerate it, or quote a floor (D-235)
- [An expectation transfers with the CHANNEL, not the subject](../feedback-an-expectation-transfers-with-the-channel-not-the-subject.md) — "the cascade will widen copper's artefact" was derived on the O₂ channel and reported on browning; it moves twice the O₂ and LESS browning, so the direction was wrong, not the size. Name the observable an expectation was derived on (D-234→D-237)
- [Pair the arms before comparing spreads](../feedback-pair-the-arms-before-comparing-spreads.md) — two ensembles over different `only=` sets draw different MEMBERS; an unpaired run reported a 1.672x widening and a 0.42x narrowing in the same table. Paired on identical draws it is 0.323 (D-240)
- [Price the refusal, not just the error](../feedback-price-the-refusal-not-just-the-error.md) — a banded Parameter can REFUSE a candidate repair outright; "wrong" is a trade-off, "does not run" is a one-way decision, and the same band usually blocks a legitimate case worth finding
- [One knob beats an unbuildable control](../feedback-one-knob-decomposition-beats-an-unbuildable-control.md) — the commensurate fixture needed a source we lack; sweeping ONLY the parameter the change resolved made the confound a constant. Report what did NOT move
- [Answering a sourcing ask can close nothing](../feedback-answering-a-sourcing-ask-can-close-nothing.md) — "one ask gates four xfails" was fully answerable and closed none; two got WORSE on their own source's must. Preregister which should move, and anchor a reimplemented harness on the old fixture first
