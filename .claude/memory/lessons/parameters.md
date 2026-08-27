---
name: lessons-parameters
description: "Where a constant's freedom lives: per-parameter bands vs joint claims, calibrated levels that decay, rate-law coordinates, and fitting vs scoring"
metadata:
  node_type: memory
  type: feedback
---

**Lessons — Parameters, bands, fits & identifiability.** Split out of `.claude/memory/MEMORY.md` on 2026-08-26; that file's
index points here by path. These rows carry **no `MEMORY.md` row of their own**, so they cost
nothing until this file is read — the same arrangement `prohibitions/` has had since D-185.

**Read this file before you write the code, not after the review.** Each row is *the trap* +
*what to do instead*; the measurement that earned it is in the linked file. If a row looks
like pedantry, open its file — every one of them cost a beat.

- [A margin is a claim about what holds it open](../feedback-a-margin-is-a-claim-about-what-holds-it-open.md) — "inert, 0.5 % margin" was incomplete: the Gay-Lussac ceiling already cleared the threshold, so the margin belonged to two byproduct-yield bands and a documented toggle closes it
- [A nominal on a band edge is not inertness](../feedback-nominal-on-a-band-edge-is-not-inertness.md) — value == an endpoint makes that edge's run the nominal run, bitwise identical by construction; a two-point screen then *invents* a straddle — classify over the band interior, and check per medium
- [Build the term that makes agreement worse first](../feedback-build-the-term-that-makes-agreement-worse-first.md) — two omitted terms of OPPOSITE sign: ship the one that costs you the headline, else the other gets credited for a cancellation (D-181: a corner that reached the measurement belonged to the error)
- [A derived yield encodes its rate law](../feedback-a-derived-yield-encodes-its-rate-law.md) — Δmeasured/Δdivisor is a yield only if production tracks that divisor; the assumption launders into `source:`, so what it IMPLIES looks measured. Cost 3 blind-alley designs defending a scaling nothing measured (D-183)
- [A shape change is a change to the pair](../feedback-a-shape-change-is-a-change-to-the-pair.md) — 82 records framed it as "the excretion shape" and never named the sink; moving one driver turned a residual into `exp(-k∫flux)`, 2.000→0.000. Measure the invariant the PAIR holds (D-189)
- [A band is per-parameter, a claim is joint](../feedback-a-band-is-per-parameter-a-claim-is-joint.md) — two notes each spanned the rival value, but it was ONE source's pair, flipping an ordering asserted as FACT — false in 2.0% of draws. Sample the joint space and COUNT; check a mutation is in-band (D-190)
- [A form too gradual for its own anchors](../feedback-a-form-can-be-too-gradual-for-its-own-anchors.md) — the cited form spans 1.51× where the source's two statements demand ~19×: no constant satisfies both, so the sweep was unanswerable. Ratio the group at each anchor BEFORE fitting (D-192)
- [A returned intermediate is a rate-law change](../feedback-a-returned-intermediate-is-a-rate-law-change.md) — feeding a quasi-steady-state node its own product gives `A/(1−s)`, so an "upper bound" was a LOWER bound; find the pole
- [Relocate the unsourced factor](../feedback-relocate-the-unsourced-factor.md) — a coefficient nothing sources constrains the algebraic FORM, not the beat; guard the property, not the trajectory
- [A pair constrains a response](../feedback-a-pair-constrains-a-response.md) — one statistic ± a species constrains a rate-free RESPONSE where one value constrains a level; broke a 4-beat "nothing adjudicates"
- [Count the anchors before adding a parameter](../feedback-count-the-anchors-before-adding-a-parameter.md) — the bar is one observable vs two free parameters; pair a refusal with the crossing value
- [A normalisation is a free parameter](../feedback-a-normalisation-is-a-free-parameter.md) — a ratio term hides its REFERENCE, so nothing fixed `pH_ref`; expressible ≠ identifiable
- [A parameter can be pinned and drawn](../feedback-a-parameter-can-be-pinned-and-drawn.md) — read at COMPILE and runtime, only the runtime half is sampled; pair against a RECOMPILED control
- [Negligible sensitivity isn't no freedom](../feedback-negligible-sensitivity-is-not-no-freedom.md) — buffering 2.1 % argued for one cheap constant, but its admissible range × a large concentration was 0.022 pH; run BOTH numbers — sensitivity says the physics is idle, span says the CHOICE isn't (D-210)
- [Fit the observable, not the consequence](../feedback-fit-the-observable-not-the-consequence.md) — the 8/8 value was REFUSED for the 7/8 one: a downstream observable sums your parameter and every unbuilt term, so fitting it hides their share in yours. Held out, its SIGN names the term (D-211)
- [A constant can be a unit definition](../feedback-a-constant-can-be-a-unit-definition.md) — two values 5.6× apart weren't rival estimates: one was BACK-COMPUTED (a residual), the other unsourced, and the answer was the unit the model's own fitted parameters live in — already on disk
- [A gap can be held open by a second anchor](../feedback-a-gap-can-be-held-open-by-a-second-anchor.md) — the knob closing a measured gap was IN BAND; a second anchor on the same knob forbade it, breaking after a fifth of the gap. Score every anchor the knob reads; make the refusal two-tier
- [A freedom can be an artefact of the frame](../feedback-a-freedom-can-be-an-artefact-of-the-frame.md) — a measured 0.0000 became 1.75 d of a 2.0 d window once the anchor's temperature was corrected; predicted 1 test would move, 6 did
- [Fit at one point, score out of sample](../feedback-fit-at-one-point-score-out-of-sample.md) — one rate fitted at 15 °C landed the two withheld temperatures within 6 %, which is what made it a test of the temperature RESPONSE and not a fit; the condition that got WORSE is the price, not a nuisance
- [A calibrated level decays when anything upstream moves](../feedback-a-calibrated-level-decays-when-anything-upstream-moves.md) — 7 constants defined as "k set to land X" stopped landing X across 3 beats (1.26-2.87×); the suite went red ONCE, because every guard was a direction or a ratio and both are scale-free
- [A constant lives in its rate law's coordinates](../feedback-a-constant-lives-in-its-rate-laws-coordinates.md) — `E_a_esters` printed 200,000 vs the observable's apparent 144,900, and a "Q10 ~ 2.6" note described the RATE CONSTANT while the beer ran 1.22; both became true when the rate law changed underneath them
- [A locked pair repairs or drifts together](../feedback-a-locked-pair-repairs-or-drifts-together.md) — same rate law ⇒ ratio fixed to 7 sig figs: history transfers with no old tree, and a repair skipping a member leaves a *guaranteed* defect (D-225 fixed 1 of 8)
- [A lump is right in total and wrong in fate](../feedback-a-lump-is-right-in-total-and-wrong-in-fate.md) — a constant back-solved to a published TOTAL absorbs contributors that LEAVE; fix by re-partition at constant total, and check the sign at both ends
- [A declared quantity can have a SECOND channel](../feedback-a-declared-quantity-can-have-a-second-channel.md) — a fit evaluated at one scenario input while a second input adds to the same quantity: a wine declaring 250 mg N/L carries 362.7, and conservation still closes (D-243)
- [A containment claim is scoped to where the ranges nest](../feedback-a-containment-claim-is-scoped-to-where-the-ranges-nest.md) — "that band contains this one, so don't sample it" was measured at ONE point and written as a property (D-243)
- [A repair that leaves the label lying](../feedback-a-repair-that-leaves-the-label-lying.md) — when two repairs both fix the consequence, the one that leaves a field's NAME meaning something else is not the faithful one; re-read the defect's own sentence
- [Two errors in the same frame cancel](../feedback-two-errors-in-the-same-frame-cancel.md) — a declaration and a release defect in the SAME frame cancel EXACTLY; repairing the visible half alone cost 6.1 %. Test the identity against the DECLARATION, then refuse the half (D-248)
