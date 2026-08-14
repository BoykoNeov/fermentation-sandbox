---
name: ascorbate-quinone-route
description: "Ascorbate on the quinone node (D-202) - BUILT and dose-gated; the published 2:1-to-1:1 signal is inexpressible in this model at any rate, and the polymerisation-band closure is ruled out"
metadata: 
  node_type: memory
  type: project
  originSessionId: 760f4220-cd88-4a64-91d5-a5002cca53b9
  modified: 2026-08-14T07:21:36.194Z
---

**Live prohibitions — ascorbate (vitamin C) at the quinone node.** Detail split out per D-185's
pattern; the status ledger points here by path. Read it when working on this subject. Every bullet
is *what it forbids* + the D-record to read for *why*. If a prohibition looks unconvincing, **go
read its D-record — do not argue past it from this file.** **Never evict an old prohibition to buy
a line.**

**BUILT at D-202, and Figure 24.12's top group is now COMPLETE.**
- `QuinoneAscorbateReduction` ships: wine-only, cascade-only, aging-gated, `{quinone, ascorbate}`.
  **Do not re-propose an ascorbate slot as unbuilt.** All four top-group members are settled — SO2
  (D-141), GSH **priced and REFUSED** (D-200, 0.32 %), H2S (D-201), ascorbate (D-202). **No member
  of that group remains available as an explanation for the branching gap.**
- **The default is 0 and it is SOURCED, not a convenience** — UWC §24.4.3.2: *"new wine has a
  negligible ascorbic acid content"*, and it is a permitted additive. **Never seed this pool the
  D-134 copper way**: dosing moves the benchmark by more than its whole margin. It enters only
  through `add_ascorbate`. This is **not** the D-45 hard-zero defect.
- `k_rel_ascorbate_quinone` = **1.0 from printed prose** (*"just as quickly"*), band [0.5, 2.0],
  **both edges CONSTRUCTED**. Never quote Figure 24.12's log plot as a `source:` (D-199).
- **There is NO unsourced coefficient here to relocate**, and that is a measured claim, not an
  omission: the reduction is two-electron on both sides, so 1:1 is fixed at both ends and
  dehydroascorbate has no second handle. **Do not "apply D-201's relocation lesson" here** — the
  law is written quinone-first on purpose.

**The source states the ascorbate/SO2 relation THREE times and they disagree — do not pick one silently.**
- Rate constant *"just as quickly"* ⇒ 1.0. Wine observation *"appears to react faster than SO2"* at
  0.34 vs 0.5 mM ⇒ ≥~1.5. Wine result *"completely prevent[s]"* ⇒ ~30 by naive arithmetic.
  **The band spans none of the top of that, DELIBERATELY** — burying the conflict inside an
  uncertainty range converts a falsifiable gap into a tuning allowance.

**THE HEADLINE: the published signal is INEXPRESSIBLE, and it is not a rate problem.**
- Table 24.3 prints a **pair**: total SO2 loss per O2 is **2:1** without ascorbic and **1:1** with.
  The signal is a **halving**. The model's response is **0.990** (1.107884 → 1.096812).
- **Never propose closing that by fitting `k_rel`.** Measured: at **100×** the printed rate — two
  decades above the source — the response is still only 0.898. The pool **depletes**; 60 mg/L is
  0.34 mmol/L against ~0.42 mmol/L of quinone through the challenge, so the dose is comparable to
  the whole flux it must intercept. **D-200's fixed-concentration arithmetic does not transfer.**
- The real reason: **the model has no room to fall.** Its *un-dosed* baseline is already **1.108**,
  at the bottom of the 2:1→1:1 gap before any ascorbate exists.

**The pair ADJUDICATES what D-199 said nothing could — and the verdict is a refusal.**
- A pair constrains a **response**, and a response needs **no absolute rate** — which is what the
  branching was blocked for. Across `k_quinone_polymerization`'s whole declared band, un-dosed →
  dosed: LOW 1.5172 → 1.3955 (0.920), nominal 1.1079 → 1.0968 (0.990), HIGH 0.9738 → 0.9736
  (0.9998). **No position reproduces both printed cells**, and the response *weakens* as
  polymerisation speeds up. **D-141's polymerisation-band closure is RULED OUT.**
- **The nucleophile-constants closure is NOT ruled out and was NOT tested** — `k_so2_oxidation` is
  shared with the **peroxide** node, so scaling it moves both halves of D-72's split at once and
  the arm would not answer the question it appears to. Separating them is a structural beat.
- **This corrects D-200's attribution**: the benchmark's in-band pass needs BOTH ascorbate's absence
  AND the polymerisation constant's position — at its band's HIGH edge an un-dosed wine reads
  0.9738, below Miao's floor with no ascorbate at all. **The statistic cannot separate them.**

**Two further effects are NOT built, on structure, and the corpus WAS searched.**
- The **pro-oxidant limb** (ascorbate + O2 → H2O2) cannot be a Process at all: H2O2 is
  quasi-steady-state, so it is a change to the shared branch-node **rate law** (D-198's class).
- The **iron-cycling acceleration** (*"white wines containing ascorbic acid … more comparable to red
  wines"*) is **mechanically distinct** and is **NOT** built by adding `ascorbate` to
  `_PHENOLIC_REDUCTANT_POOLS` — that path mints a quinone per O2 and ascorbate-driven turnover makes
  none.
- **Never call these "blocked on sourcing" without repeating the search**: 24 texts, **17 mention
  ascorbate, 203 mentions, no rate constant** for the O2 limb. The search also found an independent
  second source for the built route — *Handbook of Enology* Vol. 1: ascorbic acid *"only limits
  browning **by reducing quinones** but **does not limit oxygen consumption**"*.
- **Both omissions push the headline DOWN, same as the built route** — never credit either with
  cancelling it (D-200 verdict 2).

**The sotolon route — MEASURED and REFUSED at D-203, on IDENTIFIABILITY. Do not re-propose it as
"blocked on a rate".**
- **The source is the model's OWN keystone paper**, not a UWC §9 clause: ref [31] is **Pons *et al.*
  2010** (JAFC 58:7273), already cited by `SotolonAldolCondensation` and by `k_sotolon_aldol`'s
  `source:`. The route is its **first abstract result**. D-202's "nothing sources its rate" is true
  of the *corpus* (7/24 texts, 19 mentions, every one outside UWC sourcing 2-KB to threonine or
  methionine) but **under-reads the paper**. Never call the route thinly sourced.
- **REFUSED because one anchor cannot fix two free parameters**: the node has one observable (the
  5-20 µg/L anchor, model at 7.4461) and a build adds the yield *plus* leaves `k_sotolon_aldol` an
  author estimate. **The "it breaks the calibration" objection is MEASURED FALSE** — the calibration
  wine carries `ascorbate` = 0.0 throughout, so the route contributes exactly 0 there.
  [[feedback-count-the-anchors-before-adding-a-parameter]]
- **The item now carries a TARGET, not a shrug: ~10 % molar conversion** at 20 mg/L O2 crosses the
  8 µg/L threshold; below ~1 % it is invisible; at 5 mg/L O2 + 60 SO2 it **never** crosses.
  Ceiling (60 mg/L ascorbate, 1:1) = 34.78 mg/L α-KB = **16.97×** the model's 2.049.
- **A second precursor source is worth EXACTLY ZERO without oxygen** — bitwise, 6 s.f., across the
  whole 17× ladder, because the aldol is the *product* of two substrates and a sealed unsulfited dry
  white has acetaldehyde ≈ 0. With O2 it is 16.8-17.4×. Never price this route without naming the O2.
- Pons scopes his own measurement out **twice**: model wine at **40 °C / 6 months**, and his
  conclusion is the pathway *"is also valid in white wines with no added ascorbic acid"* because
  yeast 2-KB already reaches oxidised-wine levels. **Unblocking needs Pons 2010's FIGURES** (dose→
  sotolon at 40 °C/6 mo; and the 2-KB × acetaldehyde series, which supplies the missing 2nd anchor).
  **Closed on every host tried** — ACS 403, Wiley 402, HAL denied, OpenAlex/S2 `oa_status: closed`.
  A search summary offered "1 µg/L from 10 mg/L 2-KB + 1 mg/L acetaldehyde" — **NOT recorded as a
  value** [[feedback-transcribe-tables-not-prose]].
- If ever built, **isolability is exact and free** (mass-action, 0 substrate ⇒ bitwise 0). D-202's
  carbon cost stands: on-ledger carbon forces `ascorbate` onto the ledger + a dehydroascorbate pool.

**COUPLED — `SotolonAldolCondensation` is NOT in `_AGING_GATED_PROCESSES` and runs from t = 0
(D-203 §7). INVESTIGATED: leave it UNGATED. Do not re-propose a gate.**
- **Never quote the pre-aging SHARE (100 % / 45.3 % / 22.0 %) — it is a ratio artefact.** The
  numerator is a **CONSTANT**: 0.025302 µg/L unsulfited, 0.026757 at must-SO2 60, in *every* arm,
  = **0.32-0.33 % of the 8 µg/L threshold**. Worst case found is the **sweet** anchor at 0.281699
  (3.5 % of threshold). [[feedback-a-ratio-guard-cannot-see-a-common-factor]]
- **The D-49 "transient acetaldehyde" argument for gating is VOID** — D-49 is about *pyruvate*;
  `acetaldehyde` is D-27's green-apple transient **and the SO2-binder**, a real extracellular pool
  (peak 4.19→24.64→45.16→72.60 mg/L as must SO2 goes 0→100). Both substrates really coexist
  during fermentation, so **a gate would delete real chemistry**. The docstring's "sealed bottle,
  no living yeast" picks the **α-KB pool**, it does not switch the reaction off.
- **A guard is still owed and it is NOT a gate**: pin the **absolute** offset as a small constant,
  never the share. The counterfactual gate is a **GREEN mutation** — every assertion in
  `test_so2_protection_erodes_as_the_free_so2_is_spent` survives it though the 5 mg/L ratio moves
  0.178 → 0.145. Nothing in the suite sees this term.

D-203 receipts under `M:\claud_projects\temp\ferment\d203-ascorbate-sotolon\` — `PREREGISTER.md`,
`probe1_ceiling.py` (rate-free ceiling + the sealed null), `probe2_with_oxygen.py` (the O2 ladder
that scoped it), `probe3_gate_attribution.py` (the aging-gate attribution, disabled arm = clean 0),
`probe4_crossing_yield.py` (kills the calibration blocker, measures the crossing yields), `RECORD.md`.

D-202 receipts under `M:\claud_projects\temp\ferment\d202-ascorbate\` — `PREREGISTER.md`, `probe1`
(bridge vs D-200's 8.0915 %), `probe2` + `head-wt/` (the HEAD-plus-inert-pad control),
`probe3` (the monkeypatch that silently did nothing, kept as the worked example), `probe4`
(rate sweep + polymerisation band), `probe5` (guard envelope vs the 2.75× near-miss),
`probe6` (Table 24.3's pair across the band).

Related: [[feedback-a-pair-constrains-a-response]], [[feedback-relocate-the-unsourced-factor]],
[[feedback-a-null-result-needs-a-positive-control]], [[feedback-pin-the-band-not-the-nominal]],
[[feedback-a-ratio-guard-cannot-see-a-common-factor]], [[feedback-check-the-blocker-is-still-blocking]],
[[feedback-count-the-anchors-before-adding-a-parameter]], [[feedback-re-read-the-source-you-already-mined]],
[[feedback-a-threshold-cannot-separate-same-sign-regimes]], [[feedback-same-species-different-reaction]].
