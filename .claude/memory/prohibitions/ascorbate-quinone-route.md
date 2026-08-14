---
name: ascorbate-quinone-route
description: "Ascorbate on the quinone node (D-202) - BUILT and dose-gated; the published 2:1-to-1:1 signal is inexpressible in this model at any rate, and the polymerisation-band closure is ruled out"
metadata: 
  node_type: memory
  type: project
  originSessionId: 760f4220-cd88-4a64-91d5-a5002cca53b9
  modified: 2026-08-14T06:33:49.232Z
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

**Newly opened (D-202), and priced on its own terms, not inherited.** UWC §9: sotolon comes from
2-ketobutyric acid *"formed by either yeast metabolism **or ascorbic acid degradation**"*, and the
model already has `alpha_ketobutyrate` + `SotolonAldolCondensation` (D-107). That route deposits
**on-ledger** carbon and would force `ascorbate` onto the carbon ledger with a dehydroascorbate
product pool (the D-80 split-ledger capture). **No rate in the corpus.**

Receipts under `M:\claud_projects\temp\ferment\d202-ascorbate\` — `PREREGISTER.md`, `probe1`
(bridge vs D-200's 8.0915 %), `probe2` + `head-wt/` (the HEAD-plus-inert-pad control),
`probe3` (the monkeypatch that silently did nothing, kept as the worked example), `probe4`
(rate sweep + polymerisation band), `probe5` (guard envelope vs the 2.75× near-miss),
`probe6` (Table 24.3's pair across the band).
Related: [[feedback-a-pair-constrains-a-response]], [[feedback-relocate-the-unsourced-factor]],
[[feedback-a-null-result-needs-a-positive-control]], [[feedback-pin-the-band-not-the-nominal]],
[[feedback-a-ratio-guard-cannot-see-a-common-factor]], [[feedback-check-the-blocker-is-still-blocking]].
