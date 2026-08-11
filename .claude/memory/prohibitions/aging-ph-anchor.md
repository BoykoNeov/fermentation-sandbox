---
name: aging-ph-anchor
description: "Setting the pH a beverage AGES at (D-186) - the set_ph verb is BUILT and is cation-moving, not a pH dial"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — the aging-pH anchor.** Detail split out of
`.claude/memory/project-fermentation-sandbox.md` at D-185; that file's ledger points here by
path. Read it when working on this subject. Every bullet is *what it forbids* + the D-record
to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue
past it from this file.** **Never evict an old prohibition to buy a line.**

**The aging pH is an INPUT now (D-186). Do not re-propose it as missing.**
- **D-150's "there is no way to set an aging pH" is FALSE since D-186** — but the sentence
  survives verbatim in the copy-forward open lists of **~8 records**, alongside D-143's
  withdrawn SO2 item. **The ledger is the authority; that tail is not.** The gap it named:
  `initial_ph` anchors t=0, `Byp` then drags pH down, and a **0.3500 span arrives as 0.3045
  (87.0 %)** — *not* a constant offset (arms drift **0.0525 vs 0.0980**). **D-150's own 0.3052
  no longer reproduces** — D-182's carbonic term moved it; re-measure, never quote.
- **`set_ph` is CATION-MOVING — never describe or extend it as a pH dial.** Raising pH = carbonate
  deacidification, lowering = **cation-exchange resin**; both real, which is what makes it shippable.
  **Acidifying by ADDITION stays `add_acid`** (acid onto its own slot, cation untouched) — never
  merge the two. **D-65 anticipated the CATEGORY, in its DOSE form** ("K-tartrate *additions*");
  this ships the **TARGET** form an aging study needs — **never claim D-65 specified the API**.
  Cation moved −5.67e−3 / +1.18e−2 mol/L in the receipts.
- **`cation_charge_for_ph` must stay the inverse of `ph_of_state` TERM FOR TERM** — same per-medium
  registry (D-179), same `Byp`, same `min(evolved, C_sat(T))` (D-182). **A one-sided inverse passes
  every "pH moved the right way" test**; the guard that catches it is the **no-op round trip**
  (own pH in ⇒ same cation out, **1.24e−15**). Forward/inverse agree to **1.75e−11** over pH 2.8-4.0.
  It is a **state** anchor: `solve_cation_charge` may assume `Byp`=CO2=0 at the compile seam, this may not.
- **The reachability check CANNOT move to compile — do not "fix" it there.** The floor is a property
  of the *state* at the adjustment. It raises from the mutation, which the driver applies **BETWEEN
  segments**, so it surfaces labelled with no partial state — not a traceback from inside the RHS.
  Floor = `ph_of_state` with the cation zeroed (`solve_ph` is total, D-46); **2.1869** on the probe wine.
- **The opt-in gate's reason DIFFERS per medium — never let one justification cover both.**
  **Beer:** structural — without `initial_ph` every acid slot is 0 (D-179), so anchoring writes a
  counter-cation into an **empty** acid load. **Wine:** epistemic — the balance is populated
  (`tartaric`/`malic` seed regardless; un-anchored wine is pH 2.92, D-182), so the objection is that
  it would **manufacture** the pH information D-18 says must be an input.
- **Never write the `CO2` slot from it** (cumulative EVOLVED gas, already saturated ⇒ moves the
  dissolved term by exactly nothing) and **never claim it models cold stabilisation** — KHT removes
  tartrate too and **only the cation is booked**. **No carbon/nitrogen flow** (charge, not matter)
  and **no tier moves**; the zero is asserted **beside** a "the cation actually moved" check so it
  cannot pass vacuously. **Beer's finished pH stays a PREDICTION** (D-180) — this acts downstream of it.
- **Under an ENSEMBLE the anchor is NOMINAL-ONLY — measured, and it is the PRECEDENT.** `y0` and
  `events` are held fixed across members while each draws its own pKa set (**21** vary). Spread at
  24 members: `initial_ph` **0.1273** at t=0 vs `set_ph` **0.1292** post-event, **ratio 1.015** —
  same class, not a new weakness (a compile-time value is in no `reads` ⇒ never drawn).
  **NEVER "fix" this one anchor alone**: the two would then disagree about a member's pH.
- **The spread-ratio guard for that was DESIGNED and REJECTED — do not rebuild it.** Falsified
  before shipping: a planted pKa shift throwing the **nominal** anchor **0.29 pH** off target moved
  the ratio only **0.988 → 1.119**, green at any honest threshold. A wrong SHARED input shifts the
  **mean**, not the **spread**, so the statistic was blind to the exact defect it named — and blind
  reassuringly. **What catches it is the nominal-exactness pin, which already exists.**
  General form: **falsify a guard against the defect it names; if it stays green, ship the
  measurement WITHOUT the guard** [[feedback-a-spread-guard-misses-a-shared-input]].
- **It adds NO pH-rate dependence anywhere.** The under-response stands: real wine ~**1.62×** per pH
  unit against the model's ~**1.006×** (D-150), and the pH term on the activation node is still
  **REFUSED** on two legs. What D-186 buys is that the gap is now measurable *at a stated pH*.
  Direction check only: molecular SO2 **0.0014 vs 0.0006 mg/L** at aging pH 3.26 vs 3.61 — a nearly
  fully bound wine, so **~2.3× is a direction, NOT a magnitude claim**.
