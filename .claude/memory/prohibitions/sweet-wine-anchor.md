---
name: prohibitions-sweet-wine-anchor
description: "D-194 — the sweet-wine aging scenarios' must anchor: what shipped, why 38 and not 45, and why this beat owed a guard where D-189 did not"
metadata: 
  node_type: memory
  type: project
  originSessionId: 53cf9b2e-d59e-43b1-8918-019200d03f50
  modified: 2026-08-12T10:03:29.351Z
---

# The sweet-wine scenario anchor — D-194

**SETTLED. Do not re-propose any of this as unbuilt.** Read before touching `_SWEET_BRIX`, the
Maillard/Caramelization aging scenarios, or the sweet-wine tests' premise.

## BUILT

**`_SWEET_BRIX = 38.0`** in `tests/test_aging_scenario.py` (was 70.0). A real botrytis must —
412 g/L, ferments to the D-129 ethanol ceiling (E≈159.5) **well before** the fixed day-30
`begin_aging` breakpoint, and holds **76.7 → 76.4 g/L** across the 730-day tail.
**No physics changed** — this is a test-file constant; every parameter and Process is what D-193
left.

Plus `test_the_sweet_scenarios_age_an_ARRESTED_wine_not_a_still_fermenting_must`, which carries
the retired 70.0 as a **permanent positive control**.

## NEVER re-litigate

- **70 °Brix was not "a weak wine" — it was a must STILL FERMENTING.** `begin_aging` fired at
  **0.50 % ABV** on a must losing 4.9 g/L over days 25-30, and it fermented on through the whole
  two-year tail. D-192 §9's "they age a 1.4 % ABV wine" is corrected. Its "947 g/L raw" is not
  the model's initial `S` either — that is **880.7 g/L**.
- **The stale comment credited D-129's ethanol ceiling for arresting it.** True when written;
  D-192's osmotic brake (~76× at 880 g/L) means the ceiling is never reached. The claim existed
  **twice** (the constant's comment and the `_wine` helper's) and was checked nowhere.
- **Never anchor on the RESIDUAL.** 45 °Brix alone reproduces real Sauternes' 120-150 g/L
  residual and was **REJECTED**: it buys that by inflating the must past what botrytis must is,
  to conceal the model's **20.2 % vs real 13-14 % ABV ceiling**. `_SWEET_BRIX` **is the must's
  Brix — an input the winemaker measures**, so it takes a real value and the residual is left to
  be whatever the model predicts. **The 20.2 % ceiling is deliberately left visible, not
  compensated** — re-anchoring a scenario input does not license retuning physics.
- **The binding constraint is the FIXED day-30 breakpoint, not the brake.** D-192's "at 32-40
  °Brix the brake changes the PATH, not the DESTINATION" was measured on runs long enough to
  *arrive*, and says nothing about arriving by a fixed breakpoint. Always check `E` at the
  breakpoint and `S` flatness, never just the endpoint.
- **The brake SELF-LIFTS** below its 300 g/L threshold, so the opening rate cut costs far less
  wall-clock than its t=0 magnitude suggests. Everything 32-45 arrives. The risk is the **low
  end failing `assert S > 50.0`**, never the high end failing to arrive.
- **Do not weaken the guard's mutation arm.** One arm must be expected GREEN and one RED, or
  neither result carries information.

## The lesson that generalises

**`assert S > 50.0` passed for the WRONG REASON** — a still-fermenting must has *more* sugar than
a stuck one, so the sweetness assert witnessed the premise's collapse without noticing. All 98
tests in the file passed at 70 °Brix. **A magnitude threshold on an unqualified pool cannot tell
you which regime produced it.**

And: **the mutate-the-premise check came back GREEN, which is what LICENSED a new guard** — the
opposite verdict to D-189, where the same check went red on five existing asserts and no guard was
owed. **A count of asserts naming a thing predicts neither outcome** (D-195: pyruvate had 27,
α-KB had none, both refused without a new guard).

See [[feedback-mutate-the-premise-before-building-the-guard]],
[[feedback-verify-the-restore-between-mutation-arms]], [[feedback-pre-register-the-cheap-prediction]].
