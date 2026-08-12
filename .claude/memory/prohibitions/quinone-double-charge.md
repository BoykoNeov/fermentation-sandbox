---
name: prohibitions-quinone-double-charge
description: "D-197 — the size of the quinone double charge in the DEFAULT oxidative set: five draws, the majority of the oxygen, and why de-duplicating it naively is a regression"
metadata: 
  node_type: memory
  type: project
  originSessionId: 720cec99-de09-483a-bafb-7296d361e840
  modified: 2026-08-12T12:24:32.813Z
---

# The quinone double charge in the default set — D-197

**SETTLED. Do not re-propose any of this as unbuilt, unmeasured or blocked.** Read before
touching the direct oxidative set's O₂ bookkeeping, D-75's Strecker route, or anything citing
"the inherited quinone double-count lump".

## MEASURED, nothing built

Five tests (`tests/test_oxidative_cascade_guards.py` §Guard 9). No Process, no parameter, no RHS
change, **no YAML value** — one YAML *comment* rewritten. Suite 1744 → **1749**.

## The finding: a TRUE number scoped to the WRONG THING

D-75 conceded the lump and defended its SIZE with `k_strecker` — "~2 % of the always-on total,
well under 1 % of the O₂". **Both numbers are right.** Strecker's measured share of consumed O₂ is
**0.012–0.033 %**, and 2.0000 % is an exact identity between two constants (1.0e-5/5.0e-4), not a
property of any run. **But the lump is FIVE draws**, and `anthocyanin_fading` — identical
additive-on-top idiom, **never named in the concession** — takes **~1000× more** (30.4 %
unsulfited). Downstream total: **35.6 % unsulfited, 62.3 % at SO₂ 30, 76.1 % at SO₂ 60.** The
defence was conducted on the lump's smallest member. [[feedback-a-magnitude-defence-can-be-scoped-to-the-smallest-member]]

## NEVER re-litigate

- **The two-stage rework is NOT unbuilt — it is the cascade (D-141)**, wired, isolable, non-default
  by decision. D-75's out-of-scope call stands and is the owner's.
- **Never de-duplicate inside the parallel frame.** On Carrascón's **primary** rule (90 % consumed)
  deleting the five downstream `o2` debits with no re-baseline takes uptake **0.3840 → 0.1085
  mg O₂/L/day** against their observed **0.88–1.25** band: from 2.3× below the floor to **8.1×
  below**. **LOAD-BEARING for agreement with measured uptake** — a second reason D-75 never gave.
- **Use the 90 % limb, NEVER the 7-day cap alone.** The cap closes the window while the model's rate
  is still high (33–81 % consumed vs the reds' ~85–91 %) and flatters it to "15 % below" — the first
  draft of D-197 §4 said that and was **amended in-session**. Cleanest framing avoids rates:
  **18.6 d to consume 90 % of ~8 mg/L where eight real reds took ~7–8 d ⇒ ~2.3× too slow**, an
  **accepted deviation, recorded not tuned** (closing it moves D-72/D-73/D-74/D-81 anchors).
  The SO₂ 60 arm is **UNUSABLE** either way — free SO₂ 40.2 mg/L exceeds all eight wines (2.6–33.7).
- **The two ratio guards' margins are held open by FREE SO₂, and each asserts that premise.** At the
  fixture: free SO₂ 35.7 mg/L, sulfite alone **86.8 %**, so 0.9288 / 14.04× — but **unsulfited the
  same composition reads 35.6 % / 1.55× and BOTH thresholds fail.** Never "fix" a red there by
  loosening the threshold; the premise is what moved.
- **Σk_eff/anchor (9.04× / 25.05×) does NOT discriminate** — browning's `k_eff` contains D-132's
  *sourced* phenolic boost. Downgraded before the campaign ran; never quote it as the lump's size.
- **The partition is exactly kᵢ/Σkⱼ** (every sink first-order in `o2`), so it is level-independent;
  `closure_oxygen_ingress` is a **source**, never in the denominator.
- **No guard was OWED**: all four mutation arms were already red somewhere (67/30/5 full-suite reds,
  of which 3/1/3 are new). The guards assert the partition's **meaning**, not its first detection —
  do not cite them as having closed a hole.
- **This is NOT an opening onto the SO₂-role question** (D-72/D-137). Carrascón's "constituents
  compete with SO₂ for Fe(III)/quinones/H₂O₂" is corroboration only; the inversion's re-fix inside
  the parallel frame stays REFUSED.

Receipts: `M:\claud_projects\temp\ferment\d197-quinone-double-charge\`.
