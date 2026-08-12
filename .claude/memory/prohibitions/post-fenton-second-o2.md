---
name: prohibitions-post-fenton-second-o2
description: "D-196 — the second O2 the Fenton limb costs: what shipped, why the double-counting objection died, and why the state-dependence is load-bearing"
metadata:
  node_type: memory
  type: project
---

# The post-Fenton second oxygen — D-196

**SETTLED. Do not re-propose any of this as unbuilt or blocked.** Read before touching
`PeroxideEthanolOxidation`, the cascade's O2 bookkeeping, or anything citing "the missing
post-Fenton secondary O₂ draw".

## BUILT

`_O2_PER_ACETALDEHYDE = 1.0` in `oxidative_cascade.py`. `PeroxideEthanolOxidation.touches`
gains `"o2"`, and the limb debits `share × activation_rate` worth of O2 — the O2 the
1-hydroxyethyl radical takes. **Cascade only; the cascade is still NON-DEFAULT, so no shipped
run moves.** Suite 1742 → 1744 (two new guards).

## The item was NEVER blocked — twenty records said otherwise

D-142 → D-188 carried "the missing post-Fenton secondary O₂ draw" on their open lists —
**twenty records, D-150 through D-166 unbroken, then D-186/187/188** — and D-188 labelled it
`(D-142, unsourced)`. **It was printed all along** in *Understanding Wine Chemistry* 2nd ed.
§24.4.4.1, on disk in `_txt/` since before D-142: the 1-hydroxyethyl radical "goes on to react
with oxygen and produce acetaldehyde". Nobody re-opened the source.
[[feedback-check-the-blocker-is-still-blocking]]

## NEVER re-litigate

- **Pinned per ACETALDEHYDE, never per H₂O₂.** The source names an unquantified competitor
  (hydroxycinnamate quenching); per-acetaldehyde absorbs it, because a quenched radical yields
  neither the aldehyde nor the draw. Per-H₂O₂ would need a quench ratio nobody has published.
- **1.0 is an UPPER BOUND on the NET draw.** The hydroperoxyl limb could return an oxidising
  equivalent and push the ratio back UP. Unsourced — Danilewicz 2007's radical chain is
  *bisulfite* autoxidation, which he concludes does **not** run in wine. Named, not built, so a
  later recycling term cannot be credited with a cancellation.
- **Gate 1 (D-137) is NOT breached.** Its content is that *ground-state* molecules don't react
  with O2. A carbon-centred radical doing so is the next link of the same chain.
- **The double-counting objection was real and DIED ON A MEASUREMENT.** The instantaneous draw
  does double unsulfited (exactly 2×), but a hermetic O2 pool is **supply-limited**: cascade/direct
  budget 0.9614× → 1.0439×, so D-141's budget-reproduction claim **survives**. Beer's unexhausted
  120 d arm is where it shows: 71.3 % → 91.7 % (not shipped — beer defaults to direct).
- **The state-dependence is LOAD-BEARING, and it took TWO controls to show it.** A *rate* control
  (scale `k_o2_depletion_total`) moves the ratio 0.06 % — **vacuous**, because a rate cannot move
  a stoichiometric ratio in a supply-limited system. The **stoichiometric** control (same average
  cost as a flat per-activation surcharge) agrees to 0.21 % on the real-wine ratio and **breaks
  BOTH Danilewicz asymptotes** (0.9717 vs tol 0.01; 1.9455 vs tol 0.02) where the shipped
  share-weighted form leaves them at 0.9953/1.9926. Opposite outcome to D-195.
- **The copper "re-fit" is EXACTLY a no-op.** The new draw is copper-multiplied by construction,
  so `_COPPER_MULTIPLIED_DRAWS["cascade"]` gained it — but Ferreira's spread reads **1.3245× on
  both arms to four decimals**, because the branch fraction never touches `copper`. D-134's
  600 L/g is untouched.
- **Acetaldehyde: the production identity HOLDS** (1.003×/1.006×). The unsulfited pool falls
  3.89× on a 1.84× production change because ~87 % is consumed downstream (chiefly
  `acetaldehyde_reduction`) — the pool is a **residual**. At 60 mg/L SO₂ both move 3.6 % and
  nothing is visible. The residual amplification beyond production is **conjecture, labelled**.

## Cost recorded, not hidden

At the shared operating point the cascade reads 1.1338 → **1.1079**, still inside Miao's band.
The **dose-70 arm LEFT the band** (1.1034 → 1.0728), so "in band from ~70 upward" narrows to
"in band only at 80". The test still passes. **D-142's "would move r in the right direction" is
CONFIRMED, not corrected** — `r` is SO₂ per O₂ and D-142 wanted it DOWN, to reach wine #1's
sub-1. A `Corrects: D-142` marker was shipped and **withdrawn in the same session**.

## Still open

The hydroperoxyl recycling limb (unsourced). The activation node's pH term. The quinone
branching. Receipts: `M:\claud_projects\temp\ferment\d196-post-fenton-o2\`.
