---
name: nitrogen-storage-and-charge
description: "D-250 — the uptake surplus is INTRACELLULAR (`stored_nitrogen`, in no charge balance), the split rule, `touches_where_present`, and why the bacterial residue is CLOSED as not-a-defect"
metadata: 
  node_type: memory
  type: project
  originSessionId: c5d2411f-695e-478f-bbae-e0b90d259b55
  modified: 2026-08-28T09:23:20.762Z
---

**Live prohibitions — the intracellular nitrogen store and the charge it must not carry (D-250).**
Read before proposing anything about `stored_nitrogen`, `AssimilableNitrogenUptake`'s destination,
growth's nitrogen Monod, the D-32 swap's refund, `Process.touches_where_present`, mid-run wine pH
on a dosed must, or "MLF/Brett are starved by yeast uptake". Sibling:
`wine-nitrogen-budget.md` (D-243→D-248), which this corrects on two counts.
**If a prohibition looks unconvincing, go read D-250 — do not argue past it from here.**

- **`stored_nitrogen` is in NO charge balance and NEVER goes back into `N`.** D-248 refunded
  uptake's surplus to `N`, which `acidbase` reads at the must's mean charge per mole N (D-209) —
  so nitrogen already inside the cell went on titrating the must: `N` 300 → **436.8 mg N/L** and
  pH 3.030 → **3.216** on a 2 g/L amino-acid-dosed wine, **+0.215** against the no-uptake arm.
  Dose-scaled (0.045 at 0.5 g/L), and exactly **0.0000** undosed — D-248's isolability claim was
  never in doubt. **Nitrogen was conserved to 1e-8 the whole time: the defect is CHARGE, not
  mass.** Never restate the slot rise as nitrogen creation. It also restores
  `nitrogen_charge_excess`'s recorded invariant that no Process adds differently-charged N to `N`.

- **The draw split is PROPORTIONAL and it is ONE helper for two callers.** Growth's Monod and its
  `f_N·dX` draw read `N` + store together, split by what each holds; the D-32 swap refunds that
  same draw on the **same** split (`growth.add_assimilable_nitrogen`). **Booking the swap's refund
  wholly to `N` against a split draw sends net `dN/dt` POSITIVE and the artefact returns with the
  uptake Process innocent** — guarded, driven at a state where the store holds 95 %. Store-first
  was rejected: it puts a C⁰ kink exactly where the store empties.

- **The two CELLULAR gates read the sum, and that was forced, not chosen.** Uptake *inflated* `N`,
  so the Ehrlich fusel gate `n/(K_n+n)` ran high and H₂S de-repression `K/(K+n)` ran suppressed
  from D-248 onward — the state in which **D-248's four fusel xfails closed**. Gates on `N` alone
  collapse both. Every other `N` writer (aging, keto acids, mercaptans, oxidative cascade,
  precursor fates, the swap) is a genuine release into the medium and is untouched.

- **The whole observable footprint is pH, and that is the honest scope.** `N` and the store share
  a source and all sinks, so the SUM follows the pre-D-250 `N` exactly: **beer bit-identical
  (0.00e+00)**, wine ≤**1.9e-7** on biomass/sugar/ethanol/all eight pools/three fusels/H₂S. That
  residual is the SCHEMA, falsified not asserted — the **undosed** wine, where uptake cannot
  contribute at all, moves 1.17e-7 too (98 → 99 slots shifts BDF step selection; third instance).

- **The bacterial-nitrogen residue is CLOSED as NOT-A-DEFECT — never re-propose it.** D-248's
  `Flags: D-100` ("MLF/Brett cannot read the `N` slot uptake fills") rested on reading that slot
  as extracellular. It is intracellular, so the blindness is **correct**. Brett still loses **96 %**
  of its growth increment and MLF 46 % at the low dose, and that is the model being RIGHT: real
  yeast take essentially all the assimilable nitrogen (Crépin's 0.2 %). **The gap is a bacterial
  nitrogen source the model lacks — PEPTIDES, which yeast do not take.** Isolating the competitor
  in the MLF/Brett tests is now PERMANENT, not a stopgap. Never teach a bacterium to read the store.

- **The malate reading REVERSES across the dose — never quote the 2 g/L one alone.** At 0.5 g/L the
  starved arm leaves **4.0× more** malate at day 3 (starvation showing through); at 2 g/L the
  0.215 pH excursion lifts MLF's own pH logistic enough to over-compensate and the starved arm
  converts MORE. Two errors in one observable. Brett's confound points the OTHER way (higher pH ⇒
  less molecular SO₂ ⇒ less inhibition), which is what makes its 96 % attributable.

- **`Process.touches_where_present` is a medium-conditional DECLARED touch, not an exemption.**
  For a Process wired into BOTH media that writes a slot only one carries — growth is the one
  case. `touches` outright makes every beer ProcessSet raise on an unknown variable; declaring it
  nowhere makes every wine ProcessSet raise on a leak. An **undeclared** slot still fails strict
  mode, pinned with a leaky-Process arm. Do not widen it to a second Process without that case.

- **Named and NOT repaired:** the store is not tied to living biomass, so it does not drain on cell
  death or autolyse (conservation still closes; equally true of D-248's parking in `N`). And
  D-250 makes a **narrower reading separable for the first time** — nitrogen still OUTSIDE the
  cells, which is what Crépin sampled. `_assimilable_n_mgl` keeps counting the store, because that
  is what D-248's "40.8 % → 0.62 %" meant. On the narrower quantity **D-249's gap reads 1.71×, not
  1.59×** — and **D-249's verdict survives**: still slower than the 1.92× run containing it.
  Re-deriving D-249's headline on it is a separate beat.

Measurements: `M:\claud_projects\temp\ferment\d250-bacterial-nitrogen\` — `PREREG.md` and
`PREREG-BUILD.md` (both written before their runs), `FINDINGS.md`, `probe_size.py`,
`probe_mlf_reversal.py`, `probe_ph.py`, `probe_ph_scope.py`, `probe_dose_decomp.py`,
`probe_identity.py` + `base.json`/`new.json` (one process per tree), and `mutate.py`
(three arms, snapshot-restored and byte-checked between them).
