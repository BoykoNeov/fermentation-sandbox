---
name: quinone-branching
description: "The quinone node's consumer split (D-145 to D-199) - the named pull is on disk, and it still cannot close the branching"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — the quinone node's branching.** Detail split out per D-185's pattern; the
status ledger points here by path. Read it when working on this subject. Every bullet is *what it
forbids* + the D-record to read for *why*. If a prohibition looks unconvincing, **go read its
D-record — do not argue past it from this file.** **Never evict an old prohibition to buy a line.**

**The blocker is DEMOLISHED but the item is NOT closed (D-199).**
- **Never call Nikolantonaki & Waterhouse 2012 "the pull that would settle it" again.** Its
  **Figure 24.12 is reprinted with ACS permission in *Understanding Wine Chemistry* 2nd ed.
  p. 331**, on disk since before D-141, and rode ~50 records as a purchase blocker. Third time a
  "blocked on external sourcing" item was already on disk (D-191, D-196). Four code sites struck.
- **It does NOT settle `k_quinone_polymerization`, and the correction must not say it does.**
  The note's *"NOTHING IN HAND ADJUDICATES WHICH"* (of the two fate-gap closures) **SURVIVES the
  pull arriving** — the figure ranks **nucleophiles**, and this Process is by construction the
  fate of quinone **no nucleophile captured**. A nucleophile table can have no row for it.
  **Status is now "blocked on STRUCTURE", not "blocked on a pull."**
- **The figure ranks 2 of the 5 consumers**: `quinone_sulfonation` (SO2, letter `e`) and
  `quinone_strecker_degradation` (Met/Phe, letter `a`). It does **not** rank polymerisation,
  anthocyanin fading or ellagitannin. **Never map phloroglucinol onto `quinone_polymerization`** —
  Phl is the flavan-3-ol A-ring **addition**, a co-substrate reaction the model has **no Process
  for** (`tannin` never consumes `quinone`); that mis-map would have invented a six-decade
  contradiction against the model's largest consumer.
- **Ascorbate and glutathione are in the JOINT TOP GROUP with SO2** (letter `e`) and have **zero
  slots, zero parameters, zero Processes**. UWC's own between-nucleophile argument runs through
  exactly those two, so the branching cannot be set right without them.
- **Values off Figure 24.12 are READ FROM A LOG PLOT (±0.3 decades) and are never a `source:`.**
  The **printed** content is the significance letters plus two verbatim "orders of magnitude
  slower" sentences; that carries the ordering without the eye-reading.

**The branching IS the benchmark headline — and 1.7481 is REFUSED.**
- Scaling `k_quinone_polymerization` (a constant `QuinoneSulfonation` **never reads**) walks the
  SO2:O2 ratio **0.9738 → 1.1079 → 1.5172** across its own band, and to **1.7866** asymptotically
  (**not 2.0** — bisulfite never wins the whole H2O2 node and D-196's second O2 inflates the
  denominator). Bisulfite's share goes **1.86 % → 15.21 % → 59.69 %**.
- **`×0.01` lands at 1.7481, within 3 % of Danilewicz's 1.7 — REFUSED, do not rediscover it.**
  D-138 and D-141 both say **1.7 must EMERGE from partial quinone capture, never be fitted.**
- **The HIGH edge of the declared band reads 0.9738** — below Miao's floor, below Danilewicz's
  range, below the **1:1 blocked-quinone limit**. **This is NOT D-174's instance** (that is
  `k_o2_depletion_total`, 1.0655) — two different constants, **never merge them**. Scope it
  narrowly: live only in a **cascade ensemble**, and **no shipped test runs one**.
- The guard derives its band edges **from the parameter store, never the note** (D-154 idiom). A
  **red there means the sensitivity moved and D-199 must be re-measured** — never relax it.

**The product identity was wrong, and fixing it CORROBORATES fork D2.**
- `QuinoneSulfonation`'s sulfonate is the **MINOR** product: UWC §24.4.3.2 gives **>90 % reduction
  back to the o-diphenol at wine pH**. **`_SO2_PER_QUINONE` stays 1.0** — both routes consume one
  bisulfite per quinone — so the emergent 1:2 limit is untouched and the fix is **inert by
  construction**, not by measurement.
- **Never read this as breaking fork D2.** `OxygenActivation` never debits the o-diphenol; under
  reduction that pool is **regenerated**, so not debiting it is *right* for the major route. The
  residual sulfonate route is the part that permanently removes a phenolic and is **recorded, not
  built** — no shipped pool is its substrate (`hydroxycinnamics`/`tannin` are activation drivers).
- **The Strecker ordering is off ~6.4 decades and is INCONSEQUENTIAL** (0.0025 % of the node).
  Never promote it to a second headline. Its content is a mechanism claim: the quinone→Strecker
  route **cannot run in a well-sulfited wine**, which is UWC's own regime caveat. The model's
  Strecker law is **gated, not bilinear** in the amino acid, so the normalisation is an
  **equivalence, not an identity**.

Receipts: `M:\claud_projects\temp\ferment\d199-quinone-branching\` — `PREREGISTER.md` (**4 of 9
predictions missed**, including the one the beat turned on), `RESULTS.md`, `FIG_24_12.md`,
`probe_shares.py`, `probe_ordering.py`.
Related: [[feedback-check-the-blocker-is-still-blocking]],
[[feedback-a-named-pull-may-not-answer-the-question]],
[[feedback-name-the-field-your-predicate-read]], [[feedback-pin-the-band-not-the-nominal]].
