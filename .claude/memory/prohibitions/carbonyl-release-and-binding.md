---
name: carbonyl-release-and-binding
description: "The SO2 free/bound equilibrium's release side and the four binders' constants (D-190) - release emergent, methional blocked, an ordering corrected"
metadata: 
  node_type: memory
  type: project
  originSessionId: f829484e-0ed2-4cb0-97ec-aaac6b46ce0c
  modified: 2026-08-12T06:26:12.691Z
---

**Live prohibitions — carbonyl release and the binding constants.** Split out at D-185's
pattern; the status ledger points here by path. Read it when working on this subject. Every
bullet is *what it forbids* + the D-record to read for *why*. If a prohibition looks
unconvincing, **go read its D-record — do not argue past it from this file.** **Never evict an
old prohibition to buy a line.**

**The adduct-release route is MEASURED (D-190). Do not re-propose it as untested — D-137/D-139/
D-140 carried "deserves its own D-record" forward three times and it now has one.**
- **Release is ALREADY EMERGENT for acetaldehyde — never build it.** The equilibrium is
  **stateless** (no bound-adduct slot); release is what the re-partition does as SO2 falls.
  Bound share 99.75 % → 99.48 % over five years, leaving free **1.7×** above the
  fixed-bound-fraction counterfactual.
- **The model REPRODUCES a published two-sided literal** (UWC 2e Ch. 24 via Ch. 17: >99 % bound
  above 30 mg/L free SO2, >95 % above 2 mg/L) at **99.5-99.7 / 95.1-95.4 %**. Now pinned. **The
  BAND does not carry it**: the weak edge gives **93.3-93.6 %**, missing the 95 % floor, and the
  nominal clears it by only ~0.2 pp. **Never "fix" either side** — both are sourced.
- **The oxidised-aroma carbonyls bind NOTHING** — only acetaldehyde/pyruvate/alpha-KG/
  5-oxofructose do. So Bueno's release-vs-de-novo ordering is **not expressible** for methional
  or phenylacetaldehyde, and SO2's 11.5× methional suppression is **oxidant competition only,
  never masking**. Right direction, half the mechanism.
- **The methional adduct is BLOCKED, not refused on merit** — no locally-held Kd. UWC Table 17.2
  lists other odorants and **not** methional, and the same book calls it "weakly binding" while
  the one saturated aldehyde it does tabulate (hexanal 3.5e-6) is **strong**. **Never interpolate
  between those.** Unlock: a measured methional-bisulfite Kd.
- **Diacetyl is a sourced fifth binder and is REFUSED anyway** (Table 17.2, Kd 1.4e-4). The only
  part of that table checkable against our own provenance is the part that **disagrees**, and the
  disputed value is **the same number**, adjacent in the same column of a PDF extraction whose
  name and value columns arrive as separate blobs. Needs the rendered table or a second source.
- **The aroma readout reads TOTAL, deliberately** — switching acetaldehyde to its free share moves
  no verdict (0.49 mg/L against a ~100 mg/L threshold) and would make it the only pool of 23 whose
  OAV means something different. **Pinned as a recorded choice, not a bug.**

**The pyruvate/alpha-KG ordering is CORRECTED IN PROSE ONLY — the values are UNRESOLVED.**
- Shipped makes **alpha-KG the stronger** binder (1.4e-4 vs 5.55e-4). **All three locally-held
  secondaries reverse it** — Handbook T8.5 (3.0e-4 / 5.0e-4), Handbook Fig 8.3 "Barbe 2000"
  (2.0e-4 / 6.6e-4), UWC T17.2 (1.4e-4 / 4.9e-4). **Two of them are the sources D-130's footnote
  named as corroborating**; that footnote was never checked against their numbers.
- **NO VALUE WAS MOVED and none should be** without Burroughs & Sparks 1973 (J. Sci. Food Agric.
  24:187-198), **not held locally**. Choosing among the secondaries is a guess.
- **The bands ALREADY span both alternatives** one parameter at a time — the defect was the
  *sentence*, which asserted the ordering as fact while the overlapping bands, sampled
  independently, reverse it in **1.996 % of 200k draws**. Measured on **draws, not edges**.
- **Never quote the reversal corner's −14.5 % as its cost** — the ordering-preserving corner moves
  **+11.8 %**, so that span is the **joint band width**. The actually-cited table moves it −1.8 %.
- **A straight transposition is FORBIDDEN by pyruvate's own band low** (1.4e-4 < 1.5e-4, by 6.7 %).
  That is why the swapped-suite arm produced 461 failed / **401 errors** — a parameter **load**
  failure, not physics. Never read that red as a consequence.

Receipts: `M:\claud_projects\temp\ferment\d190-adduct-release\` — `PREREGISTER.md`, `RESULTS.md`,
`probe1..5.py`, and the control / mutated / final suite runs.
Related: [[feedback-a-band-is-per-parameter-a-claim-is-joint]],
[[feedback-pin-the-band-not-the-nominal]],
[[feedback-pair-the-red-with-an-ordering-preserving-baseline]],
[[feedback-a-notes-field-is-unchecked-storage]].
