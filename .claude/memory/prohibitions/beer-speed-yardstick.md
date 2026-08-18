---
name: beer-speed-yardstick
description: "D-218 - the published trial adjudicating beer's two speed anchors finds the brief too slow at every reading, and the survival of the 5-7 d window turns on a per-cell dry mass this repo ships two values of, 5.6x apart"
metadata: 
  node_type: memory
  type: project
  originSessionId: bd744d4f-5794-466f-81ff-deb7e49f9bff
  modified: 2026-08-18T09:13:03.318Z
---

**Live prohibitions — which yardstick beer's fermentation SPEED is calibrated against (D-218).**
Detail split out of `.claude/memory/project-fermentation-sandbox.md`; that file's ledger points
here by path and does **not** index it in `MEMORY.md` ([[feedback-measure-which-end-is-growing]]).

## What the beat settled

D-216 §11 / D-217 §8 left the owner a choice between §2.2's **5-7 d acceptance criterion** and
Tyrell's **measured extract course**, with nothing to decide it. The owner authorised a search.
**It is now decided on the literature and open only on a conversion factor.**

The adjudicating source is **Foster et al. 2022, *Front. Microbiol.* 13:747546**
(doi:10.3389/fmicb.2022.747546; full text in `kveik.xml`, Europe PMC). It is the first source the
archive has carrying the whole tuple: **12.5 °P = 1.045** hopped all-malt (within 0.003 SG of
§2.2's own wort), a **counted** pitch of 1.2 × 10⁷ cells/mL, two temperatures for the **same three
Beer 1 clade strains**, and the **same 1.010 target gravity**, which the paper takes from Parker
2008 — so the brief's **finish line** turns out to be independently sourced even though its
**duration** is not.

Corroborating: Rodríguez-Guerrero 2022 (*Foods* 11:3602, four pilot batches at **54-96 h**),
O'Brien 2024 (*Modelling* 5(1) 201, IPA primary **3-5 d**), and Reid 2021 (*Fermentation* 7(1) 13,
a fitted logistic for Speers' industrial lager series).

## Never do these

- **Never cite Foster's two temperatures as a test of this model's temperature response.** It is
  the only same-yeast two-temperature pair the archive has ever had, and it **discriminates
  nothing**: the bound is cleared at *every* printed value of `E_a_uptake` (1.59 at the low edge,
  2.21 shipped, 2.44 at the high), and the first value that fires is **90,000 J/mol, 1.43× out of
  band**. D-216 §10 / D-217 §6's rule about out-of-band mutation arms, stated about the
  *literature*. **D-217's `E_a_uptake` refusal therefore needs no re-opening.**
- **Never read Foster's "3 days" or "within 10 days" as durations.** Both are **ceilings**.
  Figure 2's panels are **timepoints** (12/48/72/120 h) with temperature on the x-axis, so
  *"after only 3 days ... (Figure 2C)"* cites the **72 h sample**; and *"within 10 days"* at 12 °C
  **is the incubation length**. An upper bound on the duration ratio needs a *lower* bound on the
  22 °C value, which the paper does not give — read at the open end the bound is 5.0, not 3.33.
- **Never treat the ratio ≈ the bare Arrhenius factor as structural.** At the shipped value they
  agree to 0.7 %, but the residual runs 1.104 → 1.007 → 0.986 across the band and **crosses 1.0 at
  ≈58,000 J/mol**. The nominal sits 2,900 J/mol short of that crossing, 9 % of the band's width.
  It is a **coincidence**, the same trap D-217 §2 named around −90,000 J/mol.
- **Never read the ratio's pitch-invariance as a pass.** It is **local** (0.007 at the nominal,
  0.10 at the low edge) and it is a **defect signature** — real ferments do show a pitch-dependent
  temperature response. The apparent *perfect* invariance was the **hourly grid's quantum**: on a
  6-minute grid the nominal spread is +0.0070, not −0.0002. [[feedback-read-a-fast-curve-on-a-fixed-grid]]
- **Never retire §2.2's 5-7 d window on this evidence.** It passes today at 6.08 d. Retiring it
  means moving `q_sugar_max`, and D-216 §4 already priced that: a rate satisfying Foster finishes
  **a day early** against both measured tails while **still missing both measured day 2s**.
  Trading an acceptance criterion the model passes for a rate that fits the *shape* worse is not a
  fidelity gain. [[feedback-fit-the-observable-not-the-consequence]]
- **Never table `q* = 1.500` as a fitted crossing.** At 18 pg/cell it is the printed band's
  **ceiling** — a saturation. Even there the model needs 3.04 d, missing Foster by 1.4 %. Cf. D-177.
- **Never revive the withdrawn 10 °C isothermal row.** Reid/Speers say *"**starting** temperature
  of 10 °C"* — a free-rising industrial lager. Running it as a set-point is **D-217 §3's failure
  mode one source over**, and it produced a headline (17.5 d, 6.9× day-2 shortfall) and a
  mechanism claim (`E_a_uptake` is the defect) that are both **refuted, not softened**.
- **Never record MDPI as a paywall.** It blocks WebFetch *and* curl at the host but reads fine
  through the browser tool — a blocked **transport**. [[feedback-paywalled-is-one-host]]

## The fork that is now the whole question

**This repo ships two per-cell dry masses, 5.6× apart, and neither is sourced to a measurement:**
**18 pg/cell** in `tests/benchmarks/test_validation_varela2004.py` and `..._palma2012.py` (*"the
standard ~18 pg/cell S. cerevisiae dry-weight figure"*, used for both wine pitches) against
**~100 pg/cell** implied by `beer_generic.yaml` (`TYRELL_SCENARIO`'s 1.0 g/L vs Tyrell's counted
9.96 × 10⁶ cells/mL). D-216 §7 sized it as *"~2× the textbook 40-60 pg"* — true of the textbook
range, and it **understates the repo's own internal disagreement by nearly 3×**. D-218 **Flags**
D-216 for it; both values still ship.

Bracketed over both readings and both ends of Foster's sampling interval, §2.2's window survives
in **exactly 1 of 8 cells** — the one needing the reading the archive itself calls anomalous
**and** the 72 h sample read as exact. At 100 pg the shipped model reaches Foster's endpoint in
**3.33 d against a published ≤3 — an 11 % miss with no knob touched.** The model is not badly
wrong about beer's speed; **the brief is.**

**Adopting 18 pg/cell is not a one-line change and D-218 did not price it:** `TYRELL_SCENARIO`'s
pitch goes 1.0 → 0.179, and D-211's `mu_max` was fitted at 1.0. A beat that picks it up inherits a
`mu_max` refit **and** a broken extract calibration.

## What the evidence does establish: the early limb, again

Not the rate scale, not the temperature response. The day-2 shortfall survives **every defensible
pitch** (closing it needs ~2.7 g/L ⇒ ~270 pg/cell, **5× textbook**) and **every rate respecting an
endpoint**. Every measured or fitted course peaks on **day 1-2** and tails; the model peaks on
**day 4**. Speers' own authors needed the **asymmetric** 5-parameter logistic on 3 of 7 industrial
datasets — a degree of freedom this engine does not have. That is the remaining defect and it is
structural. **Do not call Reid's fitted lager curve a second measured course**: a 3-parameter
logistic is symmetric and peaks at its inflection *by construction*.

## Guards (tests/test_organic_acids.py §12)

Four tests, ~18 s, falsified against **five in-band arms with a designed GREEN**. `E_a_growth` →
low edge fires three of four; `mu_max` → high edge fires two; `K_repression` → high edge fires the
fork only; `E_a_uptake` → 60,000 and `q_sugar_max` → 0.6 move nothing. **The bound guard's column
is empty by construction — that IS the finding**, and its positive control (the out-of-band
90,000 arm) lives inside the test. A first pass used `E_a_growth` → 0, which is **out of band and
died in pydantic**; third time that harness fact has cost an arm (D-216 §10, D-217 §6).

## Still open

**Source the per-cell dry mass.** Load-bearing in three places (D-216 §7, D-216 §8, D-218 §4) and
cheaper than anything else. Also unmined: de Andrés-Toro's five isothermal runs (8/12/16/20/24 °C,
every open route 403s — the *data*, separate from the fitted coefficient D-217 refused), and
Foster's Figure 2, which holds six more temperatures for these same strains.

Receipts: `M:\claud_projects\temp\ferment\d218-beer-speed-yardstick\`.
