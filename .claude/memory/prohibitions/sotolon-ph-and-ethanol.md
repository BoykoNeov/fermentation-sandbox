---
name: prohibitions-sotolon-ph-and-ethanol
description: "D-205 — Pham's pH + ethanol terms on the sotolon aldol: REFUSED on identifiability, measured not argued"
metadata: 
  node_type: memory
  type: project
  originSessionId: 00a92466-5b3f-4515-a31f-d00c098d113d
  modified: 2026-08-14T08:13:26.756Z
---

# Sotolon's pH and ethanol terms (D-205) — REFUSED, measured

**Reached by path from the ledger in [[project-fermentation-sandbox]]. No `MEMORY.md` row (D-185).**
The subject is SETTLED: do not re-propose either half as unbuilt. Detail below is *why*.

## The prose that caused this

`SotolonAldolCondensation`'s scope item (1) said Pham measures sotolon rising with decreasing pH
and decreasing ethanol, "and the model has both quantities, so this is a **real omission rather
than an inexpressible one**". True, and **ten records read it as "buildable"**. It is now marked ⚠
in the docstring. **EXPRESSIBLE IS NOT IDENTIFIABLE** — that is the whole record.

## Both halves refused, for different reasons

- **pH — refused on IDENTIFIABILITY, and it is sharper than D-203's.** Two independently motivated
  forms (acid catalysis in [H⁺]; the undissociated-α-KB fraction, pKa ~2.5) **agree with each other
  within 17 %**, so the mechanism choice is nearly free — *that* guess was wrong and is recorded as
  wrong. The free parameter is the **REFERENCE pH**: `k_sotolon_aldol` was calibrated with **no pH
  term**, so no pH is its reference, and the same form at `pH_ref` 3.4 vs 3.0 reports **two
  oxidised wines or none** (2.51× apart, matching the analytic slide to 0.2 %). No observable
  against two parameters, one of which flips a sensory verdict alone.
- **ethanol — BLOCKED on an equilibrium constant (the D-190 shape).** Mechanism IS sourced (UWC
  §9.5: acetaldehyde + ethanol → hemiacetal/1,1-diethoxyethane, which "remove acetaldehyde in its
  free form" = Pham's direction). **No Keq in any of the 24 corpus texts** — the only hits are
  brandy composition tables in mg/L. **Its crossing value is deliberately NOT computed**: it would
  need the Keq *and* a new bound-acetaldehyde pool on the carbon ledger, i.e. inventing the
  constant whose absence is the block.

## Do not re-source these

**Pham 1995 is UNREACHABLE and two hosts were checked** (HAL = Anubis challenge, AGRIS = abstract
only). The abstract is directional, full stop — no range, no coefficient, no response surface, and
**it does not say whether the medium held SO₂**. The only acid-catalysed aldol kinetics in reach is
Casale/Elrod 2007: aldehyde+aldehyde in **40-85 % sulfuric acid** on excess-acidity theory, with
acetaldehyde flagged anomalously slow. **Different reaction, different regime — never cite it here.**

## The emergent limb — measured, and NOT a substitute

pH HAS reached this rate since D-108, through `free_acetaldehyde`'s `bisulfite_fraction`. **It is a
DIFFERENT REACTION** (bisulfite competition for the carbonyl, not acid catalysis of the
condensation) — same direction only. **It does not discharge scope item (1).**

- **1.893 %** across pH 3.0-3.8 in a SEALED sulfited wine (SO₂ in EXCESS over the carbonyls, 99.3 %
  of acetaldehyde bound) — where sotolon is 0.73 % of threshold.
- **0.003 %** once O₂ is dosed: SO₂ becomes the **LIMITING REAGENT** and binds fully whatever the
  pH — and that is the arm sitting at **96.25 % of the 8 µg/L threshold**.
- **Unsulfited is pH-independent BITWISE** (the `total_so2 <= 0` short-circuit). That was the
  negative control and it is what makes every other number attributable.
- Not substrate limitation: α-KB holds ~2.11 mg/L in both arms.
- **The 1.893 % regime is one the code already calls an unphysical isolate** (zero O₂ ingress, no
  closure permeation). Real bottles walk toward the arm where the limb vanishes.

## Costs paid, don't re-pay

**NO guard was owed** — the bitwise pH-independence is exactly what a correct future build must
BREAK, so pinning it would fight the fix; falsified before declining. **D-204's four pins are
untouched** (no rate law moved). The repo has `pKa_pyruvic` but **no α-ketobutyrate pKa**, and the
keto acids sit deliberately outside `PKA_PARAM_NAMES` — so even the "sourced physical constant"
form needs a new parameter.

Receipts: `M:\claud_projects\temp\ferment\d205-sotolon-ph-ethanol\` (`PREREGISTER.md` first, then
`FINDINGS.md` + probes 1-5). Archive: D-205, corrects D-107.
