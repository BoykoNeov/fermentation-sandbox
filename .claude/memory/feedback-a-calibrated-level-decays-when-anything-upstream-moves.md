---
name: feedback-a-calibrated-level-decays-when-anything-upstream-moves
description: "A constant defined as 'k set to land X' is a claim about the MODEL, so it rots whenever anything upstream moves — and direction/ratio guards are scale-free and cannot see it"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: f44fbd81-53a8-4e00-ac37-9275b8fd1ab0
  modified: 2026-08-18T22:44:34.893Z
---

Seven beer constants carry a `conditions:` field saying **"k set to land finished X at ..."**. That is
not a fact about the constant, it is a claim about the whole model, so it decays the moment anything
upstream changes. Across D-211 → D-223 three beats moved beer's growth rate and its uptake rate, and all
seven landed levels moved with them — the five Ehrlich higher alcohols by **×2.87** and then back to
**×1.68** against four published means the model had been reproducing to <0.5 %, the two esters by
×0.79. **The 1842-test suite went red exactly ONCE in the process**, and only because an unrelated beat
had happened to pin a packaged number four commits earlier.

**Why nothing saw it:** every guard on those pools was a **direction** ("fusels rise with T"), an
**ordering** (`E_a_esters > E_a_uptake`), a **ratio** (isoamyl acetate ≈ 11 % of ethyl acetate) or a
**temperature response**. All of those are **scale-free** — invariant to a common factor — and the drift
was a common factor. A defect that multiplies everything is invisible to every test that divides.
Sibling of [[feedback-a-ratio-guard-cannot-see-a-common-factor]], one level up: there the blind guard was
one assert, here it was the entire coverage of a subsystem.

**How to apply:**
- When a parameter's provenance says it is **set to land a level**, that sentence is a specification.
  Write the test that reads it, as an equality against the **TARGET** — never a snapshot of what the
  model currently produces, which records the drift instead of catching it.
- Before shipping a change to any rate the model integrates, ask **which calibrated levels are
  downstream of it** and measure them. "What did NOT move" lists are worth writing, and D-211's and
  D-222's both omitted the same five pools.
- Suspect a **common factor** whenever a subsystem's tests are all shapes. Sweep the git history with
  one probe across many commits (worktrees + `sys.path`) — the *when* is usually one commit and it names
  the mechanism for free.
- The sharpest statement of such a drift is often a **threshold crossing**, not a mean:
  48.261 mg/L against a sourced 50 mg/L is "the model claims a defect note", which is a fidelity
  failure; "1.61× a survey mean" is only a number. [[feedback-a-margin-is-a-claim-about-what-holds-it-open]]

Instance: D-224. See [[project-fermentation-sandbox]] and `.claude/memory/prohibitions/beer-aroma-calibration-levels.md`.
