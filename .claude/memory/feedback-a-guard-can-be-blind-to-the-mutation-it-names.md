---
name: a-guard-can-be-blind-to-the-mutation-it-names
description: "An isolability test whose setup removes the confounder along with the subject tests the gate, not the wiring; run the mutation its docstring names"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1a678dd9-ee02-4677-b675-5fe9b7f13673
  modified: 2026-08-31T15:58:01.960Z
---

A test whose docstring names the defect it guards against can be **structurally incapable of
catching it** — when the setup that removes the subject also removes the confounder. Nothing goes
red, the docstring reads well, and the claim was never tested.

**Why:** at D-255, `test_no_amino_acid_dose_means_no_label_anywhere` said in prose that it guarded
against a tracer credited off the **leucine** branch (which feeds the same alcohol pool). Its arm
ran with the amino-acid dose at zero — but **with no dose there is no leucine either**, so the
mutant reads zero and passes. The arm tested the compile gate, never the wiring, from D-115 to
D-255. Running the mutation is what exposed it: the four original asserts all PASSED and only the
new arm went red.

**How to apply:**
- **Run the mutation the docstring names**, not one that merely reddens the file
  ([[feedback-run-the-mutation-the-claim-names]], [[feedback-check-the-published-test-can-fail]]).
  Then read *which assert* fired — "3 failed" is not the same as "the arm I built failed".
- The fix is an arm where the subject is removed **and everything else is kept**: full precursor
  spectrum, valine alone emptied, still a real ferment (isoamyl 0.1834 vs 0.1846 g/L). Assert the
  run is still live, or the sharp arm silently degrades into the blunt one.
- **A zeroed dose is a gate test; a zeroed species is a wiring test.** Prefer the latter, and keep
  both when the gate is also a claim ([[feedback-a-toggle-measures-nothing-on-an-exhausted-pool]]).
- Sourcing a real medium can be what makes the sharp arm *buildable* — before it there was no way
  to hold every other pool at a defensible level.
- A mutation that dies in a **schema validator** before reaching the physics is a red for the
  wrong reason and proves nothing about the assert
  ([[feedback-a-fast-gate-masks-the-slow-one-behind-it]]); build a gentler one that lands.
- Verify a restore with `git status --porcelain`, **not** a byte-diff against a `cp` snapshot —
  under autocrlf the snapshot differs on every line and the check lies
  ([[feedback-crlf-join-inflates-line-count]], [[feedback-verify-the-restore-between-mutation-arms]]).
