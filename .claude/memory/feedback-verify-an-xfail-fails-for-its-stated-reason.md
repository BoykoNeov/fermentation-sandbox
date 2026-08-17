---
name: feedback-verify-an-xfail-fails-for-its-stated-reason
description: D-215 - an xfail that fails on a typo is indistinguishable from one catching the defect; re-run under --runxfail and read what the RED actually names
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c639a01c-94f5-44c3-a91b-491833d1c5c9
  modified: 2026-08-17T18:40:08.463Z
---

**A strict xfail proves only that something failed — not that YOUR defect failed it.** At D-215 I
shipped `test_the_model_ferments_tyrells_wort_on_tyrells_schedule` as a strict xfail naming a
2.8× fermentation-speed gap. It xfailed. Green suite, correct-looking summary line. Under
`pytest --runxfail` it turned out to be failing on **`AttributeError: 'str' object has no
attribute 'spec'`** — I had called `sugar_species("beer")` when it takes a *schema*. The test had
never reached its assert, so the "expected failure" was expected for the wrong reason, and the
day the real gap closed it would have gone on xfailing forever instead of turning green.

**Why:** an xfail inverts the usual signal, so every ordinary safeguard inverts with it. A typo,
an import error, a wrong fixture and a genuine caught defect all render as the same `x`. The
whole value of the D-208 idiom — *state what is true of the source and false of the model, so a
fix turns it GREEN* — depends on the failure being the named one; otherwise it is a test that can
never pass and never complain. This is [[feedback-grep-finds-claims-not-guards]] in the test
layer: the marker is a *claim* about why it is red, and nothing checks the claim.

**How to apply:** after writing any `xfail`, re-run it with `--runxfail` and **read the assertion
message**. It must be your assert, with your numbers in it. `AssertionError: day 2: the model has
fermented 21.2% ... against Tyrell's measured 59.4%` is a verified xfail; anything that is not an
`AssertionError` from the line you wrote is a broken test wearing an xfail's clothes. Same check
applies when an xfail is *inherited* — a strict xfail that has been red for many beats is worth
re-running this way before it is cited as evidence the gap is still open. Related:
[[feedback-pair-the-red-with-an-ordering-preserving-baseline]].
