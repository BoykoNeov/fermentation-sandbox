---
name: feedback-patch-the-binding-the-code-reads
description: "`from X import name` copies by value, so patching X leaves the consumer unchanged — and every assertion downstream still passes while measuring the un-patched arm"
metadata:
  node_type: memory
  type: feedback
---

**A test that re-points a module global must patch the binding the *consumer* holds, and must
carry an assertion that fails if the patch missed.** `from carbon_routing import
LABELLED_PRECURSOR` binds a copy into `byproducts`; setting the attribute on `carbon_routing`
afterwards changes nothing the run reads.

D-256. Rollero published two tracer experiments on one medium — valine and leucine — and the
model ships the valine arm. The leucine arm was reached with no production change by
`monkeypatch.setattr(byproducts, "LABELLED_PRECURSOR", "leucine")`. Aimed one module upstream at
`carbon_routing`, the run comes back **byte-identical**, and the substantive assertion — that the
leucine branch under-attributes against the paper's floor — **still passes**, because the valine
arm is under that floor too. The test would have reported a finding about a route it never ran.

**Why:** a wrong patch target does not raise. It produces a green test measuring the wrong thing,
and the failure mode is invisible precisely when the two arms are qualitatively similar — which
is the normal case, since both arms are usually the same physics with a different label. The
substantive assertion cannot detect it; only a comparison between the arms can.

**How to apply:** assert the two arms read *differently* (`abs(a − b) > tol · b`), not merely that
each lands where expected — the anti-vacuity arm is the load-bearing half of the test and belongs
in its docstring as such. Verify it by mutation: aim the patch at the wrong module and confirm
*only* that assertion fails. Prefer reading the name through the module (`mod.NAME`) in code you
expect tests to steer, so the seam exists at all. Related:
[[feedback-an-xfail-buries-the-asserts-after-it]],
[[feedback-a-toggle-measures-nothing-on-an-exhausted-pool]] (its sibling: there the toggle ran but
the instrument could not see it; here the instrument works and the toggle never ran).
