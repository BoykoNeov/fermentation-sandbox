---
name: feedback-price-the-refusal-not-just-the-error
description: "Before calling a route merely wrong, check whether it RUNS — a banded parameter can refuse it outright, which changes a two-way decision into a one-way one"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: b530ceda-f935-4502-a0e8-cac9593dd384
  modified: 2026-08-27T12:45:52.812Z
---

When you price a candidate repair as "possible but wrong", **check that it actually executes**
before writing that down. A quantity that looks like a bad number may be an *unconstructible* one:
in this repo every `Parameter` validates its value against its own band at construction, so a
derived value outside the band raises rather than shipping a wrong result.

**Why:** D-243 priced "evaluate the yield fit at total assimilable nitrogen" as an extrapolation
reaching an unphysical 38 % cell nitrogen. D-244 ran it and found the band refuses it at a total of
444.0 mg N/L — **35 suite scenarios would simply stop compiling**. That is not the same finding.
"Wrong" is a trade-off the owner weighs against a cost; "does not run" removes the option and
leaves a one-way decision. The same probe also exposed a ceiling nobody had documented: no wine
declaring over ~444 mg N/L had EVER compiled, with an opaque error naming neither the quantity nor
its source.

**How to apply:** the check is cheap — construct the derived object, or compile one scenario at the
extreme of the range, before the decision memo is written. Extend it to the sibling question: if a
route is impossible, is something *legitimate* impossible for the same reason? A band that blocks a
bad repair usually also blocks a real case, and that case is a finding in its own right. Related:
[[feedback-a-parameter-can-be-pinned-and-drawn]], [[feedback-verify-latest-state-not-breadcrumbs]].
