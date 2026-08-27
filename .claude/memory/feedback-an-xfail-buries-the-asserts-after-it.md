---
name: feedback-an-xfail-buries-the-asserts-after-it
description: D-245 - marking a multi-assertion test xfail silently kills every assertion after the failing one; split so the mark carries ONE claim
metadata:
  node_type: memory
  type: feedback
---

**An `xfail` is a mark on a whole test, but a test stops at its first failing assert — so every
assertion after that one stops being a guard and nobody is told.** At D-244 I marked
`test_the_d104_sink_absorbs_the_d103_gate_shape_spread` xfail for its *middle* assertion (all five
alcohols inside Rollero's low band). Its **third** assertion — propanol's sink-on share is under
half its sink-off share, which is the D-112 finding's own payoff — was never reached and never run
for the whole of D-244. The same beat did it twice:
`test_the_de_novo_cap_is_inert_where_the_precursor_exhausts` died on its exhaustion *premise* one
line above the assertion that carried the content, so the measurement that decides whether D-120's
refusal survives went unmade until D-245 ran it by hand.

**Why:** an xfail looks like a *localised* admission — "this one claim is known-broken" — and reads
that way in the summary line, which prints one reason string per test. The mark is not localised at
all: it converts everything downstream of the break into dead code with a green-looking `x`. Losing
a guard is exactly what the strict-xfail idiom exists to prevent, so the failure mode is the tool
undoing its own purpose. Related: [[feedback-verify-an-xfail-fails-for-its-stated-reason]] — that
one is "the RED is not the one you named", this one is "the GREENs behind it are not running".

**How to apply:** before marking a test xfail, **count its assertions**. If more than one can fail
independently, split the test so the mark carries the failing claim alone and the rest stay live
guards — the split is a repair, not churn, and costs at most one extra fixture run. When you
*inherit* an xfail, re-run it with `--runxfail` and check which assert the RED names: any assertion
after that line has been unguarded for as long as the mark has been there, and it is worth reading
them to see what was quietly lost. Same reasoning applies to a premise assert placed above the
content assert — put the cheap premise check in its own test if it can go stale.
