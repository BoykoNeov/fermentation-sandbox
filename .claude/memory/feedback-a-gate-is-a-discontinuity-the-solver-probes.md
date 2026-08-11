---
name: feedback-a-gate-is-a-discontinuity-the-solver-probes
description: A conditional inside a solver-called function is a jump the Jacobian probe straddles; the gated version broke FEWER tests than the continuous one and was the worse option
metadata:
  type: feedback
---

Adding an `if` to anything `solve_ivp` evaluates makes the right-hand side
**discontinuous**, and BDF's `num_jac` takes differences **across the jump** — those
Jacobian entries are wrong, not merely noisy. Measured at D-182 on one 1-year run:
`charge_residual` was called 63,237 times and **34,950 of those saw the gated term
switched on**, because `num_jac` perturbs every state slot in turn and a gate testing
`> 0.0` flips when a zero slot is nudged to **1.49e-17**.

The trap is that this looks like *success*: the gated version broke **4** aging pins and
the continuous one broke **18**. Fewer reds meant a broken Jacobian nobody would look at;
more reds each named a real, continuous, measured change.

I also built the gate on a number I had not attributed. An un-anchored wine solved to
pH 2.92 and I read that as my new term's doing. It was `Byp`'s — that wine solves to 2.92
**with or without** the term (2.9242 vs 2.9217), so the term moved a fiction by 0.0025 and
fidelity was a wash on exactly the states the gate covered.

**Why:** a red count is not a quality metric, and the cheaper-looking branch can be the
defect. When fidelity is a wash between two forms, **numerical hygiene decides**, and it
decides for the continuous one — there is nothing for a probe to straddle.

**How to apply:** before adding a conditional to a function the solver calls, ask what
`num_jac` sees when it perturbs the slot the condition reads. If the answer is "the branch
flips", the gate is a defect. Instrument it — wrap the function and count how many calls
land on each side [[feedback-count-and-print-your-skips]] — rather than inferring it from
which arm has fewer failures. And run the with/without arm before attributing any number
to your own change [[feedback-pair-the-arm-with-its-baseline]]: a pre-existing pool
(`Byp`) can own the entire effect you are about to design around. Related:
[[feedback-pair-the-red-with-an-ordering-preserving-baseline]],
[[feedback-rejected-values-must-be-unreachable]].
