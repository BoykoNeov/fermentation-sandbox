---
name: feedback-a-guard-must-not-fire-on-an-improvement
description: When a record prices a cost, the test that pins it must assert the ordering that caused it, not the physical line it crossed — otherwise removing the cost turns the suite red
metadata:
  type: feedback
---

A record measures that some change is expensive; the test pins the expense by asserting the
quantity sits **past** a physical threshold. That test now requires the expense. Any later change
that reduces it — which is the outcome everyone wants — fails the suite, and the failure message,
written when the cost was the finding, tells the reader to re-price rather than to celebrate.

**Why:** D-251 priced a front-loading calibration in stored nitrogen and asserted the loaded cells
read *above* the N-replete elemental reference of 0.114 g N/g DW, with a message saying that the
calibration being cheap "is the reason D-251 left it to the owner — re-price it rather than
deleting this." When D-253 took the calibration, that message became an instruction to do the
wrong thing, and the assertion's direction meant a future Droop quota or a smaller store would go
red for storing *less*. The docstring around it already said the loading was "the footprint, not a
veto"; the assertion contradicted its own docstring, and nothing caught that because both were
green.

The same beat showed the second half: the *upper* bound (`loaded <= 0.14`) held at the shipped
capacity and fails at the capacity parameter's own band edge, where the transient reaches 0.322.
A bound satisfied at one point in a declared band is a statement about that point, not a property
of the model, and an ensemble drawing over the band will find the difference.

**How to apply:** pin the **ordering that produced the cost** — the calibrated arm against its own
control — rather than the threshold it crossed. Keep the physical line only where crossing it
would genuinely be inadmissible, and label it as a statement about the shipped value. When a cost
is later removed, the ordering guard reports it as a structural surprise worth reading; the
threshold guard reports it as a failure. And re-read a guard's **failure message** whenever its
record's verdict changes: the message is prose, it is never executed, and it is the part that ages
worst.

Related: [[feedback-a-guard-that-hardcodes-an-input-cannot-price-it]], [[feedback-price-a-transient-against-a-ceiling-of-the-same-kind]].
