---
name: feedback-right-number-wrong-condition
description: "A spec that asserts a value AND the condition it holds at is two claims; refuting the pair is not refuting the value, and collapsing a conditional source to a summary hides which half is wrong"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d9b6cdf3-1baf-48f5-bebb-63068e2c2118
  modified: 2026-08-18T11:20:43.677Z
---

When a specification pairs a **number** with the **condition it is asserted at**, that is two
claims, and evidence can refute the pair while corroborating the number. Say which half died.

At D-218/D-219 the brief's *"a beer wort reaches 1.010 in 5-7 days"* was scored against a
published trial and declared unsupported — twice, the second time explicitly *"unconditional"*.
Both readings came from the trial's **timepoint panels**, which collapse a temperature-resolved
experiment into four samples dominated by its warm arms. At D-220 the same paper's supplementary
figure turned out to hold the **course, at eight temperatures**. Resolved by temperature, the
brief's duration is **corroborated at 15 °C** (all three ale strains land at 5.0-5.9 d) and
**refuted at the 20 °C the benchmark also asserts** (2.9-3.8 d at 22 °C). The number was right.
The pairing was wrong. Two records had called the number wrong.

**Why:** a summary statistic over a conditional source silently marginalises the condition, and
the marginal can carry the opposite sign to every cell. "Nothing supports X" is then a claim
about the *projection*, not about X — and it reads as a verdict on the value, which is what gets
copied forward. The same collapse hid a second thing: with temperature resolved, the model's
error turned out to be near-constant across 12-22 °C, i.e. a **level** error, which vindicated a
standing refusal to touch the temperature term rather than re-opening it.

**How to apply:** before recording "the literature refutes X", ask what X is *conditioned on*
and whether the evidence resolves that condition or averages over it.

If it averages, say so and name the condition as unresolved rather than shipping an
unconditional verdict.

When the pair is what fails, write it as *"the number stands; it is the number TOGETHER WITH
<condition> that cannot"* — the repair differs: one edits the value, the other edits the
condition, and they are different beats with different prices.

Prefer a **bracket** over an interpolation when the resolved evidence straddles the asserted
condition: a bracket argues from measurements only, where an interpolation puts a
hard-sounding point estimate on top of a soft input.

Related: [[feedback-a-summary-statistic-is-not-the-curve]],
[[feedback-re-read-the-source-you-already-mined]], [[feedback-conceded-caveats-are-not-coverage]].
