---
name: feedback-discuss-disagreements
description: Surface and discuss disagreements with specs/handoffs before building
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
---

The user wants design disagreements **surfaced and discussed**, not silently
followed or silently overridden. They said the handoff doc is not rock-solid and
to raise anything I disagree with.

**Why:** this is a research project where the architecture decisions carry the
value; the user wants to weigh trade-offs, not discover them after the fact.

**How to apply:** when a spec/handoff conflicts with sound engineering, state the
disagreement concisely with rationale and a recommendation, and use a quick
question for the few choices that are genuinely the user's (and have strong
defaults otherwise — proceed and note them). Record resolved design decisions in
`docs/DECISIONS.md`. See [[user-boykoneov]] and [[project-fermentation-sandbox]].
