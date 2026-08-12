---
name: feedback-same-species-different-reaction
description: A source that names your species may be about a different REACTION; a negative result for one role is not a negative for another
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c5f901db-8b1f-4b2e-a63e-b0739efff44d
  modified: 2026-08-12T13:56:20.152Z
---

A passage that names the species you are asking about is not automatically about the
**reaction** you are asking about. D-197 predicted the hydroperoxyl item would close as a
"documented no" because *Understanding Wine Chemistry* says "an attempt to detect this species
was not successful". D-198 found that sentence is about HO₂· as the intermediate of the
**initiation** node `Fe(II) + O₂ → Fe(III) + HO₂·` (index p. 327-328), while the open item was
the **ethanol-limb co-product** four pages later (index p. 332). Same chemical species, two
roles, one failed detection — and it belonged to the role the model had *already* declined.

**Why:** name-matching is what a grep does, and a grep is what produces the shortlist. The
species name is the least discriminating thing in the sentence. A negative result carries a
scope — *detected in what reaction, under what conditions* — and that scope does not travel
with the name. Worse, a mis-scoped negative is self-sealing: it closes the item, so nobody
re-reads the passage. This one was one commit from being recorded as the refusal's basis.

**How to apply:** before letting a passage settle a question, name the **reaction** it is about
and check it is yours — reactants, products, and where in the chain it sits. The index is the
cheapest discriminator a book offers: two entries with different page numbers are two subjects,
however identical the words. Also check whether your own model has already taken a position on
the passage's reaction — if it has, the passage cannot be new evidence about a different limb.
And when the source *is* about your reaction but silent on the part you need, record it as
**silent, not negative**: UWC's Figure 24.13 draws the step with no co-product at all, which
answers nothing. Related: [[feedback-re-read-the-source-you-already-mined]],
[[feedback-a-scope-note-can-carry-a-mechanism-claim]], [[feedback-check-the-blocker-is-still-blocking]].
