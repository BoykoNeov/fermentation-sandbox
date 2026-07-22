---
name: feedback-verify-latest-state-not-breadcrumbs
description: "Before proposing/starting work, verify against the LATEST decision + the code — stale 'Next:' breadcrumbs in old progress notes get burned down by later decisions"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 14c4cd91-8005-4d3b-8620-c315b1e7cf0b
  modified: 2026-07-20T17:55:31.025Z
---

Before proposing a task as "open" — or starting it — **verify it against the
latest decision state and the actual code**, not against a "Next:" list written
in an older progress note or DECISIONS entry. Old backlogs are burned down by
the decisions that follow them.

**Why:** In one session (2026-07-20, at D-121) I proposed TWO already-built
tasks in a row from stale breadcrumbs. (1) "beat 1b (descriptor projection)" —
built as slice 1 (D-95) + slice 2 (D-98), and D-98 *explicitly retired* the
"Next: beat 1b slice 2" item. (2) "speciate the fusels lump" — done at **D-99**;
D-110 states "No pool in the project is lumped any more." Both were named in
plan/DECISIONS "Next:" tails that predated the decision that closed them (the
D-97 progress note and D-98's "Next:" both list "speciate fusels," and D-99 —
the very next decision — did it). Reading a single older entry's "Next:" as the
current backlog is the trap. This project's lumps are ALL gone: esters→3 (D-96),
fusels→5 (D-99), amino_acids→8 (D-100), mercaptans was a false lump = methanethiol
(D-110). See [[project-fermentation-sandbox]].

**How to apply:** When picking or proposing "what's next," treat old "Next:"
lists as leads, not truth. Confirm each candidate against (a) the most recent
DECISIONS entries and (b) the code — grep for the state slot / Process / param
it would add and check it isn't already there. Only present a task as buildable
after that check. A "Next:" tail is only reliable in the *latest* decision;
every earlier one has been partly consumed. Cross-check with
[[feedback-measure-which-side-before-building]] (measure before building).
