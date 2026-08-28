---
name: feedback-a-paper-can-print-the-same-numbers-twice-differently
description: "Results and Methods can give the same sampling landmarks different values, and the protocol constant you need may not be in the paper at all — check both sections and say which one a transcribed constant came from"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 6025169e-c573-4f1a-af77-db80fef3cfd0
  modified: 2026-08-28T07:45:16.737Z
---

One paper printed its four sampling landmarks **twice and differently**: Results said "16, 20,
40, and 110 h", Methods said "(16 h, 20 h, 28 h, and 150 h)". A constant in the repo had been
transcribed from one of them with no note that the other existed, so a later beat scoring against
it could not tell whether a 1.5× miss was really 2.2×. The same paper never stated its
**inoculum at all** (`grep -c inocul` = 0) — a load-bearing protocol constant simply absent,
while the fixture quietly carried a house default 6.25× larger than the only sourced value for
that medium (its sibling paper's).

**Why:** a transcribed constant records the value, not the section it came from, so an internal
contradiction in the source becomes invisible the moment it lands in code. And absence is worse
than contradiction: a missing protocol constant leaves no gap to notice, because a scenario input
always has *some* value — usually a convention inherited from elsewhere
[[feedback-a-limitation-can-belong-to-its-frame]].

**How to apply:** grep the full text for each landmark/protocol number rather than reading the
one section that has it, and when they disagree, transcribe the **Methods** value (it states what
was actually done) *and* record the conflict next to the constant, with which direction each
reading moves the result. Before scoring anything against a paper's must, `grep -c` its own text
for the protocol constants you are assuming — inoculum, temperature schedule, vessel — and treat
a zero-hit as an open question, not a default. Transferring one from a sibling paper is legitimate
but is a **trade**: price what it costs on the other side before adopting it (here it fixed the
timing and broke a shipped biomass anchor). Related:
[[feedback-reproduce-a-published-number-before-trusting-the-new-column]].
