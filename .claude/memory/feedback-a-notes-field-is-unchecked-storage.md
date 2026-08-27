---
name: feedback-a-notes-field-is-unchecked-storage
description: "Prose in a provenance note is pinned by nothing — three numbers were wrong in one beat, all in notes, none in a shipped value"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 5f69596c-babc-41fb-a985-7a1c89b328ef
  modified: 2026-08-12T00:03:31.076Z
---

**A `notes:`/`source:` prose field is the least-checked storage in the repo. Re-open the SOURCE
before leaning on a citation already in the file — the note is a claim of the same species as the
one you are building, not evidence for it.**

D-187 found **three wrong numbers in one sitting, every one in prose and none in a shipped value**:

- **D-136's transposition, 51 records old.** Oliveira 2013's introduction gives Godden's natural
  cork at **9.33 mg/stopper/yr** and Lopes 2005's synthetics at **1.60 mg in the first month**, in
  consecutive sentences. D-136 filed each under the other. It survived because **nothing ever reads
  a note** until a later beat needs exactly that sentence — and the number that beat needed was the
  swapped one. Correcting it destroyed a claimed corroboration (the entry now excludes Godden **on
  method**) and created a real one it had never had.
- **Two of my own, written the same day**: "~2.7 years of steady ingress" that was actually **8.9**,
  and "both of P2's caliper means fall inside the band" when one falls **3.5 % outside**.

**The same species of storage lives in tests: a docstring, and above all an xfail REASON STRING.** D-247 found three passing tests whose prose told a future reader that rescaling the availability gate "to the pool it actually gates carries propanol over the floor" — a claim that beat had just measured false. No assertion in any of them could ever have caught it: a reason string is read by humans and printed by pytest, and checked by nothing. Repair it in the beat that disproves it, not the next one.

**Why:** a shipped `value` is pinned by a test, re-derived when a band moves, and re-read whenever a
Process changes. The prose around it is written once, at the moment of most enthusiasm, and then
quoted forward as if it had been checked. Its error rate is not lower than the value's — it is
simply unmeasured, and it laundered into `source:` exactly as [[feedback-a-derived-yield-encodes-its-rate-law]] describes.

**How to apply:** when a new beat leans on a citation already in a YAML note, **fetch the source
again** rather than trusting the note ([[feedback-re-read-the-source-you-already-mined]],
[[feedback-transcribe-tables-not-prose]]). When writing new notes, pin the prose's numbers too —
and write the test against the **published literals**, never against the shipped value: that is the
one construction that can disagree with the file it is checking, and it is what caught both of
mine within minutes. See also [[feedback-pin-the-band-not-the-nominal]].
