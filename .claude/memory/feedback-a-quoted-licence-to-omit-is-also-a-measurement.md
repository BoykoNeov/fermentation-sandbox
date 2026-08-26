---
name: a-quoted-licence-to-omit-is-also-a-measurement
description: A source sentence quoted in the code as a reason to leave something out usually also SIZES what was left out - read the whole sentence
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 80611b07-f4ed-492d-8b4b-03c3866a0e21
  modified: 2026-08-26T18:40:47.513Z
---

A sentence you already quoted as **licence to omit** is usually also a **measurement of the
omission**. Read both halves before treating the subject as closed.

`M_GLUTAMIC`'s docstring has carried Peyer §5.5 verbatim since D-179: *"most amino acids
contribute only poorly to the BC in the relevant pH range … Little contribution … is made by
aspartate (pKa 3.86), glutamic acid (pKa 4.25) and histidine (pKa 6.04), which account for ca.
10 % of the total BC at a wort pH of 5.5."* Its first half is this repo's own geometry argument
for refusing malt phosphate (D-178), in the source's words. The archive kept that half, read
"little contribution" as negligible, and modelled beer's buffer as peptides alone.

The second half **names three exceptions and sizes them at 10 %**. Computed from a composition
table already in the repo, those three are 8.6 % of the model wort's local buffering — and the
back-solved peptide lump had silently absorbed their share, so the wort buffered correctly at
t=0 and went on buffering after the yeast ate them (D-239).

**Why:** "the source says X is small" is a claim about MAGNITUDE, never about FATE. A small pool
that LEAVES is not the same modelling object as a small pool that stays, and a lump fitted to a
total is right in total while wrong in both.

**How to apply:** when a docstring or a parameter note quotes a source as the reason something
is not modelled, re-read the quoted sentence in the source and ask two questions it usually
answers: *how big is the thing it dismisses*, and *does it share the fate of what absorbed it*.
Grep your own files for the quote first — it may already be transcribed
[[feedback-a-transcription-answers-more-than-its-purpose]],
[[feedback-re-read-the-source-you-already-mined]].
