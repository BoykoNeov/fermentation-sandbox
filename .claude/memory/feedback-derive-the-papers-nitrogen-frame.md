---
name: feedback-derive-the-papers-nitrogen-frame
description: A paper's headline "180 mg N/L" can be in an accounting frame your model does not share; derive the frame from the paper's own numbers and check it reproduces the printed total
metadata:
  node_type: memory
  type: feedback
---

Two papers printed "180" and "300 mg N/L of assimilable nitrogen". Neither is the number to type
into this model's `yan_mgl`. Crepin's own accounting counts **arginine at 3 N, tryptophan at 1,
histidine at 1** -- the assimilable subset, since the indole and imidazole nitrogens are never
released -- while `nitrogen_mass_fraction` counts every atom because the conservation ledger must.
Same must, 179.91 in the paper's frame and **201.11** in the model's: a 12 % gap that is not
rounding. Crepin's "180" is also **consumed** N, not the must's content.

**Why:** typing the paper's headline into a field that means something else is the same class of
error as scoring a model on a must 2.25x the source's nitrogen (D-244 section 6), just inverted --
it strands 21 mg N/L of real assimilable nitrogen and looks commensurate while being wrong.

**How to apply:** never transcribe a headline aggregate. Take the per-species numbers and
**derive** the convention: divide the paper's per-species consumed mg N/L by its per-species mM
and read off the integer N count. Then check the derived frame reproduces the printed total -- if
it lands within a fraction of a percent it is the paper's frame, not a fit. Recompute the
aggregate in YOUR model's frame from the composition, and declare that. Also worth naming when
found: one field carrying two meanings ("total N" and "assimilable N") is itself the defect.
