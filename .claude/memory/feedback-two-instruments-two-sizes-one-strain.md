---
name: two-instruments-two-sizes-one-strain
description: Two size sources for the SAME strain disagreed 4.6x in volume because they measure differently; never take a mass from one and a size from the other
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 7d706f68-7f41-4d89-b462-b170ebb9e88d
  modified: 2026-09-05T20:21:56.107Z
---

Cell **size** is not one quantity. A microscopy volume and a flow-cytometry *apparent* diameter
(forward scatter calibrated on beads) are different instruments' answers, and on the **same
laboratory strains** they disagreed by **4.63x (haploid) and 5.24x (diploid) in volume**: Klis
2014's haploid is 44 fL (4.4 µm equivalent sphere) where Okada 2023's is 7.3 µm apparent, i.e.
204 fL.

D-271 was about to use cell size to judge whether D-230's branch 1 — a brewing cell of 70.9-91.9
pg — is plausible. Through sourced constants that demand is a **188-244 fL** cell, which reads as
**2.26-2.93x a diploid** in one frame and **about one haploid** in the other. The route does not
merely fail to settle the question; it **answers both ways**, chosen by which paper you open.

**Why:** a derived quantity inherits the frame of every input, and mass-from-one-source ÷
size-from-another silently crosses two instruments. The error is invisible because both numbers
are sourced, both are correct in their own paper, and the arithmetic between them is trivial.

**How to apply:** before combining two sourced quantities about the same object, check they agree
where they **overlap** — here, the same strains. If they do not, say the route is frame-broken and
pin it with a guard rather than picking the convenient one; D-271 ships two (the frames disagree
by >4x; the demand lands in both places at once). A source from the "right" domain is no escape:
Okada's brewing strain is inside the frame that makes everything large.
[[feedback-a-mole-table-is-two-mass-frames]], [[feedback-agreement-can-be-a-frame-difference]],
[[feedback-a-units-fork-is-not-a-band]]
