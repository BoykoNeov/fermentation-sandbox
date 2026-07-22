---
name: feedback-name-guards-for-what-they-forbid
description: "Name a guard/invariant for the violation it catches, not the nearest-sounding real mechanism — a wrong label invites a fair objection and reads as modelling prejudice"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 4c5542b6-994d-42ff-9b5f-a6dbc7d14d50
  modified: 2026-07-20T09:37:51.154Z
---

When documenting a guard, assertion, or conservation check, describe **what it
forbids in the model's own terms**. Do not label it with the nearest-sounding
real-world mechanism unless the guard genuinely forbids that mechanism.

**Why:** at D-117 the `worst_c < 1.0` carbon-refund guard was documented as
forbidding **"gluconeogenesis, which fermenting yeast do not do."** The owner
objected — correctly — that gluconeogenesis *is* real in wine yeast and the model
should not run from reality. But the guard never touched it: the refund lands in
`S`, the **extracellular** sugar pool, while real gluconeogenesis makes
**intracellular** G6P for trehalose/glycogen, is glucose-repressed during active
fermentation, and is never secreted. The guard forbids sugar appearing in the must
from nowhere — a **mass-balance violation**, not a pathway. The correct name was
already in the code: the refund is a **sparing credit** (amino-acid carbon
substitutes for sugar growth's stoichiometry already charged), and `1.0` is simply
where the credit runs out — the validity boundary of a proxy, not a magic number.

**How to apply:** before writing "this forbids X", ask whether X could physically
produce the quantity the guard measures. If not, name the boundary instead ("the
sparing credit's ceiling", "more mass out than in"). A guard whose stated reason is
wrong **will eventually be argued away on the strength of that wrong reason** — and
the person arguing will be right about the science and wrong about the code, which
is the worst shape a review can take. Corollary: this project's oldest lesson
(D-96/D-102/D-108/D-109) is that *code can be right while its stated reason is
wrong* — the label on an invariant is exactly that failure mode. Related:
[[feedback-discuss-disagreements]].
