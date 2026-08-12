---
name: fined-sulfur-destination
description: "Where copper fining's bound sulfur goes (D-193) — BUILT as a 1:1 transfer into the D-135 reservoirs, with the retention-fraction variant measured and rejected"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — the destination of copper-fined sulfur.** Split out at D-185's pattern; the
status ledger points here by path. Read it when working on `add_copper`, the bound-sulfide
reservoirs, or the fining carbon ledger. Every bullet is *what it forbids* + the D-record to read
for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue past it from
this file.** **Never evict an old prohibition to buy a line.**

**The transfer is BUILT (D-193). `add_copper` moves what it binds into `bound_h2s` /
`bound_methanethiol`. Do not re-propose it as unbuilt, and never "restore" the annihilating
behaviour.**
- **1:1, WHOLE MASS, no parameter mediates it.** The sulfur is NOT scaled by
  `copper_fining_residual_fraction`. That variant was BUILT as an arm and REJECTED: 0.95 is a
  **printed lower BOUND on copper retention measured after FILTERING or RACKING** — operations
  this verb does not model (`rack` is a separate verb, and the parameter's own band note already
  says those sinks are "not part of this event"). Scaling by it invents a loss fraction and makes
  that invention the only carbon outflow. **Copper 0.95 and sulfur 1.0 differ ON PURPOSE**; a test
  fails if a future edit "harmonises" them.
- **The free pools are byte-for-byte untouched at the event** — the odour fix is unchanged and only
  the destination differs. Never let a routing change leak into the removal arithmetic.
- **D-45's carbon booking is CORRECTED, not shrunk**: the mercaptide stays dispersed, so the fining
  moves **no carbon at all** (−1.19e-23 g/L, was −1.09e-06). The ledger identity closes with a ZERO
  correction term.
- **D-135's two refusals STAND and this does not breach either.** The release **RATE** reads no
  copper (a PLS coefficient is not a stoichiometry) — `reads` unchanged. The **SPONTANEOUS**
  formation route (fermentative H₂S meeting must copper, no intervention) is still unmodelled and
  still needs the unpublished binding constant. A seeded reservoir plus an event that adds to it is
  not a copper-dependent rate law.
- **D-135's docstring parenthetical is CORRECTED**: "fining … releases more thiol from a reservoir
  fining never touched — the real, and correct, behaviour" is FALSE. Franco-Luesma's negative
  copper coefficient for bonded MeSH is about **natural wine copper levels**, not about a
  deliberate dose forming Cu(SR)₂ by the verb's own stoichiometry — the two coexist.
- **One named extrapolation, in three places:** the fined-in mercaptide is released at the
  *natural* pool's 8.1 %/yr, which nothing measures for it. **Weaker** than what it replaced
  (instant, total, permanent destruction), so do not "fix" it by reverting.
- **The 3-year numbers are scenario properties, not the finding.** Reservoir ~doubles (19.70 →
  39.19 µg/L) but free H₂S moves only +1.09 µg/L at 3 y (release is 1.9 %/yr). **The finding is
  that a fining is no longer permanent.** The fined wine still ends up far LESS reductive than an
  unfined control (26.50 vs 44.90) — that ordering is pinned; never read the routing as undoing the
  fining.
- **Beer is inert by having neither slot**, the D-191 guard idiom.

Receipts: `M:\claud_projects\temp\ferment\d193-fined-sulfur-destination\` — `PREREGISTER.md`
(wrong about the design, instructively), `RESULTS.md`, `probe1_routing.py`.
Related: [[prohibitions-residual-copper-fining]],
[[feedback-prefer-the-variant-your-guards-can-see]],
[[feedback-compute-the-clean-fix-before-adopting-it]],
[[feedback-a-scope-note-can-carry-a-mechanism-claim]].
