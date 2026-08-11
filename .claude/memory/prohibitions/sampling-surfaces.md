---
name: sampling-surfaces
description: "Which bands exist and which constraints are checked at a point (D-153 to D-162) - both archive-wide sweeps are done"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — sampling surfaces.** Detail split out of
`.claude/memory/project-fermentation-sandbox.md` at D-185; that file's ledger points here by
path. Read it when working on this subject. Every bullet is *what it forbids* + the D-record
to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue
past it from this file.** **Never evict an old prohibition to buy a line.**

**Sampled bands (D-153 → D-157) — BOTH archive-wide sweeps are DONE. Do not re-run either.** The **sampler**
surface (which bands exist) = D-153/D-156; the **assertion** surface (constraints checked at a point) = D-157.
- **THE RECURRING SHAPE, 6 instances (D-118, D-154, D-155, D-157, D-180, D-181): a constraint verified at
  a POINT where the sampler reads a BAND.** Whenever a guard or bound uses a nominal, check whether that
  quantity is itself sampled — and take the **joint** worst case over every band involved.
  **D-180 hit it while FIXING it; D-181 hit it while DOCUMENTING D-180's** — it added 3 floor bands to the
  same test and pinned the 2 pKas and 3 seeds it had just shipped. **Enumerate the drawn set from the
  registry, never by hand**, and assert the corner COUNT. It was invisible because it barely moved the
  answer (8.7-81.4 → **7.6-82.2 %**): a shape recurring 6× is how I work, not bad luck. **D-182 did NOT
  recur it** — its 3 new dims went in the SHIPPING commit (9 dims, 19683 corners); that is the fix.
- **FOUR surfaces (D-157), TWO distributions — PINNED (D-156, `tests/test_sampling_surfaces.py`): do not re-audit,
  do not "simplify".** Compile-seam distinct varying **280**, drawn **185** — both a **COMMENT in
  `test_drawability_surface.py`, not asserts ⇒ STALE GREEN: re-measure BOTH, never cite** (its "152 unaffected"
  was false from **D-179**). **NOT 279** (per-*file* sum, double-counts shared names); **structural 61, NOT
  D-157's 66** (D-159, 5 merely scenario-inert). Predicate = **declared `reads`**; **UNREGISTERED-class and
  compile-time doses are NEVER drawn**. **A shared name carries DISJOINT bands per medium** (19 of 33). `psychophysics.yaml` **UNIFORM** —
  never "fix" to triangular, **never apply the triangular mass stat**; `sensory.yaml`'s **36 NEVER sampled**.
  `SHARED_FILES` restated **on purpose** (D-108/D-109 vacuity). **A test at `x == mode` is vacuous** — go off-mode.
- **D-157's live contradiction CLOSED (D-158) — the band WON; never re-narrow 0.084 to 0.08.** Resolved
  **INTERNALLY, no fetch** (84 = Shinohara's 16.4 % E-rate; 30–80 lived ONLY in the test comment asserting it).
  Test **recomputes** all three (`abs=5e-4`; never `rel`/`round(x,3)` — pins *formatting*). Band = E-rate spread at **FIXED acetic 0.35** — the acetic SPREAD is `acetic_acid_typical`'s, D-176.
- **`reads` has TWO masters — tier propagation AND sampler scope (D-160, fixes D-159's defect).**
  `PH_SYSTEM_READS`/`SO2_BINDING_READS` (`acidbase.py`) → 19 members/5 modules. **Keep DISJOINT** (pinned) and
  **derived**, never re-listed. **`temperature_ramp_rate` stays undeclared BY DESIGN**; no tier moved.
  **RESTATEMENT DONE (D-161) — never re-run it.** Affected class = **ONE** row, a 9-name `only=` **ISOLATION** not a
  band; its "13.99/16.69/20.09 across seeds 0/1/2" is **SORTED** [[feedback-a-majority-is-not-a-direction]]. Through the shipped sampler the fix
  is **undetectable**; predicted **widening NOT observed** (10/24, p=0.54). D-159 pins consumption (`test_drawability_surface.py`).
- **Closure ordering SCOPED (D-162) — do not "fix", narrow or re-measure it.** **Three** claims: P1's **three-tier**
  grouping is its conclusion verbatim; `technical<screwcap` and `nomacorc<supremecorq` are Table-I **nominals**, not
  that sentence, and one pair carries 94 % of the chain breach. Declaration-level ONLY: one `scenario.closure` →
  the **`closure_otr` STATE slot**. **12** assertions undecided [[feedback-count-and-print-your-skips]] — **list never persisted**.
