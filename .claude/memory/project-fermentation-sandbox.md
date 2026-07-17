---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation
engine in Python (uv, scipy/numpy/pydantic). Repo:
https://github.com/BoykoNeov/fermentation-sandbox (default branch `main`).

**This file is session-boot context, not a changelog.** The full per-decision
narrative lives in the repo and is the single source of truth — do NOT copy it
here (see [[feedback-batch-end-ritual]]):

- `docs/DECISIONS.md` — every engineering decision D-1 … D-111, with the fork,
  the choice, and the reasoning. **The canonical archive.**
- `docs/ARCHITECTURE.md` — layering, package map, the core/runtime/scenario seams.
- `docs/plans/milestone-*.md` — task checklists per milestone.
- `CLAUDE.md` — prime directives (tiers, provenance, one-directional deps).

**Status (2026-07-17 — at D-111):**
- **Milestones 0, 1, 2 COMPLETE.** M1 = single-strain isothermal nitrogen-limited
  primary fermentation passing both §2.2 benchmarks (wine ~24°Brix→dry 8–14d,
  beer ~1.048→~1.010 in ~5.5d). M2 = all physics beats + post-M2 refinement.
- **Milestone 3 in progress** — the sensory/OAV + aging frontier (owner's pick at
  D-66; Tier-3 aging). Currently at **D-111**: 1171 tests passing, 16/16
  benchmarks in-run, ruff + mypy strict clean.
- **No pool in the project is lumped** (closed at D-110); the `lumped` flag and
  its machinery are kept **dormant**, not deleted — removing them would delete the
  D-66 contract along with its last instance.

**Next / open candidates** (live list; full reasoning in the latest D-record):
- **The leucine shortfall — the sharpest open item.** Model 1.12% vs Rollero's
  leucine tracer 3.4–17.3%. It is D-103's gate **shape**, so a scalar cannot fix
  it. **The lever is the 122.7% over-claim** (leucine's gate claims 90.9% of
  isoamyl while the KIC branch wants 31.8%) — the first arithmetically impossible
  statement the gate has produced.
- **D-104's inverted split** (the prize). D-111 built the KIC mechanism D-104 said
  the model lacked, **valine side only** — still untested against the inversion.
- **D-109's parsimony question** — per-species vs shared BAT1/BAT2; the prototype
  is kinetically-limited, not near-equilibrium.
- **Remaining keto-acid node routes** — KMV → isoleucine; phenylpyruvate.
- **Rollero's isoamyl *acetate* enrichment** (0–19.7% valine-labelled) — an unused
  independent check on the D-97 ATF1 coupling.
- **Closure O₂ ingress** — still load-bearing on the sotolon headline.
- **Acetaldehyde in maturation** + the 0-vs-2.7 unsulfited floor.
- **The deferred tail** (D-110 narrowed this bucket to one item and shipped it;
  the narrowing is **still unconfirmed by the owner** and the rest are NOT
  started): `oav`→`magnitude`; Pham's pH + ethanol terms; sotolon enantiomers;
  methionine's assimilation/sink + `methionol`; growth-linked excretion shape
  (D-49 option B); a peptide pool; variety-specific DMSp; sourced yeast-autolysate
  spectrum; re-anchor `f_methional`; masking (blocked on cosα); D-55's stale
  Brett-phenol prose in `chemistry.py`.

**Standing rule that outlived its bullet (D-66):** direction is the owner's call,
every time — ask before picking the next milestone (see [[feedback-discuss-disagreements]]).
