---
name: oxidation-direct-and-burst
description: "The direct oxidative set and the burst superset (D-132 to D-137, D-147, D-149 to D-152) - burst is wired and non-default"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — oxidation, direct set and burst.** Detail split out of
`.claude/memory/project-fermentation-sandbox.md` at D-185; that file's ledger points here by
path. Read it when working on this subject. Every bullet is *what it forbids* + the D-record
to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue
past it from this file.** **Never evict an old prohibition to buy a line.**

**Oxidation (D-132 → D-137, D-149 → D-152)**
- D-132's phenolic boost is **additive, never proportional**, browning-side only. D-133's `burst_antioxidant` is
  **EXCESS** over it and must **read none of** `tannin`/`anthocyanin`/`so2_total` — anti-double-count, and binds.
- `k_copper_multiplier` = **600** L/g, **§2.5 CLOSED (D-149)** — never re-open on "the printed table says copper
  is stronger" (Nguyen T3.1 → 2092, failing the rejected 2000's budget). Held by the **source**, not Ferreira's
  ceiling; re-open needs **real wine, ≥3 Cu levels**. Guard 4.
- **BOUNDED, not null (D-152):** copper-orthogonal L16 gives **k ≤ 918 L/g**, excluding 2092, 2000 and band-high
  1500 under all six arms (Guard 7). **Never rebuild on Table 2's SDs; keep condition 12; NEVER cite it as
  evidence FOR 600** (only 1 of 6 arms admits 600) — one-directional, against higher k only.
- **BAND FIXED (D-154) — was the live defect; do not re-propose it.** High **1500 → 662.8** L/g: the bound at
  `copper_typical`'s **MAXIMUM**, not the shipped-centring 918, because `copper_typical` is **itself sampled**
  and `k_bound` *decreases* in it. **Never take 918**: over 200k joint draws it still violates **5.01%** of pairs
  (1500 violates **37.6%**, so D-152's 29% **understates** it). Verify sampling claims **on draws, not edges**.
  **Never adopt D-152's printed "663"** — exact is 662.802522, so 663 ships red. Value/low
  edge/`copper_typical` **untouched**. Guard **recomputes** the bound + asserts monotonicity, never reads the
  note [[feedback-rejected-values-must-be-unreachable]]. **No Fe(III) state** — not D-134's "iron in surplus": QSS rests on a **~18× separation**.
- **A pH term on activation/`k_browning_eff` is REFUSED (D-150)** — never re-open on "Nguyen's table shows strong
  pH dependence". Strongest leg is **inseparability from copper**, corroborated by Carrasco-Quiroz's
  copper-**orthogonal** L16 (D-151) — never dismiss that as Nguyen's dosing. His is an *initial* statistic on a
  *steady* node ⇒ **the pH term's home is the burst**. Guard 5's bound is an **OBSERVED spread, NOT a ceiling**;
  re-open needs a **within-wine** steady-OCR series, ≥3 levels, measuring Cu.
- **`initial_ph` anchors t=0 only** — **no way to SET an aging pH**. First-order in `[o2]`; **no MM/`Km`**. **Do
  NOT re-attempt the Ferreira/Carrascón PLS extraction** — all three blocked (`_findings/D-134-*.md`).
- **The O2 gate is Fe(II)+O2**; SO2 is right in *size*, **INVERTED in role** (enabler, not competitor) — which is
  why **D-72's wrong mechanism yields a right-looking 1:2**; never read that ratio as confirmation. **Never
  re-fix the acetaldehyde/phenolic inversion inside the parallel frame — it provably cancels**; approve the
  cascade on MECHANISM only, **never tune to phenol ⇒ more MeCHO**. Danilewicz 2011 is **PARAPHRASE** — pull the
  abstract before it backs a `source:`; in-wine **mixing-limited**, atmospheric pH-independence does **NOT**
  transfer to wine.

**Burst — WIRED and NON-DEFAULT (D-147). Not dead code, not unbuilt, not refused.**
- `oxidative="direct_burst"` = direct **+** burst, a **superset**, not a mechanism. **No `cascade_burst`**
  (D-138 stays UNDETERMINED; its "transient modifier" is FLAGGED).
- **Never make it default** — moves all 31 D-140 pins and re-opens the Danilewicz direct arm, **moot: do not
  re-derive.** **The split is UNPINNED and stays so** — C1 pins only the **product**, C2 fails **structurally**.
  **Never re-solve it** (D-146's second-unmeasured-number trap). Unlock: **Ferreira 2015 per-cycle O₂ curves**.
- **Self-exhaustion is Ferreira's PROTOCOL, not the sink** (~1400× the cork flux) ⇒ a permanent **~37% tax that
  GROWS to a plateau**; two guards forbid calling it transient — **never relax to a "pool depletes" check**.
  **Seed follows the consumer**: **0.0** outside `direct_burst`; dosing `burst_antioxidant_gpl` with no consumer
  **raises** (D-45's "absent ≠ 0" **inverts**). Ferreira's **2.7× does not discriminate**, and is paraphrase.
