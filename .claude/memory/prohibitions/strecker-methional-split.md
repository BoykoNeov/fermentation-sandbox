---
name: strecker-methional-split
description: "Why f_methional cannot be derived from the amino-acid abundances (D-206), and the CompiledScenario reuse contract"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — the Strecker methional/phenylacetaldehyde split.** Detail reached BY PATH
from the ledger in `.claude/memory/project-fermentation-sandbox.md`. Every bullet is *what it
forbids* + the D-record to read for *why*. If a prohibition looks unconvincing, **go read its
D-record — do not argue past it from this file.** **Never evict an old prohibition to buy a line.**

**The split's derivation is REFUSED (D-101 → D-206) — do not re-open it, and do not re-band it.**
- **`f_methional` is the ONLY channel** from the model's methionine story to its methional output.
  The abundances **cannot reach it**: pool `dose·f_i/Σf` over gate `aa_i/(K·f_i + aa_i)` makes the
  fraction **CANCEL** — D-100's design, general to **all eight** pools (every gate **0.888889** at
  nominal, registry-enumerated). A **4× recompiled** change in the must's methionine moves aged
  methional **0.06 %**. Deriving would **INVENT** a sensitivity the gate removes by design; that is
  the ground, "reactivity = 1 is unsourced" the LESSER reason. **Re-banding to the abundance-implied
  range is refused on the SAME ground** — it imports a constraint from the silenced channel.
- **The D-101 note's "≈ 0.136" was WRONG when written** — it pairs Cabernet's methionine with the
  two-source phenylalanine midpoint. Consistent from the shipped nominals: **0.15152 vs the
  shipped 0.15, 1.01 % apart**, so its "close but NOT equal" reasoning had no gap to reason from.
  Its **"deriving it is queued as its own beat" clause is DELETED, not re-worded** (D-107's shape).
- **Priced, not shrugged:** the implied methionine:phenylalanine reactivity ratio over the three
  bands has median **1.2397**, 5–95 % **[0.5224, 2.7443]**, **21.86 %** of draws outside 2× either
  way. That is **this file's internal coherence, NOT anything the model reports** — never lead with
  it. In a sealed 5-year bottle the independent-vs-coupled ensembles differ by a **17 % median
  shift** (the band's own asymmetry, mode 0.15 → median 0.177) and a **1.59× width** change, but
  the verdict does **not** move: **89.5 % vs 90.2 %** of bottles over threshold.
- **No guard was owed and that was measured:** the band-floor mutation is **RED** (15 failures,
  incl. `test_old_oxidative_set_reproduces_its_trajectory[1y-methional]`). The property measured —
  the control arm's flatness — is what a correct future build must BREAK, so pinning it would
  fight the fix (D-205 §7's reasoning).

**`simulate_scheduled` LEAVES the Process set reconfigured — a CONTRACT, not a bug (D-206 §7).**
- A second run on one `CompiledScenario` starts with the first's 22 aging Processes live **from
  t = 0**: **+10.3 %** on aged methional, active count unchanged at 49, nothing raised. **Bracket
  reused sets with `enabled_snapshot()`/`restore_enabled()`**; `simulate_ensemble` already does,
  per member. **Do NOT "fix" it to restore** — measured, that fails **26 tests**, because the
  O₂-partition guards read the configuration in force *at the end* of a run and get an empty dict.
  It is now stated in `simulate_scheduled`'s and `CompiledScenario.run`'s docstrings; before D-206
  it lived only in a comment inside `test_ensemble.py`.
