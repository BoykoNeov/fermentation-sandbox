---
name: acetaldehyde-ladder
description: "Herzan 2020's SO2 ladder and post-ferment acetaldehyde (D-188) - measured, attributed, and NOT built"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — the acetaldehyde ladder.** Split out at D-185's pattern; the status
ledger points here by path. Read it when working on this subject. Every bullet is *what it
forbids* + the D-record to read for *why*. If a prohibition looks unconvincing, **go read its
D-record — do not argue past it from this file.** **Never evict an old prohibition to buy a line.**

**The ladder is MEASURED and the gap is PINNED (D-188). Do not re-propose it as unmeasured.**
- **D-108's "no maturation-phase source to make more" is FALSE** since D-136/D-187: a sealed
  sulfited bottle now makes acetaldehyde. D-108 grew a ⚠ at D-188. The item rode the
  copy-forward open lists **20 times after its blocker shipped** (30 archive mentions total).
- **The `0-vs-2.7` unsulfited floor is CLOSED as a separate item** — it is the closure menu's
  span, pinned two-sided: hermetic **0**, screwcap **0.4668**, technical cork 0.7043, natural
  cork 1.0131, Nomacorc 1.4550, SupremeCorq **2.2591** vs Herzan's 2.7 (0.837×). **Never quote
  one closure's number as "N % closed"** — that is the flat-model lottery, one number standing
  for a set. Nothing was fitted; D-136/D-187 closed it unaimed.
- **NEVER report a per-cell ratio off this ladder.** The model emits **one constant** (26.256 /
  26.507 / 26.604, 1.3 % apart) against a published 17.2–51.6, so *some* ratio is guaranteed to
  read well — the first draft's "1.53× high" was that artifact. The finding is **structural**:
  the model collapses a three-stage design onto the **must stage alone**. Only the ABSOLUTE
  must-column lift on Herzan's two MATCHED pairs is quantitative (+25.95/+26.41 vs +10.7/+25.7).
- **Do NOT group "sulfited must beats unsulfited must"** — false in the published data:
  (0/30/35)=25.9 beats (60/0/35)=17.2. The tank column can out-do the must column; hold the
  other two stages fixed or the claim is wrong.
- **The attribution is SPLIT and a one-sided fix cannot move the benchmark.** `k_so2_oxidation`→0
  **restores the sign** (so **D-47 protection WORKS** — do not indict it), but the lift is
  **+1.0 %** vs the paper's **+298 %**: both `AcetaldehydeProduction` terms (D-27 base, D-48
  SO2-induced) ride `fermentative_flux_shape`, **exactly 0 from dryness**. Competition owns the
  SIGN, a missing post-dryness source owns the SIZE.
- **Do NOT build the post-dryness source.** Its only magnitude anchor is the table under test ⇒ a
  fit wearing a benchmark's name. It also **cannot run as written**: viable `X` is **2.13e-4 g/L**
  by end of tank, so a yeast gate would fire on a wine with no yeast — the gap is the absent
  **lees / young-wine phase**, not one term. Unlock = a source measuring the **rate**, not an
  endpoint. Pinned as **3 strict xfails asserting DIRECTION only** — a magnitude in an xfail is a
  fit target handed to whoever closes it.
- **The tank phase is metered at the BOTTLE's closure rate** (one aging phase, not two) and the
  wine seeds **tannin/anthocyanin/ethyl_bridge = 0** ⇒ D-80 bridging inert. **Never re-derive
  this file's conclusion on a red wine.** "Maintained" = **bisect the shipped binding solver for
  the free-SO2 top-up every 14 days**, never a single dose. Interior timeline split is **swept,
  not asserted** (sulfited rows 26.25–26.99; floor 0.622–1.211).
