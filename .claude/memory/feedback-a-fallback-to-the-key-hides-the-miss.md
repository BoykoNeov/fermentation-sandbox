---
name: feedback-a-fallback-to-the-key-hides-the-miss
description: "A lookup that returns its own key when it misses makes the miss indistinguishable from a hit downstream; return a sentinel and COUNT the misses in the output"
metadata:
  node_type: memory
  type: feedback
---

**`consts.get(name, name)` is not a default — it is a silent failure that types correctly.**
When a resolver falls back to the thing it was asked to resolve, the miss comes back the same
SHAPE as a hit and nothing downstream can tell them apart.

D-275: `tools/gen_open_ledger.py` resolved an xfail's `reason=NAME` from module-level string
assignments **in the file carrying the marker**, then ended `return consts.get(node.id,
node.id)`. One of the seven rows in `docs/OPEN.md` — the derived answer to "what is
scientifically open" — takes its reason from a constant assigned in a sibling test module and
imported at the top of the marked file. That row printed
`_D245_D120_LEGS_GONE_DIRECTION_BACK_AT_D248` as its reason and `—` as its records column,
i.e. it asserted in the file's own idiom that **no decision record bears on it**. Six do
(D-245, D-120, D-206, D-246, D-248, D-254). The file's COUNTS were right the whole time —
seven open markers, correct — so every check on it stayed green: **an unreadable row is not a
missing row, and nothing was counting readability.**

**Why:** the fallback was written for the case where a reason is a literal and the name lookup
is a nicety. But a string-shaped miss survives every type check, every count, and every
"is the generated file current" test, because the file IS current — it faithfully renders the
wrong thing. The failure is also invisible from either end: the module that owns the text has
no marker, and the module with the marker has no text, so neither file looks wrong on its own.

**How to apply:** when a resolver can miss, make the miss a different KIND of value (`None`, a
sentinel) at the resolution boundary, and if the caller must still render something, **count
the misses and print the count in the artefact** — the idiom this repo already uses for
"marker-shaped lines the parser cannot read". A count that lands at 0 on the live tree binds
on the next occurrence and mutes nobody. Then measure the tripwire: the guard here restores
the single-file resolution and requires it to reproduce the identifier, because "0 unresolved"
passes just as well when the check is measuring nothing. And when a resolver reads names across
files, resolve to a **fixed point** — re-exports put a constant two hops from its use, and a
local assignment must shadow an import of the same name.

Related: [[feedback-a-guard-is-only-as-good-as-what-its-red-names]] (the mutation arm),
[[feedback-measure-the-surface-in-the-unit-that-fails-it]] (a monitor in the wrong unit reports
clean), [[feedback-name-guards-for-what-they-forbid]] (a check that fires constantly gets muted).
