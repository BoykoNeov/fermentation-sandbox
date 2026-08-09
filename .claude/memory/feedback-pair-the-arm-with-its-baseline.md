---
name: feedback-pair-the-arm-with-its-baseline
description: "Scoping a sampled run (only=/exclude=) shifts the RNG draw sequence, so the arm and its baseline are two different random ensembles — pair them on a fixed hypercube instead"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 98da795f-77eb-4137-b991-7ce6520375db
  modified: 2026-08-09T12:26:41.846Z
---

**When an arm changes *which* names are drawn, it also changes *what every other name draws*.** In
D-165 the consequential half compared a default ensemble against `exclude=[one_name]` and against a
log-scale sampler. Every arm looked dramatic — 313 %, 483 %, 364 % — and every number was noise:
removing one name from `sample_parameters`' ordered draw list shifts the RNG stream, so all 66
*remaining* parameters got different values too. The arm and its baseline were not the same experiment
with one thing changed; they were two unrelated random ensembles.

**Why:** this is the [[feedback-verify-the-restore-between-mutation-arms]] failure mode moved into a
stochastic harness. The comparison *runs*, produces plausible-looking percentages, and there is nothing
red to notice. It is worse than a crash because the numbers are usable — I nearly recorded a 15 %
displacement attribution built on it. Two independent tells caught it: the excluded parameter was one
I had just proven *inert* (so a 313 % move was impossible), and the moved outputs were ~1e-10
quantities whose ratios were meaningless anyway.

**How to apply:** pair the arms on a **fixed unit hypercube** and vary only the mapping. In this repo
that is `sampler="lhs"` (or `"sobol"`), which draws the hypercube from the seed and maps each column
through `_inverse_cdf` — same name list + same seed ⇒ byte-identical hypercube, so a triangular arm and
a log-triangular arm differ in nothing but the inverse CDF. For an *inertness* question, don't use
`exclude=` at all: use `only=[the_one_name]` so everything else is pinned, and assert the members come
back bitwise-identical to the nominal. And put a **magnitude floor** on any ratio you report, printing
what it dropped ([[feedback-count-and-print-your-skips]]) — a dry ferment's `S ≈ −5e-11` produced a
`med/nom` of −312 551×.

**The same beat, a second harness bug of the same family:** `Path.write_text` translated the repo's LF
files to CRLF, so a mutation harness's restore wrote correct *content* and still failed a SHA check.
The assertion fired, `git diff --stat` confirmed the content was intact, and the fix was binary I/O —
but a harness that had only compared content would have gone on to commit CRLF into an LF repo.
