---
name: feedback-an-out-of-band-mutation-breaks-loading-not-the-quantity
description: "A mutation that puts a parameter outside its own declared band is refused by the loader, so every test dies at import with nothing integrated — a whole-module RED that scores no arm and reads like maximal sensitivity"
metadata:
  node_type: memory
  type: feedback
---

[[feedback-verify-the-mutation-applied-not-just-the-restore]] covers a mutation that silently
**fails to apply**. This is the opposite failure: it applies perfectly and then never **runs**.

**Why:** a D-243 arm halved a must-spectrum share from 0.430 to 0.215 to shrink the nitrogen a pool
carries. That value is outside the parameter's own `[0.34, 0.55]` uncertainty band, so
`load_parameters` refused it and **all five tests in the module failed in 1.03 s** — no compile, no
integration, no quantity moved. Taken at face value the arm looks like a spectacular hit: the widest
possible red set. It scores nothing. Re-run at the band's own low edge (0.34, an in-band move) it
produced the intended single red.

The tell is the **duration and the width together**: a whole module dying in about a second is a
load or import failure, not a set of assertions disagreeing. A per-arm red *set* that suddenly
includes tests with no logical relation to the mutated quantity is the same signal
[[feedback-a-red-set-firing-wider-is-a-finding]] — but here the finding is about the harness, not
the code.

**How to apply:** when a repo validates parameters on load (bands, schema, provenance), a mutation
arm has an admissible range and stepping outside it converts a measurement into a smoke test.
Choose arm values **inside the declared band** — its edges are the strongest legal move and are
usually enough. If an out-of-band value is genuinely what the arm needs, widen the band in the same
edit and say so. And make the harness itself notice: assert the arm produced at least one *passing*
test as well as the predicted failures, or check that the run took long enough to have integrated
anything. An arm whose denominator is zero should never be reported as a hit.
