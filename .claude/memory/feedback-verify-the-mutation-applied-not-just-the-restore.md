---
name: feedback-verify-the-mutation-applied-not-just-the-restore
description: "A harness that edits by str.replace on a parsed float silently no-ops — 6.93e-5 round-trips as 6.93e-05 — so the arm loads its ORIGINAL value and reports success; assert the mutation took, not only that the restore did"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 9d970785-c90f-4a56-ae65-f4ba3014180b
  modified: 2026-08-25T18:50:52.252Z
---

[[feedback-verify-the-restore-between-mutation-arms]] covers the restore side. This is the **apply**
side, and it fails more quietly.

D-226's calibration harness widened a parameter's uncertainty band so an out-of-band arm would load,
using `seg.replace(f"low: {lo}", f"low: {new}")` where `lo` came from `float(...)` on the YAML text.
The file says `low: 6.93e-5`; Python renders that float as `6.93e-05`. **The replace matched
nothing, returned the string unchanged, raised nothing, and the arm loaded with its ORIGINAL band.**
It only surfaced because the schema then rejected the value — had the value happened to be in range,
the arm would have run to completion and reported a calibration that never happened.

**Why:** `str.replace` returning its input unchanged is indistinguishable from success. Any
round-trip through `float`/`repr` can reformat a literal — exponent zero-padding, trailing zeros,
`1e-3` vs `0.001` — and YAML, JSON and TOML all permit forms Python does not re-emit. The result is
an arm that is *green because it never ran the mutation*, which is the same failure shape as a
control that confirms itself.

**How to apply:**
- **Edit by SLICE INDEX, not by matching a re-rendered value.** Locate the token with
  `s.index("low:")` and the following comma, and rebuild the string around those offsets. The text
  you found is the text you replace.
- **Assert the edit changed the file**: `assert new_text != old_text`, per field, not per file.
  A per-file check passes when three of four fields applied.
- **Read the value back through the real loader** and assert it is the arm's value before running
  anything expensive. One line, and it converts a silent wrong-arm into a loud one.
- The same applies to any structured-text mutation: commit-message writers, YAML/JSON patchers,
  parameter overrides. If a harness reports "5 of 5 arms clean", ask what its denominator was
  [[feedback-count-and-print-your-skips]] and whether each arm's mutation was ever *observed* to
  land.
