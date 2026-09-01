---
name: feedback-a-shipped-constant-cannot-name-your-machine
description: "A path that is true on the machine that wrote the code is a break on every other one; derive it at runtime and guard the class, not the constant"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2e4a077a-2634-4e22-b911-4a46b672d8e0
  modified: 2026-09-01T15:01:22.229Z
---

A location written into code that other people run must be **derived at runtime**, never a
literal. `app/main.py` carried `REPORT_DIR` as an absolute path on drive `M:`, and the writer
creates the folder it is handed — so the console's *Write it* button raised on every machine
without that drive, in front of the user least able to read a traceback (D-264, correcting
D-261).

**Why:** it arrived through a rule that was right in its own scope. The working rule "temp and
scratch files go in one fixed folder" is about paths *I* write while working; applied to a
constant inside shipped software it ships my filesystem to the user. Scope the rule to the
files it was written for.

**How to apply:** derive from `Path.home()` (or a documented environment variable) and say the
full path back in the success message — a file saved somewhere the user cannot find is not
saved. Guard the **class**, not the one constant: `test_no_shipped_interface_code_names_an_absolute_path`
scans every file under `app/` for a drive letter or a rooted unix path, comments included, since
a comment naming a false location documents the same falsehood. Sweep with
`grep -nE '[A-Za-z]:[\\/]' app/ src/`; a grep for the one path you remember finds only that one.
Related: [[feedback-verify-latest-state-not-breadcrumbs]], [[feedback-a-doc-rots-where-it-duplicates]].
