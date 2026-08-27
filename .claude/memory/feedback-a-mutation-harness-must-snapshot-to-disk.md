---
name: feedback-a-mutation-harness-must-snapshot-to-disk
description: "A mutation arm that restores the file from an in-memory copy leaves the mutation live on disk if the run is killed — snapshot to a file and restore from that"
metadata:
  node_type: memory
  type: feedback
---

A mutation-arm harness that reads the target file into a variable, writes a mutated version, and
restores from that variable in a `finally` is **not crash-safe**. `finally` does not run when the
process is killed, so an interrupted arm leaves the mutation on disk — in a source tree you are
still editing.

**Why:** at D-240 I killed an arms run mid-arm (it was mutating the test module while I was
editing the same file). The `finally` never fired, and the arm's mutation — one deleted line from
a registry tuple — stayed live. I then edited the file **on top of** the mutated version, backed
up that version as my "known good" copy, and ran the suite twice against it. The two REDs it
produced looked like my own new guard failing and cost a diagnosis pass; the giveaway was that
both failures were on wine scenarios and my edit had touched neither.

**How to apply:** write the pristine text to a **snapshot file** before mutating, restore from the
snapshot, and verify by SHA-256 that the restored bytes match — three lines that make an
interrupted run recoverable by hand. Two habits go with it: never edit the file an arms run is
mutating (start the arms only when the module is final), and when a test goes RED just after a
killed run, `git diff` / grep the target file **before** debugging the test — a mutation still on
disk is indistinguishable from a real defect if you only read the failure message.
