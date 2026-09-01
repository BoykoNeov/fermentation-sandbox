---
name: feedback-a-teardown-carries-platform-assumptions
description: A test's cleanup carries platform assumptions nobody thinks of as part of the test, and an ignored join is the defect
metadata:
  type: feedback
---

Closing a file descriptor that another thread is blocked in `accept()` on is **undefined by POSIX**.
Windows aborts the pending accept at once; Linux may leave the thread parked in the kernel, still
holding the listening socket open — so the port goes on serving after `close()` and the final probe
is *right* to call it busy. `join(timeout=2)` then gives up silently and nothing checks the result,
so a teardown that did not happen surfaces as a probe that does not work. Measured under CPU load on
WSL: **7 failures in 60 rounds, and the accept thread was alive in exactly those 7**. Rewritten to
ask the thread to leave (socket timeout + `Event`) and assert the join: 0 in 60, both platforms.

**Why:** the docstring already said "the failure this pins is a Windows one" — the author knew the
platform mattered for what the test *pins*. The assumption that actually broke was in the cleanup,
which nobody thinks of as part of the test, so it was never platform-reviewed.

**How to apply:** a green test proves the property holds **on the operating system that ran it**.
Review setup and teardown for platform assumptions with the same suspicion as the assertions. Never
leave a `join`/`wait` whose result is discarded — assert it, or the failure is reported as the thing
you were measuring. To discriminate, lift the teardown into a standalone harness and run it under
load on the target OS; an idle box passes it. See [[feedback-a-stand-in-must-answer-like-the-real-thing]] (D-264→D-265).
