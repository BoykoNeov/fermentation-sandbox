---
name: feedback-a-stand-in-must-answer-like-the-real-thing
description: "A fixture standing in for a real service must reproduce the behaviour under test, or its RED is about the fixture; and two probes sharing an assumption are not corroboration"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 2e4a077a-2634-4e22-b911-4a46b672d8e0
  modified: 2026-09-01T15:01:34.534Z
---

Two halves of the same beat (D-264), both about the *checking* rather than the code.

**A stand-in must behave like the thing it stands in for.** The guard for "never hand the
server a busy port" used a socket that listened with a backlog of one and never accepted. The
first probe filled the backlog; the second was refused, so the same genuinely-busy port
answered *busy* then *free* — a RED that was entirely about the fixture. A real server accepts,
so the fixture must accept, in a thread.

**Two probes agreeing is not corroboration when they share the wrong assumption.** "Bind the
port and see if it raises" reported a port free while a console was demonstrably serving on it;
adding `SO_REUSEADDR` reported free again. On Windows the option belongs to the *existing*
listener (Uvicorn sets it), so a successful bind proves nothing about whether anyone is there —
only a refused *connection* does. Both probes were wrong in the same direction for one reason.

**Why:** it was caught only because a real server happened to be running on a known port while
the probe was written — i.e. the probe was run against a **known** answer rather than a hoped-for
one. Without that, both versions would have shipped green.

**How to apply:** before trusting any probe, run it against a case whose answer you already
know, in *both* directions. When two implementations agree, ask what assumption they share
before counting it as evidence. And prefer the question the caller actually asks — "is anyone
answering?" — over the one that is easiest to write. Related:
[[feedback-a-null-result-needs-a-positive-control]], [[feedback-mutate-the-premise-before-building-the-guard]],
[[feedback-a-shared-fixture-has-two-consumers]].
