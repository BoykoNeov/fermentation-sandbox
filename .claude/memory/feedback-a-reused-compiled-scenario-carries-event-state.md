---
name: feedback-a-reused-compiled-scenario-carries-event-state
description: "Integrating twice from one CompiledScenario contaminates the second run: scheduled events mutate the ProcessSet"
metadata:
  node_type: memory
  type: feedback
---

Probing D-248's effect on the aging pins, I compiled the scenario once and integrated it twice
with different parameters — the obvious way to isolate one knob. The second run reported A420
**+54 %**. Recompiling per run, the real move was **+1.6e-4**. The whole finding was an artefact.

**Why:** `simulate_scheduled` runs the scenario's events, and interventions like `begin_aging`
enable and disable Processes on the `ProcessSet`. That mutation persists on the compiled object,
so run 2 starts from run 1's *end* wiring. `ProcessSet` has `enabled_snapshot`/`restore_enabled`
precisely because this is a real hazard; a probe that does not use them inherits it silently, and
silently is the problem — the contaminated number looked like a dramatic finding worth chasing.

**How to apply:** in any probe that integrates more than once, **compile fresh inside the loop**.
It costs milliseconds. If a probe's result is dramatic and the shipped tests only moved by a
whisker, suspect the harness before the model — the tests recompile per run and the probe may not.
Related: [[feedback-compute-the-clean-fix-before-adopting-it]].
