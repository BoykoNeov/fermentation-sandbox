---
name: feedback-a-repair-that-leaves-the-label-lying
description: "When two repairs both fix the consequence, the one that leaves a field's NAME meaning something else is not the more faithful option"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: b530ceda-f935-4502-a0e8-cac9593dd384
  modified: 2026-08-27T12:46:02.065Z
---

When a defect has been stated as *"the declaration is wrong"*, a repair that corrects only the
downstream consequence is not a repair. Ask what each candidate leaves the input field **named
after**: a knob called yeast-assimilable nitrogen that in fact carries only the ammonium share is a
false label surviving its own fix, and every later reader inherits it.

**Why:** D-244 had two routes that both made the yield fit see the run's real nitrogen. Summing the
channels was cheaper (nothing re-declared, no trajectories restructured) and was briefly the
recommended one. It fails on this test alone: D-243's finding was *"a scenario says 250 and the run
carries 362.7"* — a **declaration** defect. Summing repairs the consequence and leaves `yan_mgl`
meaning "the ammonium part only", with its own docstring contradicting it. Partition was more
expensive and was correct.

**How to apply:** re-read the defect's own sentence before choosing a repair, and check which
candidate makes that sentence false. Two corroborations to look for: what the surrounding data
implies about intent (here the dose was split by a published *must* spectrum, so it speciates
nitrogen the must already has), and whether the record that found the defect pre-authorised one
route's churn — D-243's docstring told a future beat how to handle partition's reds and said
nothing about summing. Related: [[feedback-closer-to-reality-decides]].
