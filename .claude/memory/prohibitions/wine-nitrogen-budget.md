---
name: wine-nitrogen-budget
description: "D-243 — wine's nitrogen budget is AUDITED: numerator sourced, identity exact, ledger closed, the beer inadmissibility argument does NOT transfer, and the two-channel evaluation point is open on the owner's call"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — wine's nitrogen budget (D-243).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about wine YAN, `yan_mgl`, `amino_acids_gpl`, `biomass_N_fraction`,
Coleman's `Y_X/N` regression, the Varela comparison, or "the nitrogen budget is wrong". Every
bullet is *what it forbids* + the record to read for *why*. **If a prohibition looks unconvincing,
go read D-243 — do not argue past it from here.** Siblings: `seed-reads-repair.md` (D-241, whose
SUBSUMED verdict this corrects the scope of), `banded-undrawn-census.md` (D-240).

- **The D-232 item is CLOSED — never re-propose "wine's nitrogen budget has not been audited".**
  It rode the open list from D-232 through D-241 and is now done, both halves: the numerator is
  sourced (§2) and the identity is inverted against a measured crop (§3), which is exactly what
  D-230 did for beer. What remains open is ONE named design decision (below), not the audit.

- **The numerator is SOUND — do not re-check Varela's nitrogen for a proline subtraction.**
  Varela print the subtraction themselves: 380 / 65 mg/L total, 300 / 50 **assimilable**, proline
  (20.5 % wt/wt of their nitrogen mixture) excluded because it is not assimilable anaerobically.
  Quoted verbatim in D-243 §2 and in the test module docstring. The timepoint is also checked and
  benign (stationary-phase dry weights at 48 h / 72 h; the engine's biomass is flat from ~40 h),
  and they report assimilable nitrogen **completely depleted** — the identity's own assumption.

- **NEVER write the beer symmetry.** The tempting record is "wine's demand is outside the
  physiological range too, in the other direction". It is FALSE and was corrected before shipping.
  Beer's demand (20.2–26.2 % N) was above a hard ceiling; wine has **no comparable floor** —
  `wine_generic.yaml` calls 0.114 the *N-replete* reference that "drops under nitrogen limitation",
  and the shipped regression already runs f_N = 0.0362 at YAN=50. The demanded 0.0521 at YAN=300
  sits **between the engine's own two values**, so a nitrogen-budget explanation for wine is
  **ADMISSIBLE**. Do not cite D-230's inadmissibility argument on the wine side.

- **This is NOT a conservation defect — do not open it as one.** The ledger closes to 3.1e-14
  relative across all four sink arms, including autolysis and oxygen. The growth identity
  `X0 + YAN/f_N` reproduces the simulated crop to ratio **1.00000000** at both Varela levels.

- **The slope disagreement is a RE-EXPRESSION, not a finding.** Engine yield rises 2.466x from
  N=300 to N=50, Varela's 1.521x. D-56 finding 3 already prints both endpoint magnitudes
  (11.2 vs 19.3; 27.7 vs 30). Checked before writing, precisely because it is the shape a new
  headline would take. Pursuing the mechanism needs a third dataset that does not exist.

- **OPEN, and the OWNER'S CALL — the evaluation point.** `yan_mgl` seeds the `N` slot **and** is
  where Coleman's regression is evaluated (D-14); `amino_acids_gpl` seeds eight pools also on the
  nitrogen ledger (D-100). They **ADD**; they do not partition. At the suite's commonest 0.5 g/L
  dose a wine declaring 250 mg N/L carries **362.7**, and the fit stays at 250 — more nitrogen and
  a yield for a poorer must, compounding. D-32's own text says *"amino acids are part of YAN"*,
  which the seam does not implement; D-243 `Flags:` D-32 and D-14 rather than repairing.

- **Do NOT take the obvious repair without the owner.** Summing the channels into the fit leaves
  Coleman's fitted 70–350 mg N/L range at that same 0.5 g/L dose and reaches f_N = 0.379 at 2 g/L
  — cells 38 % nitrogen, three times anything physiological. It trades a wrong evaluation point
  for an extrapolated one. The other route, making the channels partition, changes what every
  dosed wine scenario in the suite means. Both are priced in D-243 §6; neither is free.

- **Do NOT tune anything against Varela, ever.** D-56 set the firewall: it is the project's only
  independent wine dataset and is a validation set only while it is never a calibration set. D-243
  fits nothing to it. Do not widen the Varela benchmark's characterized bands to absorb a change.

- **D-241 §2/§3's SUBSUMED verdict STANDS — its SCOPE was wrong.** Containment of the two Coleman
  coefficients' implied range inside the sampled bracket holds only over **YAN 66.0–324.8 mg N/L**,
  and the "2.11x wider" figure is 250 mg/L's value (7.16x at YAN=50, 1.22x at 350). The escapes
  are small (4.6 % of the low tail at 50; 1.1159x over the high edge at 350) and the verdict is
  right for the battery and every ensemble scenario here. **Do not restore a red by widening the
  bracket literals** — that bracket is what the ensemble actually draws for the constant that
  governs biomass. Guarded in both directions by
  `test_the_subsumed_verdict_is_scoped_to_the_yan_interval_that_contains_it`.

- **Do not re-run the sink census expecting the pre-registered answer.** A pre-registered "no
  non-biomass sink over 5 %" was wrong by five times: with autolysis on, **45.9 %** of the nitrogen
  ends in the amino-acid pools, because lees self-digestion returns it faster than anything
  re-assimilates it post-dryness. Pinned with its autolysis-off control.

- **Untouched and named, not oversights:** the tier half of the two channels (pool provenance
  tiers do not propagate into the biomass they become — same shape as D-241 §10's seed-tier gap),
  and the mechanism behind the flatter measured nitrogen-dependence, which is D-232 §5's
  "different organism / medium: untouched, unsourced" residue.

Measurements: `M:\claud_projects\temp\ferment\d243-wine-nitrogen-audit\` — `PREREGISTER.md`,
`FINDINGS.md`, `probe1_identity.py` … `probe4_two_channels.py`, `mutate.py`, `mutate2.py`.
