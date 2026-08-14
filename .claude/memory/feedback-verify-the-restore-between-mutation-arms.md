---
name: feedback-verify-the-restore-between-mutation-arms
description: "In a mutation matrix, prove the baseline came back between arms — a silent restore failure makes the arms you expect to be RED confirm themselves"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1289a7da-873a-4fc1-882a-f8c7f961f6e7
  modified: 2026-08-12T12:12:43.181Z
---

**Between mutation arms, verify the file actually got restored — and include at least one arm whose
expected outcome is the opposite of the others.** In D-158 the backup used
`cp $F /tmp/f.bak 2>/dev/null || cp $F <durable>`; git-bash has its own `/tmp`, so the *primary* `cp`
succeeded, the fallback never ran, the later restore `cp` failed on a nonexistent path, and arms B and
C both ran against a file still carrying arm A's mutation. All three reported RED — and RED was the
expected answer for two of them.

**Why:** a mutation matrix is graded against expectations, so any arm whose expectation is "red"
cannot distinguish a working guard from a broken harness. The corruption was caught only because one
arm's expected answer was GREEN (a consistent re-sourcing must pass) and it came back red. Without
that arm the matrix would have shipped looking complete. Same family as
[[feedback-count-and-print-your-skips]] — a harness that fails quietly returns the answer you already
believed.

**How to apply:** back up to an explicit, verified path under `M:\claud_projects\temp` (never a bare
`/tmp`, and never a `||` fallback whose primary can succeed somewhere unintended); `ls` the backup
before trusting it. Re-assert the baseline between arms — `git status --porcelain` clean, or grep the
value back — and run the untouched baseline as its own arm at both ends. Design at least one arm to be
green: for a guard that pins a derivation, the green arm is a *consistent* re-sourcing (move the source
number and the shipped number together), which is also the arm that proves the guard pins the
derivation rather than the literals. See [[feedback-mutate-the-premise-before-building-the-guard]].

**Generalises past mutation matrices to any scored harness (D-164).** The first run of a
*pre-registered repro* reported 9 of 9 arms raising, with two predictions marked AS PREDICTED — and
every arm had raised for one unrelated reason (the scenario carried no `temperature_schedule`, so
nothing ever reached the code under test). Two predictions confirmed themselves against a harness that
never ran the experiment. The three designed-GREEN control arms are the only thing that showed it.
So: **any arm whose expected verdict matches the harness's failure mode is unscored until a control
proves the harness works** — for exception-counting harnesses the failure mode is "raises", so the
control must be an arm that must NOT raise. Write the controls into the harness from the first run,
not after a suspicious result. See [[feedback-compute-the-clean-fix-before-adopting-it]].

**An arm's COLOUR is not coverage — capture which test it killed (D-165).** A harness that reads only
pytest's `-q` summary knows an arm went red, never *what* went red, so every arm→test attribution is
inferred. D-165's A2 anchored on `if distribution == "triangular":`, which occurs **twice** in
`ensemble.py` (8-space indent in `sample_parameters`, 4-space in `_inverse_cdf`); a first-occurrence
replace hit the wrong one, so the arm meant to cover the load-bearing *median* test covered a different
test, that test had **no arm at all**, and the record shipped claiming "5 of 5 arms as predicted". Red
was the expected answer, so nothing looked wrong. **How to apply:** run arms with `-rf --tb=no` and
**without `-x`** (which stops the file early and hides second kills — A1 was really killing two tests),
parse the `FAILED …::test_name` lines per arm, then print the set of tests killed by *at least one*
arm against the file's full test list and **name the uncovered ones in the record**. Assert
uniqueness of the anchor before mutating. A "6 of 6 arms" claim and a "6 of 8 tests covered" claim are
different claims; only the second is coverage.

**Verify the APPLY, not just the RESTORE — a schema rejection wears the same colour as a finding
(D-171).** Round 1 of D-171 moved each nominal across its ordering. For three arms the crossing point
lay *outside that parameter's own uncertainty band*, and `Parameter` enforces `low <= value <= high`,
so the store refused to load and pytest died at **collection**. All three reported RED. One of them
was the **known-RED control**: it "matched its prediction" having never executed the assertion it
existed to exercise, certifying a harness that was never tested. The restore was fine; the *apply* was
not, and it still produced the predicted colour. **How to apply:** pre-flight every mutated tree
through the real loader before running the suite, print the values it read back, and report a load
failure as **INVALID — never RED**. Where a single-parameter crossing is schema-illegal, move *both*
members toward each other onto their own band edges. See
[[feedback-pair-the-red-with-an-ordering-preserving-baseline]].

**The code under test can UNDO your mutation mid-run (D-197).** An arm that switched five O₂ sinks
off with `ProcessSet.disable()` before integrating came back **bitwise identical** to the shipped
arm — because `begin_aging` is an *event* that re-enables the aging set partway through the run, so
the mutation was reverted by the very thing being measured. It read as "removing 62 % of the oxygen
draw changes nothing", which is a finding-shaped result from a harness that never ran the
experiment. **How to apply:** mutate at a level the run cannot rewrite (patch the class's
`derivatives`, not the set's enabled flags), and keep a `noop` arm *designed* to be identical so
"identical to shipped" is a verified restore rather than an ambiguous null. Ask what re-configures
state mid-run — events, schedules, reconfigure hooks — before choosing where to cut.

**A green control certifies only its OWN mutation CHANNEL (D-201).** D-201 ran a proper null
control — the injected Process present at zero rate — and it came back bitwise identical to
baseline, exactly as designed. It certified nothing about the *other* arms, which mutated a
**parameter** instead of injecting a Process, and which silently did not mutate at all:
`CompiledScenario.param_values` is a **property** returning a fresh `parameters.resolve()` dict, so
`compiled.param_values["k"] = v` writes into a throwaway that is discarded on the next access. Those
arms returned removal fractions identical to the reference *to every printed digit* and read as a
flawless confirmation of the prediction under test. The tell was that they were **too** identical —
a real pool change perturbs the solver's step selection even when the predicted ratio is invariant.
**How to apply:** a control licenses one channel only, so if arms mutate through two channels
(inject a Process *and* override a parameter) each needs its own landing check. Make the check an
assertion on an *observable that must move* — "the pool moved 0.5x, LANDED" — and `raise SystemExit`
before reporting anything if it did not. Then D-197's trap fired again in the very next probe
(`pset.disable()` re-enabled by `begin_aging`), and the landing check is what caught that one too:
one habit covers a family of failures that each look like a clean result.

**A REUSED object is a third channel: run 2 is not a repeat of run 1 (D-206).** Two probes compiled
one scenario and integrated it many times. `begin_aging`'s reconfigure **enables 22 aging Processes
and nothing puts them back**, so every run after the first began with aging chemistry live from
t = 0 — **+10.3 %** on the measured output, the active-Process **count unchanged at 49**, and no
error. It is the D-197 trap inverted: there the run *undid* the mutation, here the run *keeps* a
change the next arm inherits. The arms still differed plausibly, so nothing looked wrong.
**How to apply:** when arms share any object — a compiled scenario, a Process set, a fixture — the
cheapest possible control is **two identical runs asserted equal bitwise**, placed before the
matrix and designed to be GREEN. It costs one extra run and it is the only thing that separates
"the arms differ because of my mutation" from "the arms differ because of each other". And before
"fixing" the leak, mutate it: making the restore unconditional here failed **26 tests**, because
guards deliberately read the configuration in force *at the end* of a run — the persistence was a
contract, undocumented at the call site. See [[feedback-mutate-the-premise-before-building-the-guard]].
