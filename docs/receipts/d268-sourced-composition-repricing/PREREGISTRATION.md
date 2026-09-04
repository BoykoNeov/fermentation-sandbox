# D-268 pre-registration — written BEFORE the first number

The beat: take D-267's flagged repair. `w_i`, the per-precursor yeast-protein composition every
growth-anchored measurement in D-259, D-260 and D-266 rests on, stops being D-259's stated
bracket and becomes the sourced composition (Lange & Heijnen 2001 Table IV) converted into the
free-acid frame the draw actually uses. Nothing in `src/` changes; what changes is what three
records report.

**Design decided before measuring.** The composition is *derived from the transcription already
in the suite* (`_lange_heijnen_shares("free")`), not re-typed as a second literal — so the only
hand-entered numbers stay Table IV's own 19 rows, which already carry a closure guard. The
residue half stops being a bracket, because it is now a measured point. What survives as a
bracket is the protein fraction of dry weight, `{0.40, 0.45, 0.50}`, which D-267 §3 corroborated
at all three edges. So `w_i(edge) = protein_fraction[edge] * free_acid_share[species] / 100`.

Per-precursor, the multiplier on `w_i` at each edge is the sourced share over the bracket's own
edge: **lo 1.607x (Leu), mid 1.285x, hi 1.071x** and the analogous ratios for the other four
(D-267 §5). The bracket therefore both **rises** and **narrows** — its span drops from 1.875x
(0.40x6.0 -> 0.50x9.0 for leucine) to 1.25x, the protein fraction's own span.

## Predictions, and what each one would mean if wrong

**P1 — direction.** Every growth-anchored split rises at every edge, on every fixture. A larger
`w_i` draws more precursor into the lump; there is no route by which it could fall. If any split
falls, the counterfactual Process is not being fed the weights this beat thinks it is.

**P2 — the bracket narrows.** The lo->hi span of every reported split shrinks. D-259's
"leucine spans 13.1-22.0 %" and D-260's corrected "21.3-33.7 %" both narrow and both shift up.

**P3 — D-104's 20.9 % moves further away, in the direction D-260 already put it.** D-259 read it
at the TOP of the uncorrected bracket; D-260's modifier correction put it just BELOW the
corrected bracket's low edge. Sourcing raises every edge, so it sits further below. It does not
return to the top of anything.

**P4 — the refusal survives.** The corrected bracket's high edge stays below Crepin's 77 %
leucine protein share. Predicted because the largest multiplier is 1.607x at the lo edge and
1.071x at hi, against a gap of 33.7 -> 77 (2.3x). If this fails, D-259's and D-260's refusal of
the growth-anchored sink does NOT survive its own composition being sourced, which is a much
larger result than this beat expects and must be recorded as one rather than absorbed.

**P5 — phenylalanine stays landed.** It reads 96-99.5 % against its sourced 97.5 % fate because
it is supply-limited, not weight-limited; raising `w_i` cannot push it past 100. The leucine
control (leucine must NOT also approach ~98 %) is the thing at risk here: if leucine rises far
enough that both are near 100, D-259's positive finding becomes an artefact of supply limitation
and must be withdrawn rather than re-pinned.

**P6 — the two D-260 dials part company.**
* `_D260_MATCHED_F` (0.174) is a CONTROL, defined as the shipped-form `f` that matches the
  growth-anchored split. The split moves, so it must be re-derived. This carries no
  interpretive weight — it is a calibration of a control to its own definition.
* `_D260_LAMBDA` (5.0) is a DIAL PRICED AT A NUMBER: "an over-draw of ~5x on growth's gated
  demand lands the split in Crepin's band". If a smaller over-draw now lands it, **the price of
  the numerator lever has dropped and that is a finding**, reported as one. Predicted: lambda
  falls, because the baseline it multiplies has risen. It will NOT be retuned silently to keep
  a window green.

**P7 — THE ONE THAT COULD RETRACT A HEADLINE.** D-266's joint arm (blend + growth sink) landed
Crepin's leucine split at 54.50 % *and* Rollero's leucine tracer INSIDE his measured 3.4-8.2 %,
which was that record's headline and the reason its fork had two directions instead of one. The
tracer is `(1 - split) * consumed / isoamyl` with `consumed` pinned at the must's supply and
`isoamyl` invariant under numerator-side repairs (D-260 §1). So **raising the split must lower
the tracer.** Registered prediction: the joint split rises well above 54.5 %, and the tracer
falls toward Rollero's 3.4 % floor. Whether it crosses is NOT predicted — it is the measurement.

* If it stays inside: D-266's headline survives at a different point on the same line, and the
  fork stays two-directional.
* If it crosses below: **D-266's headline is retracted by its own composition being sourced.**
  The joint repair then satisfies neither paper's leucine observable, the "closer to reality
  needs one direction, this has two" fork loses one of its directions, and the owner-gated
  build D-266 §9 reserved is priced differently than when it was offered.

This arm is measured FIRST, before the D-259 arms, because it is the only result that could
change what the record says rather than only what a number reads.

**P8 — costs nobody has to pay twice.** D-266 §5 recorded that the joint arm over-shoots the
ile/val/thr splits by 15-43 points against Crepin. Raising `w_i` raises those splits too, so the
over-shoot GROWS. Registered so it cannot be reported later as a surprise.

## What would make this beat vacuous

If the measured splits move by less than the closure control's own 2 % tolerance, the sourcing
changed nothing measurable and the repair is a bookkeeping change. The multipliers above
(1.07-1.61x on the weight) make that unlikely, but it is the vacuity check.
