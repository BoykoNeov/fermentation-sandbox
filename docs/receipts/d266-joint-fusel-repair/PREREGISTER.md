# D-266 pre-registration — the JOINT fusel repair, measured for the first time

Written 2026-09-02 **before the first probe number**, in the D-259/D-260 discipline. The
receipts folder is *in the repo* this time (`docs/receipts/`): the beat ran in a cloud session
with no `M:\` drive, and a scratchpad dies with its container.

## The question

D-257 built the temporal fusel repair — `fusel_rate_shape = flux · [f·N/(K_n+N) + (1−f)] ·
arrhenius`, `f = 0.79` fitted on Rollero's SM250, the five `k` re-anchored ×0.4033 to hold the
Wang 2024 level anchor, and the amino-acid attribution scaled by the Ehrlich share of the rate
— and **reverted** it because the shipped precursor sink rides the alcohol rate, so slowing
production stranded 72.9 % of the must's phenylalanine. D-259/D-260 measured the growth-anchored
sink that clears exactly that blocker and **refused** it because, under the shipped rate law,
Crépin's leucine split and Rollero's leucine tracer trade one-for-one along
`tracer = (1 − split) · consumed / isoamyl`, with the fixture's ~2× inflated isoamyl setting the
line.

Each refusal's ground is the other repair's defect. **Nobody has run the pair.** This beat does,
beside each single repair and the shipped model, on every sourced observable the suite already
transcribes. Nothing in `src/` changes; the arms are test-side counterfactuals.

## Arms

| arm | rate law | precursor sink |
|---|---|---|
| `shipped` | gated flux form (D-99) | `PrecursorNonEhrlichFates`, `f/(1−f)` × Ehrlich draw (D-104) |
| `blend` | D-257 blend, `f`=0.79, `k`×0.4033, attribution scaled | shipped |
| `growth` | shipped | `_GrowthAnchoredFates` mid composition, growth's modifiers attached (D-260 §5) |
| `joint` | D-257 blend | growth-anchored, mid (lo/hi as a bracket on Crépin's must only) |

## Observables (all already in the suite as sourced constants)

1. Crépin's own must, other precursor consumers off (D-259/D-260 fixture): leucine split vs
   77–86 %, leucine tracer vs Rollero 3.4–8.2 %, isoamyl total (µM), phenylalanine left at
   end (%), phenylalanine lump share vs 0.975, ile/val/thr splits vs 51/41/38, quadrature closure.
2. Rollero's own three musts (D-255 fixture, full Process set): isoamyl present at nitrogen-gate
   closure as a fraction of final vs 42–54 %; isoamyl EF total vs his 1066–1337 / 1314–1365 /
   793–1034 µM; nitrogen response SM425/SM70 vs his 0.74–0.77×; valine enrichment vs Table S1;
   leucine enrichment vs Table S2 (3.4–3.5 / 4.2–4.7 / 6.8–8.2 %).
3. The D-112 level anchor: undosed isoamyl at YAN 250 and 300, 20 °C, vs 140–205 mg/L.
4. The D-104 fixture (brix 24, YAN 250, pitch 0.25, aa 1.0 g/L, 20 °C): phenylalanine left vs
   D-257's 20.3 % (shipped) / 72.9 % (blend); worst joint carbon refund vs growth's draw (< 1.0);
   net dS/dt never positive; carbon and nitrogen closure.

## Predictions, scored below after the run

- **P1 (harness reproduction).** `blend` reproduces D-257: NT fractions 36.3 / 48.1 / 55.9 %
  (±2 points), anchor 172.3 mg/L (±3), phenylalanine left on the D-104 fixture 72.9 % (±3
  points). `growth` reproduces D-260: leucine split 27.58 % and tracer 5.900 % on Crépin's must
  (±0.5 points / ±0.1 points). **If P1 misses, nothing below is attributable.**
- **P2.** `joint` clears D-257's blocker: phenylalanine left on the D-104 fixture ≤ 10 %, and on
  Crépin's must ≤ 5 %. Mechanism: growth's protein demand strips the pool whatever the alcohol
  rate does (D-259 §6).
- **P3.** `joint` leucine split on Crépin's must rises ≥ 25 points above `growth`'s 27.6 %, to
  55–80 %. Mechanism: the pool is fully consumed either way and the blend cuts the Ehrlich draw
  ~3× (0.79 × 0.4033), so the lump absorbs the residue.
- **P4.** `joint` leucine tracer on Crépin's must falls **below** Rollero's 3.4 % floor — the
  D-260 one-knob identity holds and the split gain is paid on the tracer — unless isoamyl on this
  must falls under ~1500 µM, which P5 says it does not.
- **P5.** `joint` isoamyl on Crépin's must lands 1400–2000 µM (from 2123): the level is
  anchor-preserved at YAN 250 / 20 °C by construction, and this must is 180 mg N/L at 28 °C.
- **P6.** `joint` NT fractions on Rollero's three musts stay within 3 points of `blend`'s: the
  sink swap does not touch production.
- **P7.** `joint` nitrogen response of isoamyl across Rollero's range ≤ 1.4× (D-257's 1.30×).
- **P8.** `joint` leucine enrichment on Rollero's own musts sits below the Table S2 floor on all
  three (the shipped form is already under at 0.99 / 2.05 / 2.65 %, and the blend cuts the draw
  while production is level-preserved).
- **P9.** The level anchor holds under `joint` at both YAN levels (an undosed run carries no
  sink; identical to `blend`).
- **P10.** Under `joint` on the D-104 fixture: carbon and nitrogen close to ≤ 1e-6 relative,
  net dS/dt ≤ 0 everywhere, worst joint C refund < 1.0.

**Expected verdict, stated so it can be wrong:** the joint arm keeps D-257's temporal gains (P6,
P7), clears D-257's blocker (P2), moves the leucine split most of the way to Crépin (P3), and
pays for it on the leucine tracer (P4, P8). If so, the two refusals were NOT each other's defect
— the tracer miss is the joint arm's own — and what binds is the isoamyl level the anchor sets
against synthetic-medium totals, not either mechanism.

## Abandonment rule

If P1 misses, the harness is wrong and the beat records only that. If P3 and P4 both miss in the
*favourable* direction (split ≥ 77 % **and** tracer ≥ 3.4 %), the joint arm satisfies both
sourced targets and the beat's headline changes: report it as such, re-check the closure and the
patch bindings (an arms-differ assert) before believing it.
