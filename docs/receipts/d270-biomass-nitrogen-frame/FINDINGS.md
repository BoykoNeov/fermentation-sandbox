# D-270 — what `biomass_N_fraction` is in each medium

Receipts for the record. `probe.py` compiles the scenarios and prints every number below;
`findings.json` is its output. Nothing here integrates a model — the override is a pure function
of the scenario's declared assimilable nitrogen, evaluated at the compile boundary.

## 1. The sourced range

Two wine-yeast protein statements, each through its own stated convention (D-267 §2 transcribed
the anchors; §6 did this conversion and this record re-uses it rather than re-deriving it):

| statement | figure | convention | total cell N |
|---|---|---|---|
| *Concise Encyclopedia of Wine and Winemaking* | 40 % protein | "based on N × 6.25" | **0.0640** |
| *Concise Encyclopedia of Wine and Winemaking* | 45 % protein | "based on N × 6.25" | **0.0720** |
| *Understanding Wine Chemistry* 2nd ed, chapter note | 50 % protein | protein is ⅙ nitrogen | **0.0833** |

van Gulik & Heijnen 1995's 0.42 is the fourth anchor D-267 §2 lists and is deliberately excluded:
it is a protein mass from a chemostat balance with no stated nitrogen convention, so converting it
would mean picking one on the source's behalf.

**Sourced range: 0.0640 – 0.0833 g N / g dry weight.**

## 2. What each medium's compiled run actually uses

`probe.py`, run 2026-09-04. `tier` is the compiled parameter's own tier (1 = plausible,
0 = speculative — the D-244 hold drops it).

| arm | declared YAN mg N/L | compiled `f_N` | vs sourced | tier |
|---|---|---|---|---|
| wine yan=50 | 50.0 | 0.0362 | BELOW | plausible |
| wine yan=80 | 80.0 | 0.0403 | BELOW | plausible |
| wine yan=100 | 100.0 | 0.0433 | BELOW | plausible |
| wine yan=150 | 150.0 | 0.0519 | BELOW | plausible |
| wine yan=250 | 250.0 | **0.0745** | **inside** | plausible |
| wine yan=300 | 300.0 | 0.0892 | ABOVE | plausible |
| wine yan=330 | 330.0 | 0.0994 | ABOVE | plausible |
| wine yan=350 | 350.0 | 0.1068 | ABOVE | plausible |
| wine yan=400 | 400.0 | 0.1068 (held) | ABOVE | speculative |
| wine yan=500 | 500.0 | 0.1068 (held) | ABOVE | speculative |
| wine yan=250 + 2 g/L aa | 700.9 | 0.1068 (held) | ABOVE | speculative |
| wine yan=250 + 4 g/L aa | 1151.8 | 0.1068 (held) | ABOVE | speculative |
| wine yan=80 + 2 g/L aa | 530.9 | 0.1068 (held) | ABOVE | speculative |
| **beer yan=200** | 200.0 | **0.1140** (static, no override) | ABOVE | plausible |

Census: **4 below / 1 inside / 8 above**. Span 0.0362 → 0.1068, a factor **2.95**.

## 3. The circularity

`convert.py` corroborates Coleman's 4 × 10⁻¹¹ g gram by inverting his yield to cells per gram of
nitrogen and pricing it with an elemental composition he had no hand in:

```
Y_X/N(330)      = exp(3.50 − 3.61e-3 × 330) = 10.0613 g cell / g N
cells per g N   = 10.0613 / 4e-11           = 2.5153e11
dry mass/cell   = (1 / f_N) / 2.5153e11
```

| `f_N` fed in | source of that number | pg / cell |
|---|---|---|
| 0.1140 | Roels `CH₁.₈O₀.₅N₀.₂` (static, shipped) | **34.87** |
| 0.0994 | **the engine's own compiled value at 330 mg N/L** | **40.00** |
| 0.0833 | sourced high (UWC) | 47.71 |
| 0.0720 | sourced mid (45 % crude protein) | 55.22 |
| 0.0640 | sourced low (40 % crude protein) | 62.12 |
| 0.0800 | static band low edge | 49.70 |
| 0.1400 | static band high edge | 28.40 |

The 40.00 is exact, not a coincidence: compiled `f_N ≡ 1/Y_X/N`, so the expression collapses to
`Y_X/N · 4e-11 / Y_X/N`. **The compiled fraction carries no compositional information**, and
`convert.py` is right to use the static elemental value.

## 4. Beer's price

Beer's growth is nitrogen-limited, so the ceiling is the identity `dX = YAN / f_N`.

| `f_N` | biomass ceiling vs shipped | growth extent |
|---|---|---|
| 0.1140 (shipped) | 1.000× | 5.378× |
| 0.0833 (sourced high) | 1.368× | **7.36×** |
| 0.0720 (sourced mid) | 1.583× | **8.52×** |
| 0.0640 (sourced low) | 1.781× | **9.58×** |

Against Tyrell 2013's counted 2.918–3.483× and *The Chemistry of Beer*'s printed 4–5×. Every
sourced value puts the extent above both targets, the mildest by 1.47× on the printed one.

## 5. D-230's refusal, re-scored

The cell nitrogen Tyrell's counted crop would demand at the engine's own gram is 0.202–0.262.

| measured against | outside by |
|---|---|
| static band top 0.14 | 1.44 – 1.87× |
| shipped static 0.114 | 1.77 – 2.30× |
| sourced high 0.0833 | 2.42 – 3.14× |
| sourced low 0.0640 | 3.16 – 4.09× |

Sourcing the composition moves the demand **further** outside. The refusal strengthens.

## 6. The side effect, reported and not acted on

Re-priced through §3, the elemental route gives **47.71 – 62.12 pg/cell** rather than 34.87,
against D-219's settled 40 pg and its 28–50 pg band (low edge inside, high edge outside). D-230's
branch 1 reads Tyrell's counts as 70.9–91.9 pg/cell; the gap narrows from **2.03–2.64×** to
**1.14–1.93×**. Evidence on a residue D-232 could only widen — not a settlement, and D-219's gram
is the unit wine's fitted parameters live in, so re-opening it is that record's to re-open.
