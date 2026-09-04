# D-267 findings — `w_i` is sourced, and the source does not fit the bracket

Numbers from `convert.py` → `findings.json`. No model run, no `src/` change: this is a
transcription plus arithmetic, and the arithmetic is checked against the engine's own molar
masses rather than against a table retyped here.

## What was open

D-259 §2 stated it plainly: the yeast protein amino-acid composition every growth-anchored
number rests on "is **not in `src/`**, not in the D-104 record, and has no receipts folder".
For want of a source that beat stated its own bracket — protein 40/45/50 % of dry weight,
leucine 6.0/7.5/9.0 g per 100 g protein, and so on — and reported every edge. D-266 §9 carried
it forward unchanged: "still a test constant and still recorded nowhere."

## No pre-registration, and what would have been in it

This is a transcription beat, not a measurement, so there was no arm to blind. The one thing
that could honestly have been registered is the expectation this beat began with — that the
bracket was an author estimate no source would reach, and that the beat would end by attaching a
citation to numbers that stayed put. **That expectation was wrong twice over** (§3, §4).

## 1. Where the search went, in order, because "nothing sources it" has been wrong ten times

| where | result |
|---|---|
| `manual_sources/*/_txt` (42 files) and the 24 textbook extractions | protein *fraction* twice (§2); no residue composition |
| `manual_sources/aem_00670-07 …pdf`, `mathematics-13-01373-v2.pdf` — the two PDFs with no text extraction, so no grep had ever seen inside them | Coleman, Fish & Block 2007 (already the engine's growth anchor) and Moimenta *et al.* 2025. Neither carries a composition |
| `d255-rollero-medium/` — the folder where D-256's tracer table had already hidden unread | Table S3 is a **must** composition (% YAN per source), not a biomass one |
| BioNumbers BNID 112800 → its primary reference | **the table** |

The eleventh instance did not happen: this one really was not in the repo. Recording that is the
point of listing the sweep — the shape has been wrong ten times, and being able to say *which*
places were checked is what makes "externally sourced" a finding rather than an excuse.

## 2. The protein fraction of dry weight: three anchors, all inside D-259's 0.40–0.50

| value | source | strength |
|---|---|---|
| 40–45 % | *Concise Encyclopedia of Wine and Winemaking* — "Wine yeast normally contains 40–45% protein (based on N × 6.25)". **In this repo** | wine yeast, but a *crude-protein* conversion, not a measurement of protein (§5) |
| 42 % | van Gulik & Heijnen (1995) Biotechnol Bioeng 48(6):681–698, Table I p.682, composition after Verduyn *et al.* | the standard stoichiometric yeast biomass |
| ~50 % | *Understanding Wine Chemistry* 2nd ed., end-of-chapter note — reasoning from **a bread-yeast packet's nutrition label**. In this repo | weakest of the three; not a wine strain and not a measurement |

D-259's `{lo: 0.40, mid: 0.45, hi: 0.50}` is corroborated at all three points and is not moved.

## 3. The residue composition: Lange & Heijnen 2001, Table IV, in mol %

**Lange HC, Heijnen JJ (2001), *Biotechnol Bioeng* 75(3):334–344, Table IV p.339 — "Amino acid
composition of the protein as measured (mol %)"**, glucose-limited chemostat *S. cerevisiae*;
the paper states the relative abundance did not vary between cultures, which is what licenses
carrying it to a fermenting must. Image: `lange_heijnen_2001_tableIV.png`, transcribed verbatim
into `convert.py`. The table's own closure checks: the 19 rows sum to **100.02** mol %.

The table is a **mole** composition and the constant it answers is a **mass** one, so the
conversion is where the finding lives. Asx and Glx are taken as the acids (what the paper's acid
hydrolysis delivers); the amide reading moves the mean residue mass by 0.22 % and no verdict.

## 4. THE FINDING: two frames, and the one the model needs puts all five above the bracket

Protein mass is the sum of **residue** masses. But `_GrowthAnchoredFates` subtracts its draw from
the **free** amino-acid pool, and the engine's `MOLAR_MASS` is the free acid — leucine 131.175,
not the residue's 113.160 (checked in `convert.py`: **zero mismatches** against `M_FREE`). So a
draw that is to deliver *x* g of leucine residue must remove *x* × 131.175/113.160 g from the
pool. The two frames differ by exactly the water released in peptide-bond formation, which is
why the free-acid shares sum to **116.53** per 100 g of protein and the residue shares to 100.00.

g per 100 g protein, against D-259's lo/mid/hi:

| | bracket | residue frame | verdict | **free-acid frame** | verdict | over `hi` | over `mid` |
|---|---|---|---|---|---|---|---|
| leucine | 6.0 / 7.5 / 9.0 | 8.315 | inside | **9.639** | **above hi** | 1.071× | 1.285× |
| isoleucine | 4.0 / 5.0 / 6.0 | 6.115 | **above hi** | **7.088** | **above hi** | 1.181× | 1.418× |
| valine | 4.5 / 5.5 / 6.5 | 6.666 | **above hi** | **7.878** | **above hi** | 1.212× | 1.432× |
| threonine | 4.0 / 5.0 / 6.0 | 5.166 | inside | **6.087** | **above hi** | 1.015× | 1.217× |
| phenylalanine | 3.5 / 4.5 / 5.5 | 5.077 | inside | **5.698** | **above hi** | 1.036× | 1.266× |

This is a **units fork, not a band** — the lesson that says so was earned at D-209. One reading
ships: the free-acid frame, because it is not a judgement about the literature but about this
model's own arithmetic. The residue frame is named as the loser and is what the constant's
comment currently says.

**Consequence, stated and not acted on.** Every growth-anchored number in D-259, D-260 and D-266
was computed at a per-precursor draw **21.7–43.2 % below** what the only measured composition plus
the engine's own molar masses imply, and the mid composition those records quote as the headline
is the furthest from it. The direction is unambiguous: a larger `w_i` draws more precursor into
the lump, so D-260's 27.58 % leucine split and D-266's joint 54.50 % are **low estimates**, and
D-266 §6's fork was priced at a composition the source does not support.

Not repaired here. Moving the constant re-prices every arm of an owner-gated build (D-266 §9),
which is the owner's call and not a transcription beat's.

## 5. A frame that does NOT reconcile, reported because it was checked

"40–45 % protein **based on N × 6.25**" is a crude-protein convention: it asserts total nitrogen
of 6.4–7.2 % of dry weight. The *Understanding Wine Chemistry* note's own arithmetic (protein is
one-sixth nitrogen, protein is ~50 % of dry weight) asserts 8.3 %. The engine ships
`biomass_N_fraction` = **0.114** g N/g cell, from Roels' generic CH₁.₈O₀.₅N₀.₂ — **1.6–1.8× the
nitrogen those two wine-yeast statements imply**, and running the convention the other way turns
0.114 into a 71 % crude-protein cell, far above every anchor in §2.

Nothing is done with this. It is a second, separate frame question, it touches the shipped
nitrogen budget rather than a test constant, and one of its two ends is a bread-yeast label. It
is recorded so the next beat that reaches for either number knows the pair does not close.

## Files

* `lange_heijnen_2001_tableIV.png` — the source table, as served by BioNumbers from BNID 112800
* `convert.py` — the transcription, both frames, the engine molar-mass cross-check
* `findings.json` — its output
