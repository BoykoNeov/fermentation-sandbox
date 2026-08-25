---
name: prohibition-beer-aroma-calibration-levels
description: "D-224 through D-227 — beer's aroma calibration: the eight level constants, the growth-extent coupling, and the NINTH constant (the stripping sink) none of the first three audits could see"
metadata: 
  node_type: memory
  type: project
  originSessionId: f44fbd81-53a8-4e00-ac37-9275b8fd1ab0
  modified: 2026-08-18T22:44:16.620Z
---

# D-227 — THE STRIPPING SINK, and the NINTH constant the eight-constant audits could not see. Read this FIRST.

D-226's open item 1 is **ENTERED AND CLOSED**. `EsterVolatilization` now strips on
`fermentative_co2_rate` — the uptake Process's own CO₂ arithmetic, shared from `carbon_routing.py`
(the D-180 discipline) — not on the `X·S/(K+S)` Monod its docstring called a "CO2-evolution proxy".

## What must never be re-argued

* **`k_ester_volatil` was a NINTH drift-prone constant and every audit missed it.** Both yaml files
  said, in the unit comment, that it *"folds q_sugar_max*co2_yield*scale into one constant"* — a
  ferment-SPEED knob inside a stripping constant, so the pair is LOCKED. D-223 moved
  `q_sugar_max` 0.5 → 0.72 and the note still cited **0.5** four beats later.
  **D-225's registry rule could not catch it: `ESTER_SPECS`/`FUSEL_SPECS` enumerate POOLS, and this
  is ONE constant shared across three.** [[feedback-grep-finds-claims-not-guards]]
* **Total evolved CO₂ is speed-invariant to 6 figures; the retired driver spans 1.64×** (across
  `q_sugar_max`'s band AND the retired 0.5). It is a **WEAKER identity than D-226's** and must be
  quoted as such: the alcohols' `YAN/biomass_N_fraction` holds at ANY horizon, evolved CO₂ only
  **once the sugar is gone**. Say "speed-invariant wherever the ferment completes".
* **Beer's largest historical aroma drift: 13.91 % → 0.73 % (19×).** Band edges 4.51 % → 0.33 %.
  **The residue CHANGED OWNER — it is now `mu_max`'s (0.54 %), not `q_sugar_max`'s, and `mu_max`'s
  arms barely improve.** That was PRE-REGISTERED, not a disappointment: the CO₂ integral fixes the
  driver's SIZE, never its TIMING, and a first-order sink against a moving pool reads both.
  **Never claim exact ester invariance.**
* **WINE IS THE PROOF, and it is a DERIVATION not a fit.** Wine's re-anchor factor 2.52530 equals
  the analytic `1/(q_sugar_max·co2_yield·scale)` = 2.5252809 to **six figures** — one sugar, no
  repression, so in wine the beat is an EXACT reparameterisation and every wine level is unchanged
  to 6 figures (wine stays D-224/D-226's control). **Beer's 3.34686 is 17.7 % off its analytic
  2.84279, and THAT GAP IS the three-sugar catabolite repression + the per-species CO₂ yields
  (0.4886/0.5143/0.5235).** Same edit: bookkeeping in one medium, physics in the other.
* **Units changed with the rate law: `L/(g*h)` → `L/g`** (wine 1.26264e-2, beer 6.6937e-3). Bands
  are **TRANSPORTED by the same factor, never re-sourced** — neither edge was ever sourced and the
  transport does not make one so. `k_ester_volatil` is **still an author estimate in both media**;
  what changed is it now holds only the gas-volume/Henry prefactor (wine:beer ratio 1.886, pure
  prefactor, nothing to do with speed).
* **The driver EXCLUDES non-fermentative CO₂ ON PURPOSE and it is SIZED: 0.0018 % (beer) /
  0.635 % (wine).** Reading the whole `CO2` slot would make the sink self-referential — Ehrlich
  decarboxylation co-produces the alcohols whose esters it strips
  [[feedback-a-returned-intermediate-is-a-rate-law-change]]. Do not "complete" it.
* **One level moved: beer's isoamyl acetate −0.216 %** (it alone reads its precursor pool, D-97, so
  the k solve that holds the other two EXACTLY cannot hold it); `fruity` 1.8333 → 1.8294.
  **D-224's `approx(1.840, abs=0.01)` was 0.36 % off the value it pinned** — re-pinned to measured.
* **Every LEVEL guard is a designed GREEN under a full revert** (the re-anchor holds the calibration
  frame by construction), which is why the beat ships a guard that names the **DRIVER**
  [[feedback-prefer-the-variant-your-guards-can-see]]. Arm A predicted exactly 5 RED and got exactly
  those 5, its corner RED reporting D-226's own 0.9555.
* **Temperature is UNTOUCHED** (fusels 1.21822, esters 6.89570 vs D-226's solved 1.2182 / 6.8953),
  which is what makes the speed result attributable.

## Still open after D-227

1. **`mu_max`'s TIMING residue** — not fixable by a different driver; the gas really does arrive
   earlier in a faster beer. Needs a synthesis term that does not finish before the stream does.
2. **The growth Process's hard nitrogen cutoff** (day 0.92, 81 % of sugar left) — unchanged.
3. **`E_a_fusels`'s magnitude**, **wine's two calibration frames**, **Wang's Table 1** — unchanged.
4. **`k_ester_volatil` is unsourced and both bands are transported, not measured.**

---

# D-226 — THE COUPLING, ENTERED AND CLOSED FOR BEER. Read this before anything below.

D-224 §11 / D-225 §11's "the extent-coupling question is still the real one and is still not
entered" is **ENTERED AND ANSWERED**. Beer's aroma pools now ride growth EXTENT:
`k · mu·X · f_growth(T) · arrh(E_a)`, via `EsterSynthesisGrowthCoupled` /
`FuselAlcoholsEhrlichGrowthCoupled`, wired per medium in `media.py`.

## What must never be re-argued

* **The headline is FIDELITY, not drift.** Beer's esters were YAN-BLIND and slightly BACKWARDS
  (4× wort nitrogen moved them 1.009 → 0.983); they are now proportional (0.491 → 2.072). D-97
  found this and fixed it for isoamyl acetate ALONE. The five higher alcohols were already
  YAN-proportional under both forms — **the YAN result is about the ESTERS only**.
* **BEER ONLY, and it is sourced.** Wine's `E_a_esters` = 55,100 exists solely to make integrated
  synthesis T-flat under the FLUX form; under extent coupling flatness needs ≈ 0. And
  `wine_generic.yaml`'s own header says **no** value of it reproduces wine ester behaviour — wine's
  flux form is a **documented stand-in**. Do not "finish the job" by converting wine.
* **It is NOT full invariance and the notes say so.** 5 higher alcohols: **1.00000** at both speed
  knobs' band edges AND their retired values (no stripping sink). 3 esters: **4.51 %**, sign
  REVERSED, because `EsterVolatilization` still reads the flux. D-223's retired `q_sugar_max` 0.5
  would still move them **13.9 %** (was 20.9 %). Never say "the drift is closed" without "for five
  of eight".
* **A CONSTANT LIVES IN ITS RATE LAW'S COORDINATES.** Under the flux form the observable's apparent
  `E_a` was `E_a_esters − E_a_uptake` = 144,900, not the 200,000 printed; `E_a_fusels` printed
  70,000 with a "Q10 ~ 2.6" note while the output ran **Q10 1.22**. Both re-anchored (152,000 /
  14,100) to PRESERVE the measured spans, deliberately not to adopt the steeper reading — the
  magnitude is unsourced either way and two changes at once are unattributable.
* **D-19's ordering constraint is a property of the RATE LAW.** "`E_a` > `E_a_uptake`" is the FLUX
  form of "output rises with T"; the extent form is "`E_a` > 0". Beer's `E_a_fusels` 14,100 is
  BELOW `E_a_uptake` and that is correct. The test asserts per coupling AND asserts which coupling
  each medium is wired to.
* **`E_a_esters`'s low band edge is 87,000 and is MEASURED** (floor 84,890 at the joint stripping
  corner), not the retired algebraic 118,000+2,000. **Do not cite the stripping-ordering floor as
  the reason 152,000 sits below the cited 221-265k** — that floor constrains nothing now; what
  holds it there is the preserved span plus the lumped-fit-artifact argument.
* **Luedeking-Piret (`α·dX/dt + β·X`) is REFUSED, not overlooked**: two free parameters, one
  observable, no ester time-course on disk. Unblocked only by a published beer ester time-course.
* **PITCH is inert under BOTH forms** (1.000 across 0.25-4.0 g/L) and is not a discriminator.
  Over-pitching does suppress esters in reality; neither form captures it. Never claim it.
* **D-97's stranded-precursor objection is FALSE on today's tree** (precursor complete at
  growth-stop, ratio 1.048) — do not raise it against extent coupling again.

## Still open after D-226

1. ~~The stripping sink reads the flux SHAPE, not evolved CO₂~~ — **CLOSED at D-227, top block.**
2. **The growth Process's hard nitrogen cutoff**: growth stops day 0.92 with **81 % of the sugar
   left**, so all ester is made in day 1 and only stripped after. Inherited, not introduced.
3. **`E_a_fusels`'s magnitude** is unsourced in either coordinate system (1.22 shipped vs 2.66).
4. **Wine's two calibration frames** (below) — untouched and still open.

---

# D-225 — THE EIGHTH CONSTANT, and why D-224's guard could not see it

**Read this before trusting anything below about "the seven".** D-224 repaired seven aroma
constants and its own mechanism condemned **eight**. `k_ethyl_hexanoate` (beer) shares
`k_ethyl_acetate`'s rate law exactly and carried the whole drift unrepaired: **0.17321 mg/L
against its stated 0.22 (0.787×), OAV 0.825** where the file says it sits AT Meilgaard's
~0.21 threshold. Re-anchored **1.2e-6 → 1.524e-6 (×1.27)**, landing 0.21998.

## What must never be re-argued

* **The population is `ESTER_SPECS` + `FUSEL_SPECS`, never a list.** D-224's level guard read a
  hand-written dict of seven, so it could not see the eighth — and it stayed **GREEN** when the
  defect was put back (measured, mutation arm A). Guards are registry-driven now. **Do not add a
  ninth entry to a literal dict**; that repeats the defect a fourth time.
* **The census is closed and its denominator is printed.** 16 constants are defined by a landing
  level; 12 state it in `conditions:`, 4 stated it only in `notes:`. **No other aroma pool is
  level-defined at all** — not acetaldehyde, H₂S, DMS, diacetyl, methanethiol, the ethylphenols,
  methional, sotolon, furaneol, vanillin, guaiacol, eugenol, whiskey lactone or the aldehydes.
  Do not re-run this census; do not report a gap outside those 16.
* **All three beer esters are MECHANICALLY LOCKED** — same rate law, same `E_a_esters`, same
  stripping constant — so their ratios are exactly `k_i/k_j`, invariant to **7 significant
  figures** across both speed knobs' band edges and a 15 °C run. Consequences: one pool's history
  transfers to the others with **no old tree to check out**, and the next speed change moves all
  three by one factor, guaranteed rather than probable.
* **`k_ethyl_hexanoate`'s band is COMPUTED, not rescaled** — `[6.93e-7, 3.46e-6]`, spanning
  exactly its sourced **0.100-0.500 mg/L** ale range. Rescaling ×1.27 agrees almost perfectly
  with the note's stated "~0.05-0.6" (which is why ×1.27 is right for the NOMINAL) but its top
  reaches **0.6416 mg/L, 28 % above the sourced range** — a value the source rejects, reachable
  in a drawn ensemble. The "~0.05-0.6" gloss is the named LOSER. D-224 declined the computed form
  for isoamyl acetate only because THAT span is forked; this one's is single.
* **The OAV crossing is INTENDED and the axis does NOT move.** 0.825 → 1.048 at the pool;
  `fruity` stays owned by isoamyl acetate at 1.840 under the MAX rule, `above_threshold()` is
  `['fruity']` either way. Both halves are pinned.
* **WINE: three ester constants, TWO frames, and NO value was moved.** Wine's esters reproduce
  their stated levels at **YAN 80** (0.992-1.000×), its five Ehrlich constants at **YAN 250** —
  and the one wine ester that reads its precursor is **1.65× apart** between them. At YAN 80 the
  precursor sits at 100.84 mg/L against the **172 its own `conditions:` targets**. `sensory.yaml`'s
  "lands ~25" OAV note is a YAN-80 reading (41.85 at YAN 250). **Do not "fix" this by re-anchoring
  wine's `k_isoamyl_acetate`:** one constant cannot satisfy both frames (÷1.652 puts it at 0.456
  at YAN 80, below the Guth assay it is anchored to), 0.76 IS Guth's assay, and **wine is D-224's
  control**. The frames are now STATED in all three `conditions:` fields; the inconsistency is
  **open**, and resolving it needs a sourced reason to prefer one frame.
* **`aging.yaml`'s `ethyl_hexanoate_eq` is calibrated to WINE's finished level** as its young-wine
  reference (agrees today: 0.3939 vs stated ~0.4). A wine re-anchoring stales it silently; the
  dependency is now cross-referenced from the wine parameter that owns it.

## Still open after D-225

1. **The extent-coupling question** — unchanged from D-224 §11 and still the real one. New
   evidence: the mechanical lock means all three beer esters drift **together, exactly**, so every
   future speed change re-commits this. The registry guard catches it the day it happens; it does
   not prevent it.
2. **Wine's two calibration frames** (above).
3. **The narrowing is itself a claim** — `k_ethyl_hexanoate`'s low edge moved from 2.3× below its
   sourced ale floor to exactly on it. Justified by the band note claiming to BE the ale range,
   but no band edge here is separately sourced.

---

# Beer's aroma calibration levels (D-224) — SETTLED

Reached by path from the ledger in [[project-fermentation-sandbox]]. Full records: `docs/DECISIONS.md` D-224 and **D-225** (read the D-225 block above first — it corrects the scope of everything below).

## The verdict — do not re-open "which ester calibration is wrong"

**Neither.** D-223 §8's "one of the two is wrong and nothing here adjudicates it" is **ANSWERED**.
Isoamyl acetate is first-order in the `isoamyl_alcohol` pool (D-97 ATF1 coupling) and inherited a
**1.61×** error in its precursor; ethyl acetate's precursor is ethanol, which saturates the same
enzyme, so it did not. Put the precursor back on its own published mean and **both esters read 0.79×**
their targets — one common factor, the ferment-speed change D-223 shipped.

**Quote the 3 %, never the 0.4 %.** The 0.793 / 0.790 agreement is two disagreements cancelling: the
D-99 slacks disagree −2.35 %, the moves-since disagree +2.86 %. Defensible claim = *both moved by one
common factor 0.72-0.74, to within 3 %*. [[feedback-a-hit-can-be-two-errors-cancelling]]

## What is BUILT (7 values + 5 bands + 2 band rescalings) — never re-anchor without reading D-224

`k_propanol` 3.55e-4, `k_isobutanol` 3.41e-4, `k_active_amyl_alcohol` 3.65e-4, `k_isoamyl_alcohol`
1.13e-3, `k_2_phenylethanol` 9.12e-4, `k_ethyl_acetate` 1.39e-4, `k_isoamyl_acetate` 5.36e-4. All seven
land the level their own `conditions:` field states, **within 0.4 %**, at the calibration frame
**21 d / 20 °C / YAN 200 / pitch 1.0** — and all seven were inside their pre-existing bands, so no band
decision is hiding in a value.

**The targets are the STATED ones (20.0, 2.2, and Wang's four means), not what D-99 happened to land**
(21.32, 2.402). Absorbing that 6.6 % / 9.2 % slack is deliberate: it is what makes the two esters agree
**by construction** instead of by cancellation.

**Two figures already printed in `beer_generic.yaml` are RESTORED, not computed, by this:** the pool at
"~0.34 mM, ~87× below ATF1's Km" (shipped state: 0.547 mM / 54.4×) and "OAV ~0.6" against Meilgaard
(shipped: 0.965). That is the corroboration that 30.0 mg/L is the level the file was written for.

## The mechanism — the ONE sentence this reduces to

`FuselAlcoholsEhrlich` is **GATED on nitrogen (`N/(K_n+N)`) but draws its carbon from `S`**, so its
output is not **BOUNDED** by the nitrogen it is nominally made from: a slower-growing yeast holds the
gate open longer and makes more higher alcohol from the same YAN (every run ends at `N ≈ 0`). Therefore:

* **`mu_max` is the higher-alcohol knob** (band edges ×1.094 / ×0.774) and leaves ethyl acetate alone.
* **`q_sugar_max` is the ester knob** (×1.111 / ×0.897) and leaves the higher alcohols alone.
* **`isoamyl_acetate` reads BOTH.** That is why it looked like the disagreement.

**The separation is a property of the two KNOBS, not of the pools** — they share `S`. A 1.68× change in
the Ehrlich `k` moves packaged ethyl acetate **0.077 %** (0.28 % vs 0.19 % of sugar drawn), enough to
break D-223's corner pin. Never say the pools are decoupled.

## The drift, and that it was never scored

D-99 reproduced **all four** Wang/Frank/Steinhaus 2024 Table 1 beer means to **<0.5 %**. **D-211's
`mu_max` 0.098 → 0.034 multiplied every one by 2.87; D-222's → 0.058 halved it to 1.68. NEITHER RECORD
MENTIONS HIGHER ALCOHOLS** (grepped). Wine is the control and never moved (`mu_max` is per medium):
**all five** land inside 0.85 % (24.014/24.0, 32.997/33.0, 70.140/70.1, 170.539/172.0, 28.713/28.7),
and that is now a test. **Wine's five BANDS are also correct** (×0.3/×3, ×0.2/×4 for its own propanol)
even though 955ebbc shipped them in the same commit — so the draughting error is **beer-only**, and the
band guard is parametrized over both media so nobody re-asks it by assumption.

## The OUTPUT-level consequence — beer's solventy axis changed OWNER

Under the D-95 MAX rule: **`isoamyl_alcohol` (OAV 0.9652) → `ethyl_acetate` (0.6688)**, magnitude down
31 %. `fruity` keeps `isoamyl_acetate` and falls 21 % (2.3292 → 1.8404). **Both verdicts unchanged and
correct**: solventy below threshold, fruity above it — the banana-forward ale D-96 anchored for. Pinned,
because nothing asserted a descriptor's owner either.

## The sensory statement — this is why it is a DEFECT, not a preference

`isoamyl_alcohol` is the **only** beer pool of the five with a sourced in-matrix threshold (Meilgaard
1975, ~50 mg/L). The shipped model ran **48.261 mg/L, OAV 0.965** — within 3.6 % of claiming a solventy
note the file says a sound ale must not have — and **over the threshold (52.809) at `mu_max`'s low edge**.

**COUNT THE CROSSING JOINTLY, NEVER PER-PARAMETER** (512-member LHS, 75 sampled params, both arms):
`isoamyl_alcohol` P(> 50 mg/L) **33.20 % → 33.98 %** — unchanged inside ±2.1 pp noise, so **the repair
moves the NOMINAL (0.965 → 0.602) and NOT the tail**. `ethyl_acetate` P(> 30 mg/L) **3.52 % → 8.79 %** —
this one DID move, 2.5×, because its band was corrected to span its whole sourced range.
`isoamyl_acetate` P(> 1.2) 74.4 % → 80.9 % (crossed by design). The first draft used a per-parameter
marginal (29 % → 33 %) and got the right conclusion for the wrong reason
[[feedback-a-band-is-per-parameter-a-claim-is-joint]]. **Also: the ensemble MEDIAN higher-alcohol level
is 40.9 mg/L against a 30.1 nominal** — the ×0.3/×3 band is right-skewed, so a drawn beer is typically
a third fusel-ier than the nominal one. Never read an ensemble median as the model's answer.
[[feedback-a-margin-is-a-claim-about-what-holds-it-open]]

## The band fork — the loser is named, do not re-argue it

D-99's five Ehrlich bands are ×0.3/×3 (propanol ×0.2/×5) of a centre **2.05× below** the nominal shipped
in the same commit, so the multipliers in force were **×0.145/×1.45** and the notes were false from
birth (both halves landed at 955ebbc; checked with `git log -S`). [[feedback-pin-the-band-not-the-nominal]]

* **SHIPPED:** edges = the stated multiple of the corrected nominal (isoamyl alcohol spans 9-90 mg/L).
* **REJECTED:** rescale the shipped edges, preserving ×0.145/×1.45 — it puts Meilgaard's 50 mg/L
  **outside `k_isoamyl_alcohol`'s band entirely**, i.e. asserts no ale can be fusel-y. Nothing sources that.
* **THREE rules, not two, and all three are tests** — do not apply one where another belongs.
  `k_isoamyl_acetate` is **rescaled** with its nominal (D-97/D-99 convention: its stated span and its
  actual one agree, and its threshold is crossed at the NOMINAL by design, so no claim to protect; its
  band top moved DOWN 5.87 → 4.47 mg/L). **`k_ethyl_acetate` is COMPUTED, not rescaled**, band
  `[6.93e-5, 2.08e-4]` spanning exactly its sourced **10.00-30.02 mg/L** ale range. The first pass of
  D-224 rescaled it and put the top at 31.03 mg/L — **over the sourced range AND over Meilgaard's ~30
  mg/L threshold (OAV 1.034)** beside a note still claiming OAV < 1 at the nominal alone. That is this
  beat's own headline defect, committed inside the fix for it, and caught on review.
  **The band top now reads OAV 1.001 where it read 0.817 before D-224 — a RESTORATION, not an
  over-reach** (the source says a sound ale sits AT OR BELOW the threshold; the old band could not reach
  it only because the nominal was 21 % low). Joint corner with `q_sugar_max`'s low edge: **33.36 mg/L,
  OAV 1.112**, guarded and stated.

## The cost, and why the third ester option was refused

D-223's `Byp` funding constraint is exercised by a thinner slice: joint-corner flip fraction **5.37 % →
0.0767 %** (~1 draw in 1300; the corner still FORMS, +0.518 mg/L). **Anchoring ethyl acetate at Wang's
23.7 mg/L would take it to exactly 0.00 %** and leave that constraint invisible again — plus D-176's
independent reason (that survey mean carries a sour-beer tail). **Never move ethyl acetate to 23.7.**

**The 20 °C frame is load-bearing and now written down:** `E_a_esters` = 200 kJ/mol, so the same run at
15 °C lands ethyl acetate at **6.10 mg/L, below its own 10 mg/L floor**. 20 °C is a typical ale ferment;
§2.2's 15 °C is Foster's cool trial (D-221), a **different** frame. Duration is irrelevant past dryness.

## Guards (5 new) — the durable half

`test_the_finished_beer_lands_the_aroma_levels_its_rate_constants_are_defined_by` (equality against the
TARGET, never a snapshot), `..._each_drawn_speed_knob_moves_only_the_half_of_the_aroma_set_it_is_coupled_to`
(at the drawn band EDGES, **READ from the file** — both bands were rebuilt inside the last three records,
so a literal would silently become an interior point), `..._isoamyl_alcohol_stays_below_its_only_sourced_threshold_across_the_growth_band`,
and the three band rules. **The seven targets are asserted against the `conditions:` sentence that
specifies them**, not transcribed — and that clause had to be ADDED to the two ester entries, which
carried their level in `notes:` only. Falsified in 3 arms, designed GREEN in each, restore verified by
SHA-256. **Two arms missed their predictions and both misses paid.** Arm A: reverting the Ehrlich `k`
broke D-223's corner pin (the pools share `S`), correcting "leaves ethyl acetate EXACTLY alone" → 0.03 %.
Arm B: the isoamyl-acetate band test stayed GREEN, because **a MULTIPLIER guard is invariant to a joint
rescale of value and band by construction** — so of the two ester band rules the COMPUTED-edge one is
strictly the stronger guard. Not fixed: a computed span for isoamyl acetate needs a stated range, and
its note's "~0.7-4" disagrees with its own `source:` field's "~0.5-3". **Unresolved fork, not opened.**

## OPEN — named, not licensed

**The coupling itself.** Beer's aroma pools are coupled to **biomass-hours**; both cited mechanisms are
**extent**-coupled — de Andrés-Toro 1998 forms ethyl acetate as `Y_EA·mu_x·X_A` (biomass FORMED,
nitrogen-limited ⇒ invariant to `q_sugar_max`), and the Ehrlich pathway's substrate is amino acids
(⇒ invariant to `mu_max`). **Under an extent-coupled rate law none of the seven levels would have moved
at D-211/D-222/D-223 and D-224 would not exist.** NOT licence to build it: the flux coupling is
D-19/D-21/D-96/D-97/D-99's, carries the sourced temperature ORDERING, and reaches wine too. Its own beat.

Also open: **Wang Table 1's beer means for the two ACETATE esters were never sought** (paper not on
disk); `q_sugar_max`/`mu_max` are NOT re-visited — they are fitted to measured courses and the aroma
levels are the consequence.
