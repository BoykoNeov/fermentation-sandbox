# Milestone 3 — Tier-3 (speculative frontier): sensory/OAV + aging chemistry

> Status: **beat 1a built (D-67); aging axis opened (D-68); first aging Process built (D-69);
> aging-phase scenario wiring landed (D-70); oxidative aging axis opened (D-71); SO₂ scavenging
> — the first O₂ sink — built (D-72); O₂ sub-axis reworked for always-on sinks (D-73);
> `PhenolicBrowning` — the first always-on sink, and the first VISIBLE oxidative product (the
> `A420` browning index) — built (D-74); `StreckerDegradation` — the Strecker aldehydes
> (methional/phenylacetaldehyde), a wine-only substrate-gated O₂ sink — built (D-75).**
> The OAV sensory readout ships (`fermentation.sensory` + `sensory.yaml` + `tests/test_sensory_oav.py`).
> **D-69 built `EsterHydrolysis` (`core/kinetics/aging.py` + shared `aging.yaml` + `tests/test_aging.py`):
> net decay of the `esters` pool toward a lower equilibrium floor, Arrhenius warmer-ages-faster, released
> carbon split 5:2 → `fusels` + `Byp` (advisor-settled crux; carbon closes to machine precision).
> **D-70 wired it into the scenario pipeline (753 tests green): `EsterHydrolysis` into both media
> (disabled at compile — aging is inherently post-ferment), a `begin_aging` intervention verb (the
> `pitch_mlf` reconfigure pattern MINUS the state mutation) that enables it over a post-fermentation
> segment whose span is set by `duration_days`, and `aging.yaml` into `compile.py`'s `shared_files`. The
> advisor reframe: the load-bearing call was the deferred "what runs during aging" — settled Stance A
> (leave the ferment set on) because every producer of `esters`/`fusels`/`Byp` is fermentative-flux-gated
> and quiescent at dryness, so the aging signal is unconfounded. §7 slow-phase = the segment restart (no
> new machinery). An un-aged run is byte-for-byte the pre-aging core (isolable).**
> **D-71 opened the OXIDATIVE aging axis: `OxidativeAcetaldehyde` (dissolved O₂ oxidises ethanol →
> acetaldehyde, the 'sherry'/oxidised note the D-67 lens already reads) on a new carbon-free `o2` state
> pool + an `add_oxygen` dosing verb + 3 `aging.yaml` params (768 tests). Advisor crux: O₂, not ethanol,
> is the rate-limiter — first-order in the finite `o2` pool ⇒ SATURATING acetaldehyde, not the unbounded
> ethanol-first alternative; the owner chose the dissolved-O₂ pool (Approach B) over the smaller
> ethanol-first form, so `o2` is the shared substrate the whole future oxidative sub-axis (browning /
> Strecker / SO₂ consumption) will draw down. Carbon closes machine-precision (E→acetaldehyde, the D-27
> reduction reversed; O₂ off every ledger); reductive aging (`begin_aging` with no `add_oxygen`) is
> byte-for-byte the ester-only case.** D-72 built `SulfiteOxidation` (the first O₂ sink, wine-only,
> substrate-gated: "SO₂ protects until exhausted"); D-73 reworked the O₂ sub-axis so `k_ethanol_oxidation`
> is a *share* not the total, letting *always-on* sinks compose; D-74 built `PhenolicBrowning` — the first
> always-on sink and the first VISIBLE oxidative product: medium-agnostic, dominant O₂ share (diverts O₂
> from — and suppresses — oxidative acetaldehyde), accumulating a new `A420` browning-index state slot (an
> optical absorbance, off every ledger; the `iso_alpha` pattern, not the D-67 post-hoc OAV). **D-75 built
> `StreckerDegradation` — the O₂/amino-acid Strecker aldehydes methional (cooked-potato off-note) +
> phenylacetaldehyde (honey), the first aging Process to add aroma pools the D-67 lens did not read (two new
> wine slots + thresholds). CRUX (owner fork, superseding the D-71→D-74 forward-guess): gating the O₂ draw on
> `amino_acids` makes it DOUBLY substrate-gated (o2 AND amino_acids), like `SulfiteOxidation`, so it ADDS ON
> TOP with NO re-baseline of `k_ethanol_oxidation` — the "reduce k_ethanol again" note is retired. Carbon +
> nitrogen close (the D-45 mercaptan draw+deaminate idiom + a CO₂ decarb term); wine-only, silent without both
> substrates. Owner also chose TWO pools over one lump (opposite sensory valence).**
> **D-76 closed the emergent SUR-LIE → Strecker pathway with NO new physics: lees autolysis (`YeastAutolysis`,
> D-34) and `StreckerDegradation` (D-75) simply COMPOSE — opting into `autolysis_rate_per_h` + `add_oxygen` +
> `begin_aging` with no `amino_acids_gpl` dose lets dead lees self-digest, refilling `amino_acids` from the
> wine's own dead yeast, which the Strecker route degrades (both aldehydes emerge from the physically-real
> nitrogen source — the DIRECTIONAL result; the absolute level runs ~8× the D-75 dosed literature anchor, an
> order-of-magnitude figure not a prediction, a recorded open item in D-76). Owner fork (chose A: verify +
> document + test) over B (re-gate autolysis to
> the aging phase); a discriminating measurement decided it — the pre-dryness active-ferment release is ~15 mg/L
> (bounded-small), the ~385 mg/L breakpoint pool is legit post-dryness sur-lie, so autolysis-from-t0 needs no
> re-gating. Carbon + nitrogen close on the new triple-draw compose (autolysis refills while Strecker + mercaptan
> draw). Test-only change (helper kwarg + 3 scenario tests); zero production code.
>
> **D-77 built `OakExtraction` — the barrel/chip aroma-extractive axis, the FIRST NON-oxidative aging Process and a
> SEPARATE axis (draws no O₂). Four wood extractives — `whiskey_lactone` (coconut, light-toast), `vanillin` (vanilla,
> medium peak), `guaiacol` (smoky, heavy — the OAK guaiacol, distinct from Brett 4-EG), `eugenol` (clove, heavy) —
> diffuse in and rise toward a per-compound saturation ceiling: `d(C_i)/dt = k_oak·f(T)·max(0, ceiling_i − C_i)`, the
> inverse of EsterHydrolysis's decay-to-floor. `E_a_oak_extraction` deliberately WEAK (20 kJ/mol — diffusion, not
> reaction). The ceilings are SET-AND-HOLD off-ledger STATE slots the new `add_oak {oak_gpl, toast}` verb writes
> (`oak_gpl × toast-specific yield` from `oak.yaml`), enabled by `begin_aging`. KEY DESIGN TURN: the advisor's
> first-pass "mint the ceilings as provenance-backed params in the verb" recipe was overturned by primary-source
> evidence (verbs can't inject into the compiled ParameterSet; `param_update` is plain floats absent from
> `param_tiers`, and `begin_aging`-before-`add_oak` would KeyError mid-integrate) — a 2nd advisor pass agreed and
> dropped it; the state-slot dose (the `cation_charge`/`add_oxygen` idiom) meets every goal (provenance in the
> yields, D-1 moot since all oak pools floor SPECULATIVE, no KeyError window). OFF EVERY LEDGER (exogenous
> wood-derived, the `iso_alpha` precedent — cleanest aging Process, no `chemistry.py` change). Owner forks: 4
> compounds (+eugenol) and WINE-ONLY wired (like Strecker — medium-agnostic physics ≠ agnostic wiring). Explicit
> `ceiling ≤ 0` undershoot guard (the floor is 0, so `max()` alone fabricates extract). 8 wine-only slots, 838
> tests.**
> **D-78 (latest) built `EllagitanninOxidation` — the oak-tannin O₂-scavenging sink, the BRIDGE from the oak
> extractive axis to the O₂ sub-axis D-77 kept separate (854 tests, +16). Ellagitannin is oak's hydrolysable tannin:
> `OakExtraction` extracts a 5th pool (identical diffusion), and a new `EllagitanninOxidation` draws its share of the
> shared `o2` budget (bilinear `[o2]·[ellagitannin]`, the `SulfiteOxidation` form) and CONSUMES it. THE SPINE is oak
> PROTECTION: an oaked+oxygenated wine browns LESS (lower A420) and makes LESS oxidative acetaldehyde than an un-oaked
> wine at the same O₂ dose (the D-72 SO₂-protection threshold with a RENEWABLE buffer — the wood re-supplies tannin
> below the ceiling, so protection is sustained for months-to-years, unlike SO₂'s finite pool). SUBSTRATE-GATED on
> `[ellagitannin]` ⇒ adds ON TOP with NO re-baseline (the `k_ethanol+k_browning=5e-4` anchor untouched) even though
> it's a DOMINANT sink when present — proving the substrate-gated/always-on distinction, not the magnitude, is
> load-bearing. MASS-based consumption yield `y_ellag_per_o2` (g/g, no fake molar mass for the lumped macromolecule);
> off every ledger (both slots), so nothing conserved moves. ASTRINGENCY is a TASTE readout `analysis.astringency_series`
> (mg/L ellagitannin, IBU-exact, reads no threshold — the iso_alpha exclusion), NOT the OAV odor lens and NOT an
> A420-style slot (it tracks the current pool). Softening = ONE contributor (oxidative consumption); tannin–anthocyanin
> polymerization is the separate deferred beat (the `tannin` namespace left free). Toast ordering: ellagitannin DECLINES
> with toast (thermolabile). 1 advisor pass (endorsed all 4 forks). Next: beat 1b (descriptor projection), non-oxidative
> Maillard Strecker, tannin–anthocyanin polymerization (now unblocked), barrel fill-number, or barrel-beer oak — all
> deferred.**
>
> **PROGRESS (D-79, 2026-07-13): `TanninAnthocyaninCondensation` built — the DOMINANT red-wine astringency-softening +
> colour-stabilization mechanism, the EIGHTH aging Process, the SECOND non-oxidative one, and a THIRD separate axis on
> GRAPE pools (868 tests, +14). Free monomeric `anthocyanin` + condensed grape `tannin` (two GRAPE must inputs,
> `anthocyanin_gpl`/`tannin_gpl`, default 0 ⇒ white wine) condense (bilinear `k·f(T)·[anthocyanin]·[tannin]`,
> reaction-scale E_a) into a soft, SO₂/pH-STABLE polymeric pigment — the young-purple → aged-brick-red evolution.
> OAK- AND O₂-INDEPENDENT (the correctness crux): grape condensed `tannin` ≠ oak hydrolysable `ellagitannin`, and the
> Process draws NO o2 — a steel-tank red still polymerizes (reusing ellagitannin would wrongly require `add_oak`); it's
> the grape `tannin` the D-78 note left the namespace free for. DOUBLY substrate-gated ⇒ white/no-tannin wine
> byte-for-byte inert; off every ledger (grape-derived, no `chemistry.py` change). The polymeric pigment is a POST-HOC
> readout (`polymeric_pigment_series` = anthocyanin₀ − anthocyanin), NOT a slot — the A420 discriminator applied
> (anthocyanin's single fate ⇒ reconstructible), keeping v1 to 2 slots. Readouts: `astringency_series` now =
> (tannin + ellagitannin)×1000 (softens as tannin condenses; polymeric excluded = soft); new `color_series` counts
> free anthocyanin + polymeric pigment (colour RETAINED — the monomeric→polymeric shift, not vanishing; bleaching
> deferred). ACETALDEHYDE BRIDGE DEFERRED — the advisor caught it as a conservation trap (acetaldehyde is ON the carbon
> ledger; an off-ledger pigment consuming it fails `assert_conserved`) — it's the explicit named next beat. 1 advisor
> pass (before writing; adjusted 2 leanings, 1 conservation-breaking). Next: the acetaldehyde-bridged route, tannin
> self-polymerization, SO₂/pH bleaching (promotes pigment to a slot), beat 1b, Maillard Strecker, or barrel fill-number.**
>
> **PROGRESS (D-80, 2026-07-13): `AcetaldehydeBridgedCondensation` built — the acetaldehyde-bridged (ethylidene) route,
> the NINTH aging Process, the THIRD non-oxidative one, the SPLIT-LEDGER beat, and the first link from the oxidative
> sub-axis to red-wine colour (887 tests). Trilinear `k·f(T)·[free acetaldehyde]·[anthocyanin]·[tannin]`; a new ON-ledger
> `ethyl_bridge` slot captures the acetaldehyde carbon (weighted at cf(ethylidene)) so carbon closes non-trivially.
> HONEST-FRAMING: v1 delivered the carbon + O₂→pigment MECHANISM, NOT a colour behaviour change — `color_series` was
> O₂-invariant (superseded at D-81). 1 pre-work + 1 done-call advisor pass.**
>
> **PROGRESS (D-81, 2026-07-13): the SO₂/pH anthocyanin-BLEACHING beat — `AnthocyaninFading` built + polymeric pigment
> PROMOTED to a slot, so `color_series` now GENUINELY DECLINES (899 tests, two commits). User chose "Both (C)" at the
> design fork = reversible masking readout (A) + irreversible fade sink (B), a two-beat split; D-81 delivers B, D-82 (the
> masking readout) is the COMMITTED second half, still owed. The fade is O₂-COUPLED (bilinear `k·f(T)·[o2]·[anthocyanin]`
> on the shared o2 pool, anthocyanin → colourless `faded_anthocyanin`), so SO₂ colour-protection is EMERGENT (SO₂
> scavenges o2 via D-72, leaving less to fade) — advisor caught that the reversible-masking alternative would make colour
> RISE, not decline. HONEST O₂-GATING: colour declines UNDER O₂ exposure; an anaerobic red still holds flat via
> condensation. Three-slot identity anthocyanin + polymeric_pigment + faded_anthocyanin ≡ anthocyanin₀ (by construction).
> SUPERSEDES D-80's O₂-invariance framing. 1 pre-work fork-resolution + 1 done-call advisor pass. Next: D-82 (masking
> readout, owed), O₂-independent thermal fade, tannin self-polymerization, beat 1b, Maillard Strecker, barrel fill-number.**
>
> **PROGRESS (D-82, 2026-07-13): the masking readout — `observed_color_series` built, DELIVERING D-81's COMMITTED beat A
> (906 tests). A pure READOUT (no state slot, no fate): the Somers reversible SO₂/pH bleaching. Free monomeric anthocyanin
> is masked by a coloured fraction `χ = 1/(1 + K_h/h + K·[HSO₃⁻])` (`acidbase.anthocyanin_coloured_fraction`), polymeric
> pigment counted FULL (SO₂/pH-resistant); `observed = χ·anthocyanin·1000 + pigment·1000`. Advisor's load-bearing catch:
> the COMPETITIVE single denominator, NOT a product of two fractions — carbinol + bisulfite adduct are parallel drains of
> the flavylium pool, so a product form carries a spurious cross-term (bisulfite bleaching the colourless carbinol, ~4×
> error at pH 3.4/20 mg/L). Reads FREE bisulfite (`bisulfite_so2_at_ph`, after carbonyl binding), so reversibility is
> emergent (SO₂ bound/oxidised → mask lifts). OPPOSITE SO₂-sign to D-81 (here SO₂ MASKS; D-81 SO₂ PROTECTS) — different
> series, both real, comment cross-refs guard it. Two params (`pKa_flavylium_hydration` plausible, `K_anthocyanin_bisulfite`
> speculative). 1 pre-work + 1 done-call advisor pass. THE COLOUR AXIS'S "Both (C)" REQUEST IS NOW COMPLETE (B=D-81 fade +
> A=D-82 mask). Next: O₂-independent thermal fade, tannin self-polymerization / tannin-ethyl-tannin, beat 1b, Maillard
> Strecker, barrel fill-number, barrel-beer oak.**
>
> **PROGRESS (D-83/D-84/D-85, 2026-07-13): a 3-Process batch closing two deferred aging axes (942 tests). ONE pre-work
> advisor pass designed the batch and flagged the sole conservation risk; the off-ledger beats landed first, the
> ledger-touching one last. (1) `ThermalAnthocyaninFade` (D-83, 11th aging Process) — the O₂-INDEPENDENT thermal fade:
> `r = k·f(T)·[anthocyanin]` (first-order, no o2, no yield) → the same colourless `faded_anthocyanin` slot. The MIRROR of
> D-81: touching no o2, SO₂ gives NO protection, so a sealed/sulfited/anaerobic red now fades (only cold storage slows it) —
> RETIRES D-81's "anaerobic sealed red holds its colour" note (3 existing tests updated). (2) `TanninSelfPolymerization`
> (D-84, 12th) — the DIRECT tannin–tannin softener: `r = k·f(T)·[tannin]²` (bimolecular self-reaction, pure off-ledger sink,
> no destination slot per the D-79/D-80 tannin-is-a-pure-sink precedent). Softens astringency WITHOUT anthocyanin (a white
> now softens) → RETIRES the `astringency_series` "one-directional-per-pool" note. (3) `TanninEthylTanninCondensation`
> (D-85, 13th) — the acetaldehyde-BRIDGED tannin–ethyl–tannin route: `r = k·f(T)·[free acet]·[tannin]²`, the D-84 form + the
> D-80 acetaldehyde factor. Reuses D-80's split-ledger carbon capture (acetaldehyde → shared `ethyl_bridge`, own
> `y_acetaldehyde_per_tannin` — one acetaldehyde per TWO flavanols), so `total_carbon` closes to machine precision
> non-trivially; deposits NO pigment (colourless tannin–tannin polymer — the D-80 colour difference). The tannin–tannin axis
> is now built. Next: beat 1b (descriptor projection), non-oxidative Maillard Strecker, barrel fill-number, barrel-beer
> oak.**
>
> **PROGRESS (D-86, 2026-07-14): barrel-beer oak (945 tests). The oak axis (D-77 aroma four + D-78 ellagitannin), wine-only
> since D-77, is now wired into BEER too — bourbon-barrel stouts, oak-aged/foeder sours, whiskey-barrel ales. ONE advisor
> pass + ONE owner scope fork (FULL axis incl. ellagitannin, not aroma-only). Principle: extraction is a WOOD property
> (the `oak.yaml` yields are matrix-independent, transfer unchanged), only PERCEPTION is matrix-specific (4 new
> `threshold_<compound>_beer`, set below the wine values for beer's lower-ethanol/less-masking matrix). Wine stays
> BYTE-FOR-BYTE: the 10 oak `VarSpec`s factored into `core.media._oak_specs()`, called at the SAME position in `wine_schema`
> + appended in `beer_schema` (the `iso_alpha` precedent — NOT `_common_specs`, which would shift wine indices). Both oak
> Processes (`OakExtraction` + `EllagitanninOxidation`) wired into both media (always medium-agnostic in logic; `o2` already
> shared). `add_oak` needed NO logic change (the `whiskey_lactone`-slot guard auto-relaxed once beer carries the slot);
> `oav.py` moved the 4 oak compounds to a shared `_OAK` tuple; `astringency_series` guarded the wine-only grape `tannin`
> slot (beer astringency = oak ellagitannin alone). The GRAPE colour axis (D-79..D-85) stays wine-only (grape chemistry).
> Three wine-only enumeration/rejection tests FLIPPED to both-media (expectation changes, not weakenings); +beer end-to-end
> coverage incl. the D-78 protection spine on beer (oaked+oxygenated beer browns less). Next: beat 1b, non-oxidative
> Maillard Strecker, barrel fill-number depletion.**
>
> **PROGRESS (D-87, 2026-07-14): non-oxidative THERMAL Strecker route — `MaillardStrecker`, the tenth aging Process, the
> O₂-INDEPENDENT thermal mirror of `StreckerDegradation` (D-75) — the beat D-75 deferred (964 tests). Residual SUGAR + HEAT
> (α-dicarbonyls, NO O₂) degrade amino acids to the sweet-wine/Madeira aldehyde suite, so a SEALED sweet wine ages
> thermally where the O₂-only D-75 route is silent. TWO advisor passes + TWO owner scope forks (FULL scope both): (1) the
> aldehyde suite — added FOUR new wine-only pools (2-/3-methylbutanal, 2-methylpropanal, sotolon) beyond the two shared with
> D-75, not the lean 2-pool v1; (2) build the sugar→melanoidin thermal browning too — split into its own D-88
> (`Caramelization`) per one-Process-per-decision. Crux: `S` is a read-only DRIVER (not consumed here — the aldehyde carbon
> is the AMINO ACID's, drawn from `amino_acids`; booking a sugar draw would break `total_carbon` AND undercount real sugar
> loss, so it is FORCED not convenient). Sotolon is NOT a Strecker aldehyde (a furanone) → NO CO₂ term (the flag is
> load-bearing: closure holds for any CO₂ attribution, so a mis-key passes silently — the D-75 fidelity lesson; levels
> anchored to literature: sotolon ~5–20 µg/L Sauternes). Additive with D-75 over the shared `amino_acids` limiting reagent
> (the o2-sharing pattern on amino acids). Isolability on the `amino_acids` HARD gate (S is a soft driver). E_a = 100 kJ/mol
> (> the oxidative ~50: the sourced ORDERING — Maillard out-accelerates oxidation with temperature). New `thermal.yaml`.
> Next: D-88 (`Caramelization`, the first aging Process to consume core `S` — carries the `begin_aging` golden re-baseline),
> beat 1b, barrel fill-number.**
>
> **PROGRESS (D-88, 2026-07-14): non-oxidative sugar-only THERMAL browning — `Caramelization`, the eleventh aging Process,
> the O₂-INDEPENDENT thermal mirror of `PhenolicBrowning` (D-74) and the browning half of the D-87 thermal axis (978 tests).
> Residual SUGAR browns to melanoidin by HEAT (no O₂), raising the SAME A420 index D-74 accumulates — so a sealed sweet wine
> still darkens (Sauternes amber, Madeira/baked deep colour). The FIRST aging Process to consume core `S`: its carbon lands
> in a new on-ledger `melanoidin` carbon-park (caramelan stand-in C12H18O9, the debris/glucan precedent), so `total_carbon`
> closes exactly (release at the sugar fraction, redeposit at melanoidin's — the EsterHydrolysis split). CARAMELIZATION not
> Maillard (sugar-only, nitrogen-free — N-incorporating melanoidin deferred). GOLDEN AUDIT: every standard aging run
> ferments to dryness (S ≈ 0) before `begin_aging`, so it is byte-for-byte inert on them (the S ≤ 0 guard) — churn confined
> to sweet-wine runs. D-83-style supersession: "un-oxygenated aging = byte-for-byte ester-only" now holds only for DRY wines
> (a sealed sweet wine browns thermally). Wine-only v1 (a bundling choice, not physics — beer thermal browning deferred, the
> D-86 pattern). E_a = 100 kJ/mol (matches D-87 — the thermal axis). The non-oxidative thermal axis (D-87 aroma + D-88
> browning) is now complete. Next: beat 1b (descriptor projection), barrel fill-number, N-incorporating Maillard melanoidin
> / beer thermal browning.**
>
> **PROGRESS (D-89, 2026-07-14): amino-acid-incorporating THERMAL browning — `MaillardBrowning`, the twelfth aging Process,
> the N-incorporating browning branch D-88 `Caramelization` explicitly deferred (992 tests). Where D-88 browns SUGAR ALONE to
> nitrogen-free caramelan, TRUE Maillard condenses residual SUGAR + AMINO ACID → a brown polymer that RETAINS the amino-acid
> nitrogen (what makes a melanoidin *nitrogenous*), raising the SAME A420. Consumes both core `S` and `amino_acids` with no
> O₂ into a new on-ledger N-bearing `maillard_melanoidin` pool (glucose–glycine stand-in C8H12O5N, molar C:N ≈ 8:1) — the
> FIRST aging Process on the nitrogen ledger (first non-biomass, non-arginine species on `total_nitrogen`). ONE advisor pass,
> owner scope fork (FULL build, N-fate = "closest to reality" = ALL-N-retained, because D-87 already owns the deaminating
> branch). Dual carbon+nitrogen closure by SIZING both draws to the melanoidin formed (`r_m = r_sugar·c(sugar)/(c_m −
> n_m·c(arg)/n(arg))`); the denominator sign is the silent trap (>0 or it creates sugar with no test catching it — a metadata
> test pins it). The 3 thermal branches now split cleanly (sugar-only D-88 / N-retaining D-89 / N-releasing volatile D-87),
> summed over shared `S`/`amino_acids` by ProcessSet. `k_maillard_browning` calibrated (5e-8, nitrogen-limited/minor) so the
> shared-`amino_acids` competition doesn't erase the diagnostic sotolon — the D-74-suppresses-acetaldehyde precedent on
> amino_acids. Wine-only v1. The four-way interaction test became FIVE-WAY. Next: beat 1b (descriptor projection), barrel
> fill-number, beer thermal browning (the D-86 oak-to-beer pattern for the whole thermal axis).**
>
> **PROGRESS (D-90, 2026-07-14): beer thermal browning — `Caramelization` (D-88) extended to BEER, the D-86 oak-to-beer
> pattern (993 tests). A warm-stored / long-aged beer's residual dextrins (unfermented maltose/maltotriose) caramelize to
> `melanoidin` by HEAT with no O₂, raising the SAME shared `A420` browning index — the sealed-sweet-wine mechanism, now on
> beer. Scope is FORCED, not chosen: of the 3 thermal branches only D-88 is sugar-only; D-87/D-89 read `amino_acids`, untracked
> in beer (D-32), so they genuinely can't follow — beer thermal browning is caramelization ALONE (the inverse of D-86: physics
> is medium-agnostic, but reagent TRACKING is the wine-only wall). ONE advisor pass, no owner fork. The one real correctness
> pin (advisor's must-fix): PER-COMPONENT clamp (`s_clamped = y[S].clip(min=0)`), NOT `max(sum,0)` — on beer's 3-slot `S` a
> solver undershoot leaving one component negative would flip the apportioned debit POSITIVE and silently CREATE sugar (carbon
> closes either sign — the D-89-denominator trap family). Vectorized carbon-exact transfer: each sugar caramelizes at its OWN
> carbon fraction (glucose C6 / maltose C12 / maltotriose C18 differ), `carbon_released = r·Σ(sᵢ/s_total)·c(sugarᵢ)`. Wine
> reduces to the single-slot D-88 form BYTE-FOR-BYTE (every wine trajectory unchanged); `thermal.yaml` params transfer
> unchanged (`A420` is an absorbance, not a matrix-specific threshold). Isolability asymmetry (honest, documented): a beer
> finishes at S ≈ 5e-11 (positive, guard doesn't fire) not wine's exact ≤ 0, so a reductive beer browns a NEGLIGIBLE trace
> (A420 ≈ 4e-8 vs 0.27 O₂-dosed) — NOT byte-for-byte inert; the "reductive aging = ester-only" claim now holds only for truly
> DRY beverages. Load-bearing new test: a residual wort browns with per-component carbon closure to 1e-18; integrated
> high-residual big-stout beer closes `total_carbon` across the multi-slot vector over an aging year. Wiring: `melanoidin`
> appended to `beer_schema`, `_CARAMELIZATION_PROCESSES` in beer, gating free (name-guarded `_AGING_GATED_PROCESSES`). The
> thermal axis is now media-split: sugar-only browning BOTH media (D-88/D-90), N-routes wine-only (D-87/D-89). Next: beat 1b
> (descriptor projection), barrel fill-number, a beer-specific per-melanoidin A420 yield (refinement).**
>
> **PROGRESS (D-91, 2026-07-14): barrel fill-number depletion — a reused barrel extracts LESS, an `add_oak` DOSE input
> (997 tests). The long-deferred fill-number beat (forward-noted since D-77/D-86): a barrel is a DEPLETING oak source, so a
> 2nd-/3rd-/4th-fill barrel imparts progressively less wood character than a fresh first-fill one — the signature lever of
> barrel-aged BEER programs (first-fill bourbon barrel → imperial stout, then the neutralised barrel → a sour). Design fork
> (advisor): an ACROSS-FILL DOSE INPUT (`fill_number` scales the saturation ceiling AT DOSE TIME, barrel history known when
> the oak is charged), NOT a within-fill dynamic RESERVOIR `OakExtraction` draws down (the mechanistic model — new state /
> Process / conservation-adjacent — deferred); matches D-77's "the ceiling is set at the dose". PURELY DOSE-LEVEL: NO new
> Process, state slot, or schema change. `add_oak` gains an OPTIONAL `fill_number` (int ≥ 1, default 1); each of the 5
> ceilings scaled by `oak_fill_retention ** (fill_number − 1)` before the `+=` write; `OakExtraction` untouched.
> `fill_number = 1` (fresh, default) → `r**0 = 1.0` EXACTLY ⇒ every pre-D-91 wine + beer trajectory BYTE-FOR-BYTE unchanged
> (the whole existing oak suite already pins the first-fill case). New `oak_fill_retention = 0.5` (speculative, banded
> 0.3–0.7): sourced ORDERING = barrels go effectively NEUTRAL by ~4th–5th fill; one shared retention across all 5 extractives
> (per-compound retention deferred). Off every ledger (wood-derived, iso_alpha precedent). Validated int-valued ≥ 1 (zeroth /
> fractional fill rejected loudly at compile — brewers count first/second/third). SCOPE: oak-EXTRACTABLE depletion only — a
> first-fill ex-bourbon barrel's residual-SPIRIT soak-back (leached from the spirit, not the wood) is SEPARATE / deferred.
> +4 tests (byte-for-byte first fill, geometric 1:r:r³ discount across all 5 extractives, the motivating BEER end-to-end
> first-vs-fourth-fill stout reads lower oak OAVs + astringency, validation). ONE advisor pass (green-lit + sharpened six
> points; a second done-call caught stale "deferred fill-number" docstrings, fixed). Next: beat 1b (descriptor projection),
> a beer-specific per-melanoidin A420 yield, the deferred finite-reservoir / per-compound-retention refinements,
> bourbon-barrel spirit soak-back.**
>
> **PROGRESS (D-92/D-93/D-94, 2026-07-14): the BOURBON-BARREL SOAK-BACK trio — an ex-spirit barrel donates more than wood
> (1009 tests). D-91 scoped itself to oak-*extractable* depletion and flagged the residual-SPIRIT soak-back as separate; these
> three build it. All are `add_oak` DOSE-level (an optional categorical `spirit`, v1 `bourbon`) — NO new Process, state slot,
> or schema change until D-94's pool. `spirit` absent ⇒ byte-for-byte unchanged. (1) **D-92 — ETHANOL** (~1% ABV gain, the
> "barrel-aged stout gains ABV" effect): the advisor's BLOCKER was "is a carbon-bearing dose free?" — traced to
> `schedule.py:231`, the scheduler books every dose's `new_y − current_y` as an `ExternalFlow` automatically, so a DISCRETE
> dose closes the run-wide ledger where a within-Process CONTINUOUS leach would CREATE carbon within-segment. That is *why*
> it is a dose, not a Process. Own steeper `spirit_soak_retention=0.2` (first-fill is the term of art), decoupled from
> `oak_gpl`. (2) **D-93 — AROMA CONGENERS** (vanillin/whiskey_lactone/guaiacol): a CEILING BUMP drawn in GRADUALLY by
> `OakExtraction`, NOT a bolus — the load-bearing proof is that a bolus into the pool is ERASED by the extraction gate
> (`gap = ceiling − conc` ⇒ `max(X, C_wood)`, not the sum), so bumping the ceiling is the ONLY additive-with-wood form. The
> D-92 asymmetry (ethanol=bolus, aroma=gradual) is a STRENGTH: the LEDGER splits them — ethanol is ON it (forced to a dose),
> aroma ceilings are OFF it (gradual is both available AND more faithful). (3) **D-94 — CARAMEL**: `furaneol` (HDMF), a FIFTH
> oak extractive, both a toast-RISING wood yield and a bourbon ceiling bump. The feared collision with D-88's
> caramelization/`A420` axis is DISSOLVED, not relocated: `melanoidin` is caramelization's on-ledger COLOUR body, `furaneol`
> the off-ledger volatile AROMA of the same chemistry — so it cannot perturb D-88's carbon closure (an executable test pins
> it). Done-call catch: SIZE BY OAV BAND, not mass — furaneol's potency (threshold ~30× below vanillin's) meant a
> mass-matched bump read caramel ~7× more forward than vanilla, invisible to a bare `OAV>1` test. The genuinely deferred beat
> is caramel aroma from the beverage's OWN thermal caramelization (that WOULD be on-ledger).**
>
> **PROGRESS (D-95, 2026-07-15): beat 1b slice 1 — DESCRIPTOR-SPACE PROJECTION, the last unbuilt piece of this milestone's
> OPENING beat, deferred at D-66 and carried as a standing "Next:" in every entry from D-67 through D-94 while the aging axis
> grew 13 Processes underneath it (1027 tests, +18). New `sensory/descriptors.py`: wine's 19 / beer's 10 aroma pools project
> onto 14 / 9 descriptor axes. NO state, Process, ledger entry, YAML, or parameters. THE SLICE LINE IS THE ADDITIVITY SEAM
> (the advisor's crux): projection is inherently many-to-many (`malty` ← three aldehydes; `smoky` ← oak guaiacol + Brett 4-EG),
> but the layer BELOW already refused to aggregate — `SensoryProfile` never sums OAVs because additivity is contested (D-67).
> So the rule is NOT a free choice: summing would silently reintroduce what the layer beneath rejected. Hence the MAX rule —
> each descriptor reads its LOUDEST contributor and names it (`dominant`). *We never assume additivity, at any layer.* Slice 2
> (deferred, D-96) = weights/compression/masking + the params those need. HONEST FRAMING (the D-80 precedent): under max a
> descriptor clears iff one of its pools does, so `above_threshold()` is a pure REGROUPING of beat 1a's flags — slice 1 adds
> vocabulary + attribution, NOT new above-threshold information. Membership is STRUCTURE not parameters (binary ⇒ no weights
> ⇒ no YAML — weights are precisely what makes slice 2 need a provenance file); the axis set is DERIVED per medium, so beer
> can never report `barnyard` by construction. `DescriptorProjector` Protocol = the §4.2 swappable seam, PROVEN by a test that
> swaps in an alternative rule. Owner forks: ~12 many-to-many axes (over ~7 coarse / ~18 near-1:1); max-rule-v1. Done-call
> advisor catches: (1) my "D-94 FORCES the caramel/curry split" justification was internally inconsistent — D-94 governed the
> COMPOUND layer, while collapsing distinct compounds is this layer's JOB, and the same vocabulary merges guaiacol + 4-EG
> whose distinctness is flagged just as loudly; the split stands as a JUDGEMENT, recorded as such. (2) The `lumped` honesty
> flag was being DROPPED at the layer boundary — it now propagates from the dominant contributor (`sulfidic` = clean h2s +
> lumped mercaptans is the live case). Next: beat 1b SLICE 2 (the perceptual speculation), a beer-specific per-melanoidin
> A420 yield, the on-ledger thermal-caramelization aroma co-product.**
>
> **PROGRESS (D-97, 2026-07-15): the ATF1 PRECURSOR COUPLING — the banana ester becomes YAN-responsive (1039 tests, +8). The
> first of D-96's four named deferred refinements, and the only one with a NEW OBSERVABLE. `isoamyl_acetate` synthesis is now
> FIRST-ORDER in its precursor alcohol, the `fusels` pool. Before D-97 the banana was YAN-BLIND — flat at 0.759/0.758/0.756
> mg/L across YAN 40/80/250 while the `fusels` pool that IS its precursor swung 2.9x. Not a coarse reading: a MISSING
> DEPENDENCY. THE ADVISOR REVERSED ITS OWN STEER AND A PROBE IS WHY: its first crux was "gate on the N-gated DRIVER, not the
> accumulating pool" (structurally airtight — a monotone pool's gate is flat or back-loading), but `N` empties at DAY 2 with
> ~75% of the sugar unfermented, so a driver-coupled banana STOPS DEAD at day 2 while 51 mg/L of substrate sits there and the
> flux still supplies acetyl-CoA. Coupling to the RATE conflates "the precursor is being made" with "the ester is being made";
> the alcohol PERSISTS. Its no-op objection also dissolved: it reasoned about flatness IN TIME, but the discriminating signal
> is ACROSS-YAN. Lesson: a five-line probe settled in seconds what neither of us could settle by reasoning. THE FORM WAS A
> LOOKUP, NOT A JUDGEMENT: Fujii 1998 (AEM 64:4076-4078) gives ATF1's Km for isoamyl alcohol ~29.8 mM AND states the mechanism
> outright ("a major rate-limiting factor ... is the amount of isoamyl alcohol") — one paper sources BOTH that the coupling
> exists and what form it takes. The pool runs ~51x (wine) / ~63x (beer) BELOW Km ⇒ the [S]<<Km limit ⇒ linear is the MEASURED
> regime, not the convenient choice. WHY THE MEASURED Km IS NOT A PARAMETER (the advisor's dual trap — "omitting a Km that
> matters is the same sin as fitting one; the discipline is SOURCE THE FACT, not fewer params"): in this limit only the RATIO
> Vmax/Km is IDENTIFIABLE — scale Km 10x, refit Vmax, byte-identical trajectory ⇒ it would be a sourced-looking constant no
> output could validate. `k_isoamyl_acetate` IS that ratio; the Km's honest home is the PROVENANCE, justifying the FORM. THE
> ASYMMETRY IS DERIVED: ethyl acetate is ATF1 too, but its ethanol precursor (~2 M) SATURATES the enzyme ⇒ zeroth-order ⇒ no
> term. One enzyme, one rate law, two limits, decided by concentrations. READ NEVER DEBITED: carbon still from `S`, so
> `touches`/`total_carbon` untouched; the C5+C2=C7 inverse of D-69's 5:2 is the tempting elegance, DEFERRED and pinned
> meanwhile. OUTCOME PINNED AS A RATIO (D-96's lesson applied FIRST, not caught at the end): `k` is re-anchored so the finished
> value is unchanged ⇒ any absolute-band test would pass on the OLD model; the teeth assert the ester swing TRACKS the fusel
> swing on a real run. Verified: unwiring `precursor_pool` fails it + 3 others. NO NEW PARAMETER; `k` 5.0e-6→1.0e-4 (wine,
> x20) / 1.2e-5→3.05e-4 (beer, x25.4) — the factors DIFFER because each medium's fusel level does. EMERGENT, untuned: the
> sensory payoff lands in BEER, not wine — under D-95's MAX rule wine's `fruity` is apple-dominated (no YAN dependence) so it
> barely moves (79.0→78.6) while the banana's own OAV swings 14.7→42.1 (real but MASKED, honestly so); in beer the fruity
> DOMINANT LABEL FLIPS with nitrogen (YAN 100 ⇒ ethyl_hexanoate, 200+ ⇒ isoamyl_acetate) — a low-nitrogen wort's fruity stops
> being banana-led. Next: the isoamyl-acetate carbon re-route (the D-69 5:2 inverse); per-ester dH/E_a (BLOCKED on sourcing —
> an author estimate would LOWER fidelity, so it may stay deferred rather than be built for completeness); speciating the
> `fusels` lump (the D-96 pattern one pool over, which would retire this beat's inherited lump caveat); beat 1b slice 2.**
>
> **PROGRESS (D-99, 2026-07-16): the lumped `fusels` pool SPLIT into FIVE single-molecule higher alcohols (1070 tests, +31;
> benchmarks 16/16; ruff+mypy clean). The second of D-96's four deferred refinements, and the D-96 ester split one pool over:
> `propanol`/`isobutanol`/`active_amyl_alcohol`/`isoamyl_alcohol`/`2_phenylethanol`, each carbon-weighted by its own molecule,
> each produced by one shared Ehrlich shape times its OWN independently-anchored `k` (Wang/Frank/Steinhaus 2024 meta-analytic
> means; n=486-684 wine studies per compound). Unlike D-95/D-98 it is ON BY DEFAULT — chemistry, not a sensory projection.
> THE OWNER OVERRODE MY "four species" AND WAS RIGHT (3rd time a "can't source it" call was overruled by looking): active amyl
> has its OWN row over n=128 wine / n=64 beer studies — the "amyl alcohols" bundling is a GC COLUMN ARTIFACT (coelution), split
> in aroma research. MY RECOMMENDED CUT WAS AIMED AT THE WRONG MOLECULE: propanol, not active amyl, is the weak one (omitted
> from the meta-analysis in both media; beer propanol is the beat's one author estimate, flagged in-file). THE ~3.8x RISE
> (wine 86→328 mg/L) IS FORCED, NOT CHOSEN: the old lump sat below even the sum of the five species' LOW ends; each `k` was set
> by solving for a literature mean recorded BEFORE any k was picked (wine hit 327.9 vs 327.8 target first try). PAYOFF IS
> SUBTRACTIVE FIRST: the lump asserted EVERY higher alcohol smells like isoamyl (only ~52% is, by mass), reading 2-phenylethanol
> — which is ROSE, not solvent — at the wrong potency; the split removes that false claim even for the two chemistry-only pools
> (propanol/active amyl: no usable-matrix threshold ⇒ no OAV ⇒ honest silence). THEN ADDITIVE: 2-phenylethanol ~28.7 mg/L vs
> Guth's ~10 ⇒ OAV ~2.9 ⇒ WINE GAINS A FLORAL AXIS that read exactly 0 before (its only member was the oxidation-only D-75
> phenylacetaldehyde). Beer correctly gains none (no beer threshold; Meilgaard says 2-PE is sub-threshold in beer) — the mirror
> of D-97. THE D-96 RULE AND THE ONLY PAYOFF ARE THE SAME ACT — the INVERSE of D-98's "sourceable XOR consequential": a probe
> showed ratio-splitting the current lump (D-96-forbidden) leaves 2-PE at OAV 0.93 and delivers nothing, while honest anchoring
> (D-96-required) lights the axis. D-98 PREDICTED BOTH SIDE-EFFECTS VERBATIM: its caveat (iv) lumped exponent is retired
> (single-molecule now; 2 new aliphatic-alcohol exponents ordered by Cain's MEASURED chain-length trend but kept author-estimated
> with OVERLAPPING bands, so D-98's tripwire stays green), and the ATF1 Km comparison now reads 3-methylbutan-1-ol specifically.
> D-97's RATIO LESSON APPLIED PROSPECTIVELY: speciation doubled the banana's precursor, so `k_isoamyl_acetate` was re-anchored by
> an EXACT factor (x0.502 wine / x1.397 beer — opposite directions, because wine's isoamyl pool doubled while beer's FELL) to
> hold the finished ester and stop `fruity` flipping on no evidence. I NEARLY SHIPPED A PROVENANCE NOTE CLAIMING the re-anchor
> while the code didn't do it — caught by checking the descriptor outcome against D-97's recorded value. THE DONE-CALL CATCH
> BECAME THE SECOND FINDING: the 3.8x rise broke THREE tests in unrelated subsystems (Maillard, Brett growth, MLF growth) — all
> via the D-33 reroute draining the lumped `amino_acids` pool to ~0. A pre-D-99 emulation proved D-99 didn't CREATE this — the
> reroute already ate ~96.5% of the pool and sotolon passed at an 18% margin the old lump's under-production propped up. TWO
> SPECIATED-SCALE CONSUMERS CAN'T SHARE ONE LUMPED SUBSTRATE ⇒ `amino_acids` is THE NEXT LUMP ⇒ D-100. Resolved in-scope: the
> speculative reroute is ISOLATED out of the three affected tests (each tests a different Process; the reroute stays fully
> covered in test_fusel_reroute.py) and the pathology is PINNED as a D-100 tripwire. `mercaptans` is now the LAST lump in the
> project. HONEST LIMIT STATED: all five share one N-gate + one E_a ⇒ fixed SPECTRUM (not fixed composition); dynamic spectrum
> needs per-species E_a/gates = unsourced = the real D-98 trap = deferred. Next: D-100 (speciate `amino_acids`); the isoamyl-
> acetate carbon re-route (D-69 5:2 inverse); per-species E_a/gates (blocked); sourced beer propanol + in-matrix beer thresholds
> (paywalled); per-ester dH/E_a (blocked); masking (blocked on cosα); further fruity esters; the oav→magnitude rename.**
> Milestone 1 (Tier-1 validated core) and Milestone 2 (Tier-2
> plausible mechanisms) are closed — the §2.2 benchmark trio is green and §3.3
> "additives with clear mechanisms" completed at D-65 (717 tests). This plan opens
> **Tier-3**, the handoff's §4 frontier: "real chemistry, but integrating it into a
> trustworthy prediction is *not solved science*." Everything here is `speculative`,
> isolated, and clearly labelled — it must never perturb the validated core or its
> tests (prime directive #3).
>
> The two opening calls are recorded in **DECISIONS D-66**: (1) build the
> **sensory/OAV readout layer first**, aging chemistry second — *inverting* the
> handoff §6-step-5 order ("aging then sensory"); (2) handle the lumped aroma pools
> with a **representative-compound threshold per lump** (owner call, over the
> single-compound-only alternative).

> **PROGRESS (D-100, 2026-07-16): the lumped `amino_acids` pool SPLIT into EIGHT single-molecule pools — the LAST shared
> substrate (1094 tests, +24; benchmarks 16/16; ruff+mypy clean). The lump D-99 exposed, resolved: six PRECURSORS
> (leucine/isoleucine/valine/threonine/phenylalanine/methionine — each the molecule whose skeleton actually becomes a product)
> plus two IDENTITY-AGNOSTIC pools (`amino_acids`, kept as the arginine slot, + `amino_acids_generic` as glutamine). Ten
> consumers rewired. THE DECOUPLING: the D-33 reroute no longer touches `amino_acids` AT ALL — precursor and identity-agnostic
> consumers address DISJOINT pools, so fusel production CANNOT starve MLF/Brett/yeast-growth (structurally impossible, not
> tuned away; arginine ends 0.0101 with the reroute ON vs 0.0132 OFF, was ~1e-5 vs ~0.2). THE GATE RULE, one uniform rule with
> ZERO new parameters: each consumer gates on its own substrate with its existing K scaled by that substrate's must-spectrum
> share (`K_i = K·f_i`) — derived from `K_amino_acids` + the sourced spectrum, NOT a Michaelis constant and labelled as not
> being one (per-species Michaelis constants = eight unsourced numbers = the D-98 trap, declined). It REDUCES: at spectrum
> composition every gate equals the pre-split lumped gate, which is WHY every closed-form suite still asserts its OLD numbers
> — nothing was re-baselined. **The advisor caught me conflating two claims inside that property**: subset-agreement is
> STRUCTURAL, but the match to the old VALUE is contingent on Σf = 1 (verified by hand at F=1.05: subsets still agree, lumped
> match breaks) ⇒ ensembles sampling the bands shift the baseline slightly. TWO LUMPS RETIRED, not restated: D-33's ~4×
> nitrogen over-release (4-N arginine deaminated for a 1-N job) and D-45's mercaptan drawn from ARGININE — a molecule with NO
> SULFUR. THREE EMERGENT FINDINGS: (i) the anabolic/catabolic split is now EMERGENT — leucine supplies only ~7% of the isoamyl
> carbon and its own gate throttles the rest onto sugar, so **D-19's apologised-for sugar stand-in is now the CORRECT book for
> de-novo synthesis**; (ii) aging precursors are dominantly AUTOLYSIS-sourced (phenylalanine 1.6e-5 → 44.8 mg/L with lees,
> phenylacetaldehyde 0.016 → 1631 µg/L, restoring the literature phenyl-dominant ordering) — the published sur-lie mechanism
> as a CONSEQUENCE, not a rule; (iii) `f_methional` was a lump for precursor ABUNDANCE (its own comment says so) — NOT
> re-anchored here, that would be the D-99 anti-tuning violation. **THE MAGNITUDE OF (ii) IS NOT VALIDATED and is recorded as
> such** (advisor's call, the one that matters at this project's bar): the reroute's catabolic fraction is a ~0.5 lumped guess
> vs a literature 20-50%, so the model drains precursors HARDER than reality and "silent without lees" is what that fraction
> gives, not a prediction — the `<1e-9` bound is a TRIPWIRE on current behaviour, and moving it when the fraction is bounded
> is the SIGNAL, not a regression. **[D-103 CORRECTS THIS SENTENCE: both numbers were wrong. Measured exactly the fraction is
> 0.192 at this dose (never ~0.5); the "20-50%" was UNCITED and the sourced band is LOWER (Rollero 2017: isoamyl 2-8%,
> isobutanol 5-15%) ⇒ 20-50% would have ACQUITTED the model. The conclusion "drains harder than reality" is TRUE but every
> number under it was wrong, and the defect is the gate's per-species SHAPE, not a fraction to bound ⇒ NO knob, and the
> `<1e-9` tripwire correctly does NOT move. See D-103.]**
> **[D-104 CORRECTS TWO MORE OF THESE CLAIMS. (a) Finding (i) — "the anabolic/catabolic split is now EMERGENT" — is RETIRED
> for the protein split: measured against Crépin 2017 the emergent (demand-driven) split is MONOTONICALLY INVERTED (model
> leucine 20.9% → protein vs a measured 77–86%; model order leu<ile<val<thr vs Crépin's thr<val<ile<leu, exactly reversed),
> and no biomass composition repairs it — the model's Ehrlich demand pulls the AMINO ACID where reality pulls a mostly-de-novo
> KETO ACID. An emergent WRONG answer is worse than a sourced static one, so the split is now sourced and STATIC. What
> survives of (i): each alcohol's catabolic FRACTION is still emergent, being (1−f) × a supply/demand ratio. (b) "THE TRIPWIRE
> FLIPPED: sotolon OAV 0 → 3.22" was NEVER A FIX — sotolon rode on ~6.6 mg/L of leftover threonine that survived only because
> propanol's demand was too small to consume it. D-104's anabolic sink consumes it properly and sotolon died AGAIN (its third
> death, always through threonine) ⇒ re-based on 2-ketobutyrate (`de_novo=True`), which is the actual fix. See D-104.]**
>
> **[D-107 SUPERSEDES THAT LAST CLAIM — `de_novo=True` was not "the actual fix", it was the right diagnosis wearing the
> only vocabulary available. Sotolon is not a Strecker product AT ALL: it is a purely chemical aldol of
> α-ketobutyrate + acetaldehyde (Pham 1995), so it did not belong in `MaillardStrecker`, and its TWO exception flags
> (`decarboxylates=False`, `de_novo=True`) were the model saying so. D-107 built the **2-ketobutyrate pool** — the
> keto-acid NODE, whose producer (the D-45 mercaptan, discarding its C4 and under-drawing methionine 5×) and consumer
> (sotolon, inventing it from sugar and over-drawing threonine 1.5×) had BOTH been sitting in the tree since D-105 —
> moved sotolon to `SotolonAldolCondensation`, and **both flags then had one value across five rows and were deleted**.
> `_KNOWN_NON_STOICHIOMETRIC` is EMPTY. The threonine-vs-sotolon "REAL chemistry over one molecule" competition asserted
> below is now **GONE rather than wrong**: threonine is sotolon's *grand*parent, and the genuine competition is over
> α-ketobutyrate (propanol IS 2-KB decarboxylated) — inexpressible until the fusels are re-based on their keto acids,
> which is the rest of this milestone. See D-107.]**
>
> **[D-109 CONFIRMS THAT "GONE" AND REJECTS THE FIX IT WAS READ AS PRESCRIBING — the fusel side, measured before
> building. (1) **"GONE" is CORRECT, and D-109 opened intending to correct it.** `AlphaKetobutyrateExcretion`'s rate is
> flux-only; threonine's gate moves the **carbon source**, never the **rate** — `d(α-KB)/dt` is **bit-identical** at
> threonine 67 mg/L and at 0. So propanol cannot starve sotolon through threonine, *structurally and on purpose* (gating
> the rate on threonine would kill sotolon in a threonine-free wine — the D-104 canary). A probe measuring sotolon moving
> 0.42% with `k_propanol` was read as "the competition is present but small"; it is the **sugar ledger**, and **only the
> mutation test caught it** — deleting the threonine draw left the assertion passing. *"The sentence and the assertion are
> not the same claim"*, fifth occurrence, first that would have shipped as a correction to something **true**.
> (2) **"Re-based on their keto acids" must NOT be read as the excreted pool.** Sotolon's α-KB is the *extracellular
> residual* (D-49's test, the reason D-107 chose it); propanol is made intracellularly by living yeast mid-ferment = D-49's
> *flux intermediate* ⇒ **the same test rejects that pool for propanol**. It is also **infeasible**: propanol's molar demand
> is **2.79× the total α-KB the pool ever excretes** (42.6× the residual) ⇒ it would starve propanol *and* collapse
> sotolon's substrate — while looking exactly like the promised competition arriving. It would also zero a core both-media
> fusel when the keto-acids toggle off (PD#3).
> (3) **The payoff is REAL — the same 2.79× says so**: propanol's 2-KB demand is **~2.8× the entire excretion flux**, so
> an honestly partitioned node *would* couple propanol and sotolon materially. **Relocated, not dissolved.** (The RATIO,
> deliberately not "propanol is the node's dominant sink": 2-KB's *committed* route is isoleucine biosynthesis via **KMV**
> — named in this very scoping list — which the model does not carry, so propanol's share of TOTAL 2-KB synthesis is
> unmeasured. Both propanol and excretion are overflow off that route.)
> (4) **THE SCOPING RESULT: the fusel-side node is a FLUX PARTITION, not a pool — D-49's own physics.** The intracellular
> keto acid is a vanishing pool carrying an enormous flux ⇒ quasi-steady ⇒ `synthesis == Σ consumption` ⇒ a **partition**,
> not a state variable. **Two nodes, differing in kind**: the excreted keto-acids are pools (they persist, they bind SO₂);
> the fusel node is a partition of the **sourcing layer** (`FuselAminoAcidReroute`/`PrecursorNonEhrlichFates`) —
> **never the producer**, which is what keeps `FuselAlcoholsEhrlich` byte-for-byte when the keto-acids toggle off.
> **It is a MILESTONE, not a beat**: the real work is breaking `FUSEL_SPECS`' one-alcohol-one-precursor assumption
> (isoamyl ← leucine **and** valine: Crépin's 23% via KIC = D-104's named missing route); the prize is D-104's inverted
> split (**open** — D-104 measured that near-equilibrium transamination does *not* un-invert); the **parsimony question is
> open both ways** (per-species rates ≈ the five `f_i`, but **BAT1/BAT2 are shared** across leu/ile/val ⇒ a shared-BAT
> partition could be *cheaper* — prototype + source before deciding) **[D-116 CLOSES THIS: shared-BAT is a parsimony
> LOSS, and the trade was decided before any prototype ran. The transaminase rate IS sourceable (Koonthongkaew 2022
> Table 2, Km + kcat for both paralogs, six substrates), and D-109's structural premise is CONFIRMED — the paralogs
> agree to 1.2% on leucine. But the `f_i` it would replace are not cheap constants: **four of the five are sourced
> Crépin tracer measurements**, and the one speculative entry (phenylalanine) is **not a BAT substrate**, so the
> mechanism can never reach the number that most wants replacing. The trade is 3 sourced in-matrix constants → 6
> in-vitro ones + an unsourced protein flux + an unsourced [E]: more invented numbers at a worse tier. **The D-98 trap
> does not dissolve, it RELOCATES to [E]** — kcat is not a flux — so D-113's owner gate reattaches to [E] + the
> de-novo-KIC and decarboxylase fluxes, not to the transaminase rate. Per-species sourced `f_i` stand. See D-116.]**;
>
> **[D-118 BUILDS THE DE-NOVO PHENYLPYRUVATE ROUTE D-117 NAMED AS THE UNLOCK — and relocates the defect D-117
> diagnosed. D-117 said the model "has no de-novo route"; it HAS one (the D-100 depletion gate), and `byproducts.py`
> has claimed since D-100 that the anabolic/catabolic split therefore "falls out of the must spectrum". **Measured:
> true in shape, ~11x wrong for this molecule** — the model sourced **18.9%** of its 2-phenylethanol from consumed
> phenylalanine against a derived **~1.7%**. The one-line diagnosis: **the gate encodes AVAILABILITY, and the question
> was PROVENANCE.** Phenylalanine *is* available; it is simply not what most 2-PE is made from, and no half-saturation
> re-tune expresses that. The stoichiometric tell the gate structurally cannot see: the must carries ~0.17 mM Phe and
> the wine makes ~0.24 mM 2-PE — **more alcohol than precursor, in moles**, so full sourcing was impossible, not merely
> generous. SHIPPED: `DE_NOVO_FUSEL_ROUTES`, a sourcing-layer registry capping an alcohol's amino-acid branch at
> `1 - f_de_novo` (producer and undosed isolability untouched), and **`f_non_ehrlich_phenylalanine` now ships at its
> measured 0.975** with a real band 0.531 -> 0.975 (joint C refund 1.125x -> **0.584x**). Two traps caught by probing
> rather than reasoning: (1) a **half-remembered source** — Rollero 2017 labelled only leu/val and never measured 2-PE,
> so it could not speak to this (D-103 as a *mis*-citation); (2) a **constant-share objection aimed at a scenario that
> does not exist** — the breaching `amino_acids_gpl=1.0` is only **28 mg/L Phe** (Phe is 2.8% of the spectrum), i.e.
> essentially Minebois's own must, not the heavy dose that would make feedback inhibition bite. `f_de_novo` is a
> **consistency-closure, NOT a second source** (it and `f` are Minebois's single 2.5% against two denominators), tier
> capped at `f`'s, and its band is clamped at the **analytic breach point 0.971** — a MODEL limit, not an evidence
> spread — which a test recomputes rather than trusts. The designed-to-fail
> `test_the_sourced_lump_breaks_the_carbon_refund_guard` was **inverted, not deleted**, gaining a counterfactual that
> sets `f_de_novo = 0` and pins that the identical 0.975 breaches again. Carry-forwards recorded not tuned: realised
> share under-shoots (1.00% vs ~1.73%, guard-safe), the static share ignores feedback inhibition, and the de-novo
> decarboxylation CO2 is uncharged (widening a gap present since D-19, ~3e-5, both media). 1184->1185. See D-118.]**;
> and a **tracer stoichiometry trap** is already
> visible (isoamyl's valine-derived carbon is {3,4,5}/5 depending which carbons the two decarboxylations remove — **both
> close the ledger**, so the atom assignment must be sourced, not assumed). See D-109.]**
>
> **[D-111 BUILDS THAT NODE (fusel side, valine branch) — and the reason D-104 gave for needing it does not survive.
> `SecondaryFuselRoute` breaks `FUSEL_SPECS`' one-alcohol-one-precursor assumption: isoamyl ← leucine **and** valine
> (Crépin's 23% via KIC), as an algebraic partition in the **sourcing layer, never the producer** — exactly D-109's
> constraint, so an undosed run stays byte-for-byte the core. PAYOFF: valine-derived isoamyl **0% → 1.74%** vs Rollero's
> **2.1–7.5%** — untuned, and **predicted at 1.80% from the sourced shares before the code existed** (the clamp accounts
> for 1.80→1.74). (1) **D-109's {3,4,5}/5 tracer trap DISSOLVES rather than being solved** — it reaches neither the draw
> (the 3C sugar refund falls out of the *net* balance, not atom identity) nor the validation (Rollero defines enrichment
> as *"the fraction of labelled **molecule** with respect to its total production"* — a molecule fraction). Checking
> whether a number was load-bearing beat sourcing it. (2) **D-109's "the 23% is circular as an input" is OVERTURNED, by
> different denominators**: the parameter's is *consumed valine*, the validation's is *total isoamyl*. (3) **D-104's
> DIAGNOSIS is a commensurability error** — its "must under-count" cites a **valine** route against isoamyl's shortfall,
> which is **leucine**-derived (model 1.12% vs Rollero's *leucine* tracer 3.4–17.3%). A valine route cannot close a
> leucine gap; the conclusion (build it) survives, the reason doesn't. That gap is **D-103's gate SHAPE**. (4) **⭐ THE
> CLAMP EXPOSED THE GATE**: at a realistic must (aa=1.0 g/L) leucine's D-100 gate claims **90.9%** of isoamyl while the
> KIC branch wants **31.8%** — two sourced claims summing to **122.7%**. *One claim can be wrong silently; two cannot* —
> adding a second claimant turned a soft calibration miss into an arithmetic contradiction. (5) **⭐ FOUR "INDEPENDENT
> CONFIRMATIONS" PASSED ON THE WRONG ROWS**: Tables S1/S2 are EMF vector images, my first parser keyed rows on draw
> order and was shifted by two, and I reported isobutyric/isovaleric acid as isobutanol/isoamyl — *decisively*. All four
> checks shared the mapping assumption, so they were **one check run four times**; the one true tell (isobutyl acetate
> enriched where isobutanol wasn't — chemically impossible) I rationalised as noise. D-103's *"fetching a source is not
> reading it"*, one layer deeper: **extracting a source is not reading it either.** Fixed by parsing glyph
> **coordinates**. **[CORRECTS THE D-100 BLOCK BELOW: its "verified noise not drift (tracks atol ~1:1)" went STALE and an
> `atol=1e-8` rode on it — D-111's `f_valine` drop cuts the sink's N refund, sharpening N's approach to its true zero,
> and the undershoot grew to **16×** the solver atol. Re-measured, still noise: it COLLAPSES as the solver tightens
> (1e-9→-1.600e-8, 1e-11→-4.010e-12, 1e-13→-1.084e-14), trajectory unchanged — a real over-draw would converge to a
> NONZERO negative. Bound → 5e-8 with the measured table replacing the claim.]** See D-111.]**
>
> **[D-112 MEASURES D-111's "sharpest open item" (the leucine→isoamyl shortfall, 1.12% vs Rollero 3.4–17.3%) BEFORE
> building and RETIRES it — documentation + tests only, no production change (1171 → 1175, +4). The keto-acid node fixes
> **none** of the three things the gap turns out to be. (1) **D-103's gate-shape SPREAD is already absorbed by the D-104
> sink**: sink OFF the catabolic shares span 6% (isoamyl) → 67% (propanol), D-103's 11× spread; sink ON the `(1−f)`
> multiplier compresses ALL FIVE into 1.1–16.4%, roughly Rollero's uniform "2–15%" band. Four of five are reasonable;
> isoamyl is the lone outlier and it is UNDER, not over — so "minor alcohols wildly over-attributed" is stale. (2)
> **Isoamyl sits on a `(1−f)` MASS-CONSERVATION ceiling** = `(1−f)·leucine_C/isoamyl_C = 0.185×6.04% = 1.117%`; Crépin's
> f prices in every non-isoamyl leucine fate, so the node reallocates *how* leucine reaches isoamyl but not *how much* —
> the sharper statement is "no sourcing-layer change on this architecture can", not "a scalar can't". The gate cap (the
> intuitive fix, and the ADVISOR's prediction) is **INERT** (1.12→1.13%, leucine still exhausts): leucine (32 mg/L) is
> too scarce to persist under any draw rate — empirically refuting the hand-prediction, the probe as arbiter. Even f=0
> tops out at 7.0%, below Rollero's 17.3%. (3) **Most of the residual gap is an incommensurate DENOMINATOR**: leucine
> supply MATCHES Rollero (Table S3: 1.3% YAN → ~30 mg/L SM250 ≈ model's 32), but the probe must's `amino_acids_gpl=1.0`
> dose inflates isoamyl to 307 mg/L (Rollero's is 10–120 mg/L, EMF-column-verified). Since leucine consumed is fixed, the
> share ∝ 1/isoamyl, so at Rollero's isoamyl the share is ~2.9% ≈ Rollero's 3.4% floor — the 3–15× "shortfall" collapses
> to ~1× once the denominator matches. Plus a raw-enrichment-vs-net-carbon mismatch (U-¹³C picks up leucine↔KIC exchange
> a net model can't; flagged not over-parsed, per D-111's mis-mapping lesson). (4) **The isoamyl-MAGNITUDE dig
> (owner-requested)**: `k_isoamyl_alcohol` is correctly calibrated — 171.8/172.1 mg/L at 250/300 mgN/L with no aa dose,
> bang on the Wang-2024 172 anchor — so the ~2× over-production is the aa dose's deamination-N sustaining the fusel gate
> (the documented monotone-in-N branch), not a mis-set k. Recorded, not tuned. The keto-acid node keeps its REAL
> motivation — D-104's inverted anabolic split (leu<ile<val<thr vs Crépin's reverse), touched only on the valine side at
> D-111 — and loses the false one (closing the leucine gap). New `tests/test_fusel_catabolic_shape.py` pins all three
> facts; D-111/D-104/`precursor_fates.py` prose corrected. See D-112.]**
>
> **[D-113 MEASURES D-112's one surviving motivation — "does D-111's valine route touch D-104's inverted split?" — and
> the answer is NO, structurally. The inversion is a property of the EMERGENT demand sink (leucine 20.9% vs Crépin
> 77–86%); the shipped model IMPOSES the split via static `f`, so reading it back is the D-108 vacuity trap. An
> **invariance probe** (no demand-sink reconstruction → no invented biomass spectrum, the D-98 trap) settles it: (1)
> **leucine's Ehrlich branch is BIT-INVARIANT** under the route toggle — `ehrlich_draws` clamps the valine branch to
> *headroom above* leucine's full claim, so the route relieves leucine of **0%** of its isoamyl demand (D-111 Finding 5's
> 122.7% clamp cut the *KIC* branch, not leucine's 90.9%). (2) **Total biomass invariant** → leucine's emergent share is
> route-invariant for ANY `w_leu`. (3) The route moves NO concentration at all (both toggles: even the full D-111 route+`f`
> change) — every precursor exhausts, so it is pure carbon-**attribution**; valine's protein-proxy 0.85→0.62 is imposed
> bookkeeping on the *least*-inverted species (45.8 vs 41), never leucine. Advisor's "~1.74% relief" corrected to 0% by
> the code, surfaced not switched. **Un-inverting leucine needs a de-novo-KIC *leucine relief* (not a valine drain) +
> kinetically-limited transamination — an unsourced build, the owner's call.** Doc+test only, 1175→1176. See D-113.]**
>
> **[D-115 BUILDS what D-114 parked and finds the blocker was a category error. D-114 read "provenance to ride in
> state" as a **D-1**-level decision (tier/uncertainty are *derived*, never inside the state floats) and surfaced the
> item as the owner's call. But a **13C isotopologue concentration is not metadata about a value — it is a conserved
> extensive quantity in g/L** that flows and integrates like any other pool, so D-1's "don't wreck the hot loop with
> per-element objects" argument never reached it. D-1 untouched; `state.py`'s extension of it to the word *provenance*
> was over-broad in exactly one place and is corrected. TWO PARTS, one build (shipping the re-route alone would have
> killed D-114's "structurally zero" framing without delivering the observable): (1) **the D-69 5:2-INVERSE RE-ROUTE**
> — `EsterSynthesis` sources isoamyl acetate's C5 off `isoamyl_alcohol` and only C2 off `S`, reading the *same* ratio
> constants as the hydrolysis (moved `aging.py` → `carbon_routing.py`: a second consumer of a private ratio is the
> D-26/D-106 drift setup). The debit is self-limiting with NO clamp — the rate is first-order in the pool it debits,
> so a clamp would break the ester's mass balance rather than prevent anything. (2) **TWO TRACER SLOTS** in g/L of the
> labelled molecule (Rollero's enrichment is a *molecule* fraction, D-111 Finding 3), carbon-ledger weight **zero**
> because each is a sub-quantity of a pool already weighted. **THE ONE-SLOT DESIGN COLLAPSED and the collapse is the
> scoping lesson**: if the alcohol's fraction is flat the slot is REDUNDANT (it reproduces D-111's number); if not, the
> ester still needs quadrature (the D-103 defect) or its own slot ⇒ INSUFFICIENT. No regime makes exactly one right,
> and crediting the hydrolysis return at the *alcohol's* fraction would have baked the answer into the RHS (D-98/D-108
> vacuity, relocated where no test could see it). **MEASURED, not assumed:** the alcohol fraction runs **26.3% → 1.84%**
> across the ester-forming window, so one-slot inheritance errs **13.2%** — D-114's "93% forms late" was TRUE but
> licensed a FALSE inference, because the early 7% draws from a pool 10–25× richer. **Mass shares are not label
> shares.** RESULT: ester **1.93%** vs Rollero ~4%, alcohol **1.84%** vs 2.1–7.5%, **ratio 1.05** — the ratio is the
> deliverable (Rollero's ester sits at its parent alcohol's enrichment; the model reproduces that), the ~2× absolute
> shortfall is INHERITED from the alcohol's known D-111 gap, and `f_valine_to_isoamyl` is **untouched** (tuning it to
> land 4% is D-104's "a missing route, NOT a value to tune" one compound downstream). Three tests retired, all
> predicted — D-114's own docstring said it would FAIL here, and its choice to assert on the CARBON rather than the
> draw list is what made that a clean red instead of a stale green. No conservation test moved. See D-115.]**
>
> **[D-108 — D-107's ⚠ "fix D-27's zero acetaldehyde" is RETIRED, and the bug was in D-107's own aldol. The ⚠ rested on
> "real dry whites hold ~30 mg/L" — a **SULFITED** figure compared against **UNSULFITED** runs; the like-for-like target
> is **2.7 mg/L** (Herzan 2020, PMC7684598 variant (0/0/0)) vs the model's 0.000 ⇒ **D-27 acquitted on the gate**, and
> the prescribed flux-link fix measured **39.999 = a 15× regression** that also compressed dry and sweet to one
> temperature-only plateau. The real defect: `SotolonAldolCondensation` read **TOTAL** acetaldehyde while its docstring
> claimed it "reads the pool the binding depletes" — but `free = total − bound` is a read-only overlay, so SO₂ came out
> **raising** sotolon (dry + 60 mg/L must SO₂: 0.025 → **5.02 µg/L**, 200× to threshold), against Pons, for whom low
> free SO₂ is the prémox RISK factor. The bisulfite adduct's carbonyl is blocked and an aldol IS an attack on that
> carbonyl — the argument `AcetaldehydeReduction` (D-47), `AcetaldehydeBridging` (D-80) and the tannin polymerization
> **already make three times in the same file**. Fixed by reading `free_acetaldehyde` behind the exact `so2_total > 0`
> guard: **zero new constants**, unsulfited **byte-for-byte** (max|diff| = 0.0 ⇒ every D-107 output unmoved), sulfited
> dry sotolon **5.020 → 0.059**. **Emergent payoff — Pons' prémox mechanism itself**: SO₂ *does* deplete via
> `SulfiteOxidation` given O₂, so protection erodes as free SO₂ fades (O₂ 0/5/20/60 ⇒ sotolon 0.059/0.121/2.113/7.639).
> **Newly load-bearing: closure O₂ ingress** — a sealed bottle admits no O₂ here, so it never goes prémox. See D-108.]**
> THE TRIPWIRE FLIPPED: sotolon OAV 0 → 3.22 with the reroute ON; the D-99 limitation test
> deleted per its own instructions and REPLACED with a positive assertion (coverage preserved) + one pinning the autolysis
> dependency. **Methionine — flagged as highest-risk (scarcest, 3 consumers) — was the SAFEST, structurally: no fusel eats
> it** (methional held at OAV 60); the casualty was PHENYLALANINE. *Fourth time a risk call landed on the wrong molecule and
> only looking settled it.* WHAT IT DOESN'T FIX, pinned: threonine still feeds propanol AND sotolon — that competition is REAL
> chemistry over one molecule and the model SHOULD show it; what's removed is the FALSE one (fusels vs arginine, no shared
> chemistry). Guarantees widened from one species to a BLEND and re-proved by BOUNDING it (any {arginine, glutamine} blend
> lies in 1.29–2.14, below biomass's 4.3 and melanoidin's 8 ⇒ D-32/D-38/D-40/D-89 all still structural, no clamp). Autolysis
> releases the MUST spectrum (autolysate is protein hydrolysate — richer in branched-chain/sulfur, poorer in arginine — so the
> error is CONSERVATIVE: under-produces rather than invents). Two `assert_nonnegative` guards moved to atol=1e-8 — they were
> asserting tighter than the solver's OWN atol (1e-9); verified noise not drift (tracks atol ~1:1 across 1e-9/1e-11/1e-13),
> and that it only began mattering now is the fix working (the retired over-release was propping up YAN). **`mercaptans` is
> now the LAST lump in the project.** Next: BOUND the reroute's catabolic fraction (the highest-value follow-up this beat
> creates); a sourced yeast-autolysate spectrum; per-species uptake preference (NCR) + per-species E_a/gates (blocked — D-98
> trap); re-anchor `f_methional` now abundance is explicit; the isoamyl-acetate carbon re-route (D-69's 5:2 inverse);
> speciate `mercaptans`.**

## Build order (dependency-ordered; handoff §6 step 5, re-sequenced per D-66)

```
  sensory / OAV readout layer      ← FIRST beat — COMPLETE at D-95 (slice 2 deferred)
        │  (pure readout over compounds already tracked; zero core risk)
        │  1a. OAV ratio (sourced thresholds)      ← D-67: the honest, sourced part
        └─ 1b. descriptor-space projection         ← D-95 slice 1: vocabulary + the MAX rule
              (slice 2 — weighting/compression/masking, the perceptual math that needs
               parameters — deferred; the additivity seam is the slice line)
  aging chemistry (the "years" axis) ← subsequent beats, one Process at a time
        ├── ester formation/hydrolysis equilibria over time
        ├── oxidation (acetaldehyde/phenolic browning, Strecker)
        ├── oak extraction (vanillin, whiskey lactones — D-77; ellagitannins — D-78)
        ├── tannin–anthocyanin polymerization (red colour / astringency)   ← D-79 direct + D-80 bridged; D-81 fade + D-82 mask
        ├── anthocyanin fading: O₂-coupled (D-81) + O₂-independent thermal (D-83)   ← colour genuinely declines both ways
        ├── tannin–tannin softening: self-polymerization (D-84 direct) + tannin-ethyl-tannin (D-85 bridged)   ← softens with no anthocyanin
        └── micro-oxygenation / reductive–oxidative sulfide evolution; Maillard/sotolon
```

**Why sensory before aging** (full rationale in D-66): the sensory layer is a **pure
readout** over aroma-active compounds the model *already* tracks (esters, fusels,
diacetyl, acetaldehyde, H₂S, 4-ethylphenol, 4-ethylguaiacol, mercaptans) — so it adds
**no new ODE physics and zero risk to the validated core**. It ships first because it
then becomes the **acceptance lens for aging**: once OAVs exist, every aging Process's
effect on the aroma profile is visible immediately. Aging chemistry is heavier — new
speculative RHS Processes on a years-scale phase (phase-based integration, handoff §7
multi-scale) with scattered parameter sourcing — so it comes second, one Process at a
time behind its own tests. The handoff's "aging then sensory" order is reference, not
gospel (CLAUDE.md); the owner's own framing put sensory first too.

---

## Active beat: the sensory / OAV readout layer (handoff §4.2)

### Placement & the isolation firewall (§4.2 cardinal rule)

A **new top-layer package `fermentation.sensory`**, a sibling of `fermentation.analysis`
in the dependency graph:

```
scenario / validation  →  runtime  →  core  →  parameters / units
                    sensory  ─┘  (consumes Trajectory + thresholds; imported by NOTHING lower)
```

- It consumes a `runtime.Trajectory` (state series) plus a threshold table and returns
  OAV series. It imports the core/runtime *downward* only; **nothing in core/runtime/
  scenario imports it back** — the handoff §4.2 cardinal rule: *the sensory layer
  consumes the chemistry; the chemistry never depends on the sensory layer.*
- **Thresholds load directly into the sensory layer, NOT through the compile seam.**
  Unlike `acidbase.yaml` / `vicinal_diketones.yaml` (merged into every compiled scenario
  at `compile.py`'s `shared_files` because a *Process* reads them), **no RHS reads a
  perception threshold** — so a new `sensory.yaml` is loaded by the sensory module on its
  own, never merged into `CompiledScenario.param_values`. The chemistry never even sees
  the sensory params. This is a stronger isolation than any Tier-2 readout.
- **Tier floor.** Every OAV output tier is `Tier.combine(chemistry_input_tier,
  SPECULATIVE)` → **speculative**, *even when the input chemistry is validated*. The
  sensory mapping is itself speculative (§0.1 / Tier docstring names "sensory mapping"
  as the canonical speculative case), so it caps everything it touches. Read the input
  compound's tier via `ProcessSet.tier_of(pool, ...)` and combine with speculative.

### Definition of done (beat 1a — OAV ratio only)

1. `OAV_i(t) = concentration_i(t) / threshold_i` for each aroma-active compound, mapped
   over a trajectory (mirroring `analysis.ibu_series` / `molecular_so2_series`): an
   `oav_series(traj, thresholds, compound)` per compound plus a finished-profile view over
   the medium's active compounds at a chosen time. Dimensionless. **The aggregate reports
   per-compound OAVs and flags which exceed 1 (above-threshold) — NOT a single summed
   scalar**: summing OAVs assumes perceptual additivity, which is contested, so a summed
   number would over-claim (settle the exact aggregate shape in D-67).
2. A `sensory.yaml` provenance file with **real, sourced perception thresholds** (see
   sourcing below), each carrying value, unit, `source`, `conditions` (**the matrix** —
   see below), `uncertainty`, `tier: speculative`.
3. Unit tests: OAV monotone-increasing in its pool; **identically 0 when the pool is 0**
   (an unhopped/unspoiled/clean run has no false aroma); **tier is the speculative floor**
   even for a validated input; the reported matrix matches the medium. Plus a golden
   sanity check that a known concentration → the literature OAV (e.g. diacetyl at ~2×
   its lager threshold reads OAV ≈ 2).
4. **The Tier-1 suite and all conservation tests are byte-for-byte untouched** — the
   readout adds **no state slot, no RHS, no ledger entry** (prime directive #3, trivially:
   nothing in core changes). `pytest` / `ruff` / `mypy` green.

### Compounds and their thresholds

Aroma-active pools already in state (`core/media.py`), classified by medium and by how
cleanly OAV applies. **The set is medium-specific** (mirroring the beer-only `iso_alpha`):
`_common_specs` (both media) carries only `esters`, `fusels`, `diacetyl`, `acetaldehyde`,
`h2s`; `ethylphenols`, `ethylguaiacols`, and `mercaptans` are appended in `wine_schema`
**only**. So the **beer OAV set = the 5 common compounds**; the **wine OAV set = those 5
+ 4-EP + 4-EG + mercaptans**. A beer `oav_series` must not reach for the three wine-only
slots (they do not exist in the beer schema).

| Pool | Medium | Representative compound (threshold key) | Descriptor | Note |
|------|--------|-----------------------------------------|-----------|------|
| `diacetyl` | wine + beer | 2,3-butanedione | buttery | single molecule — clean OAV |
| `acetaldehyde` | wine + beer | acetaldehyde | green apple / bruised | single molecule — clean OAV |
| `h2s` | wine + beer | hydrogen sulfide | rotten egg | single molecule — clean OAV |
| `ethyl_acetate` / `isoamyl_acetate` / `ethyl_hexanoate` | wine + beer | each ITSELF (D-96 split the lumped `esters`) | solventy / banana / apple | single molecule — clean OAV |
| `isoamyl_alcohol` | wine + beer | isoamyl alcohol (3-methylbutan-1-ol) | solventy / fusel | single molecule since D-99 |
| `isobutanol` / `2_phenylethanol` | wine (threshold), both (chemistry) | each ITSELF (D-99 split the lumped `fusels`) | solventy / rose | single molecule; wine-thresholded, chemistry-only in beer |
| `propanol` / `active_amyl_alcohol` | wine + beer | — (no sourced threshold in any matrix) | — | single molecule; CHEMISTRY-ONLY, no OAV (D-99) |
| `ethylphenols` | **wine only** | 4-ethylphenol | horse-sweat / barnyard | single molecule — clean OAV |
| `ethylguaiacols` | **wine only** | 4-ethylguaiacol | clove / smoky | single molecule — clean OAV |
| `methanethiol` | **wine only** | methanethiol | reductive / drains | single molecule — clean OAV (was the false lump `mercaptans` until D-110) |

- **The lumped-pool call (D-66, owner-chosen) — NOW DOWN TO ZERO POOLS (D-110).** A lumped pool is a
  single g/L pool that really mixes several molecules whose thresholds span ~3 orders of
  magnitude; it is assigned the threshold of one **named representative compound** and carries
  **"assumes fixed lump composition"** loudly in that threshold's provenance `notes`. D-96
  split the `esters` lump into three single-molecule pools and D-99 split the `fusels` lump
  into five, leaving `mercaptans` as the last — and **D-110 retired it too, but by a different
  route: it was not a real lump at all.** No Process makes ethanethiol or any other thiol, so the
  pool held exactly one molecule under a plural name; the flag asserted an uncertainty the mass
  balance did not carry. Renamed to `methanethiol`, `lumped=False`, **provably zero output change**
  (the raw 84-slot state array is byte-identical, copper fining included) — the pool already *was*
  the molecule at every layer, which is why the relabel was free. **The D-66 lump-composition risk
  class is CLOSED.** The lumped-flag tests
  derive from `AromaCompound.lumped`, not a hardcoded list, so the caveat cannot linger on a
  pool that stopped being lumped — the D-66 honesty cost retires under test pressure as each
  pool is speciated. D-99 also showed *why* a self-consistent lump is dangerous even when its
  ledger closes: reading every higher alcohol at isoamyl alcohol's threshold mislabelled
  2-phenylethanol (rose) as solventy — self-consistency is not correctness.
- **`iso_alpha` / IBU is excluded** — it is a **taste** (bitterness), not an odor, and is
  already a direct mg/L→IBU readout (`analysis.ibu_series`, D-64). Do not shoehorn a
  bitterness into an odor-threshold OAV. (A future taste-intensity readout is separate.)

### Units & matrix (both load-bearing)

- **Units:** state is g/L; literature odor thresholds are µg/L–mg/L. Convert at the sensory
  boundary via `fermentation.units` (add g/L↔µg/L helpers if absent). OAV itself is
  dimensionless.
- **Matrix-specificity (a provenance requirement, not optional).** Ethanol and the wine/beer
  matrix shift most odor thresholds substantially — a wine-matrix threshold ≠ beer ≠ water/
  model solution. Every threshold's `conditions` **must record the matrix it was measured
  in**; the wine profile reads wine-matrix thresholds where they exist, beer reads beer-
  matrix, and any water/model-solution fallback is flagged as a matrix gap in `notes`.

### Parameters to source (provenance, like the D-12 sweep)

Perception (odor) thresholds, matrix-specific. Reading list (reconcile, don't transcribe):
Guth 1997 (wine aroma thresholds), Francis & Newton 2005 (wine flavour compounds review),
Meilgaard 1975 / Meilgaard et al. (beer flavour thresholds), Ferreira et al. 2000 (wine
odor-activity), plus diacetyl (~0.1 mg/L lager) and 4-EP/4-EG (~425 / ~110 µg/L red-wine
sensory) from the spoilage literature already cited in `vicinal_diketones.yaml` / the Brett
beats. Tiers: **all `speculative`** — thresholds are panel- and reference-dependent and swing
by orders of magnitude, and the OAV *mapping* is speculative regardless; carry wide
uncertainty bands.

### Approach (test-driven, mirrors the D-64 IBU readout)

1. Create `src/fermentation/sensory/` (package) with an `oav` module; add `sensory.yaml` to
   `parameters/data/` and load it standalone (a `load_thresholds()` helper, NOT via
   `compile.py`'s `shared_files`).
2. Implement `oav_series` per compound + a finished-profile aggregate; each: dimensionless,
   monotone in its pool, 0 when the pool is 0, tier = `combine(input_tier, SPECULATIVE)`.
3. Source + reconcile the thresholds; replace any placeholders; record the matrix in
   `conditions`; flag the lump-composition assumption in the three lumped thresholds.
4. Unit tests (per the DoD) + a golden OAV sanity check. Confirm the §2.2 trio +
   conservation + full suite stay byte-for-byte green (readout touches no core state).
5. Record outcomes in a DECISIONS entry (D-67) and update this plan + ARCHITECTURE.

---

## Deferred / later beats (in order)

- ~~**1b. Descriptor-space projection**~~ — **slice 1 BUILT at D-95**, **slice 2 BUILT AND CLOSED
  at D-98**. Slice 1 (`sensory/descriptors.py`) maps the OAV vector onto a 14-axis (wine) /
  9-axis (beer) descriptor vocabulary behind the `DescriptorProjector` seam (handoff §4.2: so an
  ML model trained on sensory-panel data could later replace it), max rule, zero params.

  **Slice 2 (`sensory/compression.py`, D-98) is DONE — and its result is that it may claim
  nothing.** Of D-95's four scoped items, three were never real: **matrix effects** were already
  discharged at beat 1a (thresholds are matrix-specific); **weights** are subsumed into the
  Stevens `k`, which a threshold pins (D-95 mis-cited `thermal.yaml`, whose weights are
  *production* fluxes); **masking/suppression** is BLOCKED on per-pair `cosα` coefficients that
  exist for no pair of our pools — **the one genuinely deferred item, with its unblock condition
  named**. Only **per-compound compression** was real, and it ships behind the seam, isolable,
  **not the default**.

  **The additivity objection is answered, not dodged**: compression is applied *per compound,
  before* the combination rule, and the rule is **still max**, so no additivity is assumed at any
  layer. Max is a deliberate **under-claim** — mixture perception is hypoadditive, so the truth
  lies between max and sum, in exactly the region the missing `cosα` would be needed to reach.

  **The result is a theorem: no dominance flip this layer produces can be trusted, ever.** A
  robust flip requires two pools on one axis with *disjoint* exponent bands; none exist, because
  the bands are wide *because* the exponents are author estimates. An honest band and a
  trustworthy flip from a guess are mutually exclusive — so slice 2 is informative only where it
  is redundant. **Do not re-propose this beat**; `test_a_robust_dominance_flip_is_impossible_at_
  these_bands` fails the day measured in-matrix exponents with narrow bands exist, and that
  failure is the only signal that would reopen it.
- **Aging chemistry (§4.1), one Process at a time on a slow/years phase** — ester
  formation/hydrolysis equilibria; oxidation (acetaldehyde generation, phenolic browning,
  Strecker degradation); **oak extraction** (diffusion-limited vanillin / whiskey lactones /
  ellagitannins as a function of surface-to-volume ratio, toast level, barrel age — new
  extracted pools, akin to dosed inputs); tannin–anthocyanin polymerization (red-wine colour
  and astringency evolution); micro-oxygenation / reductive–oxidative sulfide balance;
  long-aging Maillard / sotolon in oxidative styles. Each: `speculative`, isolable/togglable,
  **phase-based integration** (handoff §7 — do *not* integrate years at fermentation
  resolution), and validated only by the sensory lens built above. The honest framing stays
  "tracking a handful of key compounds with literature-based kinetics," **not** "predicting
  what five years in barrel tastes like."

## Risks

- **Validation is essentially absent at Tier-3** (handoff §4.3). Odor thresholds vary by
  matrix, panel, and reference by orders of magnitude; aging verdicts are not solved science.
  Carry wide uncertainty bands, tag `speculative`, **never** `validated`, and lean on
  directional checks (more diacetyl ⇒ higher buttery OAV) not magnitudes.
- **The §4.3 credibility firewall.** The Tier-3 temptation is to let plausible-looking
  speculation borrow the validated core's credibility. Every OAV / descriptor / aging output
  **must surface its speculative tier** in any plot, export, or report — the tier-floor rule
  above enforces this at the API.
- **Lump-composition assumption** (the accepted D-66 (a) call): a lump's OAV is only as
  meaningful as its assumed fixed composition. The risk was flagged in every lumped threshold's
  provenance and *retired as predicted* — the `esters` lump was speciated at D-96 and the
  `fusels` lump at D-99, both chemistry-layer changes motivated on their own merits (each `k`
  independently anchored, never to serve the sensory layer — §4.2). `mercaptans` is the last
  pool still carrying the assumption. D-99 also surfaced a **new lump risk one layer down**: the
  `amino_acids` pool (arginine standing in for five amino acids) cannot serve the now-speciated-
  scale consumers that draw on it (fusels, Maillard, MLF, Brett) — deferred as D-100 and pinned.
- **Aging multi-scale stiffness** (handoff §7): the years phase must use phase-based Process
  activation and an appropriate step regime; do not integrate aging at ferment resolution.
- **Thresholds sit outside the D-24 ensemble sweep** — a deliberate consequence of loading
  `sensory.yaml` standalone rather than through the compile seam: `simulate_ensemble` samples
  only compiled-scenario params, so it will *not* propagate threshold uncertainty into the OAV
  band. Defensible for a speculative readout (the OAV floor is already speculative), but state
  it explicitly in D-67 so it never later reads as an oversight, not a choice.

> **[D-119 READS MINEBOIS FIG. 6A — the unlock D-118 named — AND NOTHING SHIPS. `f_de_novo_2_phenylethanol` stays
> 0.9827, `f_non_ehrlich_phenylalanine` stays 0.975; only the caveat's standing changes. Read off the PMC deposit's
> figure image (Wiley still 402). Sc 2-PE: ~4.2 uM labelled / 109 uM unlabelled / ~113 uM total, printed I.E 3.7%.
> **The bar semantics needed deriving and the obvious reading is wrong** — the large number is the UNLABELLED segment
> (as total it reproduces isobutanol's I.E as 9.63% vs a printed 8.8%; as unlabelled, 8.78%, and isoamyl agrees).
> **The gain is the NUMERATOR, not the fraction:** 4.17/113 = 3.69% and 4.17/162 = 2.57% recover both printed numbers,
> so the 2.5% D-117/D-118 rest on is corroborated in uM independently of the algebra. **TWO WRONG INSTINCTS, both
> recorded because the errors are instructive:** (1) "the measured 0.963 breaches the 0.971 guard" — no, 0.963 has
> Minebois's ~113 uM total as its denominator, and transplanting it demands 8.7 uM from a branch that supplies ~4.2;
> **a de-novo fraction is not scale-invariant**; (2) "0.963 is a T4 snapshot that climbs to 0.9827" — refuted by Panel
> B, where the fraction is **FLAT at ~3.7%** (labelled 4 -> 4.8 uM as total 109 -> 130 uM). What actually holds 0.9827
> up is a single cross-must scaling assumption (a Wang-typical wine carries ~Minebois-like Phe while making ~2x her
> 2-PE) that one fermentation at one Phe dose cannot settle — so **the blocker MOVED rather than lifted**, from "figure
> unread" to "does the Phe flux scale with total 2-PE?". Residual risk is **guard-safe** (under-attribution => smaller
> refund), hence no re-ship; 0.963 must never enter a sampled field. All three alcohols are de-novo dominated in-study
> (2-PE 0.963, isoamyl 0.946, isobutanol 0.912), so **D-118's finding is a CLASS of error, not a phenylalanine
> peculiarity** — but those numbers inherit the same cross-must caveat and are NOT drop-in `f_de_novo` values for the
> isoamyl beat. Sharpened ask: Phe dose vs total 2-PE (Querol, aquerol@iata.csic.es). See D-119.]**
