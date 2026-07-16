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
| `mercaptans` | **wine only** | **methanethiol** (stand-in) | reductive / drains | **lumped** → representative threshold — the LAST lump (D-99) |

- **The lumped-pool call (D-66, owner-chosen) — NOW DOWN TO ONE POOL.** A lumped pool is a
  single g/L pool that really mixes several molecules whose thresholds span ~3 orders of
  magnitude; it is assigned the threshold of one **named representative compound** and carries
  **"assumes fixed lump composition"** loudly in that threshold's provenance `notes`. D-96
  split the `esters` lump into three single-molecule pools and D-99 split the `fusels` lump
  into five, so **`mercaptans` is the last lumped pool in the project**. The lumped-flag tests
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
