---
name: feedback-a-constant-can-be-a-unit-definition
description: "Two disagreeing values of a constant may not be rival estimates — one can be the DEFINITION of the unit the model's fitted parameters live in; ask which paper supplied them"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1134a524-b812-43b7-a9d8-82765ebb9d7d
  modified: 2026-08-18T10:18:21.908Z
---

When the repo ships two values of one constant, **do not open a literature search for a third.**
First ask: *what unit is the model's own fitted parameter expressed in, and which paper defined it?*

D-219. Per-cell yeast dry mass was "load-bearing in three places, 5.6× apart, source it". Both
shipped values turned out to be the wrong kind of thing:

- **18 pg** — asserted in two wine benchmarks, sourced to nothing.
- **~100 pg** — never chosen by anyone. **Back-computed** by dividing a scenario's `pitch_gpl = 1.0`
  by a counted pitch, so it is a **residual** absorbing the true cell mass **× every error in the
  model's per-gram rate**. That is why it landed 3× above any defensible cell mass — it was never
  measuring a cell. A number implied by dividing two other numbers is not an estimate of anything.

The answer was on disk since D-13: Coleman 2007, the paper wine's `Y_X/N`, `k'_d` and `mu_max` come
from, states *"assuming that each cell weighs 4 × 10⁻¹¹ g"*. He **counted** cells and weighed none,
so **every gram in that paper — hence every parameter fitted to it — is a count × 4e-11 g.** The
constant is a **unit identity, not a preference**; converting at anything else feeds the model a
number in a unit its own parameters do not use. A third literature value would have been a fourth
incomparable number.

**Why:** a search framed as "which value is right" cannot terminate when the values answer different
questions. Reframed as "what is the model's gram", it terminated in one paper already cited.

**How to apply:**
- Ask **who defines the unit** before asking who measured the quantity. Trace the constant to the
  paper that supplied the *parameters*, not to the best paper about the *substance*.
- Assert a unit identity **exactly, not to a tolerance** — an identity has no tolerance.
- Then get an **independent corroboration that cancels the assumption**: invert the source's yield
  back to what was actually measured (cells per g N) and price it with something the author had no
  hand in (elemental formula) → 34.9 vs 40, and it settles a frame the source left open ("cell
  mass", never "dry weight": read as wet it makes yeast 33 % N — absurd).
- Say **"assuming" out loud.** Coleman assumed it too, so tier is `plausible`; the `source:` field
  must not launder an assumption into a measurement ([[feedback-a-derived-yield-encodes-its-rate-law]]).
- **Derive the band from a parameter, not from prose** — here `biomass_N_fraction`'s own 0.08-0.14
  gives 28-50 pg, and both retired readings fall outside it.
- Refuse to build a band from **product-dosing conventions** (g/hL of dried yeast ↔ cells/mL): the
  same book gave 100-200 pg and 25-50 pg, and that spread *is* the convention's looseness
  ([[feedback-a-units-fork-is-not-a-band]]).
- **Never call a repair "inert" because its asserts survived.** All six Varela pins held and I wrote
  "held unchanged" — but the quantity they *characterise* widened at both nitrogen levels
  (1.910→2.024× and 2.236→2.422×), in the direction the file's own docstring says must be caught.
  Check whether the characterised number moved, not whether the assert passed. On the Palma arm I
  did exactly that and reported the cost; on the arm where the asserts happened to survive I did
  not ([[feedback-conceded-caveats-are-not-coverage]]).
