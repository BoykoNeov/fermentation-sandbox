---
name: closure-oxygen
description: "Closure oxygen — the steady OTR menu and the bottling burst (D-136, D-162, D-187); both columns of Lopes' table are now shipped"
metadata: 
  node_type: memory
  type: project
  originSessionId: 5f69596c-babc-41fb-a985-7a1c89b328ef
  modified: 2026-08-12T00:02:44.725Z
---

**Live prohibitions — closure oxygen, steady and burst.** Split out of
`.claude/memory/project-fermentation-sandbox.md` at D-185's pattern; that file's ledger points
here by path. Read it when working on this subject. Every bullet is *what it forbids* + the
D-record to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not
argue past it from this file.** **Never evict an old prohibition to buy a line.**

**BOTH columns of Lopes 2007 Table I are now shipped (D-136 steady, D-187 first month).**
- **`seal_bottle` is BUILT — the bottling burst is no longer an author-invented number.** Doses
  `bottling_burst_<closure>` from `scenario.closure`; **takes NO params** (an author-supplied
  charge is `add_oxygen`, and that split must not be merged). D-136's "left out deliberately …
  not for lack of data" is SPENT. **Never make it fire automatically** from `closure` alone:
  opt-in is what keeps every pre-D-187 scenario bit-identical, incl. D-140's 31 pins.
- **The dose is P1's first month NET of 30 days' steady ingress** — an anti-double-count worth
  **1 % (technical cork) to 32 % (SupremeCorq)**, never a rounding tidy-up. It is why the verb
  **refuses to precede `begin_aging`**: earlier, it subtracts a flux the run never paid.
- **ONE SCOPE PER BAND: subtract the SHIPPED `otr_<closure>` at value AND both edges.** Mixing in
  the steady *band's* edges makes the nominal stop being its own midpoint (a test asserts it is).
  **Screwcap is a BOUND, not a centre** — P1 prints `<500 µL`, so value = the ceiling, low = 0.0
  CONSTRUCTED. **Never "fix" it to a midpoint of [0, ceiling]**: that invents a second number,
  and it is deliberately in the nominal-on-edge census set — classify it over the band INTERIOR.
- **The burst ordering is NOT the closure ordering, and near-equality is NOT a mechanism claim.**
  The most permeable closure has the *smallest* burst (largest steady subtraction); the four
  cylindrical charges agree within 1.16-1.44 mg/bottle only because **P1's first-month ranges
  overlap**, a fact about the table. D-136's steady ordering correction is untouched.
- **Line/headspace oxygen is EXCLUDED on purpose** (P2's 1.4-1.9 mg per 375 mL is a property of
  the LINE). It belongs to `add_oxygen`. **The one exception is the screwcap**, where P1 says the
  charge IS headspace air. **Never fold a fixed bottling charge into every closure.**
- **D-136 transposed two cited numbers and D-187 fixed them** (`Corrects:`): Godden's natural cork
  is **9.33 mg/stopper/yr**, Lopes 2005's synthetics **1.60 mg first month** — each was filed under
  the other. Godden is now excluded **on method** (P1 names it: no liquid contact), not
  reconciled; the repaired 1.60 is what cross-checks the synthetic bursts. Never re-swap them.
- **OTR(T) is still BLOCKED** — needs one closure at two storage temperatures; warm-storage output
  stays a **lower bound**. **Bottle format** (magnum/half) stays inexpressible: the aging model is
  a **750 mL** bottle. **Closure-driven reduction** needs D-135's blocked de-novo route.
- **Under an ensemble every closure parameter is NOMINAL-ONLY** — they seed state or an event dose,
  not a per-RHS read, so bands do not propagate. Same class as `copper_typical`, `bound_h2s_initial`
  and D-186's `set_ph`. **Explore closure uncertainty as separate scenarios, never by sampling.**
