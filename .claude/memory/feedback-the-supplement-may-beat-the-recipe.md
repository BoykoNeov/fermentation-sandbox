---
name: feedback-the-supplement-may-beat-the-recipe
description: A paper citing a house medium recipe may still hold the MEASURED composition in the supplement its Methods points at; read that before chasing the cited paper
metadata:
  node_type: memory
  type: feedback
---

When a Methods section says "the medium was prepared as in <1990 paper>" it looks like the
composition lives in the cited paper and the current one is a dead end. **Read the supplement the
sentence points at first.** Crepin 2017 gave only "180 mg nitrogen . liter-1 ... (Data Set S1)";
Data Set S1 held *measured* initial concentrations in mM for all 20 species, mean of 14
independent fermentations, plus per-species consumption at four time points and the residuals at
end of fermentation. That is strictly better than the recipe it cites -- a recipe is nominal, this
was assayed -- and it answered a question (how much N is left at dryness) nobody had asked.

**Why:** the cited-recipe hop is one or two extra fetches and often lands on a paywall, while the
in-hand supplement is already open. Chasing the citation first wastes the budget and can end the
beat as "sourcing-blocked" when the number was one file away.

**How to apply:** grep the Methods for a Data Set / Table S / supplementary pointer beside the
number you want, and fetch that before the citation. Two mechanics that cost real time on D-246:
PMC gates supplement binaries behind a proof-of-work challenge (`cloudpmc-viewer-pow`: SHA-256 of
`challenge + nonce`, four leading zero hex digits, cookie `challenge,nonce`) which is scriptable
in ~20 lines; and **a wide PDF table's `Mean` row can wrap so its middle columns land on the line
labelled `SD`** -- transcribing the row as printed swaps means for standard deviations, silently.
Cross-check against the individual replicates above it. See [[feedback-a-generic-partition-is-not-a-defined-medium]],
which is the failure this lesson is the successful other half of.
