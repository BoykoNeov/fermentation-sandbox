---
name: feedback-re-read-the-source-you-already-mined
description: "Before declaring a beat blocked on external sourcing, re-read the papers already on disk in full — the number that unblocks it may be in a figure of one you already mined for a table"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: eda906e6-10f9-4458-bce2-c2bce96f3ec9
  modified: 2026-08-17T14:45:45.139Z
---

When a beat looks **blocked on an external number**, re-read the sources already on disk
**in full** before recording the block. A paper mined once for one thing is not spent.

**The case (D-180).** D-179 read Tyrell 2013 for its Table 1 — a literature compilation of
finished-beer acid levels — and stopped there, recording "beer has no organic-acid producer"
as open partly for want of a **wort** composition. The next beat's blocker was exactly that
missing t=0 end of a difference. It was in the *same PDF*: sections 2.4.2/3.2 report the
authors' **own trials**, whose figures carry one wort's full day-0-to-day-7 acid course for
four strains, **plus the pH and extract curves of the same ferments**. That single re-read
supplied the seeds, the yields, the sugar divisor, **and** an independent falsification
target — a matched dataset better than anything the archive had for the axis.

**The second case (D-183) sharpens it: not only other FIGURES — the INTERIORS of the figures
you already used.** D-180 read Figs 6-14 and Fig 4 for **two points each**, day 0 and day 7,
because a yield is a difference. Three beats later the *interiors* of those same archived crops
falsified **both halves** of the fix D-180 had proposed for acetic's transient: mapping the
acid course onto the extract course puts 86 % of the rise inside the first 15 % of the sugar
flux (production not flux-linked), and half the fall at **zero** flux (removal not flux-linked
either). Endpoints give you **sizes**; interiors give you **rate laws** — and a beat that only
needs sizes never notices it is assuming one [[feedback-a-derived-yield-encodes-its-rate-law]].

**The third case (D-211) is the cheapest miss of the three: another PANEL of a figure two
beats had already cropped.** Tyrell's Fig. 4 is three stacked panels of one trial — extract,
pH, and **total cell count**. D-180 cropped the extract panel, D-207 the pH panel, and both
saved the crops to disk. Nobody cropped the third. Meanwhile D-209 located beer's remaining pH
defect as "uptake timing", D-210 parked two candidates behind it, and a test docstring recorded
that the timing "was never calibrated against a wort FAN time course" — while the only measured
growth curve the project has, for the very ferment being scored, sat in the panel below the one
being read. It re-derived `mu_max` (2.88× too fast) with no new source at all.

**The fourth case (D-258) is cheaper than all three and it is the one to remember: the sentence
BETWEEN two sentences already quoted.** D-213 built `k_o2_uptake_beer`'s timescale argument out
of two quotes from one paragraph of *The Chemistry of Beer* — the lag phase ("several hours to
adapt … before growth begins") and the removal ("the oxygen present at the start of pitching is
rapidly used up"). Both are in the parameter's provenance. **The sentence sitting between them
is "the yeast will multiply four- or fivefold by a process of budding"** — the only printed beer
growth EXTENT in the corpus, and the quantity D-222, D-230 and D-232 then spent three records
auditing against a single experiment's cell counts, because nothing else was known to exist. The
passage was open, quoted, and cited from *both sides*.

**What that adds to the three above: the blind spot is not the source, it is the QUOTE
BOUNDARY.** Cases 1-3 are "you did not open the other figure/panel". This one is "you had the
words on screen and the extraction was scoped to the claim you arrived with". A number that
answers a question you are not currently asking reads as connective prose and gets skipped —
and it is *most* likely to be skipped in a passage you are quoting, because quoting feels like
having read.

**Why:** a paper gets mined for the one thing you went looking for. Tables are what a search
surfaces and what [[feedback-transcribe-tables-not-prose]] trains you to prefer, so a
compilation table can mask the authors' own experiment sitting a page later in figures. The
cost of re-reading a PDF already on disk is minutes; the cost of recording a false block is a
beat deferred indefinitely, and it compounds — the block gets copied into memory as fact.

**How to apply (D-258's addition first, because it is the cheapest):** when you quote a source,
**read the whole paragraph and say what every number in it is**, including the ones irrelevant to
the beat in hand — and when a later beat needs a quantity, re-read the *passages this repo already
quotes* before searching for new ones. `grep` the repo for the source's title and re-open every
passage it cites. Quoting a sentence is not evidence its neighbours were read.

Then: before writing "blocked on sourcing", list the sources already local
(`_pdf/`, `_txt/`, the `SOURCES.md` of past beats) and check each for the *missing* quantity
specifically, including **figures and methods sections**, not just tables. Ask what the
authors measured that they did not tabulate. A screening/comparison paper almost always ran
its own trials. If the number does turn up in a figure, it is a **figure read**: record the
figure number and a read tolerance, and never let it acquire table-grade provenance
downstream ([[feedback-pin-the-band-not-the-nominal]]). Related but distinct from
[[feedback-paywalled-is-one-host]], which is about other *hosts* of a source you cannot open;
this is about other *content* in a source you already have open.
