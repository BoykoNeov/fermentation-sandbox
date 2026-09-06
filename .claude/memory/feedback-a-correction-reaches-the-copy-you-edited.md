---
name: feedback-a-correction-reaches-the-copy-you-edited
description: "A retraction lands in the file you had open and leaves the other copies standing; grep the retracted WORDS across every surface, and fix the entry point first"
metadata:
  node_type: memory
  type: feedback
---

**A claim lives in more places than the one you were editing, and the copy that matters is the
one a reader lands on first — not the one you happened to open.**

D-273 retracted D-215 §3's "the three beer acid timings OPPOSE and nearly cancel". Its memory
update rewrote the prohibition file's frontmatter `description:` and its body block, and the
project memory's status paragraph. It left **two** copies standing: that file's own header,
which still said "go read D-215", and — the one that mattered — the **ledger row in the
project memory**, which is loaded at every session start and still stated the retracted claim
in full. The boot surface therefore named both the claim and its retraction, with nothing
saying which won, and a reader grepping the ledger (which is what the ledger is FOR) never
reaches the paragraph that overturns it. Two further copies of a re-measured ratio sat in a
third and fourth file.

**Why:** you correct the copy you are looking at, because that is the one that prompted the
correction. Nothing brings you back to the others — headers and ledger rows are boilerplate-
adjacent, read on the way past, and never re-read on purpose. Structure makes it worse: a
correction is *per-clause*, so a sibling file citing the same record for a different clause is
genuinely fine, and that makes "which citations are now wrong" un-derivable from the record
number alone (measured at D-274: a citation-level linter fires 547 times whole-file, 155
header-only, nearly all false).

**How to apply:** after writing a retraction, **grep the retracted WORDS**, not the record
number — the distinctive phrase and every number it carried — across the project memory,
`prohibitions/`, `lessons/` and the tests. Then fix in this order: the **entry point** first
(the ledger row, the index, the description — whatever a reader hits before opening anything),
then the body. For every other hit, read it and decide explicitly: restating the retracted
claim means rewrite; citing the same record for an untouched clause means leave it, and say in
the record that you read it and left it. Quoting a re-measured number in an argument that does
not rest on it means annotate the number and leave the argument.

Related: [[feedback-a-doc-rots-where-it-duplicates]] (ask which surface OWNS a claim),
[[feedback-verify-latest-state-not-breadcrumbs]] (the reader's side of the same failure),
[[feedback-a-cap-being-written-to-cannot-be-raised]] (the reason given for leaving the header
was a cap, and the cap was not the obstacle).
