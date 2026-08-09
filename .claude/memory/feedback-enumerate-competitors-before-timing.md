---
name: feedback-enumerate-competitors-before-timing
description: Any wall-clock number on this box is invalid until competing project suites are enumerated; other agent sessions routinely run 26-worker pytest at Normal priority
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 76b1a70a-57a3-4582-be49-6cf1bfb7f66b
  modified: 2026-08-09T14:15:44.483Z
---

This machine runs several Claude sessions at once, and they run each other's test
suites. A measurement taken without checking for competitors is not a measurement.
Observed in one session: the Fermentation suite took **119 s** under light background
load and **363 s** for the identical command while `space-station` ran a 26-worker
pytest at Normal priority — a 3× swing caused entirely by neighbours.

**Why:** The cost is not just a wrong number, it is a wrong *conclusion*. The 363 s run
was launched to compute a perfect-balance floor, and the floor it produced was
meaningless. Worse, the first run got labelled "alone" in a committed doc when it had a
background census finishing inside it — a provenance claim that was simply false, in a
repo whose prime directive is provenance. Contention is also **non-stationary**:
competitors start and stop mid-run, so it cannot be normalised away by taking ratios or
by running an A/B back-to-back. Two 6-minute runs 6 minutes apart are not comparable.

**How to apply:** Before any timing run, and again before quoting the result, enumerate
competitors:

```powershell
Get-CimInstance Win32_Process -Filter "Name='python.exe' OR Name='uv.exe'" |
  Where-Object { $_.CommandLine -notlike '*<thisproject>*' }
```

If the count is nonzero, either wait or label the number with the competing load that
produced it — never with "alone", "baseline", or "clean". Prefer
**contention-immune evidence** where it exists: counts, not durations. The scheduler
question in that session was settled for free by counting fixture *instantiations* in a
`--durations=0` file already on disk, after two expensive timing runs had settled
nothing. See [[feedback-a-majority-is-not-a-direction]] for labelling rows with the run
that produced them, and [[feedback-pin-tolerance-vs-solver-tolerance]] for the related
habit of measuring the noise floor before trusting a difference.
