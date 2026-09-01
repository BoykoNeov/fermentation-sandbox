# The Fermentation Console

A local app for describing a batch, running it, and seeing how far the answer can be
trusted. It sits above the engine and the engine never knows it exists.

```bash
uv sync --group ui
uv run streamlit run app/main.py
```

It opens in a browser at `http://localhost:8501`. Nothing is uploaded anywhere; the whole
thing runs on this machine.

## Why it exists

Plenty of software draws a fermentation curve. This engine can do something else: for any
point on that curve it can say which chemistry produced it, which published numbers that
chemistry used, what conditions each number was measured under, how wide a range the source
gave it, and therefore how much the curve is entitled to be believed. Every one of those
facts is already stored — a number cannot be loaded into the engine without them — so the
console's job is to make them reachable rather than to invent anything.

Hence the two things it is built around, and one it is built to avoid:

- **Confidence is drawn, not footnoted.** A line's dash pattern is how far it can be trusted,
  so that survives a screenshot. Solid is reserved for "checked against real measured data",
  which nothing in this engine has earned, and the legend says so. A scale whose top mark can
  never be reached is more honest than one quietly rescaled to fit.
- **The papers are reachable, not central.** The source trail lives in an expander at the
  foot of the run, not in a headline tab. Someone who wants a curve gets a curve; someone who
  wants to know where a number came from is two clicks away.
- **It never invents a warning-free number.** Readouts the project knows are misleading carry
  their warning onto the chart every time, and untouched quantities are reported as inert
  rather than as confirmed.

## Layout

| File | What it is |
| --- | --- |
| `main.py` | The Streamlit page: the form, the layout, the caching. The only file that imports Streamlit. |
| `render.py` | Finished run in, figures and panels out. No framework. This is what makes the app and the written report the same drawing rather than two that drift. |
| `readouts.py` | What may be drawn, in what units, with what warning attached. |
| `runner.py` | Compile-and-run boundary, plus the uncertainty ensemble and the convergence check. |
| `provenance.py` | Curve → chemistry → numbers → papers. |
| `report.py` | One run written to a single self-contained HTML file. |
| `fidelity.py` | The three separate things people call "fidelity". |
| `library.py` | Starter batches, and the form metadata read from the engine's own tables. |

## Three things that will bite an editor

**Do not cache a compiled scenario and run it twice.** Running one leaves its switches
flipped (decision D-206), so the second run starts with the first run's chemistry already on
from t = 0 — silently, with no error, and with a completely believable wrong curve. Under a
framework that re-executes its script on every widget change, that is the *default* outcome
unless you guard it. `runner.run_once` compiles fresh inside the cached call every time, and
`tests/test_app.py::test_each_run_starts_from_a_clean_scenario` is what stops it regressing.

**Never test a confidence mark for truthiness.** `Tier.SPECULATIVE` is the enum's zero and is
therefore falsy, so `if tier:` reports the least confident tier as no tier at all. Every check
in this package is `is None`.

**An untouched quantity reads as `validated` if you take the engine at face value.** The tier
combine returns the top tier for an empty list, so anything no chemistry writes comes back
wearing the one word nothing here has earned. `Readout.inert` and `RunResult.touched_variables`
exist for that reason alone.

## How this fits the repo's rules

- **Dependencies.** Streamlit and Plotly live in the `ui` dependency group, so a plain
  `uv sync` still installs a four-package research library. `pytest.importorskip` guards the
  tests that need them.
- **Types.** `app/` is type-checked under the same strict `mypy` settings as `src/` — which is
  why the engine now ships a `py.typed` marker, so its own types are visible from outside the
  package. Streamlit, Plotly and pandas ship no type information and are handled with an
  `ignore_missing_imports` override, exactly as `scipy` already was. Turning this on found two
  real bugs, including one crash on a line that had never run.
- **Lint.** `ruff check .` covers `app/`. `render.py` is the one file exempt from the
  "no `dict()` calls" rule, because Plotly's figure spec is nested keyword objects and its own
  documentation is written that way.
- **Dependency direction.** `app` imports the engine; the engine imports nothing from `app`.

## Theme

The console pins a single light theme in `.streamlit/config.toml` rather than following the
viewer's system setting. The chart palette is drawn from the subject — garnet, must gold,
bottle glass, tank steel, oak, hop, copper, lees — and those are mid-to-dark hues chosen to
read against a warm near-white. Confidence is carried here by line weight and dash pattern,
which is precisely what a thin dotted line loses when the ground flips to near-black.
