# The Fermentation Console

A local app for describing a batch, running it, and seeing how far the answer can be
trusted. It sits above the engine and the engine never knows it exists.

Double-click `start-console.cmd` (Windows) or run `./start-console.sh` (macOS, Linux). Either
one installs what the interface needs, finds a port nothing else is using, and opens the page.
The equivalent typed by hand is:

```bash
uv sync --group ui
uv run streamlit run app/main.py
```

It opens in a browser at `http://localhost:8501`. Nothing is uploaded anywhere; the whole
thing runs on this machine. That includes the framework's own usage telemetry, which
`.streamlit/config.toml` turns off, so the sentence holds however the console was started.

A written report lands in `Fermentation Console` under your `Documents` folder — or under your
home directory where there is no `Documents` — and `FERMENTATION_CONSOLE_REPORTS` moves it. The
location used to be a constant naming a drive that exists on one machine in the world, which
made *Write it* raise on every other one.

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
| `start_console.py` | What the double-clickable wrappers call: dependency check, free port, browser. |

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

The console offers both grounds. `.streamlit/config.toml` carries a `[theme.light]` and a
`[theme.dark]` block, so **Settings ▸ Appearance** under the ⋮ menu switches the whole page,
and the sidebar's *Ink* control decides which palette the charts are drawn in — following
the page by default. It used to pin light, on the argument that the palette is mid-to-dark
hues chosen against a warm near-white and that confidence is carried by line *weight* and
dash pattern, which is what a thin dotted line loses on a bad ground. That argument was
right about the constraint and wrong about the remedy: the fix is a second palette designed
against near-black, not a ground nobody may change. `render.LIGHT` and `render.DARK` are the
same eight subject hues chosen twice, and `TIER_STYLE` is deliberately *not* per-theme — a
ground may change a hue, never the meaning of a dash.

**Switching the page does not redraw the charts on its own.** Streamlit 1.62 applies the new
page CSS without re-running the script, so `st.context.theme.type` is never re-read and the
figures keep the ink they were drawn with — measured in both directions: the page goes light
while the traces stay `#e8798f`, and a manual re-run flips them to `#7c1d3f` immediately. The
server cannot detect this, because the script is not running when it happens. So the charts
catch up at the next interaction, the *Ink* control's help text says so, and pinning Light or
Dark there is itself an interaction and therefore takes effect at once.

## Vertical scale

The sidebar's *Log scale up the side* applies to every chart on the page and to the written
report. A log axis needs two guards that Plotly does not give you, both in `render`:

- **Floored off the data, not the dust.** A run ends at ~1e-8 g/L of sugar. Autoscaled, the
  log axis would span from that dust to 200 g/L and squash the whole ferment into the top
  decade. Each axis is floored `LOG_DECADES` (six) below its own peak.
- **Clamped, not dropped.** Plotly silently drops non-positive points, so a sugar line that
  reaches zero appears to *stop* days before the run ends. Values under the floor are drawn
  along it and the true value rides in `customdata`, so the hover never lies. Which lines
  that affected is printed, by name, in `render.log_scale_panel`.

## The chart with nothing on it

A white must has no anthocyanin and no tannin, so the colour chart is four flat lines and
reads as a broken page. `render.flat_group_panel` says so in words, triggered on the engine's
own notion — every readout on the chart is *inert* — rather than on how the pixels came out,
so a chart where three lines are flat and the fourth is the whole story gets no note. The
sentence about what would fill it lives on the `Group` as `when_empty`, next to its blurb.
Two batches reach the same empty chart for opposite reasons (a white has no pigment; a red
has pigment but no aging step to switch the chemistry on), which is why the note reports
*flat at zero* separately from *flat at the value it started from*.
