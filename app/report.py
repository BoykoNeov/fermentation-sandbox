"""One finished run, written to a single self-contained HTML file.

The same figures and panels the live app shows, from the same functions in
:mod:`app.render` — this module only decides where they land on a page. That is the whole
point of keeping the rendering framework-free: the report is not a second implementation
that will drift, it is the first one called with a file at the end of it.

The output is one file with the charts and their JavaScript inlined, so it opens with no
server and survives being emailed, attached to a decision record, or opened in two years.
"""

from __future__ import annotations

import datetime as dt
import html
from collections.abc import Sequence
from pathlib import Path

import plotly.graph_objects as go

from app import provenance, readouts, render
from app.runner import ConvergenceCheck, RunResult
from fermentation.core.tiers import Tier

_CSS = """
:root {
  color-scheme: light dark;
  --ink: #1c1a17;
  --ink-soft: #55504a;
  --paper: #ffffff;
  --panel: #f4f1ec;
  --rule: #ddd8d0;
  --garnet: #7c1d3f;
  --warn-bg: #fdf1e7;
  --warn-edge: #c8862b;
  --stop-bg: #fbeaea;
  --stop-edge: #a33b1f;
}
@media (prefers-color-scheme: dark) {
  :root {
    --ink: #ece7df;
    --ink-soft: #a49d93;
    --paper: #161513;
    --panel: #201e1b;
    --rule: #37342f;
    --garnet: #d4718f;
    --warn-bg: #2a2118;
    --warn-edge: #c8862b;
    --stop-bg: #2c1c19;
    --stop-edge: #d2603f;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--paper);
  color: var(--ink);
  font: 15px/1.6 "Iowan Old Style", "Charter", Georgia, "Times New Roman", serif;
}
main { max-width: 1080px; margin: 0 auto; padding: 48px 28px 96px; }
h1 { font-size: 30px; line-height: 1.15; margin: 0 0 6px; text-wrap: balance; }
h2 {
  font-size: 12px; letter-spacing: 0.12em; text-transform: uppercase;
  color: var(--ink-soft); margin: 56px 0 14px;
  border-top: 1px solid var(--rule); padding-top: 14px;
  font-family: ui-sans-serif, system-ui, "Segoe UI", sans-serif;
}
h3 { font-size: 17px; margin: 26px 0 6px; }
p { margin: 0 0 12px; max-width: 68ch; }
.sub { color: var(--ink-soft); margin-bottom: 28px; }
.meta { font-family: ui-monospace, "Cascadia Mono", Consolas, monospace; font-size: 12.5px; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; }
.card { background: var(--panel); border: 1px solid var(--rule); border-radius: 3px;
  padding: 12px 14px; }
.card .k { font-size: 11px; letter-spacing: 0.08em; text-transform: uppercase;
  color: var(--ink-soft); font-family: ui-sans-serif, system-ui, sans-serif; }
.card .v { font-size: 24px; font-variant-numeric: tabular-nums; margin: 4px 0 2px; }
.card .u { font-size: 12px; color: var(--ink-soft); }
.panel { border-left: 3px solid var(--rule); background: var(--panel); padding: 12px 16px;
  margin: 14px 0; border-radius: 0 3px 3px 0; }
.panel.caveat { border-left-color: var(--warn-edge); background: var(--warn-bg); }
.panel.warning { border-left-color: var(--stop-edge); background: var(--stop-bg); }
.panel b { display: block; margin-bottom: 3px; }
.tier { display: inline-block; font-family: ui-sans-serif, system-ui, sans-serif;
  font-size: 10.5px;
  letter-spacing: 0.06em; text-transform: uppercase; padding: 2px 7px; border-radius: 2px;
  border: 1px solid currentColor; }
.tier.validated { color: var(--ink-soft); border-style: dashed; opacity: 0.72; }
.tier.plausible { color: #2b5d7d; }
.tier.speculative { color: var(--garnet); }
.tier.inert { color: var(--ink-soft); opacity: 0.7; }
table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
th, td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--rule);
  vertical-align: top; }
th { font-family: ui-sans-serif, system-ui, sans-serif; font-size: 11px; letter-spacing: 0.08em;
  text-transform: uppercase; color: var(--ink-soft); }
td.num { font-variant-numeric: tabular-nums; text-align: right; }
.scroll { overflow-x: auto; }
.fig { margin: 20px 0 8px; }
footer { margin-top: 64px; padding-top: 16px; border-top: 1px solid var(--rule);
  color: var(--ink-soft); font-size: 12.5px; }
"""


def _esc(text: object) -> str:
    return html.escape(str(text))


def _tier_badge(tier: Tier | None) -> str:
    word = render.tier_word(tier)
    return f'<span class="tier {word}">{word}</span>'


def _fig_html(fig: go.Figure, first: bool) -> str:
    html_fragment: str = fig.to_html(
        full_html=False,
        include_plotlyjs=first,
        config={"displaylogo": False, "responsive": True},
    )
    return html_fragment


def build(
    result: RunResult,
    *,
    convergence: ConvergenceCheck | None = None,
    explain: Sequence[str] = (),
) -> str:
    """Render one run as a complete HTML document."""
    sc = result.scenario
    parts: list[str] = []
    a = parts.append

    a(f"<h1>{_esc(sc.name)}</h1>")
    a(
        f'<p class="sub">{_esc(sc.medium)} &middot; {sc.duration_days:g} days &middot; '
        f"{len(sc.interventions)} intervention(s) &middot; computed at "
        f"{_esc(result.fidelity.label)} &middot; "
        f"{dt.datetime.now().strftime('%Y-%m-%d %H:%M')}</p>"
    )

    # -- what it finished at -------------------------------------------------------------
    a("<h2>Where it ended up</h2>")
    a('<div class="cards">')
    for label, value, unit, tier in readouts.summary(result):
        a(
            f'<div class="card"><div class="k">{_esc(label)}</div>'
            f'<div class="v">{_esc(value)}</div>'
            f'<div class="u">{_esc(unit)} &nbsp; {_tier_badge(tier)}</div></div>'
        )
    a("</div>")

    # -- the honesty block ---------------------------------------------------------------
    for panel in render.honesty_panels(result):
        a(f'<div class="panel {panel.kind}"><b>{_esc(panel.title)}</b>{_esc(panel.body)}</div>')
    if convergence is not None:
        a(
            f'<div class="panel"><b>Would a more careful calculation have changed it?</b>'
            "Worked out again at the "
            f"'{_esc(convergence.tighter.precision)}' setting, the biggest change in any "
            f"headline quantity was {convergence.worst:.2e}, in "
            f"{_esc(convergence.worst_variable)}. {_esc(convergence.verdict)}</div>"
        )

    # -- the charts ----------------------------------------------------------------------
    a("<h2>How it got there</h2>")
    first = True
    all_drawn: list[readouts.Readout] = []
    for group, live in readouts.groups_for(result):
        fig = render.series_figure(result, group, live)
        a(f'<div class="fig">{_fig_html(fig, first)}</div>')
        first = False
        if group.blurb:
            a(f"<p>{_esc(group.blurb)}</p>")
        all_drawn.extend(live)
    for panel in render.caveat_panels(all_drawn):
        a(f'<div class="panel caveat"><b>{_esc(panel.title)}</b>{_esc(panel.body)}</div>')

    # -- the provenance walk -------------------------------------------------------------
    if explain:
        a("<h2>Where these numbers come from</h2>")
        for variable in explain:
            if variable not in result.schema:
                continue
            a(f"<h3>{_esc(variable)}</h3>")
            a(f"<p>{_esc(provenance.why_this_tier(result, variable))}</p>")
            mechs = provenance.mechanisms_for(result, variable)
            cards = provenance.constants_for(result, mechs)
            a('<div class="scroll"><table><thead><tr>')
            a(
                "<th>Number</th><th>Value</th><th>Could be anywhere from</th>"
                "<th>Confidence</th><th>Source, and what it was measured under</th>"
                "</tr></thead><tbody>"
            )
            for c in cards[:40]:
                doi = f'<br><span class="meta">doi:{_esc(c.doi)}</span>' if c.doi else ""
                a(
                    f'<tr><td class="meta">{_esc(c.name)}</td>'
                    f'<td class="num">{c.value:,.4g} {_esc(c.unit)}</td>'
                    f'<td class="num">{c.low:,.4g} &ndash; {c.high:,.4g}</td>'
                    f"<td>{_tier_badge(c.tier)}</td>"
                    f"<td>{_esc(c.source)}<br><span class='meta'>{_esc(c.conditions)}</span>"
                    f"{doi}</td></tr>"
                )
            a("</tbody></table></div>")
            if len(cards) > 40:
                a(f"<p>{len(cards) - 40} further numbers not listed.</p>")

    # -- the scenario, as run -------------------------------------------------------------
    a("<h2>The batch, exactly as the model received it</h2>")
    a(f'<div class="scroll"><pre class="meta">{_esc(sc.model_dump_json(indent=2))}</pre></div>')

    # -- the whole final state -------------------------------------------------------------
    a("<h2>Everything the model tracks, at the end</h2>")
    a('<div class="scroll"><table><thead><tr><th>What</th><th>Value at the end</th>')
    a("<th>Confidence</th></tr></thead><tbody>")
    for name, value, word in render.final_state_rows(result):
        a(
            f'<tr><td class="meta">{_esc(name)}</td><td class="num">{_esc(value)}</td>'
            f'<td><span class="tier {word}">{word}</span></td></tr>'
        )
    a("</tbody></table></div>")

    a(
        "<footer>Written by the Fermentation Console, from the fermentation-sandbox engine. "
        f"Settings: {_esc(result.fidelity.method)} solver, tolerances "
        f"{result.fidelity.rtol:g}/{result.fidelity.atol:g}, {result.fidelity.points} stored "
        f"points, &lsquo;{_esc(result.fidelity.oxidative)}&rsquo; oxygen chemistry. Set up in "
        f"{result.compile_seconds:.2f} s, worked out in {result.wall_seconds:.2f} s.</footer>"
    )

    body = "\n".join(parts)
    return (
        "<!doctype html>\n<html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>{_esc(sc.name)}</title><style>{_CSS}</style></head>"
        f"<body><main>{body}</main></body></html>"
    )


def write(
    result: RunResult,
    path: str | Path,
    *,
    convergence: ConvergenceCheck | None = None,
    explain: Sequence[str] = (),
) -> Path:
    """Write the report and return the path it landed on."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build(result, convergence=convergence, explain=explain), encoding="utf-8")
    return out
