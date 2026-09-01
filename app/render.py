"""Finished run in, figures and panels out. No framework, by design.

Every function here takes a :class:`~app.runner.RunResult` (or an
:class:`~app.runner.EnsembleResult`) and returns a Plotly figure or a plain data structure.
Nothing in this module imports Streamlit, and nothing in it knows whether it is drawing into
a live page or a written file. That is the whole reason it exists as its own module: the app
and the report are then the *same* rendering, called twice, rather than two drawings that
drift apart.

**Confidence is drawn, not annotated.** A line's dash pattern is its tier - solid for
validated, dashed for plausible, dotted for speculative - so how much a curve can be trusted
is visible at a glance and cannot be lost when someone screenshots the chart without the
caption. Solid is currently unreachable: nothing in this engine holds ``validated``, because
that tier is reserved for agreement with independent measured time-series the project does
not have. The style is defined anyway, and the legend says it is a standard not yet met. A
scale whose top mark can never be earned is more honest than one quietly rescaled to fit.

**A palette is designed against a ground, so there are two of them.** :data:`LIGHT` and
:data:`DARK` are the same eight hues taken from the subject, chosen twice - once to read on
a warm near-white and once on a warm near-black. Only the ink changes: :data:`TIER_STYLE` is
shared, so the dash pattern and the line weight that carry confidence are the same drawing
under either ground. That is the invariant a theme is not allowed to break.

**A log axis is floored, not truncated.** Half the quantities here legitimately start at
zero, and several end at solver dust of order 1e-8. Left alone, a log axis autoscales itself
to that dust and Plotly silently drops every non-positive point, so a line appears to *stop*
partway through the run. Instead each axis is floored :data:`LOG_DECADES` decades below its
own peak, values under the floor are drawn along it, and the hover text still reports the
true value out of ``customdata``. :func:`log_scale_panel` says out loud which lines that
affected.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import plotly.graph_objects as go

from app.readouts import Group, Readout, format_value
from app.runner import EnsembleResult, RunResult
from fermentation.analysis import SpreadAttribution
from fermentation.core.state import FloatArray
from fermentation.core.tiers import Tier


@dataclass(frozen=True)
class Theme:
    """One ground, and the ink chosen to read against it.

    The palette is categorical and taken from the subject rather than from a chart library:
    garnet, must gold, bottle glass, tank steel, oak, hop leaf, copper, lees. The dark set is
    those same eight lifted until a 1.0 px line still reads - contrast at the *thinnest*
    weight is the constraint, because weight is a confidence channel here and the whole
    argument for pinning one ground was that a thin dotted line is what a bad ground loses.
    """

    name: str
    palette: tuple[str, ...]
    #: Drawn for a quantity no mechanism in the run writes. Visible, and subordinate.
    inert: str
    grid: str
    #: Axis titles, tick labels, legend, annotations.
    axis: str
    #: Chart titles, and the nominal line of a band chart.
    ink: str
    event: str
    band_fill: str
    band_line: str
    bar_down: str
    bar_up: str


LIGHT = Theme(
    name="light",
    palette=(
        "#7c1d3f",
        "#c8862b",
        "#3f6b57",
        "#2b5d7d",
        "#8a5a3b",
        "#6b6f4e",
        "#a33b1f",
        "#4a4458",
    ),
    inert="#9a958c",
    grid="rgba(0,0,0,0.08)",
    axis="#55504a",
    ink="#1c1a17",
    event="rgba(0,0,0,0.28)",
    band_fill="rgba(124,29,63,0.14)",
    band_line="#7c1d3f",
    bar_down="#7c1d3f",
    bar_up="#2b5d7d",
)

DARK = Theme(
    name="dark",
    palette=(
        "#e8798f",
        "#e6b463",
        "#6fc39b",
        "#7bb9e2",
        "#d69f76",
        "#c0c68e",
        "#f48a62",
        "#b5a9d2",
    ),
    inert="#948c81",
    grid="rgba(255,255,255,0.11)",
    axis="#b8b1a6",
    ink="#ece7df",
    event="rgba(255,255,255,0.36)",
    band_fill="rgba(232,121,143,0.20)",
    band_line="#e8798f",
    bar_down="#e8798f",
    bar_up="#7bb9e2",
)

#: Both grounds, by the name the page and the report address them with.
THEMES: Mapping[str, Theme] = {LIGHT.name: LIGHT, DARK.name: DARK}

#: Tier -> (dash pattern, line width). The visual channel that carries confidence. Shared by
#: every theme on purpose: a ground may change the hue, never the meaning.
TIER_STYLE: Mapping[Tier | None, tuple[str, float]] = {
    Tier.VALIDATED: ("solid", 2.6),
    Tier.PLAUSIBLE: ("dash", 2.2),
    Tier.SPECULATIVE: ("dot", 1.8),
    # Deliberately NOT solid. Solid is reserved for "checked against real measured data",
    # which nothing has earned, and a page where the only solid lines are the ones nothing
    # touches would teach exactly the wrong reading of the chart.
    None: ("longdash", 1.0),
}

TIER_WORD: Mapping[Tier | None, str] = {
    Tier.VALIDATED: "validated",
    Tier.PLAUSIBLE: "plausible",
    Tier.SPECULATIVE: "speculative",
    None: "inert",
}

#: What each mark actually means, written for someone meeting it for the first time.
TIER_MEANING: Mapping[Tier | None, str] = {
    Tier.VALIDATED: (
        "Checked against real measurements from a real fermentation. Nothing in this model "
        "has reached that bar yet. The mark exists so you can see it is missing."
    ),
    Tier.PLAUSIBLE: (
        "Every number behind it comes from published research, but the result itself has "
        "never been compared against a real fermentation."
    ),
    Tier.SPECULATIVE: (
        "At least one number behind it is the author's own estimate rather than a "
        "measurement. Trust the shape of the line; do not quote the value."
    ),
    None: (
        "Nothing in this run affects it. It sat at its starting value from beginning to end "
        "- which is not the same as being confirmed correct."
    ),
}

#: How far below a chart's own peak a log axis is allowed to reach. Six decades is far under
#: anything a fermentation reading means, and far over the solver's dust - which is what an
#: unfloored log axis would otherwise scale itself to.
LOG_DECADES = 6.0

#: A run that ends at 1e-8 g/L of sugar has finished, not gone negative. Below this a value
#: is dust and is treated as zero - the same threshold :func:`app.readouts.format_value`
#: applies to the final-value table, kept to one number rather than two.
DUST = 1e-6


def tier_word(tier: Tier | None) -> str:
    """Tier as a word. Takes ``None`` for inert; never test a tier for truthiness."""
    return TIER_WORD[tier]


def _axis_title(readouts: Sequence[Readout], exclude: Sequence[str] = ()) -> str:
    units = []
    for r in readouts:
        if r.key in exclude:
            continue
        if r.unit not in units:
            units.append(r.unit)
    return " / ".join(units)


def _stack(values: Sequence[FloatArray]) -> FloatArray:
    """Every finite number drawn on one axis, in one array."""
    parts = [np.asarray(v, dtype=float).ravel() for v in values if np.asarray(v).size]
    if not parts:
        return np.zeros(0, dtype=float)
    stacked = np.concatenate(parts)
    return np.asarray(stacked[np.isfinite(stacked)], dtype=float)


def log_floor(values: Sequence[FloatArray]) -> float | None:
    """The lowest value a log axis carrying these series should show.

    ``None`` when nothing on the axis is positive: there is no honest log axis for a set of
    lines that never leaves zero, so the caller keeps a straight one rather than drawing an
    empty decade ruler over nothing.
    """
    finite = _stack(values)
    if finite.size == 0:
        return None
    hi = float(finite.max())
    if hi <= 0.0:
        return None
    return float(hi * 10.0**-LOG_DECADES)


def _axis_spec(values: Sequence[FloatArray], *, log_y: bool) -> dict[str, Any]:
    """Type and range for one y-axis.

    Two things Plotly's own autoscale gets wrong for this data, both fixed here:

    * On a log axis it would scale to the solver's dust, leaving the actual ferment as a
      wiggle in the top decade. The floor comes from the axis's own peak instead.
    * For a set of lines all flat at zero it picks -1 to 1, which for a concentration draws
      a range of negative values that cannot exist. Those axes get 0 to 1 instead, and
      :func:`flat_group_panel` says why the chart is empty.

    The order matters. An axis with nothing positive on it cannot be a log axis, and it is
    *also* the flat-at-zero case — so refusing the log scale has to fall through to the
    straight guards rather than return, or turning the log scale on would hand back exactly
    the negative half-axis the second fix exists to remove.
    """
    finite = _stack(values)
    if finite.size == 0:
        return {}
    lo, hi = float(finite.min()), float(finite.max())
    floor = log_floor(values) if log_y else None
    if floor is not None:
        bottom, top = float(np.log10(max(lo, floor))), float(np.log10(hi))
        if top - bottom < 0.5:  # a nearly flat line would otherwise get a hairline axis
            bottom, top = top - 0.5, top + 0.25
        pad = 0.04 * (top - bottom)
        return {"type": "log", "range": [bottom - pad, top + pad]}
    if hi - lo <= DUST and lo >= -DUST:
        return {"range": [0.0, max(hi * 2.0, 1.0)]}
    return {}


def _base_layout(
    fig: go.Figure,
    title: str,
    y_title: str,
    y2_title: str = "",
    *,
    theme: Theme = LIGHT,
    y_spec: dict[str, Any] | None = None,
    y2_spec: dict[str, Any] | None = None,
) -> None:
    fig.update_layout(
        title=dict(text=title, x=0.0, xanchor="left", font=dict(size=17, color=theme.ink)),
        font=dict(color=theme.axis),
        xaxis=dict(
            title="Days since pitch",
            showgrid=True,
            gridcolor=theme.grid,
            zerolinecolor=theme.grid,
        ),
        yaxis=dict(
            title=y_title,
            showgrid=True,
            gridcolor=theme.grid,
            zerolinecolor=theme.grid,
            **(y_spec or {}),
        ),
        legend=dict(orientation="h", yanchor="top", y=-0.24, xanchor="left", x=0.0),
        margin=dict(l=60, r=60, t=52, b=124),
        hovermode="x unified",
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        height=440,
    )
    if y2_title:
        fig.update_layout(
            yaxis2=dict(
                title=y2_title,
                overlaying="y",
                side="right",
                showgrid=False,
                **(y2_spec or {}),
            )
        )


def series_figure(
    result: RunResult,
    group: Group,
    readouts: Sequence[Readout],
    *,
    theme: Theme = LIGHT,
    log_y: bool = False,
) -> go.Figure:
    """One chart: every readout in a group, drawn with its tier as its dash pattern."""
    fig = go.Figure()
    days = result.days
    values = {r.key: r.series(result) for r in readouts}
    left = [r for r in readouts if r.key not in group.secondary]
    right = [r for r in readouts if r.key in group.secondary]
    on_left = left or list(readouts)
    floors: dict[str, float | None] = {
        "y": log_floor([values[r.key] for r in on_left]) if log_y else None,
        "y2": log_floor([values[r.key] for r in right]) if (log_y and right) else None,
    }

    for i, r in enumerate(readouts):
        tier = r.tier(result)
        dash, width = TIER_STYLE[tier]
        color = theme.inert if tier is None else theme.palette[i % len(theme.palette)]
        axis = "y2" if r.key in group.secondary else "y"
        raw = values[r.key]
        floor = floors[axis]
        # Drawn at the floor rather than dropped, so a line that reaches zero does not appear
        # to stop partway through the run. customdata keeps the hover honest about the value.
        drawn = np.maximum(raw, floor) if floor is not None else raw
        fig.add_trace(
            go.Scatter(
                x=days,
                y=drawn,
                customdata=raw,
                name=f"{r.label} ({r.unit})",
                mode="lines",
                line=dict(color=color, dash=dash, width=width),
                yaxis=axis,
                hovertemplate=(
                    f"<b>{r.label}</b><br>day %{{x:.2f}}<br>%{{customdata:.4g}} {r.unit}"
                    f"<br><i>{tier_word(tier)}</i><extra></extra>"
                ),
            )
        )
    _base_layout(
        fig,
        group.title,
        _axis_title(on_left),
        _axis_title(right) if right else "",
        theme=theme,
        y_spec=_axis_spec([values[r.key] for r in on_left], log_y=log_y),
        y2_spec=_axis_spec([values[r.key] for r in right], log_y=log_y) if right else None,
    )
    _mark_events(fig, result, theme)
    return fig


def _mark_events(fig: go.Figure, result: RunResult, theme: Theme = LIGHT) -> None:
    """Vertical rules where an intervention fired, so a jump has a visible cause."""
    from app.library import VERB_SPECS

    for iv in result.scenario.interventions:
        if iv.day > float(result.days[-1]):
            continue
        spec = VERB_SPECS.get(iv.action)
        fig.add_vline(
            x=iv.day,
            line=dict(color=theme.event, width=1, dash="dot"),
            annotation_text=spec.label if spec else iv.action,
            annotation_position="top",
            annotation_font_size=10,
            annotation_font_color=theme.axis,
        )


def band_figure(
    ens: EnsembleResult,
    readout: Readout,
    nominal: RunResult,
    *,
    theme: Theme = LIGHT,
    log_y: bool = False,
) -> go.Figure:
    """The run itself, plus the range the re-runs covered between the 5th and 95th."""
    fig = go.Figure()
    days = ens.days
    variable = readout.sources[0]
    band = ens.ensemble.band(variable)
    scale = _band_scale(readout, nominal)

    high = np.asarray(np.asarray(band.high) * scale, dtype=float)
    low = np.asarray(np.asarray(band.low) * scale, dtype=float)
    mid = np.asarray(np.asarray(band.median) * scale, dtype=float)
    nom = readout.series(nominal)
    floor = log_floor([high, mid, nom]) if log_y else None

    def draw(series: FloatArray) -> FloatArray:
        return np.maximum(series, floor) if floor is not None else series

    fig.add_trace(
        go.Scatter(
            x=np.concatenate([days, days[::-1]]),
            y=np.concatenate([draw(high), draw(low)[::-1]]),
            fill="toself",
            fillcolor=theme.band_fill,
            line=dict(color="rgba(0,0,0,0)"),
            name="middle 90% of the re-runs",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=days,
            y=draw(mid),
            customdata=mid,
            name="middle of the re-runs",
            mode="lines",
            line=dict(color=theme.band_line, dash="dash", width=1.8),
            hovertemplate=f"%{{customdata:.4g}} {readout.unit}<extra>middle</extra>",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=days,
            y=draw(nom),
            customdata=nom,
            name="using the published values",
            mode="lines",
            line=dict(color=theme.ink, width=2.4),
            hovertemplate=f"%{{customdata:.4g}} {readout.unit}<extra>published</extra>",
        )
    )
    _base_layout(
        fig,
        f"{readout.label} - how far the published ranges alone move it",
        readout.unit,
        theme=theme,
        y_spec=_axis_spec([high, low, mid, nom], log_y=log_y),
    )
    return fig


def _band_scale(readout: Readout, nominal: RunResult) -> float:
    """Unit factor between the raw state variable and the readout's displayed unit.

    Recovered by comparing the readout's own series against the raw state it is built on,
    rather than by keeping a second copy of the conversion. A readout whose relationship to
    its state variable is not a plain scaling (pH, alcohol by volume) has no band drawn -
    :func:`bandable` refuses those rather than scaling something that is not scalable.
    """
    raw = np.atleast_2d(np.asarray(nominal.traj.series(readout.sources[0]), dtype=float))[0]
    shown = readout.series(nominal)
    nonzero = np.abs(raw) > 0
    if not nonzero.any():
        return 1.0
    return float(np.median(shown[nonzero] / raw[nonzero]))


def bandable(readout: Readout, result: RunResult) -> bool:
    """Can a band be drawn honestly for this readout?

    Only for a readout that is one tracked quantity times a fixed factor. pH is worked out
    from all the acids at once, so a range computed on any single one of them and then
    relabelled "pH" would be a different quantity wearing pH's name. Better no band than a
    wrong one.
    """
    if len(readout.sources) != 1 or readout.sources[0] not in result.schema:
        return False
    raw = np.atleast_2d(np.asarray(result.traj.series(readout.sources[0]), dtype=float))
    if raw.shape[0] != 1:
        return False
    shown = readout.series(result)
    nonzero = np.abs(raw[0]) > 1e-12
    if nonzero.sum() < 3:
        return False
    ratios = shown[nonzero] / raw[0][nonzero]
    return bool(np.allclose(ratios, ratios[0], rtol=1e-6))


def compare_figure(
    results: Sequence[tuple[str, RunResult]],
    readout: Readout,
    *,
    theme: Theme = LIGHT,
    log_y: bool = False,
) -> go.Figure:
    """One quantity drawn across several saved runs, so the difference is the subject."""
    fig = go.Figure()
    drawable = [(label, r) for label, r in results if readout.available(r)]
    series = [readout.series(r) for _, r in drawable]
    floor = log_floor(series) if log_y else None
    for i, ((label, result), raw) in enumerate(zip(drawable, series, strict=True)):
        tier = readout.tier(result)
        dash, width = TIER_STYLE[tier]
        fig.add_trace(
            go.Scatter(
                x=result.days,
                y=np.maximum(raw, floor) if floor is not None else raw,
                customdata=raw,
                name=label,
                mode="lines",
                line=dict(color=theme.palette[i % len(theme.palette)], dash=dash, width=width),
                hovertemplate=f"%{{customdata:.4g}} {readout.unit}<extra>{label}</extra>",
            )
        )
    _base_layout(
        fig,
        f"{readout.label} - across runs",
        readout.unit,
        theme=theme,
        y_spec=_axis_spec(series, log_y=log_y),
    )
    return fig


def spread_figure(
    attribution: SpreadAttribution, top: int = 12, *, theme: Theme = LIGHT
) -> go.Figure:
    """Which of the model's numbers the uncertainty comes from, biggest share first."""
    ranked = attribution.ranked()[:top]
    names = [n for n, _ in ranked][::-1]
    shares = [v for _, v in ranked][::-1]
    signed = attribution.per_param_signed
    colors = [theme.bar_down if signed.get(n, 0.0) < 0 else theme.bar_up for n in names]
    fig = go.Figure(
        go.Bar(
            x=shares,
            y=names,
            orientation="h",
            marker=dict(color=colors),
            hovertemplate="%{y}<br>%{x:.1%} of the uncertainty<extra></extra>",
        )
    )
    fig.update_layout(
        title=dict(
            text=(
                f"What makes {attribution.variable} uncertain "
                f"({attribution.n_members} re-runs; {attribution.unexplained:.0%} of the "
                "movement this simple fit cannot account for)"
            ),
            x=0.0,
            xanchor="left",
            font=dict(size=17, color=theme.ink),
        ),
        font=dict(color=theme.axis),
        xaxis=dict(title="share of the uncertainty", tickformat=".0%", gridcolor=theme.grid),
        margin=dict(l=10, r=20, t=70, b=48),
        height=max(320, 30 * len(names) + 120),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


@dataclass(frozen=True)
class Panel:
    """A block of text the app and the report both display, with a severity."""

    kind: str  # "note" | "warning" | "caveat"
    title: str
    body: str


def _has(n: int) -> str:
    """ "has"/"have" for a count, so a generated sentence does not read as broken English."""
    return "has" if n == 1 else "have"


def honesty_panels(result: RunResult) -> list[Panel]:
    """Everything that has to be said about this run before anyone quotes a number off it."""
    panels: list[Panel] = []
    if not result.traj.success:
        panels.append(
            Panel(
                "warning",
                "The calculation did not finish",
                f"{result.traj.message} Nothing below describes a completed fermentation - "
                "the run stopped partway and the charts end where it gave up.",
            )
        )
    # Counted only over quantities something in the run actually changes. A quantity nothing
    # touches comes back from the engine reading "validated" - the top mark is what you get
    # from combining an empty list - so counting those would let a run full of untouched
    # slots claim a confidence it never earned.
    touched = result.touched_variables()
    tiers = [t for n, t in result.traj.tier_map.items() if n in touched]
    counts = {t: sum(1 for x in tiers if x is t) for t in Tier}
    validated_names = sorted(
        n for n, t in result.traj.tier_map.items() if n in touched and t is Tier.VALIDATED
    )
    estimated, published, checked = (
        counts[Tier.SPECULATIVE],
        counts[Tier.PLAUSIBLE],
        counts[Tier.VALIDATED],
    )
    body = (
        f"Of the {len(tiers)} quantities this run works out: {estimated} rest on at least "
        f"one number the author estimated rather than measured; {published} come entirely "
        f"from published research but have never been compared with a real fermentation; "
        f"{checked} {_has(checked)} been checked against real measured data"
    )
    if validated_names:
        body += (
            " (namely "
            + ", ".join(validated_names)
            + ", which you set yourself rather than the model working out)"
        )
    body += (
        ". That last standard is the one this project is aiming at and has not met for any "
        "chemistry yet, so read the charts below as an informed model rather than as a "
        "measurement."
    )
    panels.append(Panel("note", "How much of this can be trusted", body))

    switched = [m for m in result.mechanisms if m.switched_on_mid_run]
    if switched:
        names = ", ".join(m.name for m in switched[:6])
        more = f" and {len(switched) - 6} more" if len(switched) > 6 else ""
        panels.append(
            Panel(
                "note",
                f"{len(switched)} chemical steps only started partway through the run",
                f"Something you scheduled switched them on: {names}{more}. If a line is flat "
                "at first and then starts moving, that is usually why - the chemistry behind "
                "it was not running yet.",
            )
        )
    return panels


def caveat_panels(readouts: Sequence[Readout]) -> list[Panel]:
    """The warnings attached to these readouts, with duplicates removed.

    They are shown beside the chart, never tucked into a footnote. The whole reason a warning
    is attached to a readout rather than written into a document is that it then cannot be
    separated from the line it is about.
    """
    seen: set[str] = set()
    out: list[Panel] = []
    for r in readouts:
        if r.caveat and r.caveat not in seen:
            seen.add(r.caveat)
            out.append(Panel("caveat", r.label, r.caveat))
    return out


def flat_group_panel(result: RunResult, group: Group, readouts: Sequence[Readout]) -> Panel | None:
    """A note for a chart on which nothing in this run moves anything.

    A white wine carries no anthocyanin and no tannin, so the colour chart is four lines
    sitting on zero, and it reads as a broken page rather than as an answer. The trigger is
    the engine's own notion rather than a look at the pixels: every readout drawn here is
    *inert*, meaning no mechanism in this run writes to any variable behind it. That is a
    fact about the batch, so the note says which batch would make the chart move - the group
    supplies that sentence, because it is presentation and belongs beside the group's blurb.

    ``None`` as soon as one readout is driven: a chart where three lines are flat and the
    fourth is the whole story is not an empty chart, and a note there would be noise.
    """
    live = list(readouts)
    if not live or any(r.tier(result) is not None for r in live):
        return None
    at_zero = all(bool(np.all(np.abs(r.series(result)) <= DUST)) for r in live)
    where = "flat at zero" if at_zero else "flat at the value it started from"
    body = (
        f"No chemistry in this run writes to any of these, so every line here is {where}. "
        "That is not a failed calculation and it is not a measured zero - the chemistry that "
        "would move them is simply not part of this batch."
    )
    if group.when_empty:
        body += " " + group.when_empty
    return Panel("note", f"Nothing in this batch moves {group.title.lower()}", body)


def log_scale_panel(result: RunResult, readouts: Sequence[Readout]) -> Panel | None:
    """What a log scale is not showing, named rather than left to be discovered.

    Returned only when something actually reaches zero. A panel that appeared every time the
    log scale was on would train the reader to skip it, and it is wanted on the run where a
    line really does spend half the chart pinned to the floor.
    """
    seen: set[str] = set()
    hidden: list[str] = []
    for r in readouts:
        if r.label in seen:
            continue
        seen.add(r.label)
        series = r.series(result)
        if series.size and float(np.nanmin(series)) <= DUST:
            hidden.append(r.label)
    if not hidden:
        return None
    names = ", ".join(hidden[:8]) + (f", and {len(hidden) - 8} more" if len(hidden) > 8 else "")
    return Panel(
        "caveat",
        "What the log scale is not showing you",
        "A log scale has no zero, and these reach zero - or the solver's dust around it - at "
        f"some point in the run: {names}. Those stretches are drawn along the bottom of their "
        "axis rather than dropped, so no line appears to stop partway through. But the floor "
        "of a log axis here means 'under anything this chart can show', which is "
        f"{LOG_DECADES:.0f} decades below that chart's own peak and not a reading. Hovering "
        "still reports the true value.",
    )


def final_state_rows(result: RunResult, limit: int = 0) -> list[tuple[str, str, str]]:
    """Every quantity's final value as ``(name, value, confidence word)``.

    The raw view underneath the charts: everything the model tracks, so nothing is invisible
    just because nobody wrote a chart for it.
    """
    touched = result.touched_variables()
    rows: list[tuple[str, str, str]] = []
    final = result.traj.final()
    for name in sorted(result.schema.names):
        value = final[name]
        tier = result.traj.tier_map.get(name) if name in touched else None
        if isinstance(value, np.ndarray):
            text = ", ".join(format_value(float(v)) for v in np.atleast_1d(value))
        else:
            text = format_value(float(value))
        rows.append((name, text, tier_word(tier)))
    return rows[:limit] if limit else rows
