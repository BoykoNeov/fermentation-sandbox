"""Finished run in, figures and panels out. No framework, by design.

Every function here takes a :class:`~app.runner.RunResult` (or an
:class:`~app.runner.EnsembleResult`) and returns a Plotly figure or a plain data structure.
Nothing in this module imports Streamlit, and nothing in it knows whether it is drawing into
a live page or a written file. That is the whole reason it exists as its own module: the app
and the report are then the *same* rendering, called twice, rather than two drawings that
drift apart.

**Confidence is drawn, not annotated.** A line's dash pattern is its tier — solid for
validated, dashed for plausible, dotted for speculative — so how much a curve can be trusted
is visible at a glance and cannot be lost when someone screenshots the chart without the
caption. Solid is currently unreachable: nothing in this engine holds ``validated``, because
that tier is reserved for agreement with independent measured time-series the project does
not have. The style is defined anyway, and the legend says it is a standard not yet met. A
scale whose top mark can never be earned is more honest than one quietly rescaled to fit.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go

from app.readouts import Group, Readout, format_value
from app.runner import EnsembleResult, RunResult
from fermentation.analysis import SpreadAttribution
from fermentation.core.tiers import Tier

#: Categorical palette, taken from the subject rather than from a chart library: garnet,
#: must gold, bottle glass, tank steel, oak, hop leaf, copper, lees.
PALETTE: tuple[str, ...] = (
    "#7c1d3f",
    "#c8862b",
    "#3f6b57",
    "#2b5d7d",
    "#8a5a3b",
    "#6b6f4e",
    "#a33b1f",
    "#4a4458",
)

INERT_COLOR = "#9a958c"

#: Tier → (dash pattern, line width). The visual channel that carries confidence.
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


def _base_layout(fig: go.Figure, title: str, y_title: str, y2_title: str = "") -> None:
    fig.update_layout(
        title=dict(text=title, x=0.0, xanchor="left", font=dict(size=17)),
        xaxis=dict(title="Days since pitch", showgrid=True, gridcolor="rgba(0,0,0,0.06)"),
        yaxis=dict(title=y_title, showgrid=True, gridcolor="rgba(0,0,0,0.06)"),
        legend=dict(orientation="h", yanchor="top", y=-0.24, xanchor="left", x=0.0),
        margin=dict(l=60, r=60, t=52, b=124),
        hovermode="x unified",
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        height=440,
    )
    if y2_title:
        fig.update_layout(yaxis2=dict(title=y2_title, overlaying="y", side="right", showgrid=False))


def series_figure(result: RunResult, group: Group, readouts: Sequence[Readout]) -> go.Figure:
    """One chart: every readout in a group, drawn with its tier as its dash pattern."""
    fig = go.Figure()
    days = result.days
    for i, r in enumerate(readouts):
        tier = r.tier(result)
        dash, width = TIER_STYLE[tier]
        color = INERT_COLOR if tier is None else PALETTE[i % len(PALETTE)]
        on_right = r.key in group.secondary
        fig.add_trace(
            go.Scatter(
                x=days,
                y=r.series(result),
                name=f"{r.label} ({r.unit})",
                mode="lines",
                line=dict(color=color, dash=dash, width=width),
                yaxis="y2" if on_right else "y",
                hovertemplate=(
                    f"<b>{r.label}</b><br>day %{{x:.2f}}<br>%{{y:.4g}} {r.unit}"
                    f"<br><i>{tier_word(tier)}</i><extra></extra>"
                ),
            )
        )
    left_axis = [r for r in readouts if r.key not in group.secondary]
    right_axis = [r for r in readouts if r.key in group.secondary]
    _base_layout(
        fig,
        group.title,
        _axis_title(left_axis or readouts),
        _axis_title(right_axis) if right_axis else "",
    )
    _mark_events(fig, result)
    return fig


def _mark_events(fig: go.Figure, result: RunResult) -> None:
    """Vertical rules where an intervention fired, so a jump has a visible cause."""
    from app.library import VERB_SPECS

    for iv in result.scenario.interventions:
        if iv.day > float(result.days[-1]):
            continue
        spec = VERB_SPECS.get(iv.action)
        fig.add_vline(
            x=iv.day,
            line=dict(color="rgba(0,0,0,0.28)", width=1, dash="dot"),
            annotation_text=spec.label if spec else iv.action,
            annotation_position="top",
            annotation_font_size=10,
        )


def band_figure(ens: EnsembleResult, readout: Readout, nominal: RunResult) -> go.Figure:
    """The run itself, plus the range the re-runs covered between the 5th and 95th."""
    fig = go.Figure()
    days = ens.days
    variable = readout.sources[0]
    band = ens.ensemble.band(variable)
    scale = _band_scale(readout, nominal)

    fig.add_trace(
        go.Scatter(
            x=np.concatenate([days, days[::-1]]),
            y=np.concatenate([np.asarray(band.high) * scale, (np.asarray(band.low) * scale)[::-1]]),
            fill="toself",
            fillcolor="rgba(124,29,63,0.14)",
            line=dict(color="rgba(0,0,0,0)"),
            name="middle 90% of the re-runs",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=days,
            y=np.asarray(band.median) * scale,
            name="middle of the re-runs",
            mode="lines",
            line=dict(color="#7c1d3f", dash="dash", width=1.8),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=days,
            y=readout.series(nominal),
            name="using the published values",
            mode="lines",
            line=dict(color="#1c1a17", width=2.4),
        )
    )
    _base_layout(fig, f"{readout.label} - how far the published ranges alone move it", readout.unit)
    return fig


def _band_scale(readout: Readout, nominal: RunResult) -> float:
    """Unit factor between the raw state variable and the readout's displayed unit.

    Recovered by comparing the readout's own series against the raw state it is built on,
    rather than by keeping a second copy of the conversion. A readout whose relationship to
    its state variable is not a plain scaling (pH, alcohol by volume) has no band drawn —
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


def compare_figure(results: Sequence[tuple[str, RunResult]], readout: Readout) -> go.Figure:
    """One quantity drawn across several saved runs, so the difference is the subject."""
    fig = go.Figure()
    for i, (label, result) in enumerate(results):
        if not readout.available(result):
            continue
        tier = readout.tier(result)
        dash, width = TIER_STYLE[tier]
        fig.add_trace(
            go.Scatter(
                x=result.days,
                y=readout.series(result),
                name=label,
                mode="lines",
                line=dict(color=PALETTE[i % len(PALETTE)], dash=dash, width=width),
            )
        )
    _base_layout(fig, f"{readout.label} — across runs", readout.unit)
    return fig


def spread_figure(attribution: SpreadAttribution, top: int = 12) -> go.Figure:
    """Which of the model's numbers the uncertainty comes from, biggest share first."""
    ranked = attribution.ranked()[:top]
    names = [n for n, _ in ranked][::-1]
    shares = [v for _, v in ranked][::-1]
    signed = attribution.per_param_signed
    colors = ["#7c1d3f" if signed.get(n, 0.0) < 0 else "#2b5d7d" for n in names]
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
            font=dict(size=17),
        ),
        xaxis=dict(title="share of the uncertainty", tickformat=".0%"),
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
