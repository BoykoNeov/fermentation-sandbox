"""The Fermentation Console — describe a batch, run it, and see how far to trust the result.

Start it with::

    uv sync --group ui
    uv run streamlit run app/main.py

Everything drawn on the page comes from :mod:`app.render`, which knows nothing about
Streamlit; this file is the form and the layout around it. That division matters because
Streamlit re-runs this entire script on every widget change, and the only safe thing to do
under that model is keep the expensive, stateful work behind a cache and the drawing pure.

**Two rules this file lives by.**

*Set the scenario up fresh every time, and cache the finished result.* Running a prepared
scenario leaves its switches flipped (decision D-206), so running the same prepared object
twice starts the second run with the first one's chemistry already switched on — silently,
with no error, and with a plausible-looking wrong answer. The cache therefore holds finished,
inert results keyed on the batch and the settings, never a prepared scenario.

*Never test a confidence mark for truthiness.* The weakest mark is the enum's zero and is
therefore falsy, so ``if tier:`` quietly reports "least confident" as "no confidence mark at
all". Every check here is ``is None``.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Streamlit runs this file as a script, so the repo root is not on the path — only app/ is.
if str(Path(__file__).resolve().parents[1]) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from app import provenance, readouts, render, report  # noqa: E402
from app.fidelity import (  # noqa: E402
    METHODS,
    OXIDATIVE_BLURB,
    POINT_CHOICES,
    PRECISION_BLURB,
    PRECISION_ORDER,
    PRECISION_PRESETS,
    Fidelity,
)
from app.library import (  # noqa: E402
    CLOSURES,
    DEFAULT_INPUTS,
    MEDIA,
    PRIMARY_INPUTS,
    STARTERS,
    VERB_SPECS,
    VERBS,
    allowed_initial_keys,
    gate_problems,
    input_label,
    input_unit,
)
from app.runner import (  # noqa: E402
    EnsembleResult,
    RunResult,
    check_convergence,
    run_once,
    run_uncertainty,
    varying_constants,
)
from fermentation.analysis import attribute_spread  # noqa: E402
from fermentation.scenario import Intervention, Scenario, TemperaturePoint  # noqa: E402

REPORT_DIR = Path(r"M:\claud_projects\temp\ferm-ui")

#: Sampling strategies, named for what they do rather than for their acronyms.
SAMPLERS: dict[str, str] = {
    "lhs": "spread evenly over the ranges",
    "mc": "plain random draws",
    "sobol": "low-discrepancy (best coverage, wants a power-of-two count)",
}

#: How each oxygen chemistry reads on the control itself. The longer explanation is the
#: blurb underneath.
OXIDATIVE_LABEL: dict[str, str] = {
    "direct": "straight to the reactions (default)",
    "cascade": "activated by iron first",
    "direct_burst": "default, plus an early antioxidant",
}

st.set_page_config(
    page_title="Fermentation Console",
    page_icon="🍷",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
<style>
  .stMainBlockContainer { padding-top: 2.2rem; max-width: 1500px; }
  h1, h2, h3 { letter-spacing: -0.01em; }
  div[data-testid="stMetricValue"] { font-variant-numeric: tabular-nums; }
  .tierbadge { display:inline-block; font-size:10.5px; letter-spacing:.06em;
    text-transform:uppercase; padding:1px 7px; border-radius:2px; border:1px solid currentColor; }
  .tierbadge.validated { opacity:.55; border-style:dashed; }
  .tierbadge.plausible { color:#2b5d7d; }
  .tierbadge.speculative { color:#a03050; }
  .tierbadge.inert { opacity:.5; }
</style>
""",
    unsafe_allow_html=True,
)


# -- cached computation --------------------------------------------------------------------
#
# Keyed on the batch as JSON and the settings object. Do NOT be tempted to cache a prepared
# scenario here and run it twice: the second run would begin with the first one's chemistry
# already switched on, with no error and a believable wrong curve.


@st.cache_data(show_spinner=False, max_entries=24)
def cached_run(scenario_json: str, fidelity: Fidelity) -> RunResult:
    return run_once(Scenario.model_validate_json(scenario_json), fidelity)


@st.cache_data(show_spinner=False, max_entries=8)
def cached_ensemble(
    scenario_json: str,
    fidelity: Fidelity,
    n_members: int,
    seed: int,
    sampler: str,
    only: tuple[str, ...],
) -> EnsembleResult:
    return run_uncertainty(
        Scenario.model_validate_json(scenario_json),
        fidelity,
        n_members=n_members,
        seed=seed,
        sampler=sampler,
        only=list(only) or None,
    )


def badge(tier: object) -> str:
    word = render.tier_word(tier)  # type: ignore[arg-type]
    return f'<span class="tierbadge {word}">{word}</span>'


# -- session state -------------------------------------------------------------------------

if "interventions" not in st.session_state:
    st.session_state.interventions = []
if "saved" not in st.session_state:
    st.session_state.saved = {}
if "loaded_starter" not in st.session_state:
    st.session_state.loaded_starter = None


def load_starter(name: str) -> None:
    sc = STARTERS[name]
    st.session_state.loaded_starter = name
    st.session_state.medium = sc.medium
    st.session_state.run_name = sc.name
    st.session_state.duration = sc.duration_days
    st.session_state.closure = sc.closure or "(none)"
    st.session_state.initial = dict(sc.initial)
    st.session_state.temps = [
        {"day": p.day, "celsius": p.celsius} for p in sc.temperature_schedule
    ] or [{"day": 0.0, "celsius": 20.0}]
    st.session_state.interventions = [
        {"day": iv.day, "action": iv.action, "params": dict(iv.params)} for iv in sc.interventions
    ]


if st.session_state.loaded_starter is None:
    load_starter(next(iter(STARTERS)))


# -- sidebar: the batch --------------------------------------------------------------------

with st.sidebar:
    st.markdown("### The batch")

    starter = st.selectbox(
        "Start from",
        list(STARTERS),
        index=list(STARTERS).index(st.session_state.loaded_starter)
        if st.session_state.loaded_starter in STARTERS
        else 0,
        help="Loads a complete batch you can then change. Picking a different one starts over.",
    )
    if starter != st.session_state.loaded_starter:
        load_starter(starter)
        st.rerun()

    name = st.text_input("Call this run", value=st.session_state.get("run_name", "Untitled"))
    medium = st.selectbox(
        "Wine or beer", MEDIA, index=MEDIA.index(st.session_state.get("medium", "wine"))
    )
    if medium != st.session_state.get("medium"):
        st.session_state.medium = medium
        st.session_state.initial = dict(DEFAULT_INPUTS[medium])
        st.rerun()

    duration = st.number_input(
        "Follow it for (days)",
        min_value=0.5,
        max_value=3650.0,
        value=float(st.session_state.get("duration", 14.0)),
        step=1.0,
    )

    st.markdown("### What is in it")
    initial: dict[str, float] = {}
    stored = st.session_state.get("initial", {})
    for key in PRIMARY_INPUTS[medium]:
        unit = input_unit(key)
        initial[key] = st.number_input(
            f"{input_label(key)}" + (f" [{unit}]" if unit else ""),
            value=float(stored.get(key, DEFAULT_INPUTS[medium].get(key, 0.0))),
            step=0.1,
            format="%.4g",
            key=f"init_{key}",
        )

    extra_keys = [k for k in allowed_initial_keys(medium) if k not in PRIMARY_INPUTS[medium]]
    with st.expander(f"Everything else you can put in ({len(extra_keys)} more)"):
        st.caption(
            "The full list the model accepts for this kind of batch, read from the engine "
            "itself so it can never fall behind. Left at zero, an ingredient simply is not "
            "there."
        )
        chosen = st.multiselect(
            "Add an ingredient or a starting condition",
            extra_keys,
            default=[k for k in extra_keys if k in stored],
            format_func=lambda k: f"{input_label(k)}  ({k})",
        )
        for key in chosen:
            unit = input_unit(key)
            initial[key] = st.number_input(
                f"{input_label(key)}" + (f" [{unit}]" if unit else ""),
                value=float(stored.get(key, 0.0)),
                step=0.1,
                format="%.4g",
                key=f"extra_{key}",
            )

    st.markdown("### Temperature")
    st.caption("Straight lines between the points you give. A single row holds it steady.")
    temps_df = st.data_editor(
        pd.DataFrame(st.session_state.get("temps", [{"day": 0.0, "celsius": 20.0}])),
        num_rows="dynamic",
        width="stretch",
        column_config={
            "day": st.column_config.NumberColumn("Day", min_value=0.0, step=1.0),
            "celsius": st.column_config.NumberColumn("°C", step=0.5),
        },
        key="temp_editor",
    )

    st.markdown("### Things you do to it")
    st.caption("Additions, rackings, pitchings, bottling — each on the day it happens.")
    for i, iv in enumerate(list(st.session_state.interventions)):
        spec = VERB_SPECS.get(iv["action"])
        label = spec.label if spec else iv["action"]
        args = ", ".join(f"{k}={v}" for k, v in iv["params"].items())
        cols = st.columns([5, 1])
        cols[0].write(f"**day {iv['day']:g}** — {label}" + (f"  \n`{args}`" if args else ""))
        if cols[1].button("✕", key=f"del_{i}", help="Remove this one"):
            st.session_state.interventions.pop(i)
            st.rerun()

    with st.expander("Add one"):
        available = [v for v in VERBS if medium in VERB_SPECS.get(v, VERB_SPECS["rack"]).media]
        verb = st.selectbox(
            "What you do",
            available,
            format_func=lambda v: VERB_SPECS[v].label if v in VERB_SPECS else v,
        )
        spec = VERB_SPECS.get(verb)
        if spec and spec.note:
            st.caption(spec.note)
        day = st.number_input("On day", min_value=0.0, value=1.0, step=0.5, key="new_day")
        params: dict[str, float | str] = {}
        if spec:
            for pname, plabel, punit, pdefault in spec.numbers:
                params[pname] = st.number_input(
                    f"{plabel}" + (f" [{punit}]" if punit else ""),
                    value=float(pdefault),
                    step=0.01,
                    format="%.4g",
                    key=f"np_{verb}_{pname}",
                )
            for pname, plabel, pchoices in spec.choices:
                params[pname] = st.selectbox(plabel, pchoices, key=f"sp_{verb}_{pname}")
        if st.button("Add it", width="stretch"):
            st.session_state.interventions.append(
                {"day": float(day), "action": verb, "params": params}
            )
            st.rerun()

    closure_options = ["(none)", *CLOSURES]
    closure = st.selectbox(
        "Bottle closure (wine, for aging)",
        closure_options,
        index=closure_options.index(st.session_state.get("closure", "(none)"))
        if st.session_state.get("closure", "(none)") in closure_options
        else 0,
        disabled=medium != "wine",
        help="How much oxygen creeps in past the seal. Each closure carries its own measured "
        "rate, which is why this is a list to choose from rather than a number to invent.",
    )

    # -- settings: three separate things, never one slider ---------------------------------
    st.markdown("---")
    st.markdown("### How hard to work at it")
    st.caption(
        "Three different things get called accuracy. They are kept apart here because only "
        "the last one changes what is actually being modelled."
    )

    precision = st.select_slider(
        "How carefully to do the maths",
        options=[*PRECISION_ORDER, "custom"],
        value=st.session_state.get("precision", "standard"),
        key="precision",
        help="This makes the calculation more exact, not more realistic. It answers one "
        "question: has the answer stopped moving?",
    )
    st.caption(PRECISION_BLURB[precision])

    if precision == "custom":
        method = st.selectbox("Method", METHODS, index=0)
        if method == "RK45":
            st.warning(
                "RK45 is not built for a problem where very fast and very slow chemistry run "
                "side by side, and fermentation always is. Expect it to crawl or give up. It "
                "is here to diagnose, not to use."
            )
        rtol = st.number_input("Relative tolerance", value=1e-6, format="%.1e")
        atol = st.number_input("Absolute tolerance", value=1e-9, format="%.1e")
        cap_step = st.checkbox(
            "Limit how big a step it may take",
            value=False,
            help="The method already chooses its own step sizes, so a limit almost never "
            "changes the answer and always costs time. Here for the rare case where it does.",
        )
        max_step = (
            st.number_input("Largest step (hours)", value=1.0, min_value=1e-3)
            if cap_step
            else float("inf")
        )
        fidelity_kwargs = {
            "precision": "custom",
            "method": method,
            "rtol": rtol,
            "atol": atol,
            "max_step": max_step,
        }
    else:
        rtol, atol = PRECISION_PRESETS[precision]
        fidelity_kwargs = {"precision": precision, "method": "BDF", "rtol": rtol, "atol": atol}

    points = st.select_slider(
        "How many points to keep",
        options=list(POINT_CHOICES),
        value=st.session_state.get("points", 200),
        key="points",
        help="How much of the run gets stored for the charts. It does not change the answer "
        "at all, only how smooth the lines look. It is shared with the uncertainty band, "
        "which can only combine its re-runs if they all land on the same points.",
    )

    with st.expander("Which chemistry to include — the setting that changes the answer"):
        st.caption(
            "The two above change how hard the machine works. This one changes what is being "
            "modelled, so it is the only one that can move the result towards or away from "
            "reality."
        )
        oxidative = st.radio(
            "How oxygen gets used up",
            list(OXIDATIVE_BLURB),
            index=0,
            format_func=lambda k: OXIDATIVE_LABEL.get(k, k),
        )
        st.caption(OXIDATIVE_BLURB[oxidative])
        strict = st.checkbox(
            "Double-check the model against itself",
            value=False,
            help="Checks at every step that no piece of chemistry changes something it never "
            "declared it would. Catches a real class of bug, at a real cost in speed.",
        )

    fidelity = Fidelity(points=points, oxidative=oxidative, strict=strict, **fidelity_kwargs)  # type: ignore[arg-type]


# -- assemble the batch --------------------------------------------------------------------

st.session_state.run_name = name
st.session_state.duration = duration
st.session_state.initial = initial
st.session_state.temps = temps_df.to_dict("records")
st.session_state.closure = closure

temp_points = [
    TemperaturePoint(day=float(row["day"]), celsius=float(row["celsius"]))
    for row in temps_df.to_dict("records")
    if pd.notna(row.get("day")) and pd.notna(row.get("celsius"))
]

scenario_error: str | None = None
scenario: Scenario | None = None
try:
    scenario = Scenario(
        name=name or "Untitled",
        medium=medium,
        initial={k: float(v) for k, v in initial.items() if v},
        temperature_schedule=temp_points or [TemperaturePoint(day=0.0, celsius=20.0)],
        interventions=[
            Intervention(day=iv["day"], action=iv["action"], params=iv["params"])
            for iv in st.session_state.interventions
        ],
        closure=None if closure == "(none)" else closure,
        duration_days=float(duration),
    )
except Exception as exc:  # surfaced rather than crashed
    scenario_error = str(exc)

st.title("Fermentation Console")

if scenario_error or scenario is None:
    st.error(f"This batch does not add up yet.\n\n{scenario_error}")
    st.stop()

for problem in gate_problems(scenario):
    st.warning(problem)

result: RunResult | None = None
run_error: str | None = None
with st.spinner("Working it out…"):
    try:
        result = cached_run(scenario.model_dump_json(), fidelity)
    except Exception as exc:
        run_error = str(exc)

if run_error or result is None:
    st.error(
        "The model would not accept this batch. Its own explanation is below — it is usually "
        "specific about which number is out of range and which way to move it."
    )
    st.code(run_error or "", language="text")
    st.stop()

st.caption(
    f"{scenario.medium} · {scenario.duration_days:g} days · worked out in "
    f"{result.compile_seconds + result.wall_seconds:.2f} s at the "
    f"'{fidelity.precision}' setting"
)

tab_run, tab_spread, tab_compare, tab_data = st.tabs(
    ["The run", "How uncertain is it", "Compare runs", "Numbers and write-up"]
)


# -- tab: the run --------------------------------------------------------------------------

with tab_run:
    rows = readouts.summary(result)
    cols = st.columns(len(rows) or 1)
    for col, (label, value, unit, tier) in zip(cols, rows, strict=False):
        with col:
            st.metric(label, value, help=render.TIER_MEANING[tier])
            st.markdown(f"{unit} &nbsp; {badge(tier)}", unsafe_allow_html=True)

    for panel in render.honesty_panels(result):
        (st.error if panel.kind == "warning" else st.info)(f"**{panel.title}** — {panel.body}")

    with st.expander("Would a more careful calculation have given a different answer?"):
        st.write(
            "Choosing a setting is a claim that it was careful enough. This checks the claim "
            "instead of asserting it: the same batch is worked out again one step stricter, "
            "and the biggest disagreement between the two is reported. It costs one extra run."
        )
        if st.button("Check it"):
            with st.spinner("Working it out again, more strictly…"):
                cc = check_convergence(result, readouts.headline_variables(result))
            if cc is None:
                st.info("Already at the strictest setting, so there is nothing to compare to.")
            else:
                st.session_state.convergence = cc
        cc = st.session_state.get("convergence")
        if cc is not None:
            st.metric(
                f"Biggest change, in {cc.worst_variable}",
                f"{cc.worst:.2e}",
                help=f"Against the '{cc.tighter.precision}' setting, which took "
                f"{cc.seconds:.1f} s.",
            )
            (st.success if cc.converged else st.warning)(cc.verdict)

    st.markdown(
        "**How to read the lines.** The style of each one says how far it can be trusted: "
        "**dashes** where every number behind it is published research, **dots** where at "
        "least one is an estimate, a **thin pale line** where nothing in this run affects it "
        "at all. Solid would mean checked against real measured data — nothing here has "
        "earned that, so you will not see one."
    )

    drawn: list[readouts.Readout] = []
    groups = readouts.groups_for(result)
    for i in range(0, len(groups), 2):
        cols = st.columns(2)
        for col, (group, live) in zip(cols, groups[i : i + 2], strict=False):
            with col:
                st.plotly_chart(
                    render.series_figure(result, group, live),
                    width="stretch",
                    key=f"fig_{group.title}",
                )
                if group.blurb:
                    st.caption(group.blurb)
                for panel in render.caveat_panels(live):
                    st.warning(f"**{panel.title}** — {panel.body}")
                drawn.extend(live)

    # -- the source trail: reachable, deliberately not the headline -------------------------
    st.markdown("---")
    with st.expander("Where any of these numbers actually come from"):
        st.write(
            "Every number in the model had to arrive with a source, the conditions it was "
            "measured under, and a range, before it could be loaded at all. So any line above "
            "can be followed back: which chemistry produced it, which published numbers that "
            "chemistry used, and which paper each of those came from."
        )
        touched = sorted(result.touched_variables())
        default = "E" if "E" in touched else (touched[0] if touched else "")
        variable = st.selectbox(
            "Follow back",
            touched,
            index=touched.index(default) if default in touched else 0,
            help="These are the model's own short names, because that is what each piece of "
            "chemistry declares it changes. E is ethanol, S sugar, X yeast, N nitrogen.",
        )

        st.info(provenance.why_this_tier(result, variable))

        mechs = provenance.mechanisms_for(result, variable)
        cards = provenance.constants_for(result, mechs)
        caps = provenance.limiting(cards)

        left, right = st.columns([1, 2])
        with left:
            st.markdown("**The chemistry involved**")
            for m in mechs:
                mark = " · started partway through" if m.switched_on_mid_run else ""
                kind = "scales the rate" if m.kind == "modifier" else "changes it directly"
                st.markdown(
                    f"{m.name} {badge(m.tier)}<br>"
                    f"<span style='opacity:.7;font-size:12px'>{kind} · uses {len(m.reads)} "
                    f"numbers{mark}</span>",
                    unsafe_allow_html=True,
                )
        with right:
            st.markdown(f"**The numbers it uses ({len(cards)})**")
            if caps:
                st.caption(
                    f"{len(caps)} of them are the weakest link, and are what hold this "
                    "result's confidence down. They are listed first."
                )
            only_caps = st.checkbox("Show only the weakest links", value=False)
            shown = caps if only_caps else cards
            for c in shown[:60]:
                with st.expander(f"{c.name} = {c.value:,.4g} {c.unit}   ({c.tier.label})"):
                    st.markdown(f"**From** — {c.source}")
                    st.markdown(f"**Measured under** — {c.conditions}")
                    if c.doi:
                        st.markdown(f"**DOI** — `{c.doi}`")
                    span = c.span_fraction
                    st.markdown(
                        f"**Could be anywhere from** — {c.low:,.4g} to {c.high:,.4g} {c.unit}"
                        + (f"  (about ±{span * 50:.0f}%)" if span is not None else "")
                    )
                    if c.range_note:
                        st.caption(c.range_note)
                    if c.notes:
                        st.caption(c.notes)
                    st.caption("Used by: " + ", ".join(c.read_by))
            if len(shown) > 60:
                st.caption(f"{len(shown) - 60} more not listed.")

        st.markdown("---")
        est, total = provenance.estimate_census(result)
        census = provenance.tier_census(result)
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Numbers loaded for this run", total)
        c2.metric("Of those, author's estimates", est)
        c3.metric("From published research", census[list(census)[1]])
        c4.metric("Measured directly", census[list(census)[2]])
        st.caption(
            "A single number can be solidly measured — a physical property someone looked up "
            "— while no result of the model is, because a result is only as good as the "
            "weakest thing that went into it."
        )


# -- tab: uncertainty ----------------------------------------------------------------------

with tab_spread:
    st.subheader("How much of this is just the numbers being uncertain?")
    st.write(
        "Every number in the model comes with a range: the paper it came from reported a "
        "spread, or the author's estimate has bounds. Running the batch many times over, each "
        "time picking every number somewhere inside its own range, shows how far the curve "
        "could move on that account alone. It is the uncertainty you have before even asking "
        "whether the model is right."
    )

    scope_all = varying_constants(result)
    scope_var = st.selectbox(
        "Vary the numbers behind",
        ["everything in this run", *sorted(result.touched_variables())],
        help="Narrowing to one quantity varies only the numbers its own chemistry uses. That "
        "is a different question — how far these numbers move it, rather than all of them — "
        "and it needs far fewer re-runs before they can be ranked.",
    )
    if scope_var == "everything in this run":
        only: tuple[str, ...] = ()
        n_varying = len(scope_all)
    else:
        picked_mechs = provenance.mechanisms_for(result, scope_var)
        names = {n for m in picked_mechs for n in m.reads}
        only = tuple(varying_constants(result, only=sorted(names)))
        n_varying = len(only)

    c1, c2, c3 = st.columns(3)
    n_members = c1.number_input("How many re-runs", min_value=4, max_value=400, value=48, step=8)
    seed = c2.number_input(
        "Starting seed",
        min_value=0,
        value=0,
        step=1,
        help="The same seed gives the same draws, so a band you show someone else can be "
        "reproduced exactly.",
    )
    sampler = c3.selectbox(
        "How to pick the values", list(SAMPLERS), format_func=lambda k: SAMPLERS[k]
    )

    per_run = result.wall_seconds + result.compile_seconds * 0.1
    st.caption(
        f"{n_varying} numbers would vary. Roughly {per_run * n_members:.0f} seconds for "
        f"{n_members} re-runs. The band itself works at any size; ranking *which* number "
        f"matters most needs more re-runs than there are varying numbers, so at least "
        f"{n_varying + 1} here."
    )

    if st.button("Run them", type="primary"):
        with st.spinner(f"Running it {n_members} times…"):
            try:
                st.session_state.ens = cached_ensemble(
                    scenario.model_dump_json(), fidelity, int(n_members), int(seed), sampler, only
                )
            except Exception as exc:
                st.error(str(exc))

    ens = st.session_state.get("ens")
    if ens is not None and ens.scenario == scenario:
        for w in ens.warnings:
            st.warning(w)
        st.caption(
            f"{ens.ensemble.n_succeeded} of {ens.n_requested} re-runs finished, in "
            f"{ens.wall_seconds:.0f} s. {len(ens.ensemble.sampled_names)} numbers varied."
        )

        can_band = [
            r
            for r in readouts.READOUTS
            if scenario.medium in r.media and r.available(result) and render.bandable(r, result)
        ]
        pick = st.multiselect(
            "Draw a band for",
            [r.key for r in can_band],
            default=[k for k in ("ethanol", "sugar", "biomass") if k in {r.key for r in can_band}],
            format_func=lambda k: readouts.BY_KEY[k].label,
            help="Only quantities the model tracks directly. pH is worked out from all the "
            "acids together, so a range drawn on any one of them would not be a range on pH.",
        )
        for key in pick:
            st.plotly_chart(
                render.band_figure(ens, readouts.BY_KEY[key], result),
                width="stretch",
                key=f"band_{key}",
            )

        st.markdown("#### Which numbers the uncertainty comes from")
        st.caption(
            "If one number dominates this list, finding a better source for it is the single "
            "most useful thing anyone could do to the model."
        )
        target = st.selectbox("For", sorted(result.touched_variables()), key="attr_var")
        method = st.radio(
            "Fit",
            ["src", "srrc"],
            horizontal=True,
            format_func=lambda m: {
                "src": "straight-line fit",
                "srrc": "rank-based fit (for curved responses)",
            }[m],
        )
        try:
            att = attribute_spread(ens.ensemble, target, ens.param_tiers, method=method)
        except Exception as exc:
            st.warning(str(exc))
        else:
            st.plotly_chart(render.spread_figure(att), width="stretch", key="spreadfig")
            st.caption(
                f"This simple fit accounts for {att.r_squared:.0%} of the movement; the other "
                f"{att.unexplained:.0%} comes from numbers interacting with each other, which "
                "it cannot split apart. Blue means turning that number up turns the result "
                "up; red means the opposite."
            )
            by_tier = ", ".join(
                f"{t.label}: {v:.0%}"
                for t, v in sorted(att.per_tier.items(), key=lambda kv: -kv[1])
            )
            st.caption("Grouped by how well sourced the number is — " + by_tier)
    elif ens is not None:
        st.info("The batch has changed since that band was worked out. Run it again.")


# -- tab: compare --------------------------------------------------------------------------

with tab_compare:
    st.subheader("The same question asked of several batches")
    c1, c2 = st.columns([3, 1])
    save_as = c1.text_input("Keep this run as", value=scenario.name, key="save_as")
    if c2.button("Keep it", width="stretch"):
        st.session_state.saved[save_as] = result
        st.rerun()

    if not st.session_state.saved:
        st.info(
            "Nothing kept yet. Keep this run, change something in the sidebar, keep that one "
            "too, and the difference gets drawn here."
        )
    else:
        picked = st.multiselect(
            "Runs", list(st.session_state.saved), default=list(st.session_state.saved)
        )
        keys = [
            r.key for r in readouts.READOUTS if scenario.medium in r.media and r.available(result)
        ]
        chosen = st.multiselect(
            "Draw",
            keys,
            default=[k for k in ("sugar", "ethanol", "ph") if k in keys],
            format_func=lambda k: readouts.BY_KEY[k].label,
        )
        pairs = [(nm, st.session_state.saved[nm]) for nm in picked]
        for key in chosen:
            st.plotly_chart(
                render.compare_figure(pairs, readouts.BY_KEY[key]),
                width="stretch",
                key=f"cmp_{key}",
            )

        st.markdown("#### Where each one ended up")
        table: dict[str, dict[str, str]] = {}
        for nm, res in pairs:
            table[nm] = {lab: f"{val} {unit}" for lab, val, unit, _ in readouts.summary(res)}
        st.dataframe(pd.DataFrame(table), width="stretch")

        if st.button("Forget them all"):
            st.session_state.saved = {}
            st.rerun()


# -- tab: numbers and write-up -------------------------------------------------------------

with tab_data:
    st.subheader("Everything the model tracks, at the end of the run")
    st.caption(
        "The raw view underneath the charts. 'Inert' means nothing in this run affected it — "
        "it stayed at its starting value, which is not the same as being confirmed correct."
    )
    st.dataframe(
        pd.DataFrame(
            render.final_state_rows(result), columns=["what", "value at the end", "confidence"]
        ),
        width="stretch",
        height=420,
    )

    st.subheader("The batch, exactly as the model received it")
    st.caption("Copy this to reproduce the run elsewhere, or to paste into a record.")
    st.code(scenario.model_dump_json(indent=2), language="json")

    st.subheader("Write it up")
    st.write(
        "The same charts and the same source tables as this page, in one file that opens in "
        "any browser with nothing else installed. It is produced by the same code that draws "
        "this screen, so the two cannot end up disagreeing."
    )
    explain_vars = st.multiselect(
        "Include the source tables for",
        sorted(result.touched_variables()),
        default=[v for v in ("E", "S") if v in result.touched_variables()],
    )
    if st.button("Write it"):
        out = REPORT_DIR / f"{scenario.name.replace(' ', '-').lower()}.html"
        path = report.write(
            result,
            out,
            convergence=st.session_state.get("convergence"),
            explain=tuple(explain_vars),
        )
        st.success(f"Saved to {path}")
        st.download_button(
            "Download it", data=path.read_bytes(), file_name=path.name, mime="text/html"
        )
