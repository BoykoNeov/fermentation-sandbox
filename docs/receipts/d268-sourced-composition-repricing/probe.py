"""D-268 probe: re-price every growth-anchored arm at the SOURCED composition.

Run from the repo root.  Prints a JSON blob to stdout; progress to stderr.
The composition is swapped by re-pointing ``_D259_RESIDUE_SHARE_OF_PROTEIN`` at a
degenerate bracket (lo = mid = hi = the sourced free-acid share), which is exactly the
semantics the shipped repair will have: the residue half stops being a bracket, the
protein-fraction half stays one.
"""

from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.getcwd())

import tests.test_precursor_fates as tpf
import tests.test_fusel_joint_repair as tjr

#: The SHIPPED (post-D-268) weights: sourced composition x protein fraction.
NEW = tpf._d259_weights


def OLD(edge: str) -> dict:
    """D-259's stated bracket, reconstructed from the superseded literal."""
    return {
        species: tpf._D259_PROTEIN_FRACTION_OF_DRY_WEIGHT[edge] * shares[edge] / 100.0
        for species, shares in tpf._D259_BRACKET_SUPERSEDED_AT_D268.items()
    }


def use(fn):
    """Both modules read ``_d259_weights``; the joint-repair module imported it BY VALUE."""
    tpf._d259_weights = fn
    tjr._d259_weights = fn


def say(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


out: dict = {"weights": {}}
for name, fn in (("old", OLD), ("new", NEW)):
    out["weights"][name] = {e: fn(e) for e in tpf._D259_EDGES}

# ---- P7 FIRST: the joint arm on Crepin's must -------------------------------------------
say("P7: crepin joint/growth arms, new composition")
use(NEW)
out["crepin_new"] = {
    "growth": tjr._crepin_arm(False, True),
    "joint": tjr._crepin_arm(True, True),
    "joint_lo": tjr._crepin_arm(True, True, edge="lo"),
    "joint_hi": tjr._crepin_arm(True, True, edge="hi"),
}
say(f"  joint split={out['crepin_new']['joint']['split_pct']['leucine']:.2f} "
    f"tracer={out['crepin_new']['joint']['leucine_tracer_pct']:.3f}")
print(json.dumps(out, indent=1, default=float), flush=True)

say("P7 control: same arms, old composition (reproduction of D-266)")
use(OLD)
out["crepin_old"] = {
    "growth": tjr._crepin_arm(False, True),
    "joint": tjr._crepin_arm(True, True),
}
say(f"  joint split={out['crepin_old']['joint']['split_pct']['leucine']:.2f} "
    f"tracer={out['crepin_old']['joint']['leucine_tracer_pct']:.3f}")

# ---- D-259's own arms --------------------------------------------------------------------
say("D-259 uncorrected growth-anchored splits, new composition")
use(NEW)
out["d259_new"] = {e: tpf._d259_growth_anchored_split(e) for e in tpf._D259_EDGES}
say("D-259 uncorrected growth-anchored splits, old composition")
use(OLD)
out["d259_old"] = {e: tpf._d259_growth_anchored_split(e) for e in tpf._D259_EDGES}
print(json.dumps(out, indent=1, default=float), flush=True)

# ---- D-260 arms ---------------------------------------------------------------------------
say("D-260 corrected arms, new composition")
use(NEW)
WEIGHTS = NEW
d260: dict = {}
for e in tpf._D259_EDGES:
    d260[f"growth_mods_{e}"] = tpf._d260_arm(
        sink=tpf._D260GrowthAnchoredFates(WEIGHTS(e)), attach_growth_modifiers=True
    )
    say(f"  growth_mods_{e} split={d260[f'growth_mods_{e}']['split_pct']:.2f}")
d260["growth_nomods_mid"] = tpf._d260_arm(
    sink=tpf._D260GrowthAnchoredFates(WEIGHTS("mid")), attach_growth_modifiers=False
)
matched_f = round(d260["growth_nomods_mid"]["split_pct"] / 100.0, 4)
say(f"  growth_nomods_mid split={d260['growth_nomods_mid']['split_pct']:.3f} -> f={matched_f}")
d260["f_matched"] = tpf._d260_arm(f_leucine=matched_f)
d260["matched_f"] = matched_f
d260["shipped"] = tpf._d260_arm()

# lambda: the over-draw that lands Crepin's band.  Swept, not assumed.
say("lambda sweep")
d260["lambda_sweep"] = {}
for lam in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
    arm = tpf._d260_arm(
        sink=tpf._D260GrowthAnchoredFates(WEIGHTS("mid"), kappa=0.01, lam=lam),
        attach_growth_modifiers=True,
        points=20001,
    )
    d260["lambda_sweep"][lam] = arm
    say(f"  lam={lam} split={arm['split_pct']:.2f} tracer={arm['tracer_pct']:.3f} "
        f"closure={arm['closure']:.4f}")
out["d260_new"] = d260
print(json.dumps(out, indent=1, default=float), flush=True)

# ---- the remaining D-266 fixtures at the new composition ---------------------------------
say("D-266 rollero arms (growth, joint), new composition")
use(NEW)
rollero: dict = {}
for arm_name in ("growth", "joint"):
    blend, growth = tjr._ARMS[arm_name]
    rollero[arm_name] = {}
    for yan in tjr._YANS:
        valine = tjr._rollero_arm(blend, growth, yan, "valine")
        leucine = tjr._rollero_arm(blend, growth, yan, "leucine")
        valine["leucine_enrichment_pct"] = leucine["enrichment_pct"]
        valine["leucine_labelled_um"] = leucine["labelled_um"]
        rollero[arm_name][yan] = valine
        say(f"  {arm_name} SM{yan:.0f} leu_enrich={valine['leucine_enrichment_pct']:.2f} "
            f"val_enrich={valine['enrichment_pct']:.2f} nt={valine['nt_fraction_pct']:.1f}")
out["rollero_new"] = rollero

say("D-266 d104 arms (growth, joint), new composition")
out["d104_new"] = {n: tjr._d104_arm(*tjr._ARMS[n]) for n in ("growth", "joint")}
say(f"  joint phe_left={out['d104_new']['joint']['phe_left_pct']:.4f} "
    f"refund={out['d104_new']['joint']['worst_joint_c_refund']:.4f}")

print(json.dumps(out, indent=1, default=float), flush=True)
say("DONE")
