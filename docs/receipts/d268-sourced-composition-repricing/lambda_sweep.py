import json, os, sys, time
sys.path.insert(0, os.getcwd())
import tests.test_precursor_fates as tpf

def OLD(edge):
    return {sp: tpf._D259_PROTEIN_FRACTION_OF_DRY_WEIGHT[edge] * sh[edge] / 100.0
            for sp, sh in tpf._D259_BRACKET_SUPERSEDED_AT_D268.items()}

out = {}
for lam in (4.0, 5.0):
    arm = tpf._d260_arm(
        sink=tpf._D260GrowthAnchoredFates(OLD("mid"), kappa=0.01, lam=lam),
        attach_growth_modifiers=True, points=20001)
    out[f"old_lam_{lam}"] = arm
    print(f"OLD lam={lam} split={arm['split_pct']:.3f} tracer={arm['tracer_pct']:.4f} "
          f"closure={arm['closure']:.4f}", file=sys.stderr, flush=True)
# the over-draw that first reaches Crepin's 77 % at the SOURCED composition
for lam in (4.5, 4.75):
    arm = tpf._d260_arm(
        sink=tpf._D260GrowthAnchoredFates(tpf._d259_weights("mid"), kappa=0.01, lam=lam),
        attach_growth_modifiers=True, points=20001)
    out[f"new_lam_{lam}"] = arm
    print(f"NEW lam={lam} split={arm['split_pct']:.3f} tracer={arm['tracer_pct']:.4f} "
          f"closure={arm['closure']:.4f}", file=sys.stderr, flush=True)
json.dump(out, open("M:/claud_projects/temp/d268/lam.json", "w"), indent=1, default=float)
print("DONE", file=sys.stderr, flush=True)
