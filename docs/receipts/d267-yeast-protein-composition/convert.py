"""D-267: Lange & Heijnen 2001 Table IV (mol %) -> the two frames a mass draw can be in.

The source table is a MOLE composition. The test constant it answers,
``_D259_RESIDUE_SHARE_OF_PROTEIN``, is declared in mass ("g residue / 100 g yeast protein"),
and the Process that reads it (``_GrowthAnchoredFates``) subtracts its product from a FREE
amino-acid pool. Those are two different masses per mole, so the conversion is done both ways
here and the fork is decided in the record, not averaged.

No src/ change and no model run: this is arithmetic on a transcribed table, plus a cross-check
that the free-acid molar masses used here are the ones the engine itself carries.

Run from the repo root:  uv run --python 3.13 python docs/receipts/d267-yeast-protein-composition/convert.py
Writes findings.json beside itself.
"""

from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from fermentation.core.chemistry import MOLAR_MASS  # noqa: E402

#: Lange & Heijnen (2001), Biotechnol Bioeng 75(3):334-344, Table IV p.339, VERBATIM.
#: "Amino acid composition of the protein as measured (mol %)". Asx = Asp + Asn,
#: Glx = Glm + Gln (the paper's own footnote). Glucose-limited chemostat, S. cerevisiae;
#: the paper states the relative abundance did not vary between cultures.
TABLE_IV_MOL_PCT = {
    "Ala": 9.77, "Arg": 3.86, "Asx": 9.28, "Cys": 0.14, "Glx": 15.48,
    "Gly": 8.89, "His": 1.93, "Ile": 5.89, "Leu": 8.01, "Lys": 6.57,
    "Met": 1.14, "Orn": 0.24, "Phe": 3.76, "Pro": 4.22, "Ser": 5.33,
    "Thr": 5.57, "Trp": 0.65, "Tyr": 1.96, "Val": 7.33,
}  # fmt: skip

#: Molar mass of the FREE amino acid, g/mol. Asx is taken as Asp and Glx as Glu -- the acid
#: form, which is what acid hydrolysis (the paper's method) actually delivers. The amide
#: reading is the loser and is reported in the output so the size of the fork is visible:
#: it is under 0.3 % on the mean residue mass and does not move any verdict.
M_FREE = {
    "Ala": 89.094, "Arg": 174.201, "Asx": 133.103, "Cys": 121.159, "Glx": 147.130,
    "Gly": 75.067, "His": 155.155, "Ile": 131.175, "Leu": 131.175, "Lys": 146.189,
    "Met": 149.208, "Orn": 132.161, "Phe": 165.192, "Pro": 115.132, "Ser": 105.093,
    "Thr": 119.119, "Trp": 204.229, "Tyr": 181.191, "Val": 117.148,
}  # fmt: skip
M_FREE_AMIDE = {**M_FREE, "Asx": 132.119, "Glx": 146.146}
M_WATER = 18.0153

#: The bracket D-259 stated for want of a source, g/100 g protein: lo / mid / hi.
D259_BRACKET = {
    "Leu": (6.0, 7.5, 9.0), "Ile": (4.0, 5.0, 6.0), "Val": (4.5, 5.5, 6.5),
    "Thr": (4.0, 5.0, 6.0), "Phe": (3.5, 4.5, 5.5),
}  # fmt: skip

#: Table IV's three-letter names for the five precursors the fusel thread draws.
PRECURSOR = {
    "Leu": "leucine", "Ile": "isoleucine", "Val": "valine",
    "Thr": "threonine", "Phe": "phenylalanine",
}  # fmt: skip


def shares(m_free: dict[str, float]) -> dict[str, dict[str, float]]:
    """Both mass frames, per 100 g of protein, from one mole composition.

    ``residue`` weights each mole by the residue it contributes to the chain; ``free`` weights
    it by the free acid consumed to contribute it. The denominator is the protein mass either
    way -- the sum of RESIDUE masses -- so ``free`` sums to well over 100 by exactly the water
    released in peptide-bond formation, which is the point.
    """
    m_res = {k: v - M_WATER for k, v in m_free.items()}
    protein = sum(TABLE_IV_MOL_PCT[k] * m_res[k] for k in TABLE_IV_MOL_PCT)
    return {
        "residue": {k: 100.0 * TABLE_IV_MOL_PCT[k] * m_res[k] / protein for k in TABLE_IV_MOL_PCT},
        "free": {k: 100.0 * TABLE_IV_MOL_PCT[k] * m_free[k] / protein for k in TABLE_IV_MOL_PCT},
        "_protein_mass_per_100_mol": {"g": protein, "mean_residue_Da": protein / 100.0},
    }


def main() -> None:
    # The engine's own molar masses must be the free acids, or the "free" frame below is not
    # the frame the draw is actually in. Checked rather than asserted in prose.
    engine = {k: MOLAR_MASS[v] for k, v in PRECURSOR.items()}
    mismatch = {k: (engine[k], M_FREE[k]) for k in PRECURSOR if abs(engine[k] - M_FREE[k]) > 0.01}

    acid = shares(M_FREE)
    amide = shares(M_FREE_AMIDE)

    verdicts = {}
    for k, (lo, mid, hi) in D259_BRACKET.items():
        verdicts[k] = {
            "bracket": [lo, mid, hi],
            "residue_frame": acid["residue"][k],
            "free_frame": acid["free"][k],
            "residue_verdict": (
                "above hi" if acid["residue"][k] > hi
                else "below lo" if acid["residue"][k] < lo
                else "inside"
            ),
            "free_verdict": (
                "above hi" if acid["free"][k] > hi
                else "below lo" if acid["free"][k] < lo
                else "inside"
            ),
            "free_over_hi_ratio": acid["free"][k] / hi,
        }

    out = {
        "source": (
            "Lange HC, Heijnen JJ (2001) Biotechnol Bioeng 75(3):334-344, Table IV p.339, "
            "'Amino acid composition of the protein as measured (mol %)'"
        ),
        "mol_pct_total": sum(TABLE_IV_MOL_PCT.values()),
        "engine_molar_mass_crosscheck": {
            "checked": sorted(PRECURSOR),
            "mismatch": mismatch,
            "note": "engine MOLAR_MASS is the FREE amino acid (leucine 131.175), not the residue",
        },
        "mean_residue_mass_Da": {
            "acid_form": acid["_protein_mass_per_100_mol"]["mean_residue_Da"],
            "amide_form": amide["_protein_mass_per_100_mol"]["mean_residue_Da"],
        },
        "sum_of_shares": {
            "residue_frame": sum(acid["residue"].values()),
            "free_frame": sum(acid["free"].values()),
        },
        "all_species_residue_frame": acid["residue"],
        "all_species_free_frame": acid["free"],
        "precursor_verdicts": verdicts,
        "amide_fork_size_on_five": {
            k: abs(acid["free"][k] - amide["free"][k]) for k in D259_BRACKET
        },
    }

    here = pathlib.Path(__file__).resolve().parent
    (here / "findings.json").write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf8")

    print(f"mol %% total: {out['mol_pct_total']:.2f} (the table's own closure)")
    print(f"engine molar-mass mismatches: {mismatch or 'none'}")
    print(f"mean residue mass: {out['mean_residue_mass_Da']['acid_form']:.3f} Da (acid form)")
    print(f"shares sum: residue {out['sum_of_shares']['residue_frame']:.2f}, "
          f"free {out['sum_of_shares']['free_frame']:.2f} g/100 g protein")
    print()
    print(f"{'':4s} {'bracket lo/mid/hi':>20s} {'residue':>9s} {'':>9s} {'free':>9s} {'':>9s}")
    for k, v in verdicts.items():
        lo, mid, hi = v["bracket"]
        print(f"{k:4s} {lo:6.1f}{mid:7.1f}{hi:7.1f} {v['residue_frame']:9.3f} "
              f"{v['residue_verdict']:>9s} {v['free_frame']:9.3f} {v['free_verdict']:>9s}")


if __name__ == "__main__":
    main()
