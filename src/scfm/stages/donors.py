"""Donor mixing of the cell states, and whether it relates to their survival signal; called by ladder.

Per representation, over the states kept as signatures: how many are dominated by one donor
(top-donor share >= SINGLE_DONOR_FRAC), the median top-donor share and donor entropy, and the Spearman
correlation between a state's above-floor margin (both readings) and its top-donor share. Each state
enters once, through its top-k = DONOR_MIXING_TOPK signature. Then the best above-floor margin among
multi-donor states (top-donor share < MULTI_DONOR_FRAC), compared across representations.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C

BASE = "hvg_pca"


def states_with_margins(comp: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    s = scores[scores["topk"] == C.DONOR_MIXING_TOPK].copy()
    s["cluster"] = s["signature"].str.extract(r"^c(\d+)_k\d+$", expand=False).astype(int)
    k = comp[comp["kept"].astype(bool)].copy()
    k["resolution"] = k["resolution"].astype(float)
    s["resolution"] = s["resolution"].astype(float)
    cols = ["model", "resolution", "cluster", "signature", "above_floor", "above_floor_or", "direction"]
    d = k.merge(s[cols], on=["model", "resolution", "cluster"], how="left", validate="one_to_one")
    if d["signature"].isna().any():
        raise ValueError("a kept state has no top-k signature in scores.csv")
    return d


def _best(g: pd.DataFrame, col: str) -> dict:
    if g.empty:
        return {"value": None, "state": None, "top_cell_type": None, "top_donor_frac": None}
    r = g.loc[g[col].idxmax()]
    return {"value": float(r[col]), "state": f"{r['resolution']}/{r['signature']}",
            "top_cell_type": r["top_cell_type"], "top_donor_frac": float(r["top_donor_frac"])}


def _describe(g: pd.DataFrame) -> dict:
    single = g["top_donor_frac"] >= C.SINGLE_DONOR_FRAC
    return {"n_states": len(g), "n_single_donor": int(single.sum()),
            "share_single_donor": float(single.mean()),
            "n_multi_donor": int((g["top_donor_frac"] < C.MULTI_DONOR_FRAC).sum()),
            "median_top_donor_frac": float(g["top_donor_frac"].median()),
            "median_donor_entropy": float(g["donor_entropy"].median()),
            "median_n_donors": float(g["n_donors"].median())}


def donor_mixing(comp_csv: Path, scores: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    from scipy.stats import fisher_exact, spearmanr
    d = states_with_margins(pd.read_csv(comp_csv), scores)
    models = list(dict.fromkeys(d["model"]))
    rows = []
    for m in models:
        for res, g in [("all", d[d["model"] == m])] + [
                (str(r), g) for r, g in d[d["model"] == m].groupby("resolution")]:
            rec = {"model": m, "resolution": res, **_describe(g)}
            b = d[d["model"] == BASE] if res == "all" else d[(d["model"] == BASE)
                                                             & (d["resolution"] == float(res))]
            if m != BASE:
                sb = b["top_donor_frac"] >= C.SINGLE_DONOR_FRAC
                rec["fisher_p_single_donor_vs_baseline"] = float(fisher_exact(
                    [[rec["n_single_donor"], rec["n_states"] - rec["n_single_donor"]],
                     [int(sb.sum()), int((~sb).sum())]])[1])
            for col, tag in (("above_floor", ""), ("above_floor_or", "_or")):
                if res == "all":
                    rho, p = spearmanr(g[col], g["top_donor_frac"])
                    rec[f"spearman_rho_margin{tag}_vs_top_donor"] = float(rho)
                    rec[f"spearman_p{tag}"] = float(p)
                multi = g[g["top_donor_frac"] < C.MULTI_DONOR_FRAC]
                single = g[g["top_donor_frac"] >= C.SINGLE_DONOR_FRAC]
                bm, bs, ba = _best(multi, col), _best(single, col), _best(g, col)
                rec[f"best_margin{tag}_multi_donor"] = bm["value"]
                rec[f"best_margin{tag}_multi_donor_state"] = bm["state"]
                rec[f"best_margin{tag}_multi_donor_cell_type"] = bm["top_cell_type"]
                rec[f"best_margin{tag}_single_donor"] = bs["value"]
                rec[f"best_margin{tag}_all"] = ba["value"]
                rec[f"best_margin{tag}_all_top_donor_frac"] = ba["top_donor_frac"]
                rec[f"median_margin{tag}_single_donor"] = float(single[col].median()) if len(single) else None
                rec[f"median_margin{tag}_other"] = float(g.loc[g["top_donor_frac"] < C.SINGLE_DONOR_FRAC,
                                                               col].median())
            rows.append(rec)
    df = pd.DataFrame(rows)
    pooled = df[df["resolution"] == "all"].set_index("model")
    fms = [m for m in models if m != BASE]
    summary = {
        "definitions": {
            "states": f"clusters kept as signatures (>= {C.MIN_CELLS_PER_STATE} cells), one row each; 'all' "
                      "pools the three resolutions, whose states overlap in cells (so pooled counts and the "
                      "pooled Fisher p are descriptive)",
            "single_donor": f"top-donor share >= {C.SINGLE_DONOR_FRAC}",
            "multi_donor": f"top-donor share < {C.MULTI_DONOR_FRAC}",
            "margin": f"above-floor margin of the state's top-k = {C.DONOR_MIXING_TOPK} signature (one per "
                      "state, so a state is not counted three times); _or = oriented reading",
            "spearman": "Spearman rho between margin and top-donor share over all kept states of the "
                        "representation (three resolutions pooled)",
        },
        "per_representation": pooled.reset_index().replace({np.nan: None}).to_dict("records"),
        "checks": {
            f"{m}_fewer_single_donor_share_than_baseline": bool(
                pooled.loc[m, "share_single_donor"] < pooled.loc[BASE, "share_single_donor"]) for m in fms},
    }
    for tag in ("", "_or"):
        for m in fms:
            col = f"best_margin{tag}_multi_donor"
            summary["checks"][f"{m}_best_multi_donor_margin{tag}_above_baseline"] = bool(
                pooled.loc[m, col] > pooled.loc[BASE, col])
    return df, summary


def summary_lines(df: pd.DataFrame) -> list[str]:
    lines = ["", f"Donor mixing of the states (single donor: top-donor share >= {C.SINGLE_DONOR_FRAC}; "
                 f"margins of the k = {C.DONOR_MIXING_TOPK} signature)"]
    for _, r in df[df["resolution"] == "all"].iterrows():
        lines.append(f"  {r['model']:<11} {r['n_single_donor']}/{r['n_states']} single-donor "
                     f"({r['share_single_donor']:.0%}), median top-donor {r['median_top_donor_frac']:.2f}, "
                     f"entropy {r['median_donor_entropy']:.2f}; rho(margin, top donor) "
                     f"{r['spearman_rho_margin_vs_top_donor']:+.2f} (p {r['spearman_p']:.2g}), oriented "
                     f"{r['spearman_rho_margin_or_vs_top_donor']:+.2f} (p {r['spearman_p_or']:.2g}); best "
                     f"multi-donor margin {r['best_margin_multi_donor']:+.3f} / oriented "
                     f"{r['best_margin_or_multi_donor']:+.3f}")
    return lines
