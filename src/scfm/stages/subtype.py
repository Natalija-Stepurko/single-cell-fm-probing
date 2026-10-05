"""Subtype association of the picks — outcome-free; called by the ladder stage.

For every pick (primary, sensitivity), family maximum and the proliferation reference: the score's
distribution by subtype, Kruskal-Wallis H with a rank-based effect size, and the one-vs-rest AUC of
the subtype group the score tracks most. TCGA uses the PAM50 call (Normal kept as a level); METABRIC
repeats it with CLAUDIN_SUBTYPE (no outcomes read, aggregate statistics only). The survival
association given PAM50 is not recomputed here: it is the age + stage + PAM50 rows of added_value.csv.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.stats import auc_higher_in, kruskal_eta2

DEFINITIONS = {
    "score": "mean z of the signature genes (z per gene across the cohort's patients); higher = the "
             "programme is more expressed, whatever its survival direction",
    "eta2_h": "(H - k + 1) / (n - k), rank-based eta squared of the Kruskal-Wallis H",
    "epsilon2": "H / (n - 1)",
    "auc_<group>": "P(score of a group member > score of a non-member), ties 1/2; below 0.5 = lower "
                   "in the group",
    "tracked_group": "the group whose one-vs-rest AUC is furthest from 0.5",
    "auc_tracked_tcga_group": "METABRIC rows: the AUC for the group tracked in TCGA",
}


def _groups(levels) -> dict[str, list[str]]:
    return {g: m for g, m in C.SUBTYPE_GROUPS.items() if all(x in levels for x in m)}


def association(score: pd.Series, subtype: pd.Series) -> dict:
    d = pd.DataFrame({"s": score, "t": subtype}).dropna()
    levels = sorted(d["t"].unique())
    rec = kruskal_eta2(d["s"].to_numpy(), d["t"].to_numpy())
    for lv in levels:
        x = d.loc[d["t"] == lv, "s"]
        rec[f"n_{lv}"] = len(x)
        rec[f"median_{lv}"] = float(x.median())
        rec[f"q25_{lv}"], rec[f"q75_{lv}"] = (float(v) for v in x.quantile([0.25, 0.75]))
    aucs = {g: auc_higher_in(d["s"], d["t"].isin(m)) for g, m in _groups(levels).items()}
    for g, a in aucs.items():
        rec[f"auc_{g}"] = a
    g = max(aucs, key=lambda k: abs(aucs[k] - 0.5))
    rec.update(tracked_group=g, auc_tracked=aucs[g], auc_tracked_oriented=max(aucs[g], 1 - aucs[g]),
               higher_in_tracked_group=bool(aucs[g] > 0.5))
    return rec


def _items_by_signature(items: list[dict]) -> list[dict]:
    out: dict[tuple, dict] = {}
    for it in items:
        k = (it["model"], it["signature_key"])
        if k not in out:
            out[k] = {**it, "analyses": []}
        if it["analysis"] not in out[k]["analyses"]:
            out[k]["analyses"].append(it["analysis"])
    return list(out.values())


def _survival_given_pam50(av: pd.DataFrame | None, it: dict) -> dict:
    if av is None:
        return {}
    r = av[(av["base"] == "age+stage+PAM50") & (av["model"] == it["model"])
           & (av["signature"] == it["signature"])
           & (av["resolution"].astype(str) == str(it.get("resolution")))]
    if r.empty:
        return {}
    r = r.iloc[0]
    return {"hr_per_sd_given_pam50": float(r["hr_per_sd"]), "lrt_p_given_pam50": float(r["lrt_p_nominal"]),
            "delta_cindex_cv_fixed_given_pam50": float(r["delta_cindex_cv_fixed_mean"])}


def _metabric(metabric_dir: Path):
    from scfm.stages.replicate import load_metabric
    files = {f: metabric_dir / f for f in C.METABRIC_SHA256}
    if not all(p.exists() for p in files.values()):
        return None
    zexpr, _, clin = load_metabric(files)
    sub = clin["subtype"].astype(object)
    return zexpr, sub.where(sub.notna() & ~sub.isin(C.METABRIC_SUBTYPE_EXCLUDE), np.nan)


def subtype_association(items: list[dict], zexpr: pd.DataFrame, clin: pd.DataFrame,
                        av: pd.DataFrame | None, metabric_dir: Path | None) -> tuple[pd.DataFrame, dict]:
    sigs = _items_by_signature(items)
    rows = []
    tcga_sub = clin["subtype"]
    for it in sigs:
        genes = [g for g in it["genes"] if g in zexpr.index]
        rec = {"cohort": "TCGA", "subtype_column": C.SUBTYPE_COLUMN, "analyses": ";".join(it["analyses"]),
               "model": it["model"], "signature": it["signature"], "resolution": it.get("resolution"),
               "top_cell_type": it.get("top_cell_type"), "direction_tcga": it["direction_label"],
               "n_genes": len(it["genes"]), "n_genes_scored": len(genes),
               "coverage": len(genes) / len(it["genes"])}
        rec.update(association(zexpr.loc[genes].mean(axis=0), tcga_sub.reindex(zexpr.columns)))
        rec.update(_survival_given_pam50(av, it))
        rows.append(rec)
    note = None
    mb = _metabric(metabric_dir) if metabric_dir is not None else None
    if mb is None:
        note = f"METABRIC not run: files not found in {metabric_dir}"
    else:
        mz, msub = mb
        tracked = {(r["model"], r["resolution"], r["signature"]): r["tracked_group"] for r in rows}
        for it in sigs:
            genes = [g for g in it["genes"] if g in mz.index]
            cov = len(genes) / len(it["genes"])
            rec = {"cohort": "METABRIC", "subtype_column": "CLAUDIN_SUBTYPE",
                   "analyses": ";".join(it["analyses"]), "model": it["model"], "signature": it["signature"],
                   "resolution": it.get("resolution"), "top_cell_type": it.get("top_cell_type"),
                   "direction_tcga": it["direction_label"], "n_genes": len(it["genes"]),
                   "n_genes_scored": len(genes), "coverage": cov}
            if cov >= C.REPLICATE_MIN_COVERAGE and len(genes) >= 3:
                s = mz.loc[genes].mean(axis=0)
                rec.update(association(s, msub.reindex(mz.columns)))
                g = tracked[(it["model"], it.get("resolution"), it["signature"])]
                rec["tcga_tracked_group"] = g
                rec["auc_tracked_tcga_group"] = rec.get(f"auc_{g}", np.nan)
            rows.append(rec)
    df = pd.DataFrame(rows)
    keep = ["cohort", "analyses", "model", "resolution", "signature", "top_cell_type", "direction_tcga", "n",
            "kw_h", "kw_p", "eta2_h", "epsilon2", "tracked_group", "auc_tracked", "higher_in_tracked_group",
            "auc_tracked_tcga_group", "hr_per_sd_given_pam50", "lrt_p_given_pam50"]
    summary = {
        "definitions": DEFINITIONS,
        "levels": {"TCGA": "PAM50 call (LumA, LumB, Basal, Her2, Normal); patients without a call excluded",
                   "METABRIC": "CLAUDIN_SUBTYPE (LumA, LumB, Her2, Basal, Normal, claudin-low); NC "
                               "(not classified) and missing excluded; no outcomes read"},
        "groups": C.SUBTYPE_GROUPS,
        "survival_given_pam50": "hr_per_sd / lrt_p / delta C given age + stage + PAM50 are copied from "
                                "added_value.csv (base age+stage+PAM50); lrt_p is nominal for selected "
                                "signatures",
        "metabric_note": note,
        "rows": [{k: r[k] for k in keep if k in r and not (isinstance(r[k], float) and np.isnan(r[k]))}
                 for r in rows],
    }
    return df, summary


def _sig(r) -> str:
    return r["signature"] if pd.isna(r.get("resolution")) else f"{r['resolution']}/{r['signature']}"


def summary_lines(df: pd.DataFrame) -> list[str]:
    lines = ["", "Subtype association of the picks (outcome-free; Kruskal-Wallis, one-vs-rest AUC)"]
    for _, r in df.iterrows():
        if pd.isna(r.get("kw_h")):
            lines.append(f"  {r['cohort']:<8} {r['model']:<13} {_sig(r):<23} not scored (coverage)")
            continue
        lines.append(f"  {r['cohort']:<8} {r['model']:<13} {_sig(r):<23} {r['analyses']:<40} "
                     f"H={r['kw_h']:.1f} p={r['kw_p']:.2g} eta2_H={r['eta2_h']:.2f}  "
                     f"{r['tracked_group']} AUC={r['auc_tracked']:.2f}")
    return lines
