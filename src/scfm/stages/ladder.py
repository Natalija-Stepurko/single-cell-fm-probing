"""Stage 05 — assemble the ladder and test the four predictions.

Takes stage 04's per-signature scores and produces the one table the study is about: for each
representation, its best signature's concordance read against null, floor, baseline and
reference — with bootstrap intervals over PATIENTS, because the patient is the unit of
independence here and a bootstrap over genes or cells would be pseudo-replication.

P1  FM above floor            P2  FM − baseline margin excludes zero
P3  FM prognostic within a PAM50 subtype (the reference does not already encode it)
P4  FM signatures passing P2 are enriched outside the HVG set

The pre-registered predictions are kept exactly as coded. Beside them (each from its own random
generator, so the pre-registered numbers do not move), the statistics the design text describes:
  P1  family-wise permutation p of the family maximum, the statistic the null was built for
  P2  floor-adjusted margin with a selection-aware patient bootstrap (re-selection per resample)
  P3  within-subtype C against within-subtype outcome permutations, for every representation
  P4  only when P2 passes: hypergeometric HVG depletion and a comparison with the baseline
plus the ladder at every (resolution, top-k) setting (ladder_by_setting.csv), the value each pick
adds to age and stage (added_value.csv), and a `deviations` list in ladder.json. A non-primary
endpoint reads results/translate_<endpoint> and writes results/ladder_<endpoint>.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params
from scfm.stages.translate import endpoint_dir, subtype_dummies
from scfm.stats import (
    TIE_TOL,
    cindex_many,
    cv_cox_cindex,
    cv_folds,
    fw_pvalue,
    lrt_added,
    mc_annotate,
    orient,
    permute_within_strata,
    selection_aware_margin,
)
from scfm.survival import bootstrap_patients, cindex, score_signature, zscore_genes


def best_per_model(scores: pd.DataFrame) -> pd.DataFrame:
    """Best signature per model by above-floor margin, not raw C-index."""
    return (scores.dropna(subset=["above_floor"])
            .sort_values("above_floor", ascending=False)
            .groupby("model", as_index=False).head(1).set_index("model"))


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    return float(o)


def _key(res, sid) -> str:
    return f"{res}/{sid}"


def _family(scores: pd.DataFrame, m: str) -> pd.DataFrame:
    f = scores[scores["model"] == m].copy()
    f["key"] = [_key(r, s) for r, s in zip(f["resolution"].astype(str), f["signature"])]
    return f.reset_index(drop=True)


def family_wise(fam: pd.DataFrame, null, pick_c: float, col: str, af_col: str) -> dict | None:
    """Family maximum of `col`, every signature attaining it, and its family-wise p; also the
    p of the margin-selected pick against the same null. P1_corrected needs the family-wise p
    below alpha and every maximising signature above its own floor."""
    if null is None:
        return None
    stat = float(fam[col].max())
    arg = fam[fam[col] >= stat - TIE_TOL]
    p = fw_pvalue(stat, null)
    pp = fw_pvalue(pick_c, null)
    af = arg[af_col].astype(float).tolist()
    mc, mcp = mc_annotate(p, len(null), C.ALPHA), mc_annotate(pp, len(null), C.ALPHA)
    return {"stat": stat, "argmax": arg["key"].tolist(),
            "argmax_top_cell_type": arg["top_cell_type"].tolist(),
            "argmax_above_floor": af, "fw_p": p, "fw_p_mc_se": mc["mc_se"],
            "fw_p_borderline": mc["borderline"],
            "pick_stat": float(pick_c), "pick_fw_p": pp, "pick_fw_p_mc_se": mcp["mc_se"],
            "n_perm": len(null),
            "P1_floor_condition": "every maximising signature above its floor",
            "P1_corrected": bool(p < C.ALPHA and all(v > 0 for v in af))}


def _boot_chunk(S, T, E, idxs):
    return np.vstack([cindex_many(S[:, ix], T[ix], E[ix]) for ix in idxs])


def bootstrap_all(S, T, E, n_boot, seed):
    """C of every signature (rows of S) in each patient bootstrap resample: (n_boot, n_sig)."""
    from joblib import Parallel, delayed
    rng = np.random.default_rng(seed)
    idxs = [rng.integers(0, len(T), len(T)) for _ in range(n_boot)]
    parts = [p for p in np.array_split(np.arange(n_boot), C.N_JOBS) if len(p)]
    res = Parallel(n_jobs=C.N_JOBS, max_nbytes=None)(
        delayed(_boot_chunk)(S, T, E, [idxs[i] for i in p]) for p in parts)
    return np.vstack(res)


def p2_corrected(fams, Cb, fm, base, af_col, floor_col, oriented: bool) -> dict:
    """Floor-adjusted margin of the FM pick over the baseline pick with the selection redone in
    every resample (floors held at their full-cohort values)."""
    def prep(m):
        f = fams[m]
        cb = orient(Cb[m], f["direction"].to_numpy()[None, :]) if oriented else Cb[m]
        return cb, f[floor_col].to_numpy(float)
    cf, ff = prep(fm)
    cb, fb = prep(base)
    margins, jf, jb = selection_aware_margin(cf, ff, cb, fb)
    pick_f = int(fams[fm][af_col].idxmax())
    pick_b = int(fams[base][af_col].idxmax())
    point = float(fams[fm][af_col].max() - fams[base][af_col].max())
    lo, hi = np.percentile(margins, [2.5, 97.5])
    return {"margin": point, "ci": [float(lo), float(hi)], "P2_pass": bool(lo > 0),
            "boot_mean": float(margins.mean()), "n_boot": len(margins),
            "ci_definition": "two-sided 95% percentile interval; P2_pass is the one-sided test lower > 0, "
                             "and an upper bound near 0 is not evidence that the FM is worse",
            "pick_reselected_frac": {fm: float(np.mean(jf == pick_f)), base: float(np.mean(jb == pick_b))},
            "floors": "held at their full-cohort values"}


def p3_corrected(picks: dict, clin, Z, gidx, ev, tm, n_perm, seed) -> dict:
    """Within-subtype C of each pick (oriented by its full-cohort direction), max over subtypes
    with at least P3_MIN_EVENTS events; null from outcomes permuted within each subtype."""
    sub = clin.dropna(subset=["subtype", ev, tm])
    tested = [st for st, g in sub.groupby("subtype") if g[ev].sum() >= C.P3_MIN_EVENTS]
    sub = sub[sub["subtype"].isin(tested)]
    cols = [Z.columns.get_loc(p) for p in sub.index]
    Zs = Z.to_numpy()[:, cols]
    names = list(picks)
    S = np.vstack([Zs[[gidx[g] for g in picks[k]["genes"]]].mean(axis=0) for k in names])
    d = np.array([picks[k]["direction"] for k in names])
    T, E, strata = sub[tm].to_numpy(float), sub[ev].to_numpy(float), sub["subtype"].to_numpy()
    masks = {st: strata == st for st in tested}

    def within(t, e):
        return np.vstack([orient(cindex_many(S[:, mk], t[mk], e[mk]), d) for mk in masks.values()])

    obs = within(T, E)
    rng = np.random.default_rng(seed)
    null = np.empty((n_perm, len(names)))
    for b in range(n_perm):
        q = permute_within_strata(strata, rng)
        null[b] = within(T[q], E[q]).max(axis=0)
    out = {}
    for i, k in enumerate(names):
        stat = float(obs[:, i].max())
        p = fw_pvalue(stat, null[:, i])
        mc = mc_annotate(p, n_perm, C.ALPHA)
        out[k] = {"within_subtype_cindex": dict(zip(tested, obs[:, i].tolist())), "stat": stat,
                  "argmax_subtype": tested[int(obs[:, i].argmax())], "p": p, "P3_pass": bool(p < C.ALPHA),
                  "p_mc_se": mc["mc_se"], "borderline": mc["borderline"],
                  "null_p95": float(np.percentile(null[:, i], 95)),
                  "rule_as_coded": bool((obs[:, i] > C.P3_WITHIN_SUBTYPE_CINDEX).any())}
    return {"subtypes": tested, "n_patients": len(sub), "n_perm": n_perm,
            "note": "pick and direction fixed on the full cohort, which contains these patients",
            "picks": out}


def p4_corrected(fams, pick_row, fm, base, p2_pass: bool, universe: dict) -> dict:
    """Only when P2 passes: HVG depletion of the pick (hypergeometric lower tail: k HVGs among its
    n genes, drawn from a universe of M scored genes of which K are HVGs) and the pick's HVG
    fraction against the baseline's signatures (share of baseline signatures at or below it)."""
    from scipy.stats import hypergeom
    if not p2_pass:
        return {"P4_corrected": None, "status": "not applicable (P2 failed)"}
    k, n = int(round(pick_row["hvg_frac"] * pick_row["topk"])), int(pick_row["topk"])
    M, K = universe.get("universe_n"), universe.get("universe_n_hvg")
    hp = float(hypergeom.cdf(k, M, K, n)) if M and K else float("nan")
    bh = fams[base]["hvg_frac"].to_numpy(float)
    rank_p = float((1 + np.sum(bh <= pick_row["hvg_frac"])) / (1 + len(bh)))
    return {"P4_corrected": bool(hp < C.ALPHA and rank_p < C.ALPHA), "status": "tested",
            "hvg_in_pick": k, "topk": n, "expected": n * K / M if M else None,
            "hypergeom_depletion_p": hp, "baseline_rank_p": rank_p,
            "baseline_rank_definition": "(1 + #baseline signatures with hvg_frac <= the pick's) / "
                                        "(1 + #baseline signatures)"}


def ladder_by_setting(fams, nullmats, base) -> pd.DataFrame:
    """The ladder at every (resolution, top-k) setting, with the family-wise p of the setting's
    maximum against its own sub-family null when the permutation matrix is available."""
    rows = []
    for m, f in fams.items():
        for (res, k), g in f.groupby(["resolution", "topk"]):
            a = g.loc[g["above_floor"].idxmax()]
            o = g.loc[g["above_floor_or"].idxmax()]
            r = {"model": m, "resolution": res, "topk": k, "n_signatures": len(g),
                 "signature": a["signature"], "top_cell_type": a["top_cell_type"], "cindex": a["cindex"],
                 "floor_mean": a["floor_mean"], "above_floor": a["above_floor"],
                 "signature_or": o["signature"], "top_cell_type_or": o["top_cell_type"],
                 "direction_or": "risk" if o["direction"] == 1 else "protective",
                 "cindex_or": o["cindex_or"], "floor_or_mean": o["floor_or_mean"],
                 "above_floor_or": o["above_floor_or"]}
            nm = nullmats.get(m)
            if nm is not None:
                cols = [nm["index"][key] for key in g["key"]]
                sub = nm["cindex"][:, cols]
                r["setting_max_cindex"] = float(g["cindex"].max())
                r["setting_max_signature"] = ";".join(
                    g.loc[g["cindex"] >= r["setting_max_cindex"] - TIE_TOL, "signature"])
                r["setting_fw_p"] = fw_pvalue(r["setting_max_cindex"], sub.max(axis=1))
                r["setting_fw_p_mc_se"] = mc_annotate(r["setting_fw_p"], len(sub))["mc_se"]
                r["setting_max_cindex_or"] = float(g["cindex_or"].max())
                r["setting_max_signature_or"] = ";".join(
                    g.loc[g["cindex_or"] >= r["setting_max_cindex_or"] - TIE_TOL, "signature"])
                r["setting_fw_p_or"] = fw_pvalue(r["setting_max_cindex_or"],
                                                 np.maximum(sub, 1 - sub).max(axis=1))
                r["setting_fw_p_or_mc_se"] = mc_annotate(r["setting_fw_p_or"], len(sub))["mc_se"]
            rows.append(r)
    df = pd.DataFrame(rows)
    b = df[df["model"] == base].set_index(["resolution", "topk"])
    key = list(zip(df["resolution"], df["topk"]))
    df["above_floor_minus_baseline"] = df["above_floor"].to_numpy() - b.loc[key, "above_floor"].to_numpy()
    df["above_floor_or_minus_baseline"] = (df["above_floor_or"].to_numpy()
                                           - b.loc[key, "above_floor_or"].to_numpy())
    return df


SELECTION_RULES = {
    "primary": "largest above-floor margin (risk direction) over the family",
    "sensitivity": "largest oriented above-floor margin over the family, direction re-read",
    "family_max_primary": "largest C over the family",
    "family_max_sensitivity": "largest max(C, 1 - C) over the family",
    "reference": "pre-specified",
}


def select_in_fold(c, floor_mean, rule: str) -> int:
    """The signature a selection rule picks from training-fold C-indices `c` (floors held at their
    full-cohort values; an oriented floor of a protective signature is 1 - floor_mean)."""
    c = np.asarray(c, dtype=float)
    floor_mean = np.asarray(floor_mean, dtype=float)
    if rule == "primary":
        v = c - floor_mean
    elif rule == "sensitivity":
        risk = c >= 0.5
        v = np.where(risk, c, 1 - c) - np.where(risk, floor_mean, 1 - floor_mean)
    elif rule == "family_max_primary":
        v = c
    elif rule == "family_max_sensitivity":
        v = np.maximum(c, 1 - c)
    else:
        raise ValueError(rule)
    return int(np.nanargmax(v))


def _zscore(x):
    x = pd.Series(x)
    return (x - x.mean()) / x.std()


def _cv_lp_c(d, tm, ev, cols_of_fold, folds) -> np.ndarray:
    """Per repeat, Harrell's C of the pooled out-of-fold Cox linear predictor; the covariate list
    may differ per fold (cols_of_fold(repeat, fold, train_index) -> list of column names)."""
    from scfm.stats import cox_fit
    T, E = d[tm].to_numpy(float), d[ev].to_numpy(float)
    out = []
    for r, rep in enumerate(folds):
        lp = np.empty(len(d))
        for f, (tr, te) in enumerate(rep):
            cols = cols_of_fold(r, f, tr)
            cph = cox_fit(d.iloc[tr], tm, ev, cols)
            lp[te] = np.asarray(cph.predict_log_partial_hazard(d.iloc[te][cols]), dtype=float)
        out.append(float(cindex_many(lp[None, :], T, E)[0]))
    return np.array(out)


def added_value(items: list[dict], fams: dict, sig_genes, clin, zexpr, ev, tm) -> pd.DataFrame:
    """What each score adds to age + stage (and to age + stage + PAM50): HR per SD of the score,
    nominal likelihood-ratio p, and the cross-validated change in C on paired folds, twice: with the
    full-cohort signature fixed, and nested (the signature re-selected, and re-oriented, by its own
    rule inside every training fold). One row per role of a signature."""
    rng = np.random.default_rng(C.SEED_ADDED_VALUE)
    rows = []
    full = clin.dropna(subset=C.COVARIATES + [ev, tm])
    ref, dummies = subtype_dummies(clin.dropna(subset=C.COVARIATES + ["subtype", ev, tm]))
    roles = {}
    for it in items:
        roles.setdefault((it["model"], it["signature_key"]), []).append(it["analysis"])
    for base_name, d, base in [("age+stage", full, list(C.COVARIATES)),
                               ("age+stage+PAM50", ref, list(C.COVARIATES) + dummies)]:
        d = d.copy()
        T, E = d[tm].to_numpy(float), d[ev].to_numpy(float)
        folds = cv_folds(d[ev], C.N_CV_SPLITS, C.N_CV_REPEATS, rng)
        c_base = cv_cox_cindex(d, tm, ev, {"base": base}, folds)["base"]
        Zd = zexpr.loc[:, d.index].to_numpy(float)
        gidx = {g: i for i, g in enumerate(zexpr.index)}
        fam_S, fam_c = {}, {}

        def family_scores(m):
            if m not in fam_S:
                fam_S[m] = np.vstack([Zd[[gidx[g] for g in sig_genes(m, k)]].mean(axis=0)
                                      for k in fams[m]["key"]])
                fam_c[m] = [[cindex_many(fam_S[m][:, tr], T[tr], E[tr]) for tr, _ in rep] for rep in folds]
            return fam_S[m], fam_c[m]

        fixed, nested = {}, {}
        for it in items:
            m, key, rule = it["model"], it["signature_key"], it["analysis"]
            if (m, key) not in fixed:
                sc = score_signature(zexpr, it["genes"]).reindex(d.index)
                d["score_z"] = _zscore(sc)
                p, cph = lrt_added(d, tm, ev, base, ["score_z"])
                row = cph.summary.loc["score_z"]
                c_full = cv_cox_cindex(d, tm, ev, {"full": base + ["score_z"]}, folds)["full"]
                fixed[(m, key)] = {"hr_per_sd": float(np.exp(row["coef"])),
                                   "hr_ci_lo": float(np.exp(row["coef lower 95%"])),
                                   "hr_ci_hi": float(np.exp(row["coef upper 95%"])), "lrt_p_nominal": p,
                                   "c_full": c_full}
            fx = fixed[(m, key)]
            if rule == "reference":
                nest = {"c": fx["c_full"], "agree": np.nan, "distinct": 0}
            elif (m, rule) in nested:
                nest = nested[(m, rule)]
            else:
                S, cf = family_scores(m)
                fl = fams[m]["floor_mean"].to_numpy(float)
                picks = [[select_in_fold(cf[r][f], fl, rule) for f in range(len(rep))]
                         for r, rep in enumerate(folds)]
                for j in {j for rp in picks for j in rp}:
                    d[f"nested_{j}"] = _zscore(S[j]).to_numpy()
                c_nest = _cv_lp_c(d, tm, ev, lambda r, f, tr: base + [f"nested_{picks[r][f]}"], folds)
                d = d.drop(columns=[c for c in d.columns if c.startswith("nested_")])
                flat = [j for rp in picks for j in rp]
                full_keys = [k for k, (mm, rr) in ((it2["signature_key"], (it2["model"], it2["analysis"]))
                                                   for it2 in items) if (mm, rr) == (m, rule)]
                agree = np.mean([fams[m]["key"].iloc[j] in full_keys for j in flat])
                nest = nested[(m, rule)] = {"c": c_nest, "agree": float(agree), "distinct": len(set(flat))}
            dfx, dnest = fx["c_full"] - c_base, nest["c"] - c_base
            rows.append({"analysis": rule, "analyses": ";".join(roles[(m, key)]),
                         "model": it["model"], "signature": it["signature"],
                         "resolution": it.get("resolution"), "top_cell_type": it.get("top_cell_type"),
                         "direction_full_cohort": it["direction_label"],
                         "selected_on_outcome": rule != "reference", "selection_rule": SELECTION_RULES[rule],
                         "base": base_name, "n": len(d), "events": int(d[ev].sum()),
                         "hr_per_sd": fx["hr_per_sd"], "hr_ci_lo": fx["hr_ci_lo"], "hr_ci_hi": fx["hr_ci_hi"],
                         "lrt_p_nominal": fx["lrt_p_nominal"],
                         "cindex_base_cv": float(c_base.mean()),
                         "cindex_with_score_cv_fixed": float(fx["c_full"].mean()),
                         "delta_cindex_cv_fixed_mean": float(dfx.mean()),
                         "delta_cindex_cv_fixed_lo": float(np.percentile(dfx, 2.5)),
                         "delta_cindex_cv_fixed_hi": float(np.percentile(dfx, 97.5)),
                         "cindex_with_score_cv_nested": float(nest["c"].mean()),
                         "delta_cindex_cv_nested_mean": float(dnest.mean()),
                         "delta_cindex_cv_nested_lo": float(np.percentile(dnest, 2.5)),
                         "delta_cindex_cv_nested_hi": float(np.percentile(dnest, 97.5)),
                         "nested_pick_agreement": nest["agree"], "nested_distinct_picks": nest["distinct"]})
    return pd.DataFrame(rows)


DEVIATIONS = [
    {"id": "reference_rung",
     "design": "Reference rung: the PAM50 subtype call and a published prognostic signature.",
     "as_coded": "One in-sample Cox model of PAM50 subtype + age + stage on the 821 complete-case "
                 "patients (a smaller cohort than the signatures'); no published signature.",
     "now_reported": "Named references on one complete-case set (clinical, pam50, pam50_clinical), "
                     "each in sample and cross-validated; clinical_full (age + stage, all patients "
                     "with both); the PAM50 11-gene proliferation score with its own floor and "
                     "permutation p. The as-coded value is kept as reference_pam50."},
    {"id": "P1_statistic",
     "design": "Family-wise permutation null: the distribution of the best signature's C-index, so "
               "best-of-k is controlled at the ladder level.",
     "as_coded": "The signature with the largest above-floor margin was compared with the 95th "
                 "percentile of the null of the family maximum.",
     "now_reported": "The family maximum (C; oriented C in the sensitivity analysis) with its "
                     "family-wise p; P1_corrected = fw_p < 0.05 and the maximising signature above "
                     "its floor. The margin-selected pick's p against the same null is also given."},
    {"id": "P2_margin",
     "design": "Test against Baseline, both read as their distance above Floor; P2 passes if the "
               "bootstrap interval of the margin excludes zero.",
     "as_coded": "Raw C of the FM pick minus raw C of the baseline pick, both picks held fixed in "
                 "every bootstrap resample; floors not used.",
     "now_reported": "Above-floor margin (FM pick minus baseline pick) with a patient bootstrap that "
                     "re-selects each pick in every resample, floors held at full-cohort values; "
                     "P2_corrected = lower 2.5% bound > 0."},
    {"id": "P3_null",
     "design": "At least one FM-derived state is prognostic within a PAM50 subtype.",
     "as_coded": "Within-subtype C of the one margin-selected FM signature above 0.6 in any subtype "
                 "with at least 10 events; no null, not computed for the baseline.",
     "now_reported": "Max over subtypes (Basal, Her2, LumA, LumB) of the within-subtype C, against "
                     "1,000 within-subtype outcome permutations, for every representation including "
                     "the baseline; P3_corrected = p < 0.05. Pick and direction were fixed on the "
                     "full cohort, which contains these patients."},
    {"id": "monte_carlo_error",
     "design": "Decisions at alpha = 0.05.",
     "as_coded": "Permutation p-values from 500 or 1,000 permutations, without Monte-Carlo error.",
     "now_reported": "Every permutation p carries its Monte-Carlo SE and a borderline flag (within 2 SE "
                     "of 0.05); P3 uses 10,000 permutations and the selection-aware P2 bootstrap 1,000 "
                     "resamples (both their own streams). The family-wise null keeps its pre-registered "
                     "500 permutations."},
    {"id": "P4_gate",
     "design": "FM signatures that pass P2 are enriched for genes outside the HVG set.",
     "as_coded": "hvg_frac < 0.5 for the FM pick, evaluated whether or not P2 passed; no "
                 "enrichment test.",
     "now_reported": "Not applicable unless P2_corrected passes; when it does, a hypergeometric "
                     "HVG-depletion test of the pick against the scored-gene universe and the pick's HVG "
                     "fraction ranked among the baseline's signatures (baseline_rank_p)."},
    {"id": "per_setting_ladder",
     "design": "Cluster resolution and marker count are swept and the ladder is reported at each "
               "setting.",
     "as_coded": "One pick per representation across all nine settings, chosen by a survival-based "
                 "margin.",
     "now_reported": "ladder_by_setting.csv: the pick at every (resolution, top-k), with the "
                     "setting's family-wise p from its own sub-family null."},
    {"id": "covariate_adjustment",
     "design": "Cox models adjusted for age and stage in the primary model.",
     "as_coded": "Every ladder C-index is Harrell's C of the unadjusted score; the adjusted Cox model "
                 "(lifelines penalizer 0.01, Wald p, per unit of the mean-z score) fed only the hr and p "
                 "columns of scores.csv.",
     "now_reported": "added_value.csv: HR per SD, likelihood-ratio p and cross-validated change in C "
                     "of each pick over age + stage (n with both) and over age + stage + PAM50, one row "
                     "per role. Every pick and family maximum was selected, and its direction set, on the "
                     "outcomes of the same patients (selected_on_outcome); for those rows lrt_p_nominal is "
                     "not adjusted for selection and the fixed-signature CV change in C is optimistic. The "
                     "nested CV columns redo the selection (and orientation) inside every training fold; "
                     "only the pre-specified proliferation reference is free of selection. The lo/hi "
                     "columns span fold repeats on the same patients, not a confidence interval."},
    {"id": "selection_conditional_intervals",
     "design": "Bootstrap intervals over patients for the ladder.",
     "as_coded": "ladder.csv ci_lo/ci_hi and P2_margin_ci_as_coded bootstrap the picked signatures with "
                 "the pick held fixed, so they ignore the selection over {family_sizes} signatures "
                 "per representation and are conditional on it; P2_margin_as_coded is the bootstrap mean of "
                 "the raw C difference, not a full-cohort point estimate.",
     "now_reported": "Kept as coded; P2_point_as_coded gives the full-cohort raw C difference, and the "
                     "selection-aware P2 interval re-selects the picks in every resample."},
    {"id": "clustering_stability",
     "design": "Cell-level clustering is done once; its stability is a separate diagnostic.",
     "as_coded": "The stability diagnostic was not run.",
     "now_reported": "Not run; the states are reported at three resolutions and nine settings instead "
                     "(ladder_by_setting.csv)."},
    {"id": "reference_penalty",
     "design": "Reference rung: PAM50 subtype.",
     "as_coded": "Cox reference fits use lifelines penalizer 0.01 (a ridge of 0.01 x n on standardised "
                 "covariates).",
     "now_reported": "Each reference entry gives the penalizer, the in-sample C of the penalised fit and of "
                     "the unpenalised fit (they differ for PAM50 alone, whose LumA and Basal hazards are "
                     "nearly equal), and the cross-validated C, which is the figure to compare with."},
    {"id": "secondary_endpoint",
     "design": "PFI is the secondary endpoint; both reported.",
     "as_coded": "Only OS was analysed.",
     "now_reported": "PFI run through translate and ladder into results/translate_pfi and "
                     "results/ladder_pfi with the same statistics."},
]


AS_CODED_P2 = ("P2_margin_as_coded: bootstrap mean of raw C(FM pick) - C(baseline pick), picks fixed; "
               "P2_point_as_coded: the same difference on the full cohort; margin: the floor-adjusted "
               "full-cohort point estimate")


def corrected_analyses(scores, meta, sigs, clin, zexpr, ev, tm, best, best_or, P, P_or,
                       translate_dir: Path, states_dir: Path, n_boot: int, out: Path) -> dict:
    base = "hvg_pca"
    models = list(dict.fromkeys(scores["model"]))
    fms = [m for m in models if C.MODELS.get(m, {}).get("kind") == "fm"]
    fams = {m: _family(scores, m) for m in models}
    genes = lambda m, key: sigs[m][key.split("/")[0]][key.split("/", 1)[1]]["genes"]  # noqa: E731
    ok = clin[[tm, ev]].notna().all(axis=1)
    T, E = clin.loc[ok, tm].to_numpy(float), clin.loc[ok, ev].to_numpy(float)
    gidx = {g: i for i, g in enumerate(zexpr.index)}
    Zok = zexpr.loc[:, clin.index[ok]].to_numpy(float)
    S = {m: np.vstack([Zok[[gidx[g] for g in genes(m, k)]].mean(axis=0) for k in f["key"]])
         for m, f in fams.items()}
    for m, f in fams.items():
        assert np.allclose(cindex_many(S[m], T, E), f["cindex"], atol=1e-9), f"{m}: C does not reproduce"
    summ_path = states_dir / "summary.json"
    summ = json.load(open(summ_path)) if summ_path.exists() else {}
    universe = {k: summ.get(k) for k in ("universe_n", "universe_n_hvg")}
    nullmats = {}
    for m in fams:
        p = translate_dir / f"null_matrix_{m}.npz"
        if p.exists():
            z = np.load(p)
            nullmats[m] = {"cindex": z["cindex"], "index": {k: i for i, k in enumerate(z["signatures"])}}

    pick_key = lambda r: _key(str(r["resolution"]), r["signature"])  # noqa: E731
    n_sig = [len(f) for f in fams.values()]
    sizes = f"{min(n_sig)}-{max(n_sig)}"
    fw = {"primary": {}, "sensitivity": {}}
    for m in models:
        nm = meta["nulls"][m]
        fw["primary"][m] = family_wise(fams[m], nm.get("best_cindex_null"), best.loc[m, "cindex"],
                                       "cindex", "above_floor")
        if best_or is not None:
            fw["sensitivity"][m] = family_wise(fams[m], nm.get("best_oriented_null"),
                                               best_or.loc[m, "cindex_or"], "cindex_or", "above_floor_or")

    print("  selection-aware bootstrap …", flush=True)
    Sall = np.vstack([S[m] for m in models])
    Cb_all = bootstrap_all(Sall, T, E, C.N_P2_BOOT, C.SEED_P2_BOOT)
    bounds = np.cumsum([0] + [len(S[m]) for m in models])
    Cb = {m: Cb_all[:, bounds[i]:bounds[i + 1]] for i, m in enumerate(models)}
    p2 = {"primary": {}, "sensitivity": {}}
    for m in fms:
        p2["primary"][m] = p2_corrected(fams, Cb, m, base, "above_floor", "floor_mean", False)
        p2["primary"][m].update(P2_margin_as_coded=P[m]["P2_margin_over_baseline"],
                                P2_margin_ci_as_coded=P[m]["P2_margin_ci"],
                                P2_point_as_coded=float(best.loc[m, "cindex"] - best.loc[base, "cindex"]),
                                as_coded_definition=AS_CODED_P2)
        if best_or is not None:
            p2["sensitivity"][m] = p2_corrected(fams, Cb, m, base, "above_floor_or", "floor_or_mean", True)
            p2["sensitivity"][m].update(
                P2_margin_as_coded=P_or[m]["P2_margin_over_baseline"],
                P2_margin_ci_as_coded=P_or[m]["P2_margin_ci"],
                P2_point_as_coded=float(best_or.loc[m, "cindex_or"] - best_or.loc[base, "cindex_or"]),
                as_coded_definition=AS_CODED_P2)

    print("  within-subtype permutations …", flush=True)
    picks = {}
    for m, r in best.iterrows():
        picks[f"primary/{m}"] = {"genes": genes(m, pick_key(r)), "direction": 1}
    if best_or is not None:
        for m, r in best_or.iterrows():
            picks[f"sensitivity/{m}"] = {"genes": genes(m, pick_key(r)), "direction": int(r["direction"])}
    p3 = p3_corrected(picks, clin.loc[ok], zexpr, gidx, ev, tm, C.N_P3_PERM, C.SEED_P3_PERM)

    p4 = {"primary": {}, "sensitivity": {}}
    for m in fms:
        p4["primary"][m] = p4_corrected(fams, best.loc[m], m, base, p2["primary"][m]["P2_pass"], universe)
        if best_or is not None:
            p4["sensitivity"][m] = p4_corrected(fams, best_or.loc[m], m, base,
                                                p2["sensitivity"][m]["P2_pass"], universe)

    verdicts = {"primary": {}, "sensitivity": {}}
    for an, preds in (("primary", P), ("sensitivity", P_or)):
        for m in fms:
            if m not in preds:
                continue
            fwm = fw[an].get(m)
            p3m = p3["picks"][f"{an}/{m}"]
            verdicts[an][m] = {
                "P1_as_coded": preds[m]["P1_above_floor"],
                "P1_corrected": None if fwm is None else fwm["P1_corrected"],
                "P1_corrected_borderline": None if fwm is None else fwm["fw_p_borderline"],
                "P2_as_coded": preds[m]["P2_pass"], "P2_corrected": p2[an][m]["P2_pass"],
                "P3_as_coded": preds[m]["P3_pass"], "P3_corrected": p3m["P3_pass"],
                "P3_corrected_borderline": p3m["borderline"],
                "P4_as_coded": preds[m]["P4_pass"], "P4_corrected": p4[an][m]["P4_corrected"]}

    by_setting = ladder_by_setting(fams, nullmats, base)
    by_setting.to_csv(out / "ladder_by_setting.csv", index=False)

    print("  added value over age and stage …", flush=True)
    items = []

    def add(analysis, m, key, direction):
        f = fams[m].set_index("key").loc[key]
        items.append({"analysis": analysis, "model": m, "signature": key.split("/", 1)[1],
                      "signature_key": key, "resolution": key.split("/")[0],
                      "top_cell_type": f["top_cell_type"],
                      "genes": genes(m, key), "direction_label": "risk" if direction == 1 else "protective"})
    for m, r in best.iterrows():
        add("primary", m, pick_key(r), 1 if r["cindex"] >= 0.5 else -1)
    if best_or is not None:
        for m, r in best_or.iterrows():
            add("sensitivity", m, pick_key(r), int(r["direction"]))
    for an in fw:
        for m, v in fw[an].items():
            for key in (v or {}).get("argmax", []):
                c = float(fams[m].set_index("key").loc[key, "cindex"])
                add(f"family_max_{an}", m, key, 1 if c >= 0.5 else -1)
    prolif = meta.get("reference", {}).get("proliferation")
    if prolif:
        items.append({"analysis": "reference", "model": "proliferation", "signature": "PAM50 proliferation",
                      "signature_key": "reference", "genes": list(prolif["genes_used"].values()),
                      "direction_label": "risk" if prolif["cindex"] >= 0.5 else "protective"})
    av = added_value(items, fams, genes, clin, zexpr, ev, tm)
    av.to_csv(out / "added_value.csv", index=False)

    return {"added_value": av, "reference": meta.get("reference", {}), "family_wise": fw, "p2_corrected": p2,
            "p3_corrected": p3, "p4_corrected": p4, "verdicts": verdicts, "universe": universe,
            "monte_carlo": {"p2_bootstrap_resamples": C.N_P2_BOOT, "p3_permutations": C.N_P3_PERM,
                            "borderline": "a permutation p within 2 Monte-Carlo SE of alpha"},
            "seeds": {"p2_bootstrap": C.SEED_P2_BOOT, "p3_permutations": C.SEED_P3_PERM,
                      "added_value_cv": C.SEED_ADDED_VALUE, "reference_cv": C.SEED_REF_CV,
                      "reference_floor": C.SEED_REF_FLOOR, "reference_perm": C.SEED_REF_PERM},
            "deviations": [{k: v.format(family_sizes=sizes) for k, v in dv.items()} for dv in DEVIATIONS]}


def corrected_summary(corr: dict, av: pd.DataFrame | None = None) -> list[str]:
    lines = ["", "Corrected statistics (beside the pre-registered ones above)"]
    ref = corr["reference"]
    for k in ("clinical", "pam50", "pam50_clinical", "clinical_full"):
        if k in ref:
            r = ref[k]
            mle = r.get("cindex_in_sample_mle", np.nan)
            lines.append(f"  reference {k:<15} n={r['n']} events={r['events']}  "
                         f"in-sample C={r['cindex_in_sample']:.3f} (ridge) {mle:.3f} (MLE)  "
                         f"CV C={r['cindex_cv_mean']:.3f} (sd over fold repeats {r['cindex_cv_sd']:.3f})")
    if "proliferation" in ref:
        r = ref["proliferation"]
        lines.append(f"  reference proliferation   C={r['cindex']:.3f} floor={r['floor_mean']:.3f} "
                     f"perm p={r['perm_p']:.3g} (MC se {r.get('perm_p_mc_se', np.nan):.3f}, "
                     f"{r['n_perm']} perms)")
    for an in ("primary", "sensitivity"):
        for m, v in corr["family_wise"][an].items():
            if v:
                lines.append(f"  {an:<11} {m:<11} family max {v['stat']:.3f} ({', '.join(v['argmax'])}) "
                             f"fw p={v['fw_p']:.3f} (MC se {v['fw_p_mc_se']:.3f}"
                             f"{', borderline' if v['fw_p_borderline'] else ''}); "
                             f"pick p={v['pick_fw_p']:.3f}; P1_corrected={v['P1_corrected']}")
        for m, v in corr["p2_corrected"][an].items():
            lines.append(f"  {an:<11} {m:<11} P2 above-floor margin {v['margin']:+.3f} "
                         f"[{v['ci'][0]:+.3f},{v['ci'][1]:+.3f}] (selection-aware, B={v['n_boot']}) "
                         f"P2_corrected={v['P2_pass']}")
    for k, v in corr["p3_corrected"]["picks"].items():
        lines.append(f"  P3 {k:<23} max within-subtype C {v['stat']:.3f} ({v['argmax_subtype']}) "
                     f"p={v['p']:.4f} (MC se {v['p_mc_se']:.4f}{', borderline' if v['borderline'] else ''})")
    if av is not None:
        lines += ["", "Added value over age + stage (selected rows: nominal LRT p; nested CV redoes the "
                      "selection in each training fold)"]
        for _, r in av[av["base"] == "age+stage"].iterrows():
            lines.append(f"  {r['analysis']:<23} {r['model']:<13} {r['signature']:<19} "
                         f"HR/SD {r['hr_per_sd']:.2f}  LRT p {r['lrt_p_nominal']:.3g}  "
                         f"dC fixed {r['delta_cindex_cv_fixed_mean']:+.4f}  "
                         f"nested {r['delta_cindex_cv_nested_mean']:+.4f}")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--endpoint", default=None, choices=list(C.ENDPOINTS),
                    help="selects the default translate directory (default: the primary endpoint); the "
                         "output directory follows the endpoint recorded in the translate output")
    ap.add_argument("--translate-dir", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--n-boot", type=int, default=C.N_BOOTSTRAP)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    args.translate_dir = args.translate_dir or str(
        C.RESULTS / endpoint_dir("translate", args.endpoint or C.PRIMARY_ENDPOINT))
    if args.dry_run:
        print(json.dumps({"n_boot": args.n_boot, "unit": "patient"}, indent=2)); return

    rng = np.random.default_rng(args.seed)
    scores = pd.read_csv(Path(args.translate_dir) / "scores.csv")
    meta = json.load(open(Path(args.translate_dir) / "nulls.json"))
    if args.endpoint is not None and args.endpoint != meta["endpoint"]:
        raise SystemExit(f"--endpoint {args.endpoint} but {args.translate_dir} holds {meta['endpoint']}")
    args.endpoint = meta["endpoint"]
    args.out_dir = args.out_dir or str(C.RESULTS / endpoint_dir("ladder", meta["endpoint"]))
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    sigs = json.load(open(Path(args.states_dir) / "signatures.json"))
    ev, tm = C.ENDPOINTS[meta["endpoint"]]
    expr = pd.read_parquet(Path(args.data_dir) / "bulk_expr.parquet")
    clin = pd.read_csv(Path(args.data_dir) / "bulk_clinical.csv", index_col=0)
    zexpr = zscore_genes(expr)

    best = best_per_model(scores)
    genes_of = lambda r: sigs[r.name][str(r["resolution"])][r["signature"]]["genes"]

    # patient bootstrap of every best signature at once, so margins are paired
    boot = {m: [] for m in best.index}
    for bs in bootstrap_patients(clin, args.n_boot, rng):
        for m, r in best.iterrows():
            boot[m].append(cindex(bs.assign(score=score_signature(zexpr, genes_of(r))),
                                  tm, ev, "score"))
    boot = {m: np.array(v) for m, v in boot.items()}
    lo, hi = lambda a: float(np.nanpercentile(a, 2.5)), lambda a: float(np.nanpercentile(a, 97.5))

    ladder = []
    for m, r in best.iterrows():
        ladder.append({"model": m, "kind": r["kind"], "signature": r["signature"],
                       "resolution": r["resolution"], "top_cell_type": r["top_cell_type"],
                       "cindex": r["cindex"],
                       "ci_lo": lo(boot[m]), "ci_hi": hi(boot[m]),
                       "floor_mean": r["floor_mean"], "above_floor": r["above_floor"],
                       "null_p95": meta["nulls"][m]["best_cindex_null_p95"],
                       "hvg_frac": r["hvg_frac"]})
    ladder = pd.DataFrame(ladder).set_index("model")
    ref = meta.get("reference", {}).get("pam50_cindex", np.nan)

    # predictions
    base = "hvg_pca"
    P = {}
    for m in [x for x in ladder.index if ladder.loc[x, "kind"] == "fm"]:
        margin = boot[m] - boot[base]
        P[m] = {
            "P1_above_floor": bool(ladder.loc[m, "above_floor"] > 0
                                   and ladder.loc[m, "cindex"] > ladder.loc[m, "null_p95"]),
            "P2_margin_over_baseline": float(np.nanmean(margin)),
            "P2_margin_ci": [lo(margin), hi(margin)],
            "P2_pass": bool(lo(margin) > 0),
            "P4_hvg_frac": float(ladder.loc[m, "hvg_frac"]),
            "P4_pass": bool(ladder.loc[m, "hvg_frac"] < C.P4_MAX_HVG_FRAC),
        }
        # P3: prognostic within each PAM50 subtype
        p3 = {}
        for st, sub in clin.dropna(subset=["subtype"]).groupby("subtype"):
            if sub[ev].sum() >= 10:
                p3[st] = cindex(sub.assign(score=score_signature(zexpr, genes_of(best.loc[m]))),
                                tm, ev, "score")
        P[m]["P3_within_subtype_cindex"] = p3
        P[m]["P3_pass"] = bool(any(v > C.P3_WITHIN_SUBTYPE_CINDEX for v in p3.values()))

    # ── sensitivity analysis: each signature read in the direction it acts ──────────────
    # added after the primary run; everything above is the pre-registered ladder, unchanged
    sens_ladder, P_or = None, {}
    if "above_floor_or" in scores.columns:
        best_or = (scores.dropna(subset=["above_floor_or"])
                   .sort_values("above_floor_or", ascending=False)
                   .groupby("model", as_index=False).head(1).set_index("model"))
        orient = lambda c, d: c if d == 1 else 1 - c      # direction fixed on the full cohort
        rng_or = np.random.default_rng(args.seed + 1)
        boot_or = {m: [] for m in best_or.index}
        for bs in bootstrap_patients(clin, args.n_boot, rng_or):
            for m, r in best_or.iterrows():
                boot_or[m].append(orient(cindex(bs.assign(score=score_signature(zexpr, genes_of(r))),
                                                tm, ev, "score"), r["direction"]))
        boot_or = {m: np.array(v) for m, v in boot_or.items()}
        rows_or = []
        for m, r in best_or.iterrows():
            rows_or.append({"model": m, "kind": r["kind"], "signature": r["signature"],
                            "resolution": r["resolution"], "top_cell_type": r["top_cell_type"],
                            "direction": "risk" if r["direction"] == 1 else "protective",
                            "cindex": r["cindex_or"], "ci_lo": lo(boot_or[m]), "ci_hi": hi(boot_or[m]),
                            "floor_mean": r["floor_or_mean"], "above_floor": r["above_floor_or"],
                            "null_p95": meta["nulls"][m]["best_oriented_null_p95"],
                            "hvg_frac": r["hvg_frac"]})
        sens_ladder = pd.DataFrame(rows_or).set_index("model")
        for m in [x for x in sens_ladder.index if sens_ladder.loc[x, "kind"] == "fm"]:
            margin = boot_or[m] - boot_or[base]
            p3 = {}
            for st, sub in clin.dropna(subset=["subtype"]).groupby("subtype"):
                if sub[ev].sum() >= 10:
                    p3[st] = orient(cindex(sub.assign(score=score_signature(zexpr, genes_of(best_or.loc[m]))),
                                           tm, ev, "score"), best_or.loc[m, "direction"])
            P_or[m] = {
                "P1_above_floor": bool(sens_ladder.loc[m, "above_floor"] > 0
                                       and sens_ladder.loc[m, "cindex"] > sens_ladder.loc[m, "null_p95"]),
                "P2_margin_over_baseline": float(np.nanmean(margin)),
                "P2_margin_ci": [lo(margin), hi(margin)],
                "P2_pass": bool(lo(margin) > 0),
                "P3_within_subtype_cindex": p3,
                "P3_pass": bool(any(v > C.P3_WITHIN_SUBTYPE_CINDEX for v in p3.values())),
                "P4_hvg_frac": float(sens_ladder.loc[m, "hvg_frac"]),
                "P4_pass": bool(sens_ladder.loc[m, "hvg_frac"] < C.P4_MAX_HVG_FRAC),
            }
        sens_ladder.to_csv(out / "ladder_sensitivity.csv")

    ladder.to_csv(out / "ladder.csv")
    corr = corrected_analyses(scores, meta, sigs, clin, zexpr, ev, tm, best,
                              best_or if sens_ladder is not None else None, P, P_or,
                              Path(args.translate_dir), Path(args.states_dir), args.n_boot, out)
    av = corr.pop("added_value")
    json.dump({"ladder": ladder.reset_index().to_dict("records"), "reference_pam50": ref,
               "predictions": P, "endpoint": meta["endpoint"],
               "sensitivity": None if sens_ladder is None else {
                   "note": "each signature read in the direction it acts; added after the primary run",
                   "ladder": sens_ladder.reset_index().to_dict("records"), "predictions": P_or},
               "n_patients": meta["n_patients"], "n_events": meta["n_events"], **corr},
              open(out / "ladder.json", "w"), indent=2, default=_json_default)

    lines = [f"Ladder — {meta['endpoint']}, {meta['n_patients']} patients, {meta['n_events']} events",
             f"reference (PAM50) C = {ref:.3f}" if np.isfinite(ref) else "reference: n/a", ""]
    for m, r in ladder.iterrows():
        lines.append(f"  {m:<11} C={r['cindex']:.3f} [{r['ci_lo']:.3f},{r['ci_hi']:.3f}]  "
                     f"floor={r['floor_mean']:.3f}  null95={r['null_p95']:.3f}  "
                     f"hvg={r['hvg_frac']:.2f}  {r['top_cell_type']}")
    lines.append("")
    for m, p in P.items():
        lines.append(f"  {m}: P1={p['P1_above_floor']} P2={p['P2_pass']} "
                     f"(margin {p['P2_margin_over_baseline']:+.3f} "
                     f"[{p['P2_margin_ci'][0]:+.3f},{p['P2_margin_ci'][1]:+.3f}]) "
                     f"P3={p['P3_pass']} P4={p['P4_pass']}")
    if sens_ladder is not None:
        lines += ["", "Sensitivity — each signature read in the direction it acts "
                      "(added after the primary run)"]
        for m, r in sens_ladder.iterrows():
            lines.append(f"  {m:<11} C={r['cindex']:.3f} [{r['ci_lo']:.3f},{r['ci_hi']:.3f}]  "
                         f"floor={r['floor_mean']:.3f}  null95={r['null_p95']:.3f}  "
                         f"hvg={r['hvg_frac']:.2f}  {r['direction']}  {r['top_cell_type']}")
        for m, p in P_or.items():
            lines.append(f"  {m}: P1={p['P1_above_floor']} P2={p['P2_pass']} "
                         f"(margin {p['P2_margin_over_baseline']:+.3f} "
                         f"[{p['P2_margin_ci'][0]:+.3f},{p['P2_margin_ci'][1]:+.3f}]) "
                         f"P3={p['P3_pass']} P4={p['P4_pass']}")
    lines += corrected_summary(corr, av)
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    record_params(out, args, extra={"predictions": P, "verdicts": corr["verdicts"]})


if __name__ == "__main__":
    main()
