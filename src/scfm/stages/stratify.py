"""Stage stratify — do a representation's cell states, used jointly, stratify patients better than
age + stage alone, and better than the baseline's states?

Patients with age and stage (and the outcome). Features per representation: the z-scored scores of
all its signatures with top-k = MARKER_TOPK_STRATIFY, at every resolution. Each model is scored by
repeated event-stratified K-fold cross-validation on one shared set of folds (so differences are
paired per repeat): the pooled out-of-fold risk of a repeat gives one Harrell's C.

  clinical            Cox, age + stage
  ridge_<model>       ridge Cox (l1_ratio 0), age + stage + states; the states' penalizer by inner
                      3-fold CV, age and stage penalised as in the clinical model
  xgb_<model>         gradient-boosted Cox (fixed hyperparameters), age + stage + states
  clinical_pam50      Cox, age + stage + PAM50 subtype, on the complete cases with a PAM50 call (own
                      folds; compared with `clinical` refitted on the same patients and folds)
  clinical_prolif     Cox, age + stage + PAM50 proliferation score

Output: results/stratify/cv_results.csv (one row per model and repeat), cv_summary.json.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params
from scfm.stats import cindex_many, cox_fit, cv_folds
from scfm.survival import score_signature, zscore_genes

BASELINE = "hvg_pca"
CLINICAL_PENALIZER = 0.01


def state_features(sigs: dict, model: str, zexpr: pd.DataFrame, patients, topk: int) -> pd.DataFrame:
    """Patients × signatures: mean-z score of every signature of `model` with this top-k, z-scored
    across `patients` (no outcome is used)."""
    cols = {}
    for res, d in sigs[model].items():
        for sid, s in d.items():
            if s["topk"] == topk:
                sc = score_signature(zexpr, s["genes"]).reindex(patients)
                if sc.notna().all() and sc.std() > 0:
                    cols[f"{model}:{res}/{sid}"] = (sc - sc.mean()) / sc.std()
    return pd.DataFrame(cols, index=patients)


def _fold_rng(r: int, f: int) -> np.random.Generator:
    return np.random.default_rng([C.SEED_STRATIFY, r, f])


def _cox_lp(d, tm, ev, cols, tr, te, penalizer=0.01):
    cph = cox_fit(d.iloc[tr], tm, ev, cols, penalizer=penalizer)   # lifelines l1_ratio 0: ridge
    return np.asarray(cph.predict_log_partial_hazard(d.iloc[te][cols]), dtype=float)


def ridge_penalties(cols, pen: float) -> np.ndarray:
    """Per-covariate lifelines penalizer: the clinical covariates keep the clinical model's 0.01, so
    with every state coefficient shrunk to zero the ridge model is the clinical model."""
    return np.array([CLINICAL_PENALIZER if c in C.COVARIATES else pen for c in cols])


def choose_penalizer(d, tm, ev, cols, tr, grid, rng) -> tuple[float, dict]:
    """Inner event-stratified 3-fold CV on the training patients: the penalizer whose pooled
    out-of-fold linear predictor has the largest Harrell's C (the smallest one on a tie)."""
    inner = d.iloc[tr]
    folds = cv_folds(inner[ev], C.STRATIFY_INNER_SPLITS, 1, rng)[0]
    T, E = inner[tm].to_numpy(float), inner[ev].to_numpy(float)
    score = {}
    for pen in grid:
        lp = np.empty(len(inner))
        for itr, ite in folds:
            lp[ite] = _cox_lp(inner, tm, ev, cols, itr, ite, penalizer=ridge_penalties(cols, pen))
        score[pen] = float(cindex_many(lp[None, :], T, E)[0])
    best = max(grid, key=lambda p: (score[p], -p))
    return best, score


def xgb_risk(X: np.ndarray, T: np.ndarray, E: np.ndarray, tr, te) -> np.ndarray:
    """Out-of-fold log-hazard from xgboost's Cox objective (label: +time if event, -time if censored)."""
    import xgboost as xgb
    y = np.where(E[tr] == 1, T[tr], -T[tr])
    m = xgb.XGBRegressor(objective="survival:cox", **C.STRATIFY_XGB, n_jobs=1)
    m.fit(X[tr], y)
    return m.predict(X[te], output_margin=True)


def _repeat(kind, d, tm, ev, cols, rep, r):
    """One repeat of one model: (pooled out-of-fold risk, per-fold chosen penalizers)."""
    lp = np.empty(len(d))
    chosen = []
    T, E = d[tm].to_numpy(float), d[ev].to_numpy(float)
    X = d[cols].to_numpy(float)
    for f, (tr, te) in enumerate(rep):
        if kind == "cox":
            lp[te] = _cox_lp(d, tm, ev, cols, tr, te)
        elif kind == "ridge":
            pen, _ = choose_penalizer(d, tm, ev, cols, tr, C.STRATIFY_PENALIZERS, _fold_rng(r, f))
            chosen.append(pen)
            lp[te] = _cox_lp(d, tm, ev, cols, tr, te, penalizer=ridge_penalties(cols, pen))
        elif kind == "xgb":
            lp[te] = xgb_risk(X, T, E, tr, te)
        else:
            raise ValueError(kind)
    return float(cindex_many(lp[None, :], T, E)[0]), chosen


def run_cv(specs: dict, d, tm, ev, folds, n_jobs: int) -> dict:
    """specs: name -> (kind, columns). Returns name -> {"cindex": per-repeat array, "penalizers"}."""
    from joblib import Parallel, delayed
    tasks = [(name, r) for name in specs for r in range(len(folds))]
    res = Parallel(n_jobs=n_jobs)(delayed(_repeat)(specs[n][0], d, tm, ev, specs[n][1], folds[r], r)
                                  for n, r in tasks)
    out = {n: {"cindex": np.empty(len(folds)), "penalizers": [None] * len(folds)} for n in specs}
    for (n, r), (c, pens) in zip(tasks, res):
        out[n]["cindex"][r] = c
        out[n]["penalizers"][r] = pens
    return out


def paired(a: np.ndarray, b: np.ndarray) -> dict:
    d = np.asarray(a) - np.asarray(b)
    return {"mean": float(d.mean()), "sd": float(d.std(ddof=1)), "min": float(d.min()),
            "max": float(d.max()), "per_repeat": d.tolist()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "stratify"))
    ap.add_argument("--endpoint", default=C.PRIMARY_ENDPOINT, choices=list(C.ENDPOINTS))
    ap.add_argument("--topk", type=int, default=C.MARKER_TOPK_STRATIFY)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    if args.dry_run:
        print(json.dumps({"endpoint": args.endpoint, "topk": args.topk, "cv": [C.N_CV_REPEATS, C.N_CV_SPLITS],
                          "penalizers": C.STRATIFY_PENALIZERS, "xgb": C.STRATIFY_XGB}, indent=2)); return
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    ev, tm = C.ENDPOINTS[args.endpoint]

    expr = pd.read_parquet(Path(args.data_dir) / "bulk_expr.parquet")
    clin = pd.read_csv(Path(args.data_dir) / "bulk_clinical.csv", index_col=0)
    sigs = json.load(open(Path(args.states_dir) / "signatures.json"))
    zexpr = zscore_genes(expr)
    d = clin.dropna(subset=C.COVARIATES + [ev, tm]).copy()
    print(f"  {len(d)} patients with age and stage, {int(d[ev].sum())} {args.endpoint} events", flush=True)

    from scfm.stages.translate import resolve_proliferation, subtype_dummies
    prolif = score_signature(zexpr, list(resolve_proliferation(zexpr.index).values())).reindex(d.index)
    d["prolif_z"] = (prolif - prolif.mean()) / prolif.std()
    feats = {}
    for m in sigs:
        F = state_features(sigs, m, zexpr, d.index, args.topk)
        feats[m] = list(F.columns)
        d = pd.concat([d, F], axis=1)
        print(f"  {m}: {F.shape[1]} state features (top-k {args.topk})", flush=True)

    rng = np.random.default_rng(C.SEED_STRATIFY)
    folds = cv_folds(d[ev], C.N_CV_SPLITS, C.N_CV_REPEATS, rng)
    cov = list(C.COVARIATES)
    specs = {"clinical": ("cox", cov), "clinical_prolif": ("cox", cov + ["prolif_z"])}
    for m in sigs:
        specs[f"ridge_{m}"] = ("ridge", cov + feats[m])
        specs[f"xgb_{m}"] = ("xgb", cov + feats[m])
    n_jobs = min(C.N_JOBS, len(specs) * len(folds))
    res = run_cv(specs, d, tm, ev, folds, n_jobs)

    ref, dummies = subtype_dummies(d.dropna(subset=["subtype"]))
    folds_ref = cv_folds(ref[ev], C.N_CV_SPLITS, C.N_CV_REPEATS, rng)
    specs_ref = {"clinical_pam50": ("cox", cov + dummies), "clinical_on_pam50_set": ("cox", cov)}
    res_ref = run_cv(specs_ref, ref, tm, ev, folds_ref, n_jobs)
    all_specs = {**specs, **specs_ref}

    def model_type(name):
        return name.split("_", 1)[0] if name.startswith(("ridge_", "xgb_")) else None

    rows, summary = [], {}
    cohort = {"n": len(d), "events": int(d[ev].sum())}
    cohort_ref = {"n": len(ref), "events": int(ref[ev].sum())}
    for name, r in list(res.items()) + list(res_ref.items()):
        on_ref = name in res_ref
        base_name = "clinical_on_pam50_set" if on_ref else "clinical"
        base = (res_ref if on_ref else res)[base_name]["cindex"]
        kind = all_specs[name][0]
        s = {"kind": kind, **(cohort_ref if on_ref else cohort),
             "n_features": len(all_specs[name][1]),
             "cindex_mean": float(r["cindex"].mean()), "cindex_sd": float(r["cindex"].std(ddof=1)),
             "cindex_per_repeat": r["cindex"].tolist(),
             "vs_clinical": (None if name == base_name
                             else {"against": base_name, **paired(r["cindex"], base)})}
        t = model_type(name)
        rep = name.split("_", 1)[1] if t else None
        if t and rep != BASELINE and f"{t}_{BASELINE}" in res:
            s["vs_baseline_same_type"] = {"against": f"{t}_{BASELINE}",
                                          **paired(r["cindex"], res[f"{t}_{BASELINE}"]["cindex"])}
        if kind == "ridge":
            s["penalizers_chosen"] = r["penalizers"]
        summary[name] = s
        for i, c in enumerate(r["cindex"]):
            rows.append({"model": name, "kind": kind, "representation": rep, "repeat": i, "cindex": c,
                         "cindex_clinical_same_folds": float(base[i]), "n": s["n"], "events": s["events"],
                         "penalizers": ";".join(map(str, r["penalizers"][i])) if r["penalizers"][i] else ""})
        print(f"  {name:<24} C = {s['cindex_mean']:.4f} (sd over repeats {s['cindex_sd']:.4f})"
              + (f"  vs {base_name} {s['vs_clinical']['mean']:+.4f}" if s["vs_clinical"] else "")
              + (f"  vs {BASELINE} {s['vs_baseline_same_type']['mean']:+.4f}"
                 if "vs_baseline_same_type" in s else ""), flush=True)

    pd.DataFrame(rows).to_csv(out / "cv_results.csv", index=False)
    n_feats = {m: len(v) for m, v in feats.items()}
    json.dump({"endpoint": args.endpoint, "topk": args.topk, "features": n_feats,
               "cv": f"{C.N_CV_REPEATS} x {C.N_CV_SPLITS}-fold, stratified by event, seed {C.SEED_STRATIFY}; "
                     "pooled out-of-fold risk, Harrell's C per repeat; all models share one set of folds",
               "cv_pam50": "clinical_pam50 and clinical_on_pam50_set: the complete cases with a PAM50 call, "
                           "their own folds (drawn next from the same generator)",
               "cindex_sd_definition": "sd over fold repeats on the same patients "
                                       "(not a confidence interval)",
               "cox": "lifelines CoxPHFitter, penalizer 0.01 for the clinical and reference models",
               "ridge": {"l1_ratio": 0.0, "penalizers": C.STRATIFY_PENALIZERS,
                         "penalised": "state features only; age and stage keep penalizer 0.01",
                         "scale": "lifelines penalizer on the mean partial log-likelihood (a ridge of "
                                  "penalizer x n on standardised covariates)",
                         "inner_cv": f"{C.STRATIFY_INNER_SPLITS}-fold, stratified by event, largest pooled "
                                     "out-of-fold C (smallest penalizer on a tie)"},
               "xgboost": {"objective": "survival:cox", **C.STRATIFY_XGB},
               "features_note": "signature scores z-scored across the cohort "
                                "(unsupervised, outside the folds)",
               "models": summary}, open(out / "cv_summary.json", "w"), indent=2)
    record_params(out, args, extra={"n": len(d), "events": int(d[ev].sum()), "features": n_feats})


if __name__ == "__main__":
    main()
