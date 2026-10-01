"""Stage 04 — reverse-translate every signature into the bulk cohort and score it.

For each signature: score every TCGA patient (mean z of the signature genes), fit an age- and
stage-adjusted Cox model, and record the concordance index of the score. Then the two controls
that make those numbers readable:

  null   the family-wise permutation null — outcomes permuted across patients, and for each
         permutation the BEST C-index across all of a representation's signatures is kept, so
         "best of k" is corrected at the level it was selected at
  floor  matched random gene sets — same size, same mean-expression profile — so a signature is
         credited only for what its *particular* genes add over any genes of that shape

Plus the reference rung, on one complete-case set (age, stage, PAM50 call, outcome): Cox models of
clinical covariates, PAM50 subtype, and both, each in sample and cross-validated; age + stage on
every patient with both (the base model for the added-value analysis); and a published
signature, the PAM50 11-gene proliferation score, scored like any signature with its own floor and
permutation p. The pre-registered rung (in-sample Cox of subtype + age + stage) is kept as
`pam50_cindex`.

Sensitivity analysis (added after the primary run; the primary ladder is unchanged): the
pre-registered C-index credits a signature only when a higher score means worse survival.
The *_or columns read each signature in the direction it acts, with a floor oriented the same
way and a family-wise null over max(C, 1 - C) built from the same permutations.
Output: results/translate/scores.csv, nulls.json (with the per-permutation family maxima),
params.json, and null_matrix_<model>.npz (permutation × signature C, not tracked). A non-primary
endpoint writes to results/translate_<endpoint>.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params
from scfm.stats import cindex_many, cv_cox_cindex, cv_folds, fw_pvalue, insample_cindex, mc_se
from scfm.survival import cindex, fit_cox, matched_random_sets, score_signature, zscore_genes


def load_bulk(data_dir: Path):
    expr = pd.read_parquet(data_dir / "bulk_expr.parquet")           # genes × patients
    clin = pd.read_csv(data_dir / "bulk_clinical.csv", index_col=0)
    return expr, clin


def score_one(zexpr, clin, genes, ev, tm):
    df = clin.copy()
    df["score"] = score_signature(zexpr, genes)
    ci = cindex(df, tm, ev, "score")
    try:
        _, row = fit_cox(df, tm, ev, ["score"] + C.COVARIATES)
        hr, p = float(np.exp(row["coef"])), float(row["p"])
    except Exception:
        hr, p = float("nan"), float("nan")
    return ci, hr, p


def _scores(Z, ix):
    return Z[ix].mean(axis=0) if len(ix) >= 3 else None


def _ci(t, e, s):
    from lifelines.utils import concordance_index
    return float(concordance_index(t, -s, e)) if e.sum() >= 5 else float("nan")


def _floor_chunk(Z, T, E, sigs_ix):
    """Per signature: C-indices of its matched random sets."""
    out = []
    for sets in sigs_ix:
        vals = [_ci(T, E, sc) for ix in sets if (sc := _scores(Z, ix)) is not None]
        out.append([v for v in vals if np.isfinite(v)])
    return out


def _null_chunk(S, T, E, seeds):
    """Per permutation seed, the C-index of every row of S against permuted outcomes; the family
    maxima (best C; best max(C, 1 - C)) are read from these rows."""
    out = []
    for sd in seeds:
        p = np.random.default_rng(int(sd)).permutation(len(T))
        out.append(np.array([_ci(T[p], E[p], S[j]) for j in range(len(S))]))
    return out


def _oriented(ci0, floor):
    """Sensitivity analysis: read the signature in the direction it acts on the full cohort.
    A protective signature (C < 0.5) is scored as 1 - C, and so is each of its floor sets,
    so the floor asks whether the signature beats random genes in its own direction."""
    if not np.isfinite(ci0):
        return {"direction": np.nan, "cindex_or": np.nan, "floor_or_mean": np.nan,
                "floor_or_p95": np.nan, "above_floor_or": np.nan}
    d = 1 if ci0 >= 0.5 else -1
    f = [v if d == 1 else 1 - v for v in floor]
    c = ci0 if d == 1 else 1 - ci0
    return {"direction": d, "cindex_or": c,
            "floor_or_mean": float(np.mean(f)) if f else np.nan,
            "floor_or_p95": float(np.percentile(f, 95)) if f else np.nan,
            "above_floor_or": c - float(np.mean(f)) if f else np.nan}


def endpoint_dir(stage: str, endpoint: str) -> str:
    """Output directory name: the primary endpoint keeps the plain name, others get a suffix."""
    return stage if endpoint == C.PRIMARY_ENDPOINT else f"{stage}_{endpoint.lower()}"


def subtype_dummies(d: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    dm = pd.get_dummies(d["subtype"], prefix="subtype", drop_first=True, dtype=float)
    return pd.concat([d, dm], axis=1), list(dm.columns)


REF_PENALIZER = 0.01
PENALTY_NOTE = ("lifelines penalizer 0.01 on the mean partial log-likelihood, i.e. a ridge of 0.01 x n "
                "on standardised covariates; cindex_in_sample_mle is the unpenalised fit. With few "
                "groups of near-equal hazard (PAM50 alone) the two can differ: compare the CV C, not the "
                "in-sample C, with single-signature C-indices")
CV_SD_NOTE = "sd over fold repeats on the same patients (not a confidence interval)"
DEFINITIONS = {
    "scores.hr, scores.p": "per unit of the mean-z score, Cox with age and stage, lifelines penalizer "
                           "0.01 (ridge), Wald p",
    "reference.proliferation.hr_adjusted, p_adjusted": "same definition as scores.hr and scores.p",
    "reference.*.cindex_cv_sd": CV_SD_NOTE,
    "*_mc_se": "Monte-Carlo standard error sqrt(p(1-p)/n) of a permutation p from n permutations",
}


def clinical_references(clin, ev, tm, endpoint) -> dict:
    """Named Cox references: in-sample C and cross-validated C (repeated event-stratified K-fold,
    pooled out-of-fold linear predictor, Harrell's C per repeat)."""
    sref = clin.dropna(subset=["age", "stage_ord", "subtype", ev, tm])
    sref, dummies = subtype_dummies(sref)
    sets = {"clinical": C.COVARIATES, "pam50": dummies, "pam50_clinical": dummies + C.COVARIATES}
    rng = np.random.default_rng(C.SEED_REF_CV)
    folds = cv_folds(sref[ev], C.N_CV_SPLITS, C.N_CV_REPEATS, rng)
    cv = cv_cox_cindex(sref, tm, ev, sets, folds)
    full = clin.dropna(subset=C.COVARIATES + [ev, tm])
    cv_full = cv_cox_cindex(full, tm, ev, {"clinical_full": C.COVARIATES},
                            cv_folds(full[ev], C.N_CV_SPLITS, C.N_CV_REPEATS, rng))
    out = {}
    for name, cols, d, cvs in [(k, v, sref, cv[k]) for k, v in sets.items()] + [
            ("clinical_full", C.COVARIATES, full, cv_full["clinical_full"])]:
        out[name] = {"covariates": cols, "n": len(d), "events": int(d[ev].sum()),
                     "cindex_in_sample": insample_cindex(d, tm, ev, cols),
                     "cindex_in_sample_mle": insample_cindex(d, tm, ev, cols, penalizer=0.0),
                     "penalizer": REF_PENALIZER, "penalty": PENALTY_NOTE,
                     "cindex_cv_mean": float(cvs.mean()), "cindex_cv_sd": float(cvs.std(ddof=1)),
                     "cindex_cv_sd_definition": CV_SD_NOTE,
                     "cindex_cv_repeats": cvs.tolist(),
                     "cv": f"{C.N_CV_REPEATS} x {C.N_CV_SPLITS}-fold, stratified by event, "
                           f"seed {C.SEED_REF_CV}"}
        print(f"  reference {name}: n={len(d)}, in-sample C = {out[name]['cindex_in_sample']:.3f} "
              f"(unpenalised {out[name]['cindex_in_sample_mle']:.3f}), "
              f"CV C = {out[name]['cindex_cv_mean']:.3f} ± {out[name]['cindex_cv_sd']:.3f}", flush=True)
    return out


def resolve_proliferation(symbols) -> dict[str, str]:
    """Published symbol -> symbol present in bulk (the published one, else its first present alias)."""
    have = set(symbols)
    out = {}
    for g, aliases in C.PROLIFERATION_GENES.items():
        hit = next((a for a in [g] + aliases if a in have), None)
        if hit is not None:
            out[g] = hit
    return out


def proliferation_reference(clin, zexpr, gene_mean, T, E, ok, ev, tm, n_floor, n_perm) -> dict:
    used = resolve_proliferation(zexpr.index)
    genes = list(used.values())
    ci0, hr, p = score_one(zexpr, clin, genes, ev, tm)
    gidx = {g: i for i, g in enumerate(zexpr.index)}
    Zok = zexpr.to_numpy(dtype=np.float64)[:, ok]
    rng_f = np.random.default_rng(C.SEED_REF_FLOOR)
    sets = matched_random_sets(genes, gene_mean, n_floor, C.EXPRESSION_BINS, rng_f)
    floor = [v for v in _floor_chunk(Zok, T, E, [[[gidx[g] for g in gs] for gs in sets]])[0]]
    sc = _scores(Zok, [gidx[g] for g in genes])
    rng_p = np.random.default_rng(C.SEED_REF_PERM)
    # scoring patient i against outcome q[i] is scoring outcome j against sc[argsort(q)][j]
    perm_scores = np.vstack([sc[np.argsort(rng_p.permutation(len(T)))] for _ in range(n_perm)])
    null = cindex_many(perm_scores, T, E)
    o = _oriented(ci0, floor)
    out = {"name": "PAM50 11-gene proliferation score (Nielsen et al. 2010; ROR-P)",
           "genes_published": list(C.PROLIFERATION_GENES), "genes_used": used,
           "n_genes": len(genes), "n": int(ok.sum()), "events": int(E.sum()),
           "cindex": ci0, "hr_adjusted": hr, "p_adjusted": p,
           "hr_adjusted_definition": DEFINITIONS["scores.hr, scores.p"],
           "floor_mean": float(np.mean(floor)), "floor_p95": float(np.percentile(floor, 95)),
           "above_floor": ci0 - float(np.mean(floor)), "n_floor": n_floor,
           "perm_p": fw_pvalue(ci0, null), "n_perm": n_perm,
           "perm_p_mc_se": mc_se(fw_pvalue(ci0, null), n_perm),
           "seeds": {"floor": C.SEED_REF_FLOOR, "perm": C.SEED_REF_PERM}, **o}
    print(f"  reference proliferation: C = {ci0:.3f}, floor {out['floor_mean']:.3f}, "
          f"perm p = {out['perm_p']:.3g}", flush=True)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--out-dir", default=None,
                    help="default results/translate, or results/translate_<endpoint> off the primary one")
    ap.add_argument("--endpoint", default=C.PRIMARY_ENDPOINT, choices=list(C.ENDPOINTS))
    ap.add_argument("--n-perm", type=int, default=C.N_PERMUTATIONS)
    ap.add_argument("--n-floor", type=int, default=C.N_FLOOR_SETS)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    if args.out_dir is None:
        args.out_dir = str(C.RESULTS / endpoint_dir("translate", args.endpoint))
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    ev, tm = C.ENDPOINTS[args.endpoint]

    if args.dry_run:
        print(json.dumps({"endpoint": args.endpoint, "covariates": C.COVARIATES,
                          "n_perm": args.n_perm, "n_floor": args.n_floor}, indent=2)); return

    rng = np.random.default_rng(args.seed)
    expr, clin = load_bulk(Path(args.data_dir))
    zexpr = zscore_genes(expr)
    gene_mean = expr.mean(axis=1)
    sigs = json.load(open(Path(args.states_dir) / "signatures.json"))
    n_pat, n_ev = len(clin), int(clin[ev].sum())
    print(f"  {n_pat:,} patients, {n_ev} {args.endpoint} events", flush=True)

    # everything below works on numpy arrays: one score matrix per model, and the two control
    # loops (floor sets, outcome permutations) fanned out over processes; each worker gets one
    # copy of the arrays per chunk, not per task
    from joblib import Parallel, delayed
    ok = clin[[tm, ev]].notna().all(axis=1).to_numpy()
    Zok = np.ascontiguousarray(zexpr.to_numpy(dtype=np.float64)[:, ok])
    gidx = {g: i for i, g in enumerate(zexpr.index)}
    T, E = clin[tm].to_numpy()[ok], clin[ev].to_numpy()[ok]
    scores_of = lambda genes: _scores(Zok, [gidx[g] for g in genes if g in gidx])

    par = Parallel(n_jobs=C.N_JOBS, max_nbytes=None)   # loky's memmap cleanup logs spurious KeyErrors
    rows, nulls = [], {}
    for model, by_res in sigs.items():
        flat = [(res, sid, s) for res, d in by_res.items() for sid, s in d.items()]
        print(f"  {model}: {len(flat)} signatures", flush=True)
        floor_sets = [matched_random_sets(s["genes"], gene_mean, args.n_floor,
                                          C.EXPRESSION_BINS, rng) for _, _, s in flat]
        ix_sets = [[[gidx[g] for g in gs if g in gidx] for gs in fs] for fs in floor_sets]
        parts = [p for p in np.array_split(np.arange(len(ix_sets)), C.N_JOBS) if len(p)]
        floors = [f for chunk in par(delayed(_floor_chunk)(Zok, T, E, [ix_sets[i] for i in part])
                                     for part in parts)
                  for f in chunk]
        scored, keys = [], []
        for (res, sid, s), floor in zip(flat, floors):
            ci0, hr, p = score_one(zexpr, clin, s["genes"], ev, tm)
            sc = scores_of(s["genes"])
            if sc is not None:
                scored.append(sc)
                keys.append(f"{res}/{sid}")
            rows.append({"model": model, "kind": C.MODELS[model]["kind"], "resolution": res,
                         "signature": sid, "topk": s["topk"], "n_cells": s["n_cells"],
                         "top_cell_type": s["top_cell_type"], "hvg_frac": s["hvg_frac"],
                         "n_genes_in_bulk": sum(g in gidx for g in s["genes"]),
                         "cindex": ci0, "hr": hr, "p": p,
                         "floor_mean": float(np.mean(floor)) if floor else np.nan,
                         "floor_p95": float(np.percentile(floor, 95)) if floor else np.nan,
                         "above_floor": ci0 - float(np.mean(floor)) if floor else np.nan,
                         **_oriented(ci0, floor)})
        # family-wise null: best C-index across this model's signatures per permutation
        S = np.vstack(scored)
        seeds = rng.integers(1 << 62, size=args.n_perm)
        cmat = np.vstack([v for chunk in par(delayed(_null_chunk)(S, T, E, c)
                                             for c in np.array_split(seeds, C.N_JOBS)) for v in chunk])
        best = [float(c.max()) for c in cmat]
        best_or = [float(np.maximum(c, 1 - c).max()) for c in cmat]
        np.savez_compressed(out / f"null_matrix_{model}.npz", cindex=cmat, signatures=np.array(keys))
        nulls[model] = {"best_cindex_null_mean": float(np.mean(best)),
                        "best_cindex_null_p95": float(np.percentile(best, 95)),
                        "best_oriented_null_mean": float(np.mean(best_or)),
                        "best_oriented_null_p95": float(np.percentile(best_or, 95)),
                        "n_perm": args.n_perm,
                        "best_cindex_null": best, "best_oriented_null": best_or}
        print(f"  {model}: family-wise null p95 = {nulls[model]['best_cindex_null_p95']:.3f}"
              f"  (either direction: {nulls[model]['best_oriented_null_p95']:.3f})", flush=True)

    # reference rung: PAM50 subtype as a categorical predictor
    ref = {}
    if clin["subtype"].notna().sum() > 50:
        d = pd.get_dummies(clin[["subtype", ev, tm] + C.COVARIATES].dropna(), columns=["subtype"],
                           drop_first=True, dtype=float)
        from lifelines import CoxPHFitter
        cph = CoxPHFitter(penalizer=0.01).fit(d, duration_col=tm, event_col=ev)
        d["score"] = cph.predict_partial_hazard(d)
        ref["pam50_cindex"] = cindex(d, tm, ev, "score")
        print(f"  reference as coded (PAM50 + age + stage, in sample): C = {ref['pam50_cindex']:.3f}",
              flush=True)
        ref.update(clinical_references(clin, ev, tm, args.endpoint))
    ref["proliferation"] = proliferation_reference(clin, zexpr, gene_mean, T, E, ok, ev, tm,
                                                   args.n_floor, C.N_REF_PERM)

    pd.DataFrame(rows).to_csv(out / "scores.csv", index=False)
    json.dump({"nulls": nulls, "reference": ref, "definitions": DEFINITIONS, "endpoint": args.endpoint,
               "n_patients": n_pat, "n_events": n_ev}, open(out / "nulls.json", "w"), indent=2)
    record_params(out, args, extra={"n_signatures": len(rows), "n_events": n_ev})


if __name__ == "__main__":
    main()
