"""Stage 04 — reverse-translate every signature into the bulk cohort and score it.

For each signature: score every TCGA patient (mean z of the signature genes), fit an age- and
stage-adjusted Cox model, and record the concordance index of the score. Then the two controls
that make those numbers readable:

  null   the family-wise permutation null — outcomes permuted across patients, and for each
         permutation the BEST C-index across all of a representation's signatures is kept, so
         "best of k" is corrected at the level it was selected at
  floor  matched random gene sets — same size, same mean-expression profile — so a signature is
         credited only for what its *particular* genes add over any genes of that shape

Plus the reference rung: the PAM50 subtype call's own C-index in the same cohort.
Output: results/translate/scores.csv, nulls.json, params.json.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import config as C
import qc_common as qc
from surv_common import cindex, fit_cox, matched_random_sets, score_signature, zscore_genes


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
    """Per permutation seed: the best C-index over all rows of S against permuted outcomes."""
    best = []
    for sd in seeds:
        p = np.random.default_rng(int(sd)).permutation(len(T))
        best.append(max(_ci(T[p], E[p], S[j]) for j in range(len(S))))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "translate"))
    ap.add_argument("--endpoint", default=C.PRIMARY_ENDPOINT, choices=list(C.ENDPOINTS))
    ap.add_argument("--n-perm", type=int, default=C.N_PERMUTATIONS)
    ap.add_argument("--n-floor", type=int, default=C.N_FLOOR_SETS)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
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
        scored = []
        for (res, sid, s), floor in zip(flat, floors):
            ci0, hr, p = score_one(zexpr, clin, s["genes"], ev, tm)
            sc = scores_of(s["genes"])
            if sc is not None:
                scored.append(sc)
            rows.append({"model": model, "kind": C.MODELS[model]["kind"], "resolution": res,
                         "signature": sid, "topk": s["topk"], "n_cells": s["n_cells"],
                         "top_cell_type": s["top_cell_type"], "hvg_frac": s["hvg_frac"],
                         "n_genes_in_bulk": sum(g in gidx for g in s["genes"]),
                         "cindex": ci0, "hr": hr, "p": p,
                         "floor_mean": float(np.mean(floor)) if floor else np.nan,
                         "floor_p95": float(np.percentile(floor, 95)) if floor else np.nan,
                         "above_floor": ci0 - float(np.mean(floor)) if floor else np.nan})
        # family-wise null: best C-index across this model's signatures per permutation
        S = np.vstack(scored)
        seeds = rng.integers(1 << 62, size=args.n_perm)
        best = [v for chunk in par(delayed(_null_chunk)(S, T, E, c)
                                   for c in np.array_split(seeds, C.N_JOBS)) for v in chunk]
        nulls[model] = {"best_cindex_null_mean": float(np.mean(best)),
                        "best_cindex_null_p95": float(np.percentile(best, 95)),
                        "n_perm": args.n_perm}
        print(f"  {model}: family-wise null p95 = {nulls[model]['best_cindex_null_p95']:.3f}",
              flush=True)

    # reference rung: PAM50 subtype as a categorical predictor
    ref = {}
    if clin["subtype"].notna().sum() > 50:
        d = pd.get_dummies(clin[["subtype", ev, tm] + C.COVARIATES].dropna(), columns=["subtype"],
                           drop_first=True, dtype=float)
        from lifelines import CoxPHFitter
        cph = CoxPHFitter(penalizer=0.01).fit(d, duration_col=tm, event_col=ev)
        d["score"] = cph.predict_partial_hazard(d)
        ref["pam50_cindex"] = cindex(d, tm, ev, "score")
        print(f"  reference (PAM50): C = {ref['pam50_cindex']:.3f}", flush=True)

    pd.DataFrame(rows).to_csv(out / "scores.csv", index=False)
    json.dump({"nulls": nulls, "reference": ref, "endpoint": args.endpoint,
               "n_patients": n_pat, "n_events": n_ev}, open(out / "nulls.json", "w"), indent=2)
    qc.record_params(out, args, extra={"n_signatures": len(rows), "n_events": n_ev})


if __name__ == "__main__":
    main()
