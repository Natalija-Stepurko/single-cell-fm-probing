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
from surv_common import (cindex, fit_cox, matched_random_sets, permute_outcomes,
                         score_signature, zscore_genes)


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

    rows, nulls = [], {}
    for model, by_res in sigs.items():
        flat = [(res, sid, s) for res, d in by_res.items() for sid, s in d.items()]
        print(f"  {model}: {len(flat)} signatures", flush=True)
        scored = {}
        for res, sid, s in flat:
            ci, hr, p = score_one(zexpr, clin, s["genes"], ev, tm)
            floor = [cindex(clin.assign(score=score_signature(zexpr, g)), tm, ev, "score")
                     for g in matched_random_sets(s["genes"], gene_mean, args.n_floor,
                                                  C.EXPRESSION_BINS, rng)]
            floor = [f for f in floor if np.isfinite(f)]
            scored[sid] = s["genes"]
            rows.append({"model": model, "kind": C.MODELS[model]["kind"], "resolution": res,
                         "signature": sid, "topk": s["topk"], "n_cells": s["n_cells"],
                         "top_cell_type": s["top_cell_type"], "hvg_frac": s["hvg_frac"],
                         "n_genes_in_bulk": sum(g in zexpr.index for g in s["genes"]),
                         "cindex": ci, "hr": hr, "p": p,
                         "floor_mean": float(np.mean(floor)) if floor else np.nan,
                         "floor_p95": float(np.percentile(floor, 95)) if floor else np.nan,
                         "above_floor": ci - float(np.mean(floor)) if floor else np.nan})
        # family-wise null: best C-index across this model's signatures per permutation
        best = []
        for _ in range(args.n_perm):
            pc = permute_outcomes(clin, tm, ev, rng)
            best.append(max(cindex(pc.assign(score=score_signature(zexpr, g)), tm, ev, "score")
                            for g in scored.values()))
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
