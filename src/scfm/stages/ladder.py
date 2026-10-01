"""Stage 05 — assemble the ladder and test the four predictions.

Takes stage 04's per-signature scores and produces the one table the study is about: for each
representation, its best signature's concordance read against null, floor, baseline and
reference — with bootstrap intervals over PATIENTS, because the patient is the unit of
independence here and a bootstrap over genes or cells would be pseudo-replication.

P1  FM above floor            P2  FM − baseline margin excludes zero
P3  FM prognostic within a PAM50 subtype (the reference does not already encode it)
P4  FM signatures passing P2 are enriched outside the HVG set
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params
from scfm.survival import bootstrap_patients, cindex, score_signature, zscore_genes


def best_per_model(scores: pd.DataFrame) -> pd.DataFrame:
    """Best signature per model by above-floor margin, not raw C-index."""
    return (scores.dropna(subset=["above_floor"])
            .sort_values("above_floor", ascending=False)
            .groupby("model", as_index=False).head(1).set_index("model"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--translate-dir", default=str(C.RESULTS / "translate"))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "ladder"))
    ap.add_argument("--n-boot", type=int, default=C.N_BOOTSTRAP)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    if args.dry_run:
        print(json.dumps({"n_boot": args.n_boot, "unit": "patient"}, indent=2)); return

    rng = np.random.default_rng(args.seed)
    scores = pd.read_csv(Path(args.translate_dir) / "scores.csv")
    meta = json.load(open(Path(args.translate_dir) / "nulls.json"))
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
    json.dump({"ladder": ladder.reset_index().to_dict("records"), "reference_pam50": ref,
               "predictions": P, "endpoint": meta["endpoint"],
               "sensitivity": None if sens_ladder is None else {
                   "note": "each signature read in the direction it acts; added after the primary run",
                   "ladder": sens_ladder.reset_index().to_dict("records"), "predictions": P_or},
               "n_patients": meta["n_patients"], "n_events": meta["n_events"]},
              open(out / "ladder.json", "w"), indent=2, default=float)

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
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    record_params(out, args, extra={"predictions": P})


if __name__ == "__main__":
    main()
