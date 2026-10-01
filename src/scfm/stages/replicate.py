"""Stage replicate — the frozen TCGA signatures in an independent cohort (METABRIC).

Two steps, kept apart on purpose:

  --freeze            write results/replicate/frozen_signatures.json from the TCGA outputs
                      (results/ladder/ladder.json, states/signatures.json, translate/scores.csv): every
                      signature with its genes, its TCGA direction and the analyses it belongs to. No
                      selection or orientation happens after this file is written.
  --real-outcomes     score the frozen signatures in METABRIC against its outcomes
  --shuffle-outcomes  the same run with (time, OS event, DSS event) permuted jointly across patients
                      (a test of the stage; writes to results/replicate_shuffled)

The run reads only the frozen file. Per signature and endpoint (OS primary, DSS secondary): gene
coverage (scored only above REPLICATE_MIN_COVERAGE); C read in the TCGA direction; a matched-random
floor built in METABRIC and read the same way; an age-adjusted Cox model stratified by cohort with the
score per SD in the TCGA direction. A signature replicates on an endpoint if its oriented C exceeds its
METABRIC floor p95 AND its HR is above 1 (the TCGA direction) with p < 0.05. Also: P3 leads within
their named subtype (oriented C, within-subtype permutation p), and FM - baseline floor-adjusted
margins with a paired patient bootstrap. METABRIC is ODbL: only aggregate statistics are written.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params
from scfm.stats import cindex_many, fw_pvalue, orient
from scfm.survival import matched_random_sets, zscore_genes

ENDPOINTS = {"OS": "os_event", "DSS": "dss_event"}
CRITERION = ("A signature replicates on an endpoint if its oriented C exceeds its METABRIC floor p95 AND "
             "its age-adjusted, cohort-stratified HR is in the TCGA direction with p < 0.05.")


# ---- data -------------------------------------------------------------------------------------
def fetch_metabric(raw: Path) -> dict[str, Path]:
    from scfm.stages.data import fetch
    return {f: fetch(C.METABRIC_URL.format(file=f), raw / f, want=sha)
            for f, sha in C.METABRIC_SHA256.items()}


def read_clinical(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype={"COHORT": "string"}, low_memory=False)


def read_expression(path: Path, patients) -> pd.DataFrame:
    """Genes × patients log2 intensities for the expression columns in `patients`: one row per Hugo
    symbol (the duplicate with the highest mean), missing values filled with that gene's mean."""
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    raw = raw.dropna(subset=["Hugo_Symbol"]).drop(columns=["Entrez_Gene_Id"], errors="ignore")
    raw = raw.set_index("Hugo_Symbol")
    keep = set(patients)
    raw = raw[[c for c in raw.columns if c in keep]].astype(float)
    raw = raw[raw.notna().any(axis=1)]
    mean = raw.mean(axis=1)
    order = np.argsort(-mean.to_numpy(), kind="stable")
    x = raw.iloc[order]
    x = x[~x.index.duplicated(keep="first")]
    X = x.to_numpy(copy=True)
    rows, cols = np.nonzero(np.isnan(X))
    X[rows, cols] = np.nanmean(X, axis=1)[rows]
    return pd.DataFrame(X, index=x.index, columns=x.columns)


def load_metabric(files: dict[str, Path]):
    """(z-scored expression genes × patients, gene means, clinical frame indexed by patient)."""
    pat = read_clinical(files["data_clinical_patient.txt"]).set_index("PATIENT_ID")
    samp = read_clinical(files["data_clinical_sample.txt"])
    if not (samp["PATIENT_ID"] == samp["SAMPLE_ID"]).all():
        raise ValueError("METABRIC: PATIENT_ID and SAMPLE_ID differ")
    expr = read_expression(files["data_mrna_illumina_microarray.txt"], pat.index)
    patients = list(expr.columns)
    p = pat.loc[patients]
    os_status = p["OS_STATUS"].astype("string")
    vital = p["VITAL_STATUS"].astype("string")
    clin = pd.DataFrame({
        "time": pd.to_numeric(p["OS_MONTHS"], errors="coerce") * C.DAYS_PER_MONTH,
        "os_event": np.where(os_status.isna(), np.nan,
                             os_status.eq("1:DECEASED").fillna(False).astype(float)),
        # DSS: deaths of other causes are censored at death
        "dss_event": np.where(vital.isna(), np.nan,
                              vital.eq("Died of Disease").fillna(False).astype(float)),
        "age": pd.to_numeric(p["AGE_AT_DIAGNOSIS"], errors="coerce"),
        "subtype": p["CLAUDIN_SUBTYPE"].astype("string"),
        "cohort": p["COHORT"].astype("string"),
    }, index=pd.Index(patients, name="patient"))
    return zscore_genes(expr), expr.mean(axis=1), clin


def shuffle_outcomes(clin: pd.DataFrame, seed: int) -> pd.DataFrame:
    """(time, OS event, DSS event) permuted jointly across patients; everything else stays."""
    q = np.random.default_rng(seed).permutation(len(clin))
    out = clin.copy()
    cols = ["time", "os_event", "dss_event"]
    out[cols] = clin[cols].to_numpy()[q]
    return out


# ---- freeze -----------------------------------------------------------------------------------
def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def freeze(ladder_json: Path, signatures_json: Path, scores_csv: Path, nulls_json: Path) -> dict:
    """The frozen replication spec from the TCGA outputs: picks (primary, sensitivity), family-wise
    maximisers, P3 leads (the picks whose within-subtype test passed or was borderline, with their
    subtype) and the proliferation reference; one entry per (signature, direction) with every
    analysis it serves."""
    lad = json.load(open(ladder_json))
    sigs = json.load(open(signatures_json))
    sc = pd.read_csv(scores_csv)
    sc["key"] = sc["resolution"].astype(str) + "/" + sc["signature"]
    direction_or = {(m, k): int(d) for m, k, d in zip(sc["model"], sc["key"], sc["direction"])}
    entries: dict[tuple, dict] = {}

    def add(model, key, direction, analysis, **extra):
        res, sid = key.split("/", 1)
        k = (model, key, direction)
        if k not in entries:
            s = sigs[model][res][sid]
            row = sc[(sc["model"] == model) & (sc["key"] == key)].iloc[0]
            entries[k] = {"id": f"{model}:{key}" + ("" if direction == 1 else ":protective"),
                          "model": model, "resolution": res, "signature": sid, "topk": s["topk"],
                          "top_cell_type": s["top_cell_type"], "genes": s["genes"], "direction": direction,
                          "direction_label": "risk" if direction == 1 else "protective",
                          "tcga_cindex": float(row["cindex"]), "analyses": []}
        entries[k]["analyses"].append(analysis)
        entries[k].update(extra)
        return entries[k]["id"]

    picks = {"primary": {}, "sensitivity": {}}
    for r in lad["ladder"]:
        picks["primary"][r["model"]] = add(r["model"], f"{r['resolution']}/{r['signature']}", 1, "primary")
    for r in (lad.get("sensitivity") or {}).get("ladder", []):
        key = f"{r['resolution']}/{r['signature']}"
        picks["sensitivity"][r["model"]] = add(r["model"], key, direction_or[(r["model"], key)],
                                               "sensitivity")
    for an, by_model in lad.get("family_wise", {}).items():
        for m, v in by_model.items():
            for key in (v or {}).get("argmax", []):
                add(m, key, 1 if an == "primary" else direction_or[(m, key)], f"family_max_{an}")
    for name, v in lad.get("p3_corrected", {}).get("picks", {}).items():
        if v["P3_pass"] or v.get("borderline"):
            an, m = name.split("/")
            e = next(e for e in entries.values() if e["id"] == picks[an][m])
            e["analyses"].append(f"p3_lead_{an}")
            e.setdefault("p3_subtypes", []).append(v["argmax_subtype"])
    prolif = json.load(open(nulls_json))["reference"]["proliferation"]
    out = list(entries.values())
    out.append({"id": "reference:proliferation", "model": "reference", "resolution": None,
                "signature": "PAM50 proliferation", "topk": None, "top_cell_type": None,
                "genes": list(prolif["genes_used"].values()), "direction": 1, "direction_label": "risk",
                "tcga_cindex": float(prolif["cindex"]), "analyses": ["reference"]})
    base = "hvg_pca"
    comparisons = [{"analysis": an, "fm": picks[an][m], "baseline": picks[an][base]}
                   for an in ("primary", "sensitivity") for m in picks[an]
                   if m != base and base in picks[an]]
    return {"cohort": "METABRIC (cBioPortal datahub brca_metabric)", "datahub_commit": C.METABRIC_COMMIT,
            "files_sha256": C.METABRIC_SHA256, "criterion": CRITERION,
            "endpoints": {"OS": "primary", "DSS": "secondary (deaths of other causes censored at death)"},
            "min_gene_coverage": C.REPLICATE_MIN_COVERAGE,
            "n_floor_sets": C.N_REPLICATE_FLOOR, "n_perm_p3": C.N_REPLICATE_PERM,
            "n_boot_margin": C.N_REPLICATE_BOOT,
            "seeds": {"floor": C.SEED_REPLICATE_FLOOR, "p3_perm": C.SEED_REPLICATE_PERM,
                      "margin_boot": C.SEED_REPLICATE_BOOT},
            "p3_lead_rule": "picks whose within-subtype permutation test passed in TCGA (P3_corrected) or "
                            "was borderline (p within 2 Monte-Carlo SE of 0.05), tested within the subtype "
                            "where TCGA's maximum sat",
            "frozen_from": {"ladder.json": _sha(ladder_json),
                            "signatures.json": _sha(signatures_json), "scores.csv": _sha(scores_csv),
                            "nulls.json": _sha(nulls_json)},
            "signatures": out, "comparisons": comparisons}


# ---- analysis ---------------------------------------------------------------------------------
def usable(clin: pd.DataFrame, evcol: str) -> np.ndarray:
    """Patients with a positive follow-up time and a known event status for this endpoint."""
    return (clin[["time", evcol]].notna().all(axis=1) & (clin["time"] > 0)).to_numpy()


def _perm_scores(s: np.ndarray, n: int, rng) -> np.ndarray:
    """Rows: the score re-assigned under n outcome permutations (scoring patient i against outcome
    q[i] equals scoring outcome j against s[argsort(q)][j])."""
    return np.vstack([s[np.argsort(rng.permutation(len(s)))] for _ in range(n)])


def cox_stratified(df: pd.DataFrame) -> dict:
    """Age-adjusted Cox of the score (per SD, TCGA direction), stratified by cohort; unpenalised."""
    from lifelines import CoxPHFitter
    d = df[["time", "event", "score_dir_z", "age", "cohort"]].dropna()
    cph = CoxPHFitter(penalizer=0.0).fit(d, duration_col="time", event_col="event", strata=["cohort"])
    r = cph.summary.loc["score_dir_z"]
    return {"n_cox": len(d), "events_cox": int(d["event"].sum()), "hr_per_sd": float(np.exp(r["coef"])),
            "hr_ci": [float(np.exp(r["coef lower 95%"])), float(np.exp(r["coef upper 95%"]))],
            "p": float(r["p"])}


def analyse(spec: dict, zexpr: pd.DataFrame, gene_mean: pd.Series, clin: pd.DataFrame) -> dict:
    rng_floor = np.random.default_rng(C.SEED_REPLICATE_FLOOR)
    rng_perm = np.random.default_rng(C.SEED_REPLICATE_PERM)
    gidx = {g: i for i, g in enumerate(zexpr.index)}
    Z = zexpr.to_numpy(float)
    score_of = {}
    per_sig = []
    for sig in spec["signatures"]:
        present = [g for g in sig["genes"] if g in gidx]
        cov = len(present) / len(sig["genes"])
        rec = {"id": sig["id"], "model": sig["model"], "analyses": sig["analyses"],
               "n_genes": len(sig["genes"]),
               "n_genes_measured": len(present), "coverage": cov,
               "genes_missing": [g for g in sig["genes"] if g not in gidx],
               "scored": bool(cov >= spec["min_gene_coverage"] and len(present) >= 3)}
        if rec["scored"]:
            s = Z[[gidx[g] for g in present]].mean(axis=0)
            score_of[sig["id"]] = s
            sets = matched_random_sets(present, gene_mean, spec["n_floor_sets"], C.EXPRESSION_BINS, rng_floor)
            F = np.vstack([Z[[gidx[g] for g in gs]].mean(axis=0) for gs in sets])
            rec["endpoints"] = {}
            for ep, evcol in ENDPOINTS.items():
                ok = usable(clin, evcol)
                T, E = clin["time"].to_numpy(float)[ok], clin[evcol].to_numpy(float)[ok]
                c = float(cindex_many(s[ok][None, :], T, E)[0])
                fl = orient(cindex_many(F[:, ok], T, E), sig["direction"])
                c_or = orient(c, sig["direction"])
                sz = sig["direction"] * (s - s[ok].mean()) / s[ok].std(ddof=1)
                df = clin.assign(event=clin[evcol], score_dir_z=sz)[ok]
                cx = cox_stratified(df)
                r = {"n": int(ok.sum()), "events": int(E.sum()), "cindex": c, "cindex_oriented": c_or,
                     "floor_mean": float(np.mean(fl)), "floor_p95": float(np.percentile(fl, 95)),
                     "above_floor": c_or - float(np.mean(fl)), **cx}
                r["replicates"] = bool(c_or > r["floor_p95"] and cx["hr_per_sd"] > 1 and cx["p"] < C.ALPHA)
                rec["endpoints"][ep] = r
            for st in sig.get("p3_subtypes", []):
                rec.setdefault("p3", {})[st] = {}
                for ep, evcol in ENDPOINTS.items():
                    ok = usable(clin, evcol) & clin["subtype"].eq(st).fillna(False).to_numpy()
                    T, E = clin["time"].to_numpy(float)[ok], clin[evcol].to_numpy(float)[ok]
                    c_or = orient(float(cindex_many(s[ok][None, :], T, E)[0]), sig["direction"])
                    null = orient(cindex_many(_perm_scores(s[ok], spec["n_perm_p3"], rng_perm), T, E),
                                  sig["direction"])
                    rec["p3"][st][ep] = {"n": int(ok.sum()), "events": int(E.sum()), "cindex_oriented": c_or,
                                         "perm_p": fw_pvalue(c_or, null), "n_perm": spec["n_perm_p3"]}
        per_sig.append(rec)

    by_id = {r["id"]: r for r in per_sig}
    dir_of = {s["id"]: s["direction"] for s in spec["signatures"]}
    comps = []
    need = sorted({x for c in spec["comparisons"] for x in (c["fm"], c["baseline"]) if x in score_of})
    for ep, evcol in ENDPOINTS.items():
        ok = usable(clin, evcol)
        T, E = clin["time"].to_numpy(float)[ok], clin[evcol].to_numpy(float)[ok]
        rng_boot = np.random.default_rng([C.SEED_REPLICATE_BOOT, list(ENDPOINTS).index(ep)])
        idx = [rng_boot.integers(0, ok.sum(), ok.sum()) for _ in range(spec["n_boot_margin"])]
        S = np.vstack([score_of[k][ok] for k in need]) if need else np.empty((0, ok.sum()))
        Cb = np.vstack([cindex_many(S[:, ix], T[ix], E[ix]) for ix in idx]) if need else None
        for c in spec["comparisons"]:
            row = {"analysis": c["analysis"], "endpoint": ep, "fm": c["fm"], "baseline": c["baseline"]}
            if c["fm"] in score_of and c["baseline"] in score_of:
                af = {k: by_id[k]["endpoints"][ep]["above_floor"] for k in (c["fm"], c["baseline"])}
                fl = {k: by_id[k]["endpoints"][ep]["floor_mean"] for k in (c["fm"], c["baseline"])}
                b = {k: orient(Cb[:, need.index(k)], dir_of[k]) - fl[k] for k in (c["fm"], c["baseline"])}
                m = b[c["fm"]] - b[c["baseline"]]
                lo, hi = np.percentile(m, [2.5, 97.5])
                row.update(margin=af[c["fm"]] - af[c["baseline"]], ci=[float(lo), float(hi)],
                           n_boot=len(m), floors="held at their METABRIC full-cohort values")
            else:
                row["margin"] = None
            comps.append(row)
    return {"signatures": per_sig, "comparisons": comps}


def summary_lines(res: dict) -> list[str]:
    lines = [f"METABRIC replication — outcomes: {res['outcomes']}", res["criterion"], "",
             f"{'signature':<40} {'cov':>5} {'endpoint':<4} {'C_or':>6} {'floor95':>7} {'HR/SD':>6} "
             f"{'p':>8}  replicates"]
    for r in res["signatures"]:
        if not r["scored"]:
            lines.append(f"{r['id']:<40} {r['coverage']:>5.2f}  not scored (coverage)")
            continue
        for ep, e in r["endpoints"].items():
            lines.append(f"{r['id']:<40} {r['coverage']:>5.2f} {ep:<4} {e['cindex_oriented']:>6.3f} "
                         f"{e['floor_p95']:>7.3f} {e['hr_per_sd']:>6.2f} {e['p']:>8.2g}  {e['replicates']}")
        for st, eps in r.get("p3", {}).items():
            for ep, e in eps.items():
                lines.append(f"    P3 within {st} {ep}: C_or {e['cindex_oriented']:.3f}, "
                             f"perm p {e['perm_p']:.3g}")
    lines.append("")
    for c in res["comparisons"]:
        if c["margin"] is not None:
            lines.append(f"  {c['analysis']:<11} {c['endpoint']:<4} {c['fm']} - {c['baseline']}: "
                         f"{c['margin']:+.3f} [{c['ci'][0]:+.3f},{c['ci'][1]:+.3f}]")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--freeze", action="store_true", help="write the frozen spec from the TCGA outputs")
    mode.add_argument("--real-outcomes", action="store_true", help="the replication analysis")
    mode.add_argument("--shuffle-outcomes", action="store_true", help="test run on permuted outcomes")
    mode.add_argument("--dry-run", action="store_true")
    ap.add_argument("--frozen", default=str(C.RESULTS / "replicate" / "frozen_signatures.json"))
    ap.add_argument("--raw-dir", default=str(C.DATA / "raw" / "metabric"))
    ap.add_argument("--out-dir", default=None,
                    help="default results/replicate, or results/replicate_shuffled with --shuffle-outcomes")
    ap.add_argument("--force", action="store_true", help="--freeze: overwrite an existing frozen spec")
    args = ap.parse_args(argv)
    frozen = Path(args.frozen)

    if args.dry_run:
        print(json.dumps({"frozen": str(frozen), "files": list(C.METABRIC_SHA256), "criterion": CRITERION},
                         indent=2)); return
    if args.freeze:
        if frozen.exists() and not args.force:
            raise SystemExit(f"{frozen} exists; the spec is frozen (pass --force to rewrite it)")
        spec = freeze(C.RESULTS / "ladder" / "ladder.json", C.RESULTS / "states" / "signatures.json",
                      C.RESULTS / "translate" / "scores.csv", C.RESULTS / "translate" / "nulls.json")
        frozen.parent.mkdir(parents=True, exist_ok=True)
        frozen.write_text(json.dumps(spec, indent=2) + "\n")
        print(f"  {len(spec['signatures'])} signatures frozen to {frozen}")
        return

    spec = json.load(open(frozen))
    out = Path(args.out_dir or C.RESULTS / ("replicate_shuffled" if args.shuffle_outcomes else "replicate"))
    out.mkdir(parents=True, exist_ok=True)
    files = fetch_metabric(Path(args.raw_dir))
    zexpr, gene_mean, clin = load_metabric(files)
    if args.shuffle_outcomes:
        clin = shuffle_outcomes(clin, C.SEED_REPLICATE_SHUFFLE)
    print(f"  {zexpr.shape[1]} patients with expression, {zexpr.shape[0]} genes", flush=True)
    res = analyse(spec, zexpr, gene_mean, clin)
    res = {"outcomes": (f"SHUFFLED (seed {C.SEED_REPLICATE_SHUFFLE}; not a result)" if args.shuffle_outcomes
                        else "real"),
           "criterion": CRITERION, "frozen_sha256": _sha(frozen),
           "n_patients_expression": int(zexpr.shape[1]), "n_genes": int(zexpr.shape[0]),
           "endpoints": {ep: {"n": int(usable(clin, c).sum()),
                              "events": int(clin.loc[usable(clin, c), c].sum())}
                         for ep, c in ENDPOINTS.items()}, **res}
    json.dump(res, open(out / "metabric.json", "w"), indent=2)
    lines = summary_lines(res)
    (out / "metabric_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    record_params(out, args, extra={"outcomes": res["outcomes"]})


if __name__ == "__main__":
    main()
