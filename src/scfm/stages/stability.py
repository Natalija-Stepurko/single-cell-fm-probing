"""Stage stability — do the cell states survive dropping donors?

For each representation and resolution, N_REPEATS times: keep a random KEEP_DONOR_FRAC of the
donors (one draw per repeat, shared by all representations), rebuild the 15-NN graph on the kept
cells' embedding and recluster with the Leiden settings of the states stage; compare with the full
clustering restricted to the kept cells (adjusted Rand index).

For every pick and family maximum of the primary-endpoint ladder: the subsample cluster that best
matches the state (maximum Jaccard index on cells), and the Jaccard index of that cluster's top-k
markers with the signature's genes. Markers are recomputed as in the states stage: Wilcoxon on the
log-normalised expression, the same gene universe and the same k; "the rest" is the kept cells in
subsample clusters of at least MIN_CELLS_PER_STATE cells. A state with no kept cells (all its donors
dropped) scores 0 on both for that repeat. Pick summaries are given over all repeats and over the
repeats that kept the state's top donor.

Reads results/ladder/ladder.json for the picks, so it runs after ladder; it is not part of 'all'
(about 30 min on 6 cores).

Outputs (results/states/): stability.csv (per representation x resolution), stability_picks.csv
(per pick), stability_repeats.csv (per pick and repeat), stability_params.json.
"""
import argparse
import json
import os
import time
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params

SEED_STABILITY = C.SEED + 20      # donor subsets; Leiden and the graph keep the states seed
N_REPEATS = 3 if C.SMOKE else 20
KEEP_DONOR_FRAC = 0.8
N_NEIGHBORS = 15
CELL_JACCARD_MIN = 0.5
N_JOBS = min(6, C.N_JOBS)


def donor_subsets(donors, n_repeats: int, frac: float, seed: int) -> list[np.ndarray]:
    """`n_repeats` random subsets of round(frac * n) donors, from one generator."""
    pool = np.array(sorted(set(map(str, donors))))
    n_keep = int(round(frac * len(pool)))
    rng = np.random.default_rng(seed)
    return [np.sort(rng.choice(pool, size=n_keep, replace=False)) for _ in range(n_repeats)]


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


def best_match(target: np.ndarray, labels: np.ndarray) -> tuple[str | None, float]:
    """The cluster in `labels` with the largest cell Jaccard index with the cells in `target`
    (a boolean mask over the same cells), and that index; (None, 0.0) if `target` is empty."""
    if not target.any():
        return None, 0.0
    cats, codes = np.unique(labels, return_inverse=True)
    inter = np.bincount(codes[target], minlength=len(cats))
    union = target.sum() + np.bincount(codes, minlength=len(cats)) - inter
    j = inter / union
    i = int(np.argmax(j))
    return str(cats[i]), float(j[i])


def stability_units(ladder: dict) -> pd.DataFrame:
    """One row per distinct (model, resolution, cluster, topk) among the picks and family maxima,
    with every role it plays."""
    roles: dict[tuple, list[str]] = {}

    def add(model, res, sig, role):
        cl, k = sig.split("_k")
        roles.setdefault((model, float(res), cl.lstrip("c"), int(k)), []).append(role)

    for an, rows in (("primary", ladder["ladder"]), ("sensitivity", ladder["sensitivity"]["ladder"])):
        for r in rows:
            add(r["model"], r["resolution"], r["signature"], f"{an} pick")
    for an, d in ladder["family_wise"].items():
        for m, v in d.items():
            for key in (v or {}).get("argmax", []):
                res, sig = key.split("/")
                add(m, res, sig, f"{an} family max")
    return pd.DataFrame([{"model": m, "resolution": r, "cluster": cl, "topk": k,
                          "signature": f"c{cl}_k{k}", "roles": "; ".join(v)}
                         for (m, r, cl, k), v in roles.items()])


def summarise(df: pd.DataFrame, suffix: str = "") -> dict:
    n = len(df)
    out = {f"n_repeats{suffix}": n,
           f"n_cells_left_mean{suffix}": float(df["n_cells_left"].mean()) if n else np.nan}
    for col in ("cell_jaccard", "marker_jaccard"):
        out[f"{col}_mean{suffix}"] = float(df[col].mean()) if n else np.nan
        out[f"{col}_sd{suffix}"] = float(df[col].std(ddof=1)) if n > 1 else np.nan
    out[f"frac_cell_jaccard_ge_{CELL_JACCARD_MIN}{suffix}"] = (
        float((df["cell_jaccard"] >= CELL_JACCARD_MIN).mean()) if n else np.nan)
    return out


@cache
def _atlas(path: str, universe_path: str | None):
    import anndata as ad
    a = ad.read_h5ad(path)
    if universe_path:
        a = a[:, a.var_names.isin(set(json.load(open(universe_path))["genes"]))]
    return ad.AnnData(X=a.X.copy(), obs=a.obs[["donor_id"]].astype(str).copy(), var=a.var[[]].copy(),
                      uns={"log1p": {"base": None}})


@cache
def _embedding(path: str) -> np.ndarray:
    return np.load(path, allow_pickle=True)["X"]


def markers(adata, key: str, groups: list[str], n_genes: int) -> dict[str, list[str]]:
    import scanpy as sc
    sc.tl.rank_genes_groups(adata, key, groups=groups, reference="rest", method="wilcoxon",
                            n_genes=n_genes)
    return {g: list(sc.get.rank_genes_groups_df(adata, group=g)["names"]) for g in groups}


def run_repeat(model: str, r: int, kept_donors, full_labels: pd.DataFrame, units: pd.DataFrame,
               sig_genes: dict, top_donor: dict, atlas_path: str, universe_path: str | None,
               emb_path: str, resolutions, min_cells: int, seed: int):
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score
    atlas = _atlas(atlas_path, universe_path)
    mask = atlas.obs["donor_id"].isin(set(kept_donors)).to_numpy()
    a = atlas[mask].copy()
    a.obsm["X_rep"] = _embedding(emb_path)[mask]
    sc.pp.neighbors(a, use_rep="X_rep", n_neighbors=N_NEIGHBORS, random_state=seed)
    ari_rows, unit_rows = [], []
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(a, resolution=res, key_added=key, random_state=seed, flavor="igraph",
                     n_iterations=2, directed=False)
        sub = a.obs[key].astype(str).to_numpy()
        full = full_labels[key].to_numpy()[mask]
        ari_rows.append({"model": model, "resolution": res, "repeat": r,
                         "n_donors_kept": len(kept_donors), "n_cells_kept": int(mask.sum()),
                         "n_clusters_full": int(len(np.unique(full))),
                         "n_clusters_sub": int(len(np.unique(sub))),
                         "ari": float(adjusted_rand_score(full, sub))})
        here = units[(units["model"] == model) & (units["resolution"] == res)]
        if here.empty:
            continue
        matched = {}
        for cl in here["cluster"].unique():
            target = full == cl
            matched[cl] = (int(target.sum()), *best_match(target, sub))
        groups = sorted({m for _, m, _ in matched.values() if m is not None})
        found = {}
        if groups:
            sizes = pd.Series(sub).value_counts()
            pool = np.isin(sub, list(set(sizes[sizes >= min_cells].index) | set(groups)))
            b = a[pool].copy()
            b.obs[key] = pd.Categorical(sub[pool])
            found = markers(b, key, groups, int(here["topk"].max()))
        for u in here.itertuples():
            n_left, m, cj = matched[u.cluster]
            mj = jaccard(found[m][:u.topk], sig_genes[(model, res, u.signature)]) if m is not None else 0.0
            unit_rows.append({"model": model, "resolution": res, "signature": u.signature, "repeat": r,
                              "n_cells_left": n_left, "matched_cluster": m,
                              "matched_size": int((sub == m).sum()) if m is not None else 0,
                              "cell_jaccard": cj, "marker_jaccard": mj,
                              "top_donor_kept": top_donor[(model, res, u.cluster)] in set(kept_donors)})
    return ari_rows, unit_rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--emb-dir", default=str(C.RESULTS / "embeddings"))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--ladder", default=str(C.RESULTS / "ladder" / "ladder.json"))
    ap.add_argument("--models", nargs="+", default=list(C.MODELS))
    ap.add_argument("--resolutions", nargs="+", type=float, default=C.LEIDEN_RESOLUTIONS)
    ap.add_argument("--repeats", type=int, default=N_REPEATS)
    ap.add_argument("--keep-frac", type=float, default=KEEP_DONOR_FRAC)
    ap.add_argument("--min-cells", type=int, default=C.MIN_CELLS_PER_STATE)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--subset-seed", type=int, default=SEED_STABILITY)
    ap.add_argument("--n-jobs", type=int, default=N_JOBS)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    out = Path(args.states_dir)
    if args.dry_run:
        print(json.dumps({"models": args.models, "resolutions": args.resolutions, "repeats": args.repeats,
                          "keep_frac": args.keep_frac, "subset_seed": args.subset_seed,
                          "reads": [args.ladder, str(out / "signatures.json")]}, indent=2)); return

    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[v] = "1"
    from joblib import Parallel, delayed

    from scfm.stages import stability  # workers import the worker by module, not from __main__
    t0 = time.time()
    atlas_path = str(Path(args.data_dir) / "atlas.h5ad")
    gm = Path(args.data_dir) / "gene_map.json"
    universe_path = str(gm) if gm.exists() else None
    atlas = _atlas(atlas_path, universe_path)
    sigs = json.load(open(out / "signatures.json"))
    comp = pd.read_csv(out / "state_composition.csv", dtype={"cluster": str})
    comp_ix = comp.set_index(["model", "resolution", "cluster"])
    units = stability_units(json.load(open(args.ladder)))
    units = units[units["model"].isin(args.models) & units["resolution"].isin(args.resolutions)]
    sig_genes = {(u.model, u.resolution, u.signature): sigs[u.model][str(u.resolution)][u.signature]["genes"]
                 for u in units.itertuples()}
    top_donor = {(m, r, cl): comp_ix.loc[(m, r, cl), "top_donor"]
                 for m, r, cl in units[["model", "resolution", "cluster"]].itertuples(index=False)}
    subsets = donor_subsets(atlas.obs["donor_id"], args.repeats, args.keep_frac, args.subset_seed)
    print(f"  {atlas.n_obs:,} cells, {atlas.obs['donor_id'].nunique()} donors, keeping "
          f"{len(subsets[0])} per repeat; {len(units)} picks; {args.repeats} repeats", flush=True)

    labels = {}
    for m in args.models:
        cells = pd.read_parquet(out / f"cells_{m}.parquet")
        assert list(cells["obs_name"]) == list(atlas.obs_names), f"{m}: cell order mismatch"
        labels[m] = cells.set_index("obs_name")[[f"leiden_{r}" for r in args.resolutions]].astype(str)
    n_donors = int(atlas.obs["donor_id"].nunique())
    del atlas
    _atlas.cache_clear()
    tasks = [(m, r) for m in args.models for r in range(args.repeats)]
    res = Parallel(n_jobs=args.n_jobs, max_nbytes=None)(
        delayed(stability.run_repeat)(m, r, subsets[r], labels[m], units, sig_genes, top_donor, atlas_path,
                            universe_path, str(Path(args.emb_dir) / f"{m}.npz"), args.resolutions,
                            args.min_cells, args.seed)
        for m, r in tasks)
    ari = pd.DataFrame([x for a, _ in res for x in a])
    reps = pd.DataFrame([x for _, u in res for x in u])

    rows = []
    for (m, r), g in ari.groupby(["model", "resolution"], sort=False):
        rows.append({"model": m, "resolution": r, "n_repeats": len(g),
                     "n_donors_kept": int(g["n_donors_kept"].iloc[0]),
                     "n_cells_kept_mean": float(g["n_cells_kept"].mean()),
                     "n_clusters_full_mean": float(g["n_clusters_full"].mean()),
                     "n_clusters_sub_mean": float(g["n_clusters_sub"].mean()),
                     "ari_mean": float(g["ari"].mean()), "ari_sd": float(g["ari"].std(ddof=1)),
                     "ari_min": float(g["ari"].min()), "ari_max": float(g["ari"].max())})
    stab = pd.DataFrame(rows)
    picks = []
    for u in units.itertuples():
        g = reps[(reps["model"] == u.model) & (reps["resolution"] == u.resolution)
                 & (reps["signature"] == u.signature)]
        c = comp_ix.loc[(u.model, u.resolution, u.cluster)]
        picks.append({"model": u.model, "resolution": u.resolution, "signature": u.signature,
                      "cluster": u.cluster, "topk": u.topk, "roles": u.roles,
                      "top_cell_type": c["top_cell_type"], "n_cells": int(c["n_cells"]),
                      "top_donor": c["top_donor"], "top_donor_frac": float(c["top_donor_frac"]),
                      "single_donor": bool(c["single_donor"]),
                      **summarise(g), **summarise(g[g["top_donor_kept"]], "_donor_kept")})
    picks = pd.DataFrame(picks)
    stab.to_csv(out / "stability.csv", index=False)
    picks.to_csv(out / "stability_picks.csv", index=False)
    reps.to_csv(out / "stability_repeats.csv", index=False)
    secs = round(time.time() - t0, 1)
    record_params(out, args, filename="stability_params.json",
                  extra={"seconds": secs, "n_donors": n_donors, "donors_kept_per_repeat": len(subsets[0]),
                         "ladder_source": args.ladder})
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(stab.round(3).to_string(index=False))
        print(picks.drop(columns=["top_donor"]).round(3).to_string(index=False))
    print(f"  {secs:.0f} s", flush=True)


if __name__ == "__main__":
    main()
