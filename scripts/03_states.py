"""Stage 03 — turn each representation's cell states into gene signatures.

For every embedding: build a neighbour graph, cluster with Leiden at each resolution in the
sweep, and rank marker genes per cluster (Wilcoxon on log-normalised expression). Each cluster
at each (resolution, top-k) becomes one candidate signature. Signature derivation is the step
most likely to dominate the answer, so the ladder is reported at every setting in the sweep.

Also records, per signature, the fraction of its genes that are HVGs — the quantity P4 needs.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import config as C
import qc_common as qc


def signatures_for(adata, X: np.ndarray, resolutions, topks, min_cells, seed, universe=None):
    import scanpy as sc
    a = adata.copy()
    a.obsm["X_rep"] = X
    sc.pp.neighbors(a, use_rep="X_rep", n_neighbors=15, random_state=seed)
    hvg = set(a.var_names[a.var["highly_variable"]])
    out = {}
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(a, resolution=res, key_added=key, random_state=seed, flavor="igraph",
                     n_iterations=2, directed=False)
        sizes = a.obs[key].value_counts()
        keep = sizes[sizes >= min_cells].index.tolist()
        if len(keep) < 2:
            continue
        sub = a[a.obs[key].isin(keep)]
        # markers are drawn only from genes the bulk cohort measures, so a top-k signature
        # is scored in bulk with all k of its genes
        sub = (sub[:, sub.var_names.isin(universe)] if universe is not None else sub).copy()
        sc.tl.rank_genes_groups(sub, key, method="wilcoxon", n_genes=max(topks))
        out[str(res)] = {}
        for cl in keep:
            names = list(sc.get.rank_genes_groups_df(sub, group=cl)["names"])
            top_ct = sub.obs.loc[sub.obs[key] == cl, "cell_type"].value_counts()
            for k in topks:
                genes = names[:k]
                out[str(res)][f"c{cl}_k{k}"] = {
                    "cluster": str(cl), "resolution": res, "topk": k, "genes": genes,
                    "n_cells": int(sizes[cl]),
                    "top_cell_type": str(top_ct.index[0]),
                    "top_cell_type_frac": float(top_ct.iloc[0] / top_ct.sum()),
                    "hvg_frac": float(np.mean([g in hvg for g in genes])),
                }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--emb-dir", default=str(C.RESULTS / "embeddings"))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--models", nargs="+", default=list(C.MODELS))
    ap.add_argument("--resolutions", nargs="+", type=float, default=C.LEIDEN_RESOLUTIONS)
    ap.add_argument("--topk", nargs="+", type=int, default=C.MARKER_TOPK)
    ap.add_argument("--min-cells", type=int, default=C.MIN_CELLS_PER_STATE)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(json.dumps({"models": args.models, "resolutions": args.resolutions,
                          "topk": args.topk, "min_cells": args.min_cells}, indent=2)); return

    import anndata as ad
    adata = ad.read_h5ad(Path(args.data_dir) / "atlas.h5ad")
    gm = Path(args.data_dir) / "gene_map.json"
    universe = set(json.load(open(gm))["genes"]) if gm.exists() else None
    print(f"  marker universe: {len(universe) if universe else adata.n_vars:,} genes", flush=True)
    sigs, summary = {}, {}
    for m in args.models:
        p = Path(args.emb_dir) / f"{m}.npz"
        if not p.exists():
            print(f"  {m}: no embedding, skipping"); continue
        z = np.load(p, allow_pickle=True)
        assert list(z["obs_names"]) == list(adata.obs_names), f"{m}: cell order mismatch"
        print(f"  {m}: clustering …", flush=True)
        sigs[m] = signatures_for(adata, z["X"], args.resolutions, args.topk, args.min_cells,
                                 args.seed, universe)
        n = sum(len(v) for v in sigs[m].values())
        summary[m] = {"n_signatures": n,
                      "n_states": {r: len({s["cluster"] for s in v.values()})
                                   for r, v in sigs[m].items()}}
        print(f"  {m}: {n} signatures across {len(sigs[m])} resolutions", flush=True)

    json.dump(sigs, open(out / "signatures.json", "w"))
    json.dump(summary, open(out / "summary.json", "w"), indent=2)
    qc.record_params(out, args, extra=summary)


if __name__ == "__main__":
    main()
