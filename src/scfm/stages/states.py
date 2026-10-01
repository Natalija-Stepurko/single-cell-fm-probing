"""Stage 03 — turn each representation's cell states into gene signatures.

For every embedding: build a neighbour graph, cluster with Leiden at each resolution in the
sweep, and rank marker genes per cluster (Wilcoxon on log-normalised expression). Each cluster
at each (resolution, top-k) becomes one candidate signature. Signature derivation is the step
most likely to dominate the answer, so the ladder is reported at every setting in the sweep.

Also records, per signature, the fraction of its genes that are HVGs — the quantity P4 needs.

Diagnostics written beside the signatures (none of them feeds back into the signatures):
  state_composition.csv  per (model, resolution, cluster): cell, donor, dataset, assay and cell-type
                         make-up, and the share of ribosomal, mitochondrial and nuclear-lncRNA genes
                         among each top-k marker list
  cells_<model>.parquet  per-cell cluster labels and a 2-D UMAP of the same neighbour graph,
                         computed after clustering (not tracked)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm.provenance import record_params

TECHNICAL = {"ribo": r"^RP[SL]\d", "mito": r"^MT-"}
NUCLEAR_LNC = {"MALAT1", "NEAT1", "XIST"}


def marker_shares(genes: list[str]) -> dict:
    g = pd.Series(genes, dtype=str)
    out = {k: float(g.str.contains(rx).mean()) for k, rx in TECHNICAL.items()}
    out["lnc"] = float(g.isin(NUCLEAR_LNC).mean())
    return out


def composition(obs: pd.DataFrame, key: str, min_cells: int) -> list[dict]:
    """Make-up of every cluster under `key`: who and what its cells come from."""
    rows = []
    ngenes = "n_genes_by_counts" if "n_genes_by_counts" in obs else "nnz"
    for cl, g in obs.groupby(key, observed=True):
        donors = g["donor_id"].astype(str).value_counts()
        p = donors.to_numpy() / donors.sum()
        ct = g["cell_type"].astype(str).value_counts()
        rows.append({
            "cluster": str(cl), "n_cells": len(g), "kept": len(g) >= min_cells,
            "n_donors": int(len(donors)), "top_donor": donors.index[0],
            "top_donor_frac": float(p[0]), "donor_entropy": float(-(p * np.log2(p)).sum()),
            "n_datasets": int(g["dataset_id"].nunique()),
            "top_dataset_frac": float(g["dataset_id"].astype(str).value_counts(normalize=True).iloc[0]),
            "top_assay_frac": float(g["assay"].astype(str).value_counts(normalize=True).iloc[0]),
            "median_genes_per_cell": float(g[ngenes].median()),
            "top_cell_type": ct.index[0], "top_cell_type_frac": float(ct.iloc[0] / ct.sum()),
            "second_cell_type": ct.index[1] if len(ct) > 1 else "",
            "second_cell_type_frac": float(ct.iloc[1] / ct.sum()) if len(ct) > 1 else 0.0,
            "single_donor": bool(p[0] >= C.SINGLE_DONOR_FRAC),
        })
    return rows


def signatures_for(adata, X: np.ndarray, resolutions, topks, min_cells, seed, universe=None,
                   details: dict | None = None):
    """Signatures per resolution. When `details` is a dict it receives the per-cell labels
    ("labels"), the composition rows ("composition") and a UMAP of the graph ("umap")."""
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
        comp = {r["cluster"]: r for r in composition(a.obs, key, min_cells)}
        for cl in keep:
            names = list(sc.get.rank_genes_groups_df(sub, group=cl)["names"])
            for k in topks:
                comp[str(cl)].update({f"{n}_frac_k{k}": v for n, v in marker_shares(names[:k]).items()})
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
        if details is not None:
            details.setdefault("composition", []).extend({"resolution": res, **r} for r in comp.values())
    if details is not None:
        details["labels"] = a.obs[[c for c in a.obs if c.startswith("leiden_")]].astype(str)
        sc.tl.umap(a, random_state=seed)
        details["umap"] = np.asarray(a.obsm["X_umap"])
    return out


def main(argv=None):
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
    args = ap.parse_args(argv)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(json.dumps({"models": args.models, "resolutions": args.resolutions,
                          "topk": args.topk, "min_cells": args.min_cells}, indent=2)); return

    import anndata as ad
    adata = ad.read_h5ad(Path(args.data_dir) / "atlas.h5ad")
    gm = Path(args.data_dir) / "gene_map.json"
    universe = set(json.load(open(gm))["genes"]) if gm.exists() else None
    print(f"  marker universe: {len(universe) if universe else adata.n_vars:,} genes", flush=True)
    in_universe = adata.var_names.isin(universe) if universe is not None else np.ones(adata.n_vars, bool)
    universe_n = int(in_universe.sum())
    universe_n_hvg = int((in_universe & adata.var["highly_variable"].to_numpy()).sum())
    sigs, summary, comp_rows = {}, {}, []
    for m in args.models:
        p = Path(args.emb_dir) / f"{m}.npz"
        if not p.exists():
            print(f"  {m}: no embedding, skipping"); continue
        z = np.load(p, allow_pickle=True)
        assert list(z["obs_names"]) == list(adata.obs_names), f"{m}: cell order mismatch"
        print(f"  {m}: clustering …", flush=True)
        det = {}
        sigs[m] = signatures_for(adata, z["X"], args.resolutions, args.topk, args.min_cells,
                                 args.seed, universe, details=det)
        comp_rows.extend({"model": m, **r} for r in det["composition"])
        cells = det["labels"].reset_index(names="obs_name")
        cells["umap_1"], cells["umap_2"] = det["umap"][:, 0], det["umap"][:, 1]
        cells.to_parquet(out / f"cells_{m}.parquet", index=False)
        n = sum(len(v) for v in sigs[m].values())
        summary[m] = {"n_signatures": n,
                      "n_states": {r: len({s["cluster"] for s in v.values()})
                                   for r, v in sigs[m].items()}}
        print(f"  {m}: {n} signatures across {len(sigs[m])} resolutions", flush=True)

    json.dump(sigs, open(out / "signatures.json", "w"))
    summary.update({"universe_n": universe_n, "universe_n_hvg": universe_n_hvg})
    json.dump(summary, open(out / "summary.json", "w"), indent=2)
    pd.DataFrame(comp_rows).to_csv(out / "state_composition.csv", index=False)
    record_params(out, args, extra=summary)


if __name__ == "__main__":
    main()
