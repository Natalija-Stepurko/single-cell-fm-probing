"""Stage 02 — one embedding per cell, per representation, in one schema.

Three arms share the output format so every downstream stage treats them identically:
  hvg_pca     the linear baseline: 2,000 HVGs → 50 PCs
  geneformer  6-layer V1 checkpoint, last-layer mean-pooled cell embedding
  scgpt       whole-human checkpoint, CLS cell embedding (run in its own environment)

Output: results/embeddings/<model>.npz with X [n_cells, d] and obs_names, plus params.json.

CPU throughput: cells are sorted by token count and packed into batches of a fixed token
budget, so padding is small; the forward pass runs in bfloat16 (AMX on this CPU); results are
written in shards so an interrupted run resumes where it stopped. `--precision fp32 --limit N`
embeds the first N cells in full precision for the bf16 agreement check.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from scfm import config as C
from scfm.geneformer import embed_geneformer
from scfm.provenance import record_params


def embed_hvg_pca(adata, n_comps: int, seed: int) -> np.ndarray:
    import scanpy as sc
    a = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=n_comps, random_state=seed)
    return a.obsm["X_pca"].astype(np.float32)


# ───────────────────────────────────────────────────────────────── scGPT
def embed_scgpt(adata, spec, precision, shard_dir):
    """Hand raw counts to scfm/scgpt_worker.py, which runs in the scGPT environment."""
    import scipy.sparse as sp
    shard_dir.mkdir(parents=True, exist_ok=True)
    inp = shard_dir / "input.npz"
    X = sp.csr_matrix(adata.layers["counts"])
    np.savez(inp, data=X.data.astype(np.float32), indices=X.indices, indptr=X.indptr,
             shape=np.array(X.shape), genes=np.array(adata.var_names, dtype=str))
    env = dict(os.environ, OMP_NUM_THREADS=str(C.TORCH_THREADS))
    env.pop("VIRTUAL_ENV", None)
    cmd = [C.SCGPT_PYTHON, str(Path(__file__).resolve().parents[1] / "scgpt_worker.py"),
           "--input", str(inp), "--shard-dir", str(shard_dir), "--precision", precision,
           "--hf-repo", spec["hf_repo"], "--revision", spec["revision"], "--max-len", str(spec["max_len"]),
           "--tokens-per-batch", str(C.TOKENS_PER_BATCH), "--shard-cells", str(C.SHARD_CELLS),
           "--threads", str(C.TORCH_THREADS), "--seed", str(C.SEED)]
    rc = subprocess.call(cmd, env=env)
    if rc != 0:
        raise RuntimeError(f"scgpt_worker exited {rc}")
    return np.load(shard_dir / "embedding.npy")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "embeddings"))
    ap.add_argument("--models", nargs="+", default=list(C.MODELS), choices=list(C.MODELS))
    ap.add_argument("--precision", default=C.EMBED_DTYPE, choices=["bf16", "fp32"])
    ap.add_argument("--limit", type=int, default=None,
                    help="embed only the first N cells (precision check); writes <model>_<prec>_N.npz")
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(json.dumps({"models": {m: C.MODELS[m] for m in args.models},
                          "precision": args.precision, "threads": C.TORCH_THREADS,
                          "tokens_per_batch": C.TOKENS_PER_BATCH,
                          "shard_cells": C.SHARD_CELLS}, indent=2)); return

    import anndata as ad
    adata = ad.read_h5ad(Path(args.data_dir) / "atlas.h5ad")
    if args.limit:
        adata = adata[: args.limit].copy()
    done, timing, failed = {}, {}, {}
    for m in args.models:
        tag = f"{m}_{args.precision}_{args.limit}" if args.limit else m
        dest = out / f"{tag}.npz"
        if dest.exists():
            print(f"  {m}: exists, skipping"); done[m] = str(dest); continue
        spec = C.MODELS[m]
        print(f"  {m}: embedding {adata.n_obs:,} cells ({args.precision}) …", flush=True)
        t0 = time.time()
        try:
            if m == "hvg_pca":
                X = embed_hvg_pca(adata, spec["n_comps"], args.seed)
            elif m == "geneformer":
                X = embed_geneformer(adata, spec, args.precision, out / f"_shards_{tag}")
            else:
                X = embed_scgpt(adata, spec, args.precision, out / f"_shards_{tag}")
        except Exception as e:                      # one arm failing must not sink the others
            print(f"  {m}: FAILED — {type(e).__name__}: {e}", flush=True)
            failed[m] = f"{type(e).__name__}: {e}"; continue
        assert np.isfinite(X).all(), f"{m}: non-finite embedding values"
        timing[m] = round(time.time() - t0, 1)
        np.savez_compressed(dest, X=X, obs_names=np.array(adata.obs_names, dtype=object))
        done[m] = str(dest)
        print(f"  {m}: [{X.shape[0]:,} × {X.shape[1]}] -> {dest.name}  "
              f"{timing[m]:.0f}s ({adata.n_obs / max(timing[m], 1e-9):.1f} cells/s)", flush=True)

    record_params(out, args, extra={"embedded": done, "seconds": timing, "failed": failed})
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
