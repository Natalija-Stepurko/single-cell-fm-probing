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

sys.path.insert(0, str(Path(__file__).parent))
import config as C
import qc_common as qc


def embed_hvg_pca(adata, n_comps: int, seed: int) -> np.ndarray:
    import scanpy as sc
    a = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=n_comps, random_state=seed)
    return a.obsm["X_pca"].astype(np.float32)


def length_batches(lengths: np.ndarray, tokens_per_batch: int):
    """Indices sorted longest-first, cut into batches whose padded size fits the budget."""
    order = np.argsort(-lengths, kind="stable")
    batches, cur, cur_max = [], [], 0
    for i in order:
        m = max(cur_max, int(lengths[i]))
        if cur and m * (len(cur) + 1) > tokens_per_batch:
            batches.append(cur); cur, m = [], int(lengths[i])
        cur.append(int(i)); cur_max = m
    if cur:
        batches.append(cur)
    return batches


def run_sharded(n_cells, lengths, dim, shard_dir: Path, forward, tokens_per_batch, shard_cells):
    """Drive `forward(batch_indices) -> [len(batch), dim]` over all cells, shard by shard."""
    shard_dir.mkdir(parents=True, exist_ok=True)
    batches = length_batches(lengths, tokens_per_batch)
    shards, cur = [], []
    for b in batches:
        cur.append(b)
        if sum(map(len, cur)) >= shard_cells:
            shards.append(cur); cur = []
    if cur:
        shards.append(cur)
    out = np.zeros((n_cells, dim), dtype=np.float32)
    t0, done = time.time(), 0
    for k, shard in enumerate(shards):
        f = shard_dir / f"shard_{k:04d}.npz"
        idx = np.concatenate([np.asarray(b) for b in shard])
        if f.exists():
            z = np.load(f)
            assert np.array_equal(z["idx"], idx), f"{f.name}: shard layout changed; clear {shard_dir}"
            out[idx] = z["X"]
        else:
            X = np.concatenate([forward(b) for b in shard])
            np.savez(f, idx=idx, X=X)
            out[idx] = X
            done += len(idx)
        el = time.time() - t0
        print(f"    shard {k + 1}/{len(shards)}  {int(idx.size)} cells  "
              f"{done / el if el > 0 and done else 0:.1f} cells/s", flush=True)
    return out


def torch_setup(precision):
    import torch
    torch.set_num_threads(C.TORCH_THREADS)
    dtype = torch.bfloat16 if precision == "bf16" else torch.float32
    return torch, dtype


# ───────────────────────────────────────────────────────────────── Geneformer
def geneformer_tokens(adata, spec):
    """Rank-value encoding exactly as geneformer.tokenizer.TranscriptomeTokenizer (V1):
    counts / n_counts * 10,000 / gene median, nonzero genes ranked descending, first 2,048."""
    import pickle
    import scipy.sparse as sp
    from huggingface_hub import hf_hub_download

    d = spec["dict_dir"]
    load = lambda n: pickle.load(open(hf_hub_download(spec["hf_repo"], f"{d}/{n}"), "rb"))
    token_dict = load("token_dictionary_gc30M.pkl")
    median_dict = load("gene_median_dictionary_gc30M.pkl")
    mapping = load("ensembl_mapping_dict_gc30M.pkl")

    ens = adata.var["ensembl_id"].astype(str).str.upper().map(mapping)
    assert ens.dropna().is_unique, "Ensembl ids collapse onto each other; sum them first"
    loc = np.where([isinstance(e, str) and e in median_dict for e in ens])[0]
    norm = np.array([median_dict[e] for e in ens.iloc[loc]])
    toks = np.array([token_dict[e] for e in ens.iloc[loc]])

    X = sp.csr_matrix(adata.layers["counts"])[:, loc].astype(np.float64)
    n_counts = adata.obs["n_counts"].to_numpy()[:, None]
    Xn = sp.csr_matrix(X.multiply(1.0 / n_counts).multiply(10_000).multiply(1.0 / norm[None, :]))
    Xn.sort_indices()
    out = []
    for i in range(Xn.shape[0]):
        row = Xn[i]
        out.append(toks[row.indices][np.argsort(-row.data)][: spec["max_len"]])
    print(f"    tokenised {len(out):,} cells over {len(loc):,} Geneformer genes", flush=True)
    return out


def embed_geneformer(adata, spec, precision, shard_dir):
    from huggingface_hub import snapshot_download
    from transformers import BertModel
    torch, dtype = torch_setup(precision)

    tokens = geneformer_tokens(adata, spec)
    path = snapshot_download(spec["hf_repo"], allow_patterns=[f"{spec['variant']}/*"])
    model = BertModel.from_pretrained(Path(path) / spec["variant"], add_pooling_layer=False,
                                      torch_dtype=dtype).eval()

    def forward(b):
        ids = torch.nn.utils.rnn.pad_sequence(
            [torch.as_tensor(tokens[i], dtype=torch.long) for i in b], batch_first=True)
        mask = (ids != 0).long()
        with torch.inference_mode():
            h = model(input_ids=ids, attention_mask=mask).last_hidden_state.float()
        m = mask.unsqueeze(-1).float()
        return ((h * m).sum(1) / m.sum(1)).numpy()

    lengths = np.array([len(t) for t in tokens])
    return run_sharded(len(tokens), lengths, model.config.hidden_size, shard_dir, forward,
                       C.TOKENS_PER_BATCH, C.SHARD_CELLS)


# ───────────────────────────────────────────────────────────────── scGPT
def embed_scgpt(adata, spec, precision, shard_dir):
    """Hand raw counts to scripts/scgpt_worker.py, which runs in the scGPT environment."""
    import scipy.sparse as sp
    shard_dir.mkdir(parents=True, exist_ok=True)
    inp = shard_dir / "input.npz"
    X = sp.csr_matrix(adata.layers["counts"])
    np.savez(inp, data=X.data.astype(np.float32), indices=X.indices, indptr=X.indptr,
             shape=np.array(X.shape), genes=np.array(adata.var_names, dtype=str))
    env = dict(os.environ, OMP_NUM_THREADS=str(C.TORCH_THREADS))
    env.pop("VIRTUAL_ENV", None)
    cmd = [C.SCGPT_PYTHON, str(Path(__file__).parent / "scgpt_worker.py"),
           "--input", str(inp), "--shard-dir", str(shard_dir), "--precision", precision,
           "--hf-repo", spec["hf_repo"], "--max-len", str(spec["max_len"]),
           "--tokens-per-batch", str(C.TOKENS_PER_BATCH), "--shard-cells", str(C.SHARD_CELLS),
           "--threads", str(C.TORCH_THREADS), "--seed", str(C.SEED)]
    rc = subprocess.call(cmd, env=env)
    if rc != 0:
        raise RuntimeError(f"scgpt_worker exited {rc}")
    return np.load(shard_dir / "embedding.npy")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "embeddings"))
    ap.add_argument("--models", nargs="+", default=list(C.MODELS), choices=list(C.MODELS))
    ap.add_argument("--precision", default=C.EMBED_DTYPE, choices=["bf16", "fp32"])
    ap.add_argument("--limit", type=int, default=None,
                    help="embed only the first N cells (precision check); writes <model>_<prec>_N.npz")
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
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

    qc.record_params(out, args, extra={"embedded": done, "seconds": timing, "failed": failed})
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
