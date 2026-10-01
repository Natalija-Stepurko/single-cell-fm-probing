"""Length-bucketed batching and sharded, resumable forward passes.

numpy-only and independent of the rest of the package: the scGPT worker imports it from its own
Python 3.11 environment.
"""
import time
from pathlib import Path

import numpy as np


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
