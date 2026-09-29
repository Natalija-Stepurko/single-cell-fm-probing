"""scGPT cell embeddings, run inside the scGPT environment (Python 3.11, torch 2.3, torchtext).

Called by 02_embed.py; not a pipeline stage. Mirrors scgpt.tasks.cell_emb.embed_data: genes
matched to the checkpoint vocabulary, per-cell expression binning, <cls> prepended, CLS output
of the encoder, L2-normalised. Differences are throughput only: cells are sorted by length and
packed to a token budget, the forward pass can run in bfloat16, and results are
written in shards so a restart resumes. scGPT bins ties and subsamples long cells at random,
so every batch is seeded from its first cell index, which makes the output independent of
shard boundaries.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import torch

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module

emb = import_module("02_embed")          # length_batches, run_sharded: numpy-only helpers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--shard-dir", required=True)
    ap.add_argument("--hf-repo", required=True)
    ap.add_argument("--precision", default="bf16", choices=["bf16", "fp32"])
    ap.add_argument("--max-len", type=int, default=1200)
    ap.add_argument("--tokens-per-batch", type=int, default=32768)
    ap.add_argument("--shard-cells", type=int, default=2000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)

    from huggingface_hub import snapshot_download
    from scgpt.data_collator import DataCollator
    from scgpt.model import TransformerModel
    from scgpt.tokenizer import GeneVocab
    from scgpt.utils import load_pretrained

    mdir = Path(snapshot_download(a.hf_repo))
    cfg = json.load(open(mdir / "args.json"))
    vocab = GeneVocab.from_file(mdir / "vocab.json")
    for s in ["<pad>", "<cls>", "<eoc>"]:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])

    z = np.load(a.input)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    genes = z["genes"].tolist()
    keep = np.array([g in vocab for g in genes])
    X = X[:, np.where(keep)[0]].tocsr()
    X.sort_indices()
    gene_ids = np.array(vocab([g for g, k in zip(genes, keep) if k]), dtype=int)
    print(f"    scGPT: matched {keep.sum():,}/{len(genes):,} genes to the vocabulary", flush=True)

    model = TransformerModel(
        ntoken=len(vocab), d_model=cfg["embsize"], nhead=cfg["nheads"], d_hid=cfg["d_hid"],
        nlayers=cfg["nlayers"], nlayers_cls=cfg["n_layers_cls"], n_cls=1, vocab=vocab,
        dropout=cfg["dropout"], pad_token=cfg["pad_token"], pad_value=cfg["pad_value"],
        do_mvc=True, do_dab=False, use_batch_labels=False, domain_spec_batchnorm=False,
        explicit_zero_prob=False, use_fast_transformer=False, pre_norm=False)
    load_pretrained(model, torch.load(mdir / "best_model.pt", map_location="cpu"), verbose=False)
    # whole-model cast: autocast would disable PyTorch's fused encoder path, which skips padding
    dtype = torch.bfloat16 if a.precision == "bf16" else torch.float32
    model.eval().to(dtype)

    collator = DataCollator(do_padding=True, pad_token_id=vocab[cfg["pad_token"]],
                            pad_value=cfg["pad_value"], do_mlm=False, do_binning=True,
                            max_length=a.max_len, sampling=True, keep_first_n_tokens=1)
    cls_id, pad_id = vocab["<cls>"], vocab[cfg["pad_token"]]

    def example(i):
        row = X[i]
        g = np.insert(gene_ids[row.indices], 0, cls_id)
        v = np.insert(row.data.astype(np.float32), 0, cfg["pad_value"])
        return {"id": i, "genes": torch.from_numpy(g).long(), "expressions": torch.from_numpy(v).float()}

    def forward(b):
        np.random.seed(a.seed + b[0]); torch.manual_seed(a.seed + b[0])
        d = collator([example(i) for i in b])
        ids = d["gene"]
        with torch.inference_mode():
            h = model._encode(ids, d["expr"].to(dtype), src_key_padding_mask=ids.eq(pad_id))
        e = h[:, 0, :].float().numpy()
        return e / np.linalg.norm(e, axis=1, keepdims=True)

    lengths = np.minimum(np.diff(X.indptr) + 1, a.max_len)
    t0 = time.time()
    E = emb.run_sharded(X.shape[0], lengths, cfg["embsize"], Path(a.shard_dir), forward,
                        a.tokens_per_batch, a.shard_cells)
    np.save(Path(a.shard_dir) / "embedding.npy", E)
    print(f"    scGPT: {X.shape[0]:,} cells in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
