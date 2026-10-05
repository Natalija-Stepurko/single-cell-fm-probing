"""Geneformer V1: rank-value tokenisation and mean-pooled last-layer cell embeddings."""
from pathlib import Path

import numpy as np

from scfm import config as C
from scfm.batching import run_sharded


def torch_setup(precision):
    import torch
    torch.set_num_threads(C.TORCH_THREADS)
    dtype = torch.bfloat16 if precision == "bf16" else torch.float32
    return torch, dtype


def geneformer_tokens(adata, spec):
    """Rank-value encoding exactly as geneformer.tokenizer.TranscriptomeTokenizer (V1):
    counts / n_counts * 10,000 / gene median, nonzero genes ranked descending, first 2,048."""
    import pickle

    import scipy.sparse as sp
    from huggingface_hub import hf_hub_download

    d = spec["dict_dir"]

    def load(n):
        path = hf_hub_download(spec["hf_repo"], f"{d}/{n}", revision=spec.get("revision"))
        with open(path, "rb") as f:
            return pickle.load(f)

    token_dict = load("token_dictionary_gc30M.pkl")
    median_dict = load("gene_median_dictionary_gc30M.pkl")
    mapping = load("ensembl_mapping_dict_gc30M.pkl")

    ens = adata.var["ensembl_id"].astype(str).str.upper().map(mapping)
    assert ens.dropna().is_unique, "Ensembl ids collapse onto each other; sum them first"
    loc = np.where([isinstance(e, str) and e in median_dict for e in ens])[0]
    norm = np.array([median_dict[e] for e in ens.iloc[loc]])
    toks = np.array([token_dict[e] for e in ens.iloc[loc]])

    # same arithmetic and dtypes as the official tokenizer, chunk by chunk (float32 counts
    # divided by n_counts and scaled, then float64 division by the gene medians), so that
    # near-tied genes rank in the same order
    X = sp.csr_matrix(adata.layers["counts"]).astype(np.float32)
    n_counts = adata.obs["n_counts"].to_numpy()[:, None]
    out = []
    for c0 in range(0, X.shape[0], 512):
        Xn = sp.csr_matrix(X[c0:c0 + 512][:, loc] / n_counts[c0:c0 + 512] * 10_000 / norm)
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
    path = snapshot_download(spec["hf_repo"], revision=spec.get("revision"),
                             allow_patterns=[f"{spec['variant']}/*"])
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
