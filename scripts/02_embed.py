"""Stage 02 — one embedding per cell, per representation, in one schema.

Three arms share the output format so every downstream stage treats them identically:
  hvg_pca     the linear baseline the sceptics defend: 2,000 HVGs → 50 PCs
  scgpt       whole-human checkpoint, last-layer cell embedding (CLS)
  geneformer  6-layer checkpoint, last-layer mean-pooled cell embedding

Output: results/embeddings/<model>.npz with X [n_cells, d] and obs_names, plus params.json.
Foundation-model loaders are isolated in their own functions so a failure to load one model
leaves the others usable; each is CPU batch-1 by design.
"""
import argparse
import json
import sys
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


def embed_geneformer(adata, variant: str, batch_size: int = 1) -> np.ndarray:
    """Rank-value tokenisation → BERT → mean-pooled last hidden state.

    Uses Geneformer's own tokenizer so the gene vocabulary and rank-value encoding match the
    checkpoint exactly; re-implementing it here would be the main source of silent error.
    """
    import torch
    from geneformer import TranscriptomeTokenizer
    from transformers import BertModel

    torch.set_num_threads(max(1, torch.get_num_threads()))
    tok = TranscriptomeTokenizer(custom_attr_name_dict=None, nproc=1)
    model = BertModel.from_pretrained(f"ctheodoris/Geneformer", subfolder=variant,
                                      output_hidden_states=False).eval()
    # Geneformer expects raw counts in adata.X and Ensembl ids in var; the tokenizer handles both
    tokens = tok.tokenize_anndata(adata)          # list of input_id lists, one per cell
    out = np.zeros((len(tokens), model.config.hidden_size), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, len(tokens), batch_size):
            batch = tokens[i:i + batch_size]
            ids = torch.nn.utils.rnn.pad_sequence([torch.tensor(t) for t in batch],
                                                  batch_first=True)
            mask = (ids != 0).long()
            h = model(input_ids=ids, attention_mask=mask).last_hidden_state
            pooled = (h * mask.unsqueeze(-1)).sum(1) / mask.sum(1, keepdim=True)
            out[i:i + len(batch)] = pooled.numpy()
    return out


def embed_scgpt(adata, checkpoint: str, batch_size: int = 1) -> np.ndarray:
    """scGPT whole-human checkpoint → CLS cell embedding, via the package's own embedder."""
    import scgpt as scg
    return scg.tasks.embed_data(adata, model_dir=checkpoint, gene_col="index",
                                batch_size=batch_size, device="cpu",
                                return_new_adata=True).obsm["X_scGPT"].astype(np.float32)


EMBEDDERS = {"hvg_pca": embed_hvg_pca, "geneformer": embed_geneformer, "scgpt": embed_scgpt}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "embeddings"))
    ap.add_argument("--models", nargs="+", default=list(C.MODELS), choices=list(C.MODELS))
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(json.dumps({m: C.MODELS[m] for m in args.models}, indent=2)); return

    import anndata as ad
    adata = ad.read_h5ad(Path(args.data_dir) / "atlas.h5ad")
    done = {}
    for m in args.models:
        dest = out / f"{m}.npz"
        if dest.exists():
            print(f"  {m}: exists, skipping"); done[m] = str(dest); continue
        spec = C.MODELS[m]
        print(f"  {m}: embedding {adata.n_obs:,} cells …", flush=True)
        try:
            if m == "hvg_pca":
                X = embed_hvg_pca(adata, spec["n_comps"], args.seed)
            elif m == "geneformer":
                X = embed_geneformer(adata, spec["variant"])
            else:
                X = embed_scgpt(adata, spec["checkpoint"])
        except Exception as e:                      # one arm failing must not sink the others
            print(f"  {m}: FAILED — {type(e).__name__}: {e}", flush=True); continue
        np.savez_compressed(dest, X=X, obs_names=np.array(adata.obs_names, dtype=object))
        done[m] = str(dest)
        print(f"  {m}: [{X.shape[0]:,} × {X.shape[1]}] -> {dest.name}", flush=True)

    qc.record_params(out, args, extra={"embedded": done})


if __name__ == "__main__":
    main()
