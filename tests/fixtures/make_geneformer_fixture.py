"""Regenerate the Geneformer tokeniser fixture with the official tokenizer.

Needs the `geneformer` package (huggingface.co/ctheodoris/Geneformer at the revision pinned in scfm.config)
and its V1 (gc30M) dictionaries. Usage: python make_geneformer_fixture.py <dir of the gc30M dictionaries>.
Writes the two files tests/test_geneformer_tokens.py reads.
"""
import json
import pickle
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from geneformer.tokenizer import TranscriptomeTokenizer

D = sys.argv[1]   # .../geneformer/gene_dictionaries_30m
with open(D + "/gene_median_dictionary_gc30M.pkl", "rb") as f:
    med = pickle.load(f)
rng=np.random.default_rng(7)
genes=sorted(med)[:15000:5]                   # 3,000 real V1 genes
genes = genes + ["ENSG99999999990", "ENSG99999999991"]  # two genes outside the vocabulary
n_cells=30
X = sp.random(n_cells, len(genes), density=0.35, random_state=7, format="csr",
              data_rvs=lambda k: rng.integers(1, 40, k)).astype(np.float32)
X[3,:]=0; X[3,5]=2; X[3,9]=2; X.eliminate_zeros()   # a cell with an exact tie
X=X.tolil()
X[4,:2600]=rng.integers(1,60,2600)            # more expressed genes than the 2,048-token input
X[5,:]=0; X[5,10]=3; X[5,11]=3; X[5,12]=7      # equal counts on genes with near-identical medians
X=X.tocsr(); X.eliminate_zeros()
# totals include genes not in the matrix
n_counts = np.asarray(X.sum(1)).ravel() + rng.integers(0, 500, n_cells)
obs = pd.DataFrame({"n_counts": n_counts.astype(np.float32)}, index=[f"c{i}" for i in range(n_cells)])
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame({"ensembl_id": genes}, index=genes))
a.write_h5ad("gf_fixture_in.h5ad")
tk=TranscriptomeTokenizer(model_version="V1", nproc=1)
cells,_,_=tk.tokenize_anndata("gf_fixture_in.h5ad")
np.savez_compressed("geneformer_tokens_fixture.npz", data=X.data, indices=X.indices, indptr=X.indptr,
                    shape=np.array(X.shape), ensembl_id=np.array(genes), n_counts=n_counts.astype(np.float32))
expected = [np.asarray(c[:tk.model_input_size]).tolist() for c in cells]
with open("geneformer_tokens_expected.json", "w") as f:
    json.dump({"model_input_size": tk.model_input_size, "tokens": expected}, f)
print(len(cells), [len(c) for c in cells][:5])
