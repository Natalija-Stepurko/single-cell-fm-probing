"""The tokeniser against tokens produced by the official geneformer TranscriptomeTokenizer (V1)
on a 30-cell fixture that includes a tied cell and genes outside the vocabulary."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from scfm import config as C

FIX = Path(__file__).parent / "fixtures"
SPEC = C.MODELS["geneformer"]
DICTS = ("token_dictionary_gc30M.pkl", "gene_median_dictionary_gc30M.pkl", "ensembl_mapping_dict_gc30M.pkl")


def _dictionaries_available() -> bool:
    try:
        from huggingface_hub import hf_hub_download
        for name in DICTS:
            hf_hub_download(SPEC["hf_repo"], f"{SPEC['dict_dir']}/{name}", revision=SPEC["revision"])
        return True
    except Exception:
        return False


@pytest.mark.hf
def test_tokens_match_official_tokenizer():
    if not _dictionaries_available():
        pytest.skip("Geneformer V1 dictionaries not available (no network and no local HF cache)")
    import anndata as ad

    from scfm.geneformer import geneformer_tokens

    z = np.load(FIX / "geneformer_tokens_fixture.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    genes = z["ensembl_id"].astype(str)
    obs = pd.DataFrame({"n_counts": z["n_counts"]}, index=[f"c{i}" for i in range(X.shape[0])])
    adata = ad.AnnData(X=X.copy(), obs=obs, var=pd.DataFrame({"ensembl_id": genes}, index=genes))
    adata.layers["counts"] = X

    expected = json.loads((FIX / "geneformer_tokens_expected.json").read_text())
    assert expected["model_input_size"] == SPEC["max_len"]
    got = geneformer_tokens(adata, SPEC)
    assert len(got) == len(expected["tokens"]) == 30
    for i, (g, e) in enumerate(zip(got, expected["tokens"], strict=True)):
        assert np.asarray(g).tolist() == e, f"cell {i}: tokens differ"
