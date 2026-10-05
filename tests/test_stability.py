import numpy as np
import pandas as pd
import pytest

from scfm.stages.stability import best_match, donor_subsets, jaccard, markers, stability_units, summarise


def test_donor_subsets_are_reproducible_draws_without_replacement():
    donors = ["d3", "d1", "d2", "d1"] + [f"e{i}" for i in range(6)]
    a = donor_subsets(donors, 5, 0.8, seed=7)
    assert len(a) == 5 and all(len(s) == 7 for s in a)          # round(0.8 * 9 distinct donors)
    assert all(len(set(s)) == len(s) and list(s) == sorted(s) for s in a)
    assert all(set(s) <= set(donors) for s in a)
    b = donor_subsets(list(reversed(donors)), 5, 0.8, seed=7)
    assert all((x == y).all() for x, y in zip(a, b))            # input order does not matter
    assert len({tuple(s) for s in a}) > 1


def test_jaccard():
    assert jaccard(["A", "B", "C"], ["B", "C", "D"]) == 0.5
    assert jaccard([], []) == 0.0


def test_best_match_takes_the_max_cell_jaccard():
    labels = np.array(["x"] * 4 + ["y"] * 6)
    target = np.zeros(10, bool)
    target[[0, 1, 2, 4]] = True                     # 3 of x's 4 cells, 1 of y's 6
    cl, j = best_match(target, labels)
    assert cl == "x" and j == pytest.approx(3 / 5)
    assert best_match(np.zeros(10, bool), labels) == (None, 0.0)


def test_stability_units_merge_roles_of_the_same_signature():
    lad = {
        "ladder": [{"model": "m", "resolution": 1.0, "signature": "c11_k50"}],
        "sensitivity": {"ladder": [{"model": "m", "resolution": 1.0, "signature": "c8_k50"}]},
        "family_wise": {"primary": {"m": {"argmax": ["1.0/c11_k50"]}},
                        "sensitivity": {"m": {"argmax": ["0.3/c14_k50", "1.0/c8_k100"]}, "b": None}},
    }
    u = stability_units(lad).set_index(["resolution", "signature"])
    assert len(u) == 4
    assert u.loc[(1.0, "c11_k50"), "roles"] == "primary pick; primary family max"
    assert u.loc[(1.0, "c8_k100"), ["cluster", "topk"]].tolist() == ["8", 100]
    assert u.loc[(0.3, "c14_k50"), "roles"] == "sensitivity family max"


def test_summarise_threshold_is_inclusive_and_empty_is_nan():
    df = pd.DataFrame({"cell_jaccard": [0.5, 0.2, 1.0, 0.0], "marker_jaccard": [0.4, 0.1, 0.9, 0.0],
                       "n_cells_left": [10, 10, 10, 0]})
    s = summarise(df)
    assert s["n_repeats"] == 4 and s["frac_cell_jaccard_ge_0.5"] == 0.5
    assert s["cell_jaccard_mean"] == pytest.approx(0.425)
    e = summarise(df.iloc[:0], "_donor_kept")
    assert e["n_repeats_donor_kept"] == 0 and np.isnan(e["cell_jaccard_mean_donor_kept"])


def test_markers_rank_a_planted_gene_first_against_the_rest():
    import anndata as ad
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, (60, 8)).astype(np.float32)
    X[:20, 3] += 5
    a = ad.AnnData(X=np.log1p(X), uns={"log1p": {"base": None}})
    a.var_names = [f"g{i}" for i in range(8)]
    a.obs["k"] = pd.Categorical(["a"] * 20 + ["b"] * 20 + ["c"] * 20)
    m = markers(a, "k", ["a"], 3)
    assert list(m) == ["a"] and m["a"][0] == "g3" and len(m["a"]) == 3
