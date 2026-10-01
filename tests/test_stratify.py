import numpy as np
import pandas as pd
import pytest

from scfm.stages.stratify import choose_penalizer, paired, run_cv, state_features, xgb_risk
from scfm.stats import cindex_many, cv_folds


def _cohort(n=240, p=6, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    t = rng.exponential(np.exp(-(0.8 * X[:, 0] - 0.5 * X[:, 1])))
    d = pd.DataFrame(X, columns=[f"x{i}" for i in range(p)])
    d["T"], d["E"] = t, (rng.random(n) < 0.7).astype(float)
    return d


def test_state_features_take_one_topk_and_zscore():
    z = pd.DataFrame(np.random.default_rng(1).normal(size=(6, 20)), index=list("ABCDEF"),
                     columns=[f"p{i}" for i in range(20)])
    sigs = {"m": {"0.3": {"c0_k50": {"topk": 50, "genes": ["A", "B", "C"]},
                          "c0_k25": {"topk": 25, "genes": ["A", "B", "C"]}},
                  "1.0": {"c1_k50": {"topk": 50, "genes": ["D", "E", "F", "ZZ"]},
                          "c2_k50": {"topk": 50, "genes": ["A", "Y1", "Y2"]}}}}   # < 3 genes measured
    F = state_features(sigs, "m", z, z.columns, 50)
    assert list(F.columns) == ["m:0.3/c0_k50", "m:1.0/c1_k50"]
    np.testing.assert_allclose(F.mean(), 0, atol=1e-12)
    np.testing.assert_allclose(F.std(), 1)


def test_choose_penalizer_is_deterministic_and_from_the_grid():
    d = _cohort()
    tr = np.arange(len(d))
    a = choose_penalizer(d, "T", "E", ["x0", "x1", "x2"], tr, [0.1, 10.0], np.random.default_rng(4))
    b = choose_penalizer(d, "T", "E", ["x0", "x1", "x2"], tr, [0.1, 10.0], np.random.default_rng(4))
    assert a == b and a[0] in (0.1, 10.0) and set(a[1]) == {0.1, 10.0}


def test_xgb_risk_learns_a_cox_signal():
    d = _cohort(400)
    X, T, E = d[["x0", "x1"]].to_numpy(), d["T"].to_numpy(), d["E"].to_numpy()
    tr, te = np.arange(300), np.arange(300, 400)
    r = xgb_risk(X, T, E, tr, te)
    assert cindex_many(r[None, :], T[te], E[te])[0] > 0.6
    np.testing.assert_array_equal(r, xgb_risk(X, T, E, tr, te))


def test_run_cv_pairs_models_on_shared_folds():
    d = _cohort(200)
    folds = cv_folds(d["E"], 3, 2, np.random.default_rng(0))
    specs = {"a": ("cox", ["x0"]), "b": ("cox", ["x0", "x1"]), "r": ("ridge", ["x0", "x1", "x2"])}
    res = run_cv(specs, d, "T", "E", folds, n_jobs=1)
    assert res["a"]["cindex"].shape == (2,)
    assert res["b"]["cindex"].mean() > res["a"]["cindex"].mean()
    assert all(len(p) == 3 for p in res["r"]["penalizers"])
    pr = paired(res["b"]["cindex"], res["a"]["cindex"])
    assert pr["per_repeat"] == pytest.approx((res["b"]["cindex"] - res["a"]["cindex"]).tolist())


def test_ridge_penalises_only_the_state_features():
    from scfm import config as C
    from scfm.stages.stratify import ridge_penalties
    cols = list(C.COVARIATES) + ["m:0.3/c0_k50", "m:1.0/c1_k50"]
    np.testing.assert_allclose(ridge_penalties(cols, 10.0), [0.01] * len(C.COVARIATES) + [10.0, 10.0])
