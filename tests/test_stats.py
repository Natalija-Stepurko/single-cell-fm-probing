import numpy as np
import pandas as pd
import pytest
from lifelines.utils import concordance_index

from scfm.stats import (
    best_above_floor,
    cindex_many,
    cv_folds,
    direction_of,
    fw_pvalue,
    orient,
    permute_within_strata,
    selection_aware_margin,
)


def test_fw_pvalue_counts_ties_as_extreme():
    null = [0.50, 0.55, 0.60, 0.60, 0.70]
    assert fw_pvalue(0.60, null) == pytest.approx((1 + 3) / 6)
    assert fw_pvalue(0.80, null) == pytest.approx(1 / 6)
    assert fw_pvalue(0.40, null) == pytest.approx(1.0)


def test_orientation_helpers():
    assert list(direction_of([0.6, 0.5, 0.4])) == [1, 1, -1]
    assert orient(0.4, -1) == pytest.approx(0.6)
    assert orient(0.4, 1) == pytest.approx(0.4)
    np.testing.assert_allclose(orient(np.array([[0.3, 0.7]]), np.array([[-1, -1]])), [[0.7, 0.3]])


def test_cindex_many_matches_lifelines_with_ties():
    rng = np.random.default_rng(0)
    T = rng.integers(1, 15, 60).astype(float)          # many tied times
    E = rng.random(60) < 0.4
    S = np.round(rng.normal(size=(5, 60)), 1)          # tied scores
    got = cindex_many(S, T, E)
    want = [concordance_index(T, -s, E) for s in S]
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_cindex_many_on_a_bootstrap_resample_with_duplicates():
    rng = np.random.default_rng(1)
    T, E, s = rng.exponential(size=40), rng.random(40) < 0.5, rng.normal(size=40)
    ix = rng.integers(0, 40, 40)
    assert cindex_many(s[None, ix], T[ix], E[ix])[0] == pytest.approx(concordance_index(T[ix], -s[ix], E[ix]))


def test_selection_aware_margin_reselects_per_resample():
    floor_fm, floor_b = np.array([0.50, 0.45]), np.array([0.40, 0.50])
    Cb_fm = np.array([[0.60, 0.50],    # fm pick 0 (+0.10)
                      [0.52, 0.60],    # fm pick 1 (+0.15)
                      [0.55, 0.55]])   # fm pick 1 (+0.10)
    Cb_b = np.array([[0.48, 0.60],     # base pick 1 (+0.10)
                     [0.55, 0.52],     # base pick 0 (+0.15)
                     [0.40, 0.45]])    # base pick 0 (+0.00)
    m, jf, jb = selection_aware_margin(Cb_fm, floor_fm, Cb_b, floor_b)
    np.testing.assert_allclose(m, [0.0, 0.0, 0.10], atol=1e-12)
    assert jf.tolist() == [0, 1, 1] and jb.tolist() == [1, 0, 0]
    af, j = best_above_floor(np.array([[np.nan, 0.6]]), np.array([0.5, 0.5]))
    assert j[0] == 1 and af[0] == pytest.approx(0.1)


def test_permute_within_strata_keeps_strata():
    strata = np.array(list("aaabbbbcc") * 3)
    a = permute_within_strata(strata, np.random.default_rng(3))
    b = permute_within_strata(strata, np.random.default_rng(3))
    assert sorted(a) == list(range(len(strata)))
    assert (strata[a] == strata).all()
    assert (a == b).all()
    assert pd.Series(strata[a]).value_counts().equals(pd.Series(strata).value_counts())
    assert not (a == np.arange(len(strata))).all()


def test_cv_folds_stratified_and_deterministic():
    ev = np.r_[np.ones(20), np.zeros(80)]
    f1 = cv_folds(ev, 5, 2, np.random.default_rng(7))
    f2 = cv_folds(ev, 5, 2, np.random.default_rng(7))
    assert len(f1) == 2 and all(len(r) == 5 for r in f1)
    for rep in f1:
        assert sorted(np.concatenate([te for _, te in rep]).tolist()) == list(range(100))
        assert all(ev[te].sum() == 4 for _, te in rep)
    assert all((a[1] == b[1]).all() for ra, rb in zip(f1, f2) for a, b in zip(ra, rb))
    assert not all((a[1] == b[1]).all() for a, b in zip(f1[0], f1[1]))


def test_fw_pvalue_counts_a_tie_lost_to_rounding():
    stat = 0.1 + 0.2                                 # 0.30000000000000004
    null = [0.7 - 0.4, 0.2]                          # 0.29999999999999993
    assert null[0] < stat
    assert fw_pvalue(stat, null) == pytest.approx(2 / 3)


def test_mc_annotate_flags_borderline_p():
    from scfm.stats import mc_annotate
    a = mc_annotate(0.048, 1000)
    assert a["mc_se"] == pytest.approx(np.sqrt(0.048 * 0.952 / 1000))
    assert a["borderline"]
    assert not mc_annotate(0.01, 10_000)["borderline"]
