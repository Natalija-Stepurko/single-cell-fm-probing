import numpy as np
import pandas as pd
import pytest
from lifelines.utils import concordance_index
from scipy.stats import kruskal, mannwhitneyu

from scfm import config as C
from scfm.stages.donors import donor_mixing, states_with_margins
from scfm.stages.ladder import family_wise, family_wise_large, fw_nulls_large
from scfm.stages.subtype import association
from scfm.stats import auc_higher_in, family_max_nulls, fw_pvalue, kruskal_eta2, permuted_cindex


def _surv(n=60, seed=0):
    rng = np.random.default_rng(seed)
    T = rng.integers(1, 40, n).astype(float)        # integer times: ties in T
    E = (rng.random(n) < 0.5).astype(float)
    S = np.round(rng.normal(size=(5, n)), 1)        # rounded scores: ties in S
    return S, T, E


def test_permuted_cindex_matches_lifelines_on_permuted_outcomes():
    S, T, E = _surv()
    rng = np.random.default_rng(1)
    perms = [rng.permutation(len(T)) for _ in range(4)]
    fast = permuted_cindex(S, T, E, perms)
    for b, p in enumerate(perms):
        for j in range(len(S)):
            assert fast[b, j] == pytest.approx(concordance_index(T[p], -S[j], E[p]), abs=1e-12)


def test_family_max_nulls_reads_both_directions():
    cmat = np.array([[0.40, 0.55], [0.30, 0.52]])
    best, best_or = family_max_nulls(cmat)
    np.testing.assert_allclose(best, [0.55, 0.52])
    np.testing.assert_allclose(best_or, [0.60, 0.70])


def test_fw_pvalue_with_a_10k_null():
    null = np.r_[np.full(9_950, 0.52), np.full(49, 0.60), 0.61]
    assert fw_pvalue(0.60, null) == pytest.approx(51 / 10_001)   # ties at 0.60 count
    assert fw_pvalue(0.62, null) == pytest.approx(1 / 10_001)


def test_fw_nulls_large_reproducible_and_saved(tmp_path):
    S, T, E = _surv()
    a = fw_nulls_large({"m": S}, T, E, 50, 7, tmp_path)
    b = fw_nulls_large({"m": S}, T, E, 50, 7, tmp_path)
    np.testing.assert_array_equal(a["nulls"]["m"]["best_cindex_null"], b["nulls"]["m"]["best_cindex_null"])
    z = np.load(tmp_path / "null10k_m.npz")
    assert len(z["best_oriented_null"]) == 50
    assert (z["best_oriented_null"] >= z["best_cindex_null"]).all()
    assert a["check"]["max_abs_diff_vs_lifelines"] < 1e-12


def test_family_wise_large_reports_changed_verdicts():
    fam = pd.DataFrame({"key": ["1/a", "1/b"], "top_cell_type": ["t", "t"], "signature": ["a", "b"],
                        "cindex": [0.60, 0.55], "above_floor": [0.05, 0.02],
                        "cindex_or": [0.60, 0.55], "above_floor_or": [0.05, 0.02]})
    best = fam.iloc[[0]].assign(model="m").set_index("model")
    small = np.r_[np.full(470, 0.5), np.full(30, 0.61)]         # p = 31/501 > 0.05
    large = np.r_[np.full(9_700, 0.5), np.full(300, 0.61)]      # p ~ 0.030
    fw500 = {an: {"m": family_wise(fam, small, 0.60, col, af)}
             for an, col, af in (("primary", "cindex", "above_floor"),
                                 ("sensitivity", "cindex_or", "above_floor_or"))}
    nl = {"nulls": {"m": {"best_cindex_null": large, "best_oriented_null": large}}, "check": {}}
    blk, dev = family_wise_large({"m": fam}, fw500, nl, best, best, 10_000)
    assert blk["primary"]["m"]["P1_corrected"] and not blk["primary"]["m"]["P1_corrected_500"]
    assert len(blk["verdicts_changed"]) == 2 and "these verdicts change" in dev["now_reported"]
    nl["nulls"]["m"] = {"best_cindex_null": np.r_[np.full(9_000, 0.5), np.full(1_000, 0.61)]} | {
        "best_oriented_null": np.r_[np.full(9_000, 0.5), np.full(1_000, 0.61)]}
    blk, dev = family_wise_large({"m": fam}, fw500, nl, best, best, 10_000)
    assert blk["verdicts_changed"] == [] and "no verdict depends on the change" in dev["now_reported"]


def test_kruskal_eta2_formula():
    rng = np.random.default_rng(0)
    v = np.r_[rng.normal(0, 1, 30), rng.normal(1, 1, 30), rng.normal(2, 1, 40)]
    g = np.repeat(["a", "b", "c"], [30, 30, 40])
    r = kruskal_eta2(v, g)
    h, p = kruskal(v[:30], v[30:60], v[60:])
    assert r["kw_h"] == pytest.approx(h) and r["kw_p"] == pytest.approx(p)
    assert r["eta2_h"] == pytest.approx((h - 3 + 1) / (100 - 3))
    assert r["epsilon2"] == pytest.approx(h / 99)
    assert 0 < r["eta2_h"] < 1


def test_auc_orientation_and_ties():
    s = np.array([1.0, 2.0, 3.0, 4.0])
    assert auc_higher_in(s, [False, False, True, True]) == 1.0      # higher in the group
    assert auc_higher_in(s, [True, True, False, False]) == 0.0      # lower in the group
    assert auc_higher_in([1.0, 1.0], [True, False]) == 0.5           # a tie counts one half
    rng = np.random.default_rng(3)
    x, g = np.round(rng.normal(size=80), 1), rng.random(80) < 0.3
    u = mannwhitneyu(x[g], x[~g]).statistic
    assert auc_higher_in(x, g) == pytest.approx(u / (g.sum() * (~g).sum()))


def test_association_tracks_the_furthest_group():
    sub = pd.Series(["Basal"] * 20 + ["LumA"] * 30 + ["LumB"] * 20 + ["Her2"] * 10 + [None] * 5)
    score = pd.Series(np.r_[np.full(20, -2.0), np.linspace(0, 1, 30), np.linspace(0.2, 1.2, 20),
                            np.zeros(10), np.full(5, 9.0)])
    r = association(score, sub)
    assert r["n"] == 80 and r["n_levels"] == 4
    assert r["tracked_group"] == "Basal" and r["auc_tracked"] == 0.0 and not r["higher_in_tracked_group"]
    assert r["auc_tracked_oriented"] == 1.0
    assert r["median_Basal"] == -2.0
    assert "auc_claudin-low" not in r and "auc_Luminal" in r


def _comp_scores():
    comp = pd.DataFrame({
        "model": ["hvg_pca"] * 3 + ["fm"] * 3, "resolution": [0.3] * 6, "cluster": [0, 1, 2] * 2,
        "n_cells": 300, "kept": True, "n_donors": [1, 5, 9, 4, 6, 9],
        "top_donor_frac": [0.95, 0.85, 0.2, 0.9, 0.3, 0.25], "donor_entropy": [0.1, 0.4, 2.0, 0.5, 1.8, 2.1],
        "top_cell_type": "t"})
    rows = []
    for _, r in comp.iterrows():
        for k in (25, 50, 100):
            af = {0: 0.08, 1: 0.05, 2: 0.01}[r["cluster"]] if k == 50 else 0.5
            rows.append({"model": r["model"], "resolution": 0.3, "signature": f"c{r['cluster']}_k{k}",
                         "topk": k, "above_floor": af, "above_floor_or": af, "direction": 1})
    return comp, pd.DataFrame(rows)


def test_donor_mixing_one_row_per_state(tmp_path):
    comp, scores = _comp_scores()
    d = states_with_margins(comp, scores)
    assert len(d) == 6 and set(d["signature"].str[-3:]) == {"k50"}
    comp.to_csv(tmp_path / "c.csv", index=False)
    df, js = donor_mixing(tmp_path / "c.csv", scores)
    allr = df[df["resolution"] == "all"].set_index("model")
    assert allr.loc["hvg_pca", "n_single_donor"] == 2 and allr.loc["fm", "n_single_donor"] == 1
    assert allr.loc["hvg_pca", "best_margin_multi_donor"] == pytest.approx(0.01)
    assert allr.loc["fm", "best_margin_multi_donor"] == pytest.approx(0.05)
    assert allr.loc["hvg_pca", "spearman_rho_margin_vs_top_donor"] == pytest.approx(1.0)
    assert js["checks"]["fm_fewer_single_donor_share_than_baseline"]
    assert js["checks"]["fm_best_multi_donor_margin_above_baseline"]
    assert C.SINGLE_DONOR_FRAC == 0.8


def test_subtype_association_keeps_resolutions_apart():
    from scfm.stages.subtype import subtype_association
    rng = np.random.default_rng(0)
    pats = [f"p{i}" for i in range(40)]
    genes = [f"g{i}" for i in range(8)]
    z = pd.DataFrame(rng.normal(size=(8, 40)), index=genes, columns=pats)
    clin = pd.DataFrame({"subtype": ["Basal", "LumA", "LumB", "Her2"] * 10}, index=pats)
    items = [{"analysis": "primary", "model": "m", "signature": "c1_k50", "signature_key": "0.3/c1_k50",
              "resolution": "0.3", "genes": genes[:4], "direction_label": "risk"},
             {"analysis": "sensitivity", "model": "m", "signature": "c1_k50", "signature_key": "1.0/c1_k50",
              "resolution": "1.0", "genes": genes[4:], "direction_label": "protective"}]
    av = pd.DataFrame({"base": "age+stage+PAM50", "model": "m", "signature": "c1_k50",
                       "resolution": ["0.3", "1.0"], "hr_per_sd": [1.5, 0.5], "lrt_p_nominal": 0.1,
                       "delta_cindex_cv_fixed_mean": 0.0})
    df, js = subtype_association(items, z, clin, av, None)
    assert len(df) == 2 and df["kw_h"].nunique() == 2
    assert df.set_index("resolution")["hr_per_sd_given_pam50"].to_dict() == {"0.3": 1.5, "1.0": 0.5}
    assert js["metabric_note"].startswith("METABRIC not run")
