import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

from scfm import config as C
from scfm.stages.ladder import (
    family_wise,
    ladder_by_setting,
    p2_corrected,
    p3_corrected,
    p4_corrected,
    select_in_fold,
)


def _fam(keys, **cols):
    f = pd.DataFrame({"key": keys, "top_cell_type": ["t"] * len(keys), **cols})
    f["signature"] = [k.split("/")[1] for k in keys]
    return f


def test_family_wise_uses_every_tied_maximiser():
    fam = _fam(["0.3/a", "0.5/b", "1.0/c"], cindex=[0.60, 0.60, 0.55], above_floor=[0.05, -0.01, 0.02])
    null = np.r_[np.full(97, 0.52), 0.60, 0.61, 0.62]
    r = family_wise(fam, null, 0.55, "cindex", "above_floor")
    assert r["argmax"] == ["0.3/a", "0.5/b"]
    assert r["argmax_above_floor"] == [0.05, -0.01]
    assert r["fw_p"] == pytest.approx(4 / 101)          # the tie at 0.60 counts
    assert not r["P1_corrected"]                         # one maximiser is below its floor
    fam.loc[1, "above_floor"] = 0.01
    assert family_wise(fam, null, 0.55, "cindex", "above_floor")["P1_corrected"]
    assert r["pick_fw_p"] == pytest.approx(4 / 101)
    assert r["fw_p_mc_se"] == pytest.approx(np.sqrt(r["fw_p"] * (1 - r["fw_p"]) / 100))


def test_p2_corrected_orients_before_reselecting():
    fams = {"fm": _fam(["1/a", "1/b"], direction=[1, -1], floor_or_mean=[0.5, 0.5],
                       above_floor_or=[0.04, 0.10]),
            "base": _fam(["1/x", "1/y"], direction=[1, 1], floor_or_mean=[0.5, 0.45],
                         above_floor_or=[0.02, 0.05])}
    Cb = {"fm": np.array([[0.54, 0.40],      # oriented 0.54, 0.60 -> pick b (+0.10)
                          [0.60, 0.45]]),    # oriented 0.60, 0.55 -> pick a (+0.10)
          "base": np.array([[0.52, 0.50],    # +0.02, +0.05 -> y
                            [0.55, 0.47]])}  # +0.05, +0.02 -> x
    r = p2_corrected(fams, Cb, "fm", "base", "above_floor_or", "floor_or_mean", oriented=True)
    assert r["margin"] == pytest.approx(0.10 - 0.05)
    assert r["pick_reselected_frac"] == {"fm": 0.5, "base": 0.5}
    np.testing.assert_allclose(r["ci"], np.percentile([0.05, 0.05], [2.5, 97.5]))
    assert r["P2_pass"]
    raw = p2_corrected(fams, Cb, "fm", "base", "above_floor_or", "floor_or_mean", oriented=False)
    assert raw["ci"][0] < r["ci"][0]                     # unoriented, b reads 0.40 - 0.5 < 0


def _strata_cohort(rng, n_per=60):
    strata = np.repeat(["A", "B"], n_per)
    T = rng.exponential(size=2 * n_per) * np.where(strata == "A", 1.0, 3.0)
    E = (rng.random(2 * n_per) < 0.6).astype(float)
    clin = pd.DataFrame({"subtype": strata, "OS": E, "OS.time": T},
                        index=[f"p{i}" for i in range(2 * n_per)])
    Z = pd.DataFrame(rng.normal(size=(5, 2 * n_per)), index=[f"g{i}" for i in range(5)], columns=clin.index)
    return clin, Z


def test_p3_corrected_is_calibrated_under_the_null():
    rng = np.random.default_rng(11)
    ps = []
    for _ in range(40):
        clin, Z = _strata_cohort(rng)
        gidx = {g: i for i, g in enumerate(Z.index)}
        out = p3_corrected({"s": {"genes": ["g0", "g1", "g2"], "direction": -1}}, clin, Z, gidx,
                           "OS", "OS.time", 49, int(rng.integers(1 << 30)))
        assert out["subtypes"] == ["A", "B"] and out["n_patients"] == 120
        ps.append(out["picks"]["s"]["p"])
    ps = np.array(ps)
    assert 0.3 < ps.mean() < 0.7
    assert (ps < 0.05).mean() <= 0.15


def test_p3_corrected_detects_a_within_stratum_signal():
    rng = np.random.default_rng(5)
    clin, Z = _strata_cohort(rng, 80)
    Z.loc["g0"] = Z.loc["g1"] = Z.loc["g2"] = -np.log(clin["OS.time"].to_numpy()) + rng.normal(0, 0.3, 160)
    gidx = {g: i for i, g in enumerate(Z.index)}
    out = p3_corrected({"s": {"genes": ["g0", "g1", "g2"], "direction": 1}}, clin, Z, gidx,
                       "OS", "OS.time", 99, 3)
    assert out["picks"]["s"]["p"] == pytest.approx(1 / 100)
    assert not out["picks"]["s"]["borderline"]


def test_p4_corrected_hypergeometric_parameterisation():
    fams = {"fm": _fam(["1/a"], hvg_frac=[0.1]), "base": _fam([f"1/b{i}" for i in range(9)],
                                                              hvg_frac=np.linspace(0.2, 1.0, 9))}
    pick = pd.Series({"hvg_frac": 0.1, "topk": 50})
    uni = {"universe_n": 8000, "universe_n_hvg": 1500}
    r = p4_corrected(fams, pick, "fm", "base", True, uni)
    assert r["hvg_in_pick"] == 5
    assert r["hypergeom_depletion_p"] == pytest.approx(hypergeom.cdf(5, 8000, 1500, 50))
    assert r["expected"] == pytest.approx(50 * 1500 / 8000)
    assert r["baseline_rank_p"] == pytest.approx(1 / 10)
    assert r["P4_corrected"] is False                    # rank p 0.1 is not below alpha
    assert p4_corrected(fams, pick, "fm", "base", False, uni)["P4_corrected"] is None


def test_select_in_fold_rules():
    c = np.array([0.56, 0.40, 0.58])
    floor = np.array([0.48, 0.50, 0.55])
    assert select_in_fold(c, floor, "primary") == 0            # +0.08, -0.10, +0.03
    assert select_in_fold(c, floor, "sensitivity") == 1        # 0.60 - 0.50 = +0.10
    assert select_in_fold(c, floor, "family_max_primary") == 2
    assert select_in_fold(c, floor, "family_max_sensitivity") == 1
    with pytest.raises(ValueError):
        select_in_fold(c, floor, "reference")


def test_ladder_by_setting_names_the_setting_maximum():
    def fam(m, cs):
        f = _fam([f"0.3/{m}{i}" for i in range(3)], cindex=cs, floor_mean=[0.5, 0.45, 0.5],
                 direction=[1 if x >= 0.5 else -1 for x in cs], floor_or_mean=[0.5, 0.45, 0.5])
        f["above_floor"] = f["cindex"] - f["floor_mean"]
        f["cindex_or"] = np.maximum(f["cindex"], 1 - f["cindex"])
        f["above_floor_or"] = f["cindex_or"] - f["floor_or_mean"]
        f["resolution"], f["topk"] = 0.3, 50
        return f
    fams = {"hvg_pca": fam("h", [0.55, 0.53, 0.40]), "scgpt": fam("s", [0.58, 0.52, 0.51])}
    null = {m: {"cindex": np.full((9, 3), 0.5), "index": {k: i for i, k in enumerate(f["key"])}}
            for m, f in fams.items()}
    df = ladder_by_setting(fams, null, "hvg_pca").set_index("model")
    assert df.loc["hvg_pca", "signature"] == "h1"              # above-floor pick: 0.53 - 0.45
    assert df.loc["hvg_pca", "setting_max_signature"] == "h0"  # largest C
    assert df.loc["hvg_pca", "setting_max_signature_or"] == "h2"
    assert df.loc["hvg_pca", "setting_fw_p"] == pytest.approx(1 / 10)
    assert df.loc["scgpt", "above_floor_minus_baseline"] == pytest.approx(0.08 - 0.08)
    assert C.ALPHA == 0.05
