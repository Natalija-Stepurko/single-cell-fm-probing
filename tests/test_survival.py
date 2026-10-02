import numpy as np
import pandas as pd
import pytest

from scfm.survival import (
    cindex,
    matched_random_sets,
    permute_outcomes,
    score_signature,
    stage_to_ordinal,
    zscore_genes,
)


def harrell_c(time, event, risk):
    """Harrell's C for a risk score: a pair is usable when the earlier time is an event (an event
    and a censoring at the same time count, event first); tied risks score one half."""
    num = den = 0.0
    for i in range(len(time)):
        if not event[i]:
            continue
        for j in range(len(time)):
            if i == j or not (time[i] < time[j] or (time[i] == time[j] and not event[j])):
                continue
            den += 1
            num += 1.0 if risk[i] > risk[j] else 0.5 if risk[i] == risk[j] else 0.0
    return num / den


@pytest.fixture
def cohort():
    return pd.DataFrame({
        "time":  [5, 5, 5, 8, 3, 10, 10, 2, 7, 12, 4, 9],
        "event": [1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0],
        "risk":  [0.9, 0.2, 0.9, 0.5, 0.1, 0.3, 0.3, 0.7, 0.4, 0.05, 0.8, 0.35],
    }, index=[f"p{i}" for i in range(12)])


def test_cindex_matches_hand_written_harrell(cohort):
    expected = harrell_c(cohort["time"].to_numpy(), cohort["event"].to_numpy(), cohort["risk"].to_numpy())
    assert cindex(cohort, "time", "event", "risk") == pytest.approx(expected)


def test_cindex_direction(cohort):
    # at a shared time the event is the earlier failure, so it gets the higher risk
    perfect = cohort.assign(risk=-cohort["time"] + 0.1 * cohort["event"])
    assert cindex(perfect, "time", "event", "risk") == pytest.approx(1.0)
    assert cindex(perfect.assign(risk=-perfect["risk"]), "time", "event", "risk") < 0.1


def test_cindex_needs_five_events(cohort):
    few = cohort.assign(event=[1, 1, 1, 1] + [0] * 8)
    assert np.isnan(cindex(few, "time", "event", "risk"))


def test_cindex_drops_missing_scores(cohort):
    with_nan = cohort.copy()
    with_nan.loc["p9", "risk"] = np.nan
    expected = cindex(cohort.drop(index="p9"), "time", "event", "risk")
    assert cindex(with_nan, "time", "event", "risk") == pytest.approx(expected)


@pytest.fixture
def expr():
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(40)]
    df = pd.DataFrame(rng.normal(5, 2, (40, 15)), index=genes, columns=[f"p{i}" for i in range(15)])
    df.loc["g0"] = 3.0
    return df


def test_zscore_genes(expr):
    z = zscore_genes(expr)
    assert np.allclose(z.drop(index="g0").mean(axis=1), 0)
    assert np.allclose(z.drop(index="g0").std(axis=1), 1)
    assert (z.loc["g0"] == 0).all()


def test_score_signature(expr):
    z = zscore_genes(expr)
    assert score_signature(z, ["g1", "g2", "missing"]).isna().all()
    s = score_signature(z, ["g1", "g2", "g3", "missing"])
    assert np.allclose(s, z.loc[["g1", "g2", "g3"]].mean(axis=0))
    assert list(s.index) == list(z.columns)


def test_matched_random_sets():
    rng = np.random.default_rng(1)
    universe = pd.Series(rng.gamma(2.0, 2.0, 300), index=[f"g{i}" for i in range(300)])
    genes = list(universe.sort_values().index[[0, 1, 150, 151, 152, 299]]) + ["not_in_universe"]
    sets = matched_random_sets(genes, universe, 25, 10, np.random.default_rng(7))
    assert len(sets) == 25

    bins = pd.qcut(universe.rank(method="first"), 10, labels=False)
    want = bins.reindex(genes).dropna().astype(int).value_counts().sort_index()
    for s in sets:
        assert len(s) == 6 and len(set(s)) == 6
        assert bins.reindex(s).value_counts().sort_index().equals(want)

    again = matched_random_sets(genes, universe, 25, 10, np.random.default_rng(7))
    assert again == sets
    assert matched_random_sets(genes, universe, 25, 10, np.random.default_rng(8)) != sets


def test_stage_to_ordinal():
    s = pd.Series(["Stage I", "Stage IIA", "stage iiib", "Stage IV", "Stage X", None, "[Not Available]"])
    out = stage_to_ordinal(s)
    assert out.iloc[:4].tolist() == [1.0, 2.0, 3.0, 4.0]
    assert out.iloc[4:].isna().all()


def test_permute_outcomes_keeps_pairs(cohort):
    out = permute_outcomes(cohort, "time", "event", np.random.default_rng(3))
    pairs = lambda d: sorted(zip(d["time"], d["event"]))
    assert pairs(out) == pairs(cohort)
    assert out["risk"].equals(cohort["risk"])
    assert list(out.index) == list(cohort.index)
    assert not out[["time", "event"]].equals(cohort[["time", "event"]])
