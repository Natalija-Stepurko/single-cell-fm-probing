"""Statistics for the corrected ladder: family-wise p-values, orientation, fast Harrell's C,
selection-aware bootstrap margins, within-strata permutations and cross-validated Cox models."""
from __future__ import annotations

import numpy as np
import pandas as pd


def fw_pvalue(stat: float, null) -> float:
    """Permutation p-value with the +1 correction: (1 + #{null >= stat}) / (1 + n_perm)."""
    null = np.asarray(null, dtype=float)
    return float((1 + np.sum(null >= stat)) / (1 + len(null)))


def direction_of(c):
    """+1 when a higher score means worse survival (C >= 0.5), -1 when it means better."""
    return np.where(np.asarray(c, dtype=float) >= 0.5, 1, -1)


def orient(c, d):
    """A C-index read in direction d: C for a risk score (d = +1), 1 - C for a protective one."""
    c = np.asarray(c, dtype=float)
    out = np.where(np.asarray(d) == 1, c, 1 - c)
    return float(out) if out.ndim == 0 else out


def comparable_pairs(T, E):
    """Harrell's usable pairs (i, j), lifelines convention: i has the event and T_j > T_i, or
    T_j == T_i with j censored."""
    T = np.asarray(T, dtype=float)
    E = np.asarray(E).astype(bool)
    di = np.flatnonzero(E)
    m = (T[None, :] > T[di, None]) | ((T[None, :] == T[di, None]) & ~E[None, :])
    ii, jj = np.nonzero(m)
    return di[ii], jj


def cindex_many(S, T, E, chunk: int = 64) -> np.ndarray:
    """Harrell's C of every row of S (risk scores, higher = worse) against one outcome; equal to
    lifelines' concordance_index(T, -s, E) for each row, tied scores counting one half."""
    S = np.atleast_2d(np.asarray(S, dtype=float))
    I, J = comparable_pairs(T, E)
    out = np.empty(len(S))
    for a in range(0, len(S), chunk):
        x, y = S[a:a + chunk][:, I], S[a:a + chunk][:, J]
        out[a:a + chunk] = ((x > y).sum(axis=1) + 0.5 * (x == y).sum(axis=1)) / len(I)
    return out


def best_above_floor(Cb: np.ndarray, floor: np.ndarray):
    """Per bootstrap row of Cb (resamples × signatures): the largest C - floor and its column."""
    af = Cb - floor[None, :]
    j = np.nanargmax(af, axis=1)
    return af[np.arange(len(af)), j], j


def selection_aware_margin(Cb_fm, floor_fm, Cb_base, floor_base):
    """Margin of the FM's best above-floor signature over the baseline's, re-selected in every
    resample; floors stay at their full-cohort values. Returns (margins, fm_picks, base_picks)."""
    af_fm, j_fm = best_above_floor(Cb_fm, floor_fm)
    af_b, j_b = best_above_floor(Cb_base, floor_base)
    return af_fm - af_b, j_fm, j_b


def permute_within_strata(strata, rng: np.random.Generator) -> np.ndarray:
    """An index permutation that only moves patients within their own stratum."""
    strata = np.asarray(strata)
    idx = np.arange(len(strata))
    out = idx.copy()
    for s in pd.unique(strata):
        pos = idx[strata == s]
        out[pos] = pos[rng.permutation(len(pos))]
    return out


def cv_folds(event, n_splits: int, n_repeats: int, rng: np.random.Generator):
    """Repeated event-stratified K-fold splits: a list (one per repeat) of (train, test) pairs."""
    from sklearn.model_selection import StratifiedKFold
    event = np.asarray(event).astype(int)
    reps = []
    for _ in range(n_repeats):
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=int(rng.integers(2**31 - 1)))
        reps.append(list(skf.split(np.zeros(len(event)), event)))
    return reps


def cox_fit(d: pd.DataFrame, duration: str, event: str, cols: list[str], penalizer: float = 0.01):
    from lifelines import CoxPHFitter
    return CoxPHFitter(penalizer=penalizer).fit(d[[duration, event] + cols], duration_col=duration,
                                                event_col=event)


def cv_cox_cindex(d: pd.DataFrame, duration: str, event: str, feature_sets: dict[str, list[str]],
                  folds) -> dict[str, np.ndarray]:
    """Per repeat, Harrell's C of the pooled out-of-fold Cox linear predictor, for each feature set
    on the same folds (so differences between sets are paired)."""
    T, E = d[duration].to_numpy(float), d[event].to_numpy(float)
    out = {k: [] for k in feature_sets}
    for rep in folds:
        for name, cols in feature_sets.items():
            lp = np.empty(len(d))
            for tr, te in rep:
                cph = cox_fit(d.iloc[tr], duration, event, cols)
                lp[te] = np.asarray(cph.predict_log_partial_hazard(d.iloc[te][cols]), dtype=float)
            out[name].append(float(cindex_many(lp[None, :], T, E)[0]))
    return {k: np.array(v) for k, v in out.items()}


def insample_cindex(d: pd.DataFrame, duration: str, event: str, cols: list[str]) -> float:
    cph = cox_fit(d, duration, event, cols)
    lp = np.asarray(cph.predict_log_partial_hazard(d[cols]), dtype=float)
    return float(cindex_many(lp[None, :], d[duration].to_numpy(float), d[event].to_numpy(float))[0])


def lrt_added(d: pd.DataFrame, duration: str, event: str, base: list[str], added: list[str]):
    """Unpenalised nested Cox models: (likelihood-ratio p of base + added vs base, full fitter)."""
    from scipy.stats import chi2
    b = cox_fit(d, duration, event, base, penalizer=0.0)
    f = cox_fit(d, duration, event, base + added, penalizer=0.0)
    stat = 2 * (f.log_likelihood_ - b.log_likelihood_)
    return float(chi2.sf(max(stat, 0.0), df=len(added))), f
