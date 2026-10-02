"""Survival-analysis helpers shared by the data, translate, ladder and report stages.

Everything here operates on pandas frames indexed by patient.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def zscore_genes(expr: pd.DataFrame) -> pd.DataFrame:
    """Per-gene z-score across patients. expr: genes × patients, log-scale."""
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0).fillna(0.0)


def score_signature(zexpr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean z-score of the signature genes present in the matrix, one value per patient.

    The simplest transparent scorer. ssGSEA is not implemented: it would add a
    tunable that the ladder would then have to sweep.
    """
    present = [g for g in genes if g in zexpr.index]
    if len(present) < 3:
        return pd.Series(np.nan, index=zexpr.columns)
    return zexpr.loc[present].mean(axis=0)


def fit_cox(df: pd.DataFrame, duration: str, event: str, covariates: list[str]):
    """Cox proportional-hazards fit. Returns (fitter, summary_row_for_first_covariate)."""
    from lifelines import CoxPHFitter
    cols = [duration, event] + covariates
    d = df[cols].dropna()
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(d, duration_col=duration, event_col=event)
    return cph, cph.summary.loc[covariates[0]]


def cindex(df: pd.DataFrame, duration: str, event: str, score: str) -> float:
    """Harrell's concordance of a risk score against right-censored outcomes."""
    from lifelines.utils import concordance_index
    d = df[[duration, event, score]].dropna()
    if d[event].sum() < 5:
        return float("nan")
    # concordance_index expects higher predicted *survival* for longer times; negate a risk score
    return float(concordance_index(d[duration], -d[score], d[event]))


def permute_outcomes(df: pd.DataFrame, duration: str, event: str, rng: np.random.Generator) -> pd.DataFrame:
    """Permute (time, event) jointly across patients: destroys the score→outcome link only."""
    idx = rng.permutation(len(df))
    out = df.copy()
    out[[duration, event]] = df[[duration, event]].to_numpy()[idx]
    return out


def matched_random_sets(genes: list[str], universe_mean: pd.Series, n_sets: int, n_bins: int,
                        rng: np.random.Generator) -> list[list[str]]:
    """Random gene sets matched to `genes` on size AND on mean-expression bin.

    Prognostic signal correlates with expression level, so a floor drawn from all genes would sit
    too low and flatter every signature. Bin the universe by mean expression, then draw from each
    bin as many genes as the signature has there.
    """
    bins = pd.qcut(universe_mean.rank(method="first"), n_bins, labels=False)
    sig_bins = bins.reindex([g for g in genes if g in bins.index]).dropna().astype(int)
    need = sig_bins.value_counts().to_dict()
    pool = {b: bins.index[bins == b].to_numpy() for b in range(n_bins)}
    out = []
    for _ in range(n_sets):
        draw = []
        for b, k in need.items():
            draw.extend(rng.choice(pool[b], size=min(k, len(pool[b])), replace=False).tolist())
        out.append(draw)
    return out


def bootstrap_patients(df: pd.DataFrame, n_boot: int, rng: np.random.Generator):
    """Yield patient-level bootstrap resamples. The patient is the unit of independence."""
    n = len(df)
    for _ in range(n_boot):
        yield df.iloc[rng.integers(0, n, n)]


def stage_to_ordinal(stage: pd.Series) -> pd.Series:
    """AJCC pathologic stage → 1..4; sub-stages collapse to their stage; unknown → NaN."""
    s = stage.astype(str).str.upper().str.extract(r"STAGE\s*(IV|III|II|I)", expand=False)
    return s.map({"I": 1, "II": 2, "III": 3, "IV": 4}).astype(float)
