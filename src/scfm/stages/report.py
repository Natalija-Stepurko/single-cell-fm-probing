"""Stage report — figures, figures.json, the candidate shortlist and the validation template.

Every number drawn is read from results/ (ladder, translate, states, stratify, replicate); the only
quantities computed here are descriptive (Kaplan-Meier curves, numbers at risk, UMAP highlights).
The family-wise p-values and null percentiles drawn are the 10,000-permutation ones (ladder.json
family_wise_10k).

  fig_ladder.png       C of each representation's pick against its floor, its family-wise null and
                       three references; pre-specified (risk direction) and sensitivity (either)
  fig_states.png       where the picks sit in each representation's UMAP, with their programme
  fig_added_value.png  what each pick adds to age + stage (and to age + stage + PAM50)
  fig_stratify.png     cross-validated C of multivariable models with and without the states
  fig_replication.png  the frozen signatures in TCGA and in METABRIC (OS, DSS, within subtype)
  fig_km.png           Kaplan-Meier by score tertile for the lead replicated signature
  fig_margin.png       foundation model minus baseline (P2, floor-adjusted), TCGA and METABRIC
  fig_subtype.png      each pick's score by PAM50 subtype, TCGA and METABRIC (outcome-free)
  fig_donor.png        donor mixing of every kept state, and its above-floor margin against it
  figures.json         what each figure shows, how to read it, and the files it read

shortlist.md           FM signatures above their floor and the family-wise null
validation.md          the three-step validation each shortlisted state would need
atlas_duplicates.json  atlas cells whose raw counts duplicate another atlas cell (needs data/atlas.h5ad)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scfm import config as C
from scfm import figstyle as S
from scfm.provenance import record_params, relpath
from scfm.survival import score_signature, zscore_genes

FIGURES = ["fig_ladder", "fig_states", "fig_added_value", "fig_stratify", "fig_replication", "fig_km",
           "fig_margin", "fig_subtype", "fig_donor"]
ANALYSIS_NAME = {"primary": "Pre-specified", "sensitivity": "Sensitivity (post hoc)"}

# Programme of each pick, from its marker genes, donor and dataset composition. Each entry names the
# pick it describes and genes that must be in its signature, so a change in the picks fails loudly.
PROGRAMMES = {
    ("primary", "hvg_pca"): {"pick": "0.3/c21_k50", "label": "single-dataset artefact", "genes": []},
    ("primary", "scgpt"): {"pick": "1.0/c11_k50", "label": "luminal epithelium (ERBB2+)",
                           "genes": ["ERBB2", "EPCAM", "KRT19", "FOXA1"]},
    ("primary", "geneformer"): {"pick": "0.3/c14_k100", "label": "basal keratin programme",
                                "genes": ["KRT5", "KRT14", "KRT17"]},
    ("sensitivity", "hvg_pca"): {"pick": "1.0/c17_k100", "label": "ER+/PR+ luminal programme",
                                 "genes": ["PGR", "GATA3", "XBP1"]},
    ("sensitivity", "scgpt"): {"pick": "1.0/c8_k50", "label": "luminal ER programme",
                               "genes": ["GATA3", "XBP1", "TFF1"]},
    ("sensitivity", "geneformer"): {"pick": "1.0/c14_k25", "label": "nuclear-RNA / stress",
                                    "genes": ["MALAT1", "NEAT1", "JUND"]},
}


# ---- inputs -----------------------------------------------------------------------------------
def load(args) -> dict:
    R = {k: Path(getattr(args, f"{k}_dir")) for k in ("states", "translate", "ladder", "stratify",
                                                        "replicate", "data")}
    files = {
        "ladder": R["ladder"] / "ladder.json", "added_value": R["ladder"] / "added_value.csv",
        "scores": R["translate"] / "scores.csv", "nulls": R["translate"] / "nulls.json",
        "signatures": R["states"] / "signatures.json", "composition": R["states"] / "state_composition.csv",
        "stratify": R["stratify"] / "cv_summary.json", "stratify_rows": R["stratify"] / "cv_results.csv",
        "metabric": R["replicate"] / "metabric.json", "frozen": R["replicate"] / "frozen_signatures.json",
        "clinical": R["data"] / "bulk_clinical.csv", "expression": R["data"] / "bulk_expr.parquet",
        "subtype": R["ladder"] / "subtype_association.csv", "donor_mixing": R["ladder"] / "donor_mixing.json",
    }
    for m in S.MODELS:
        files[f"cells_{m}"] = R["states"] / f"cells_{m}.parquet"
    D = {"files": files}
    D["L"] = json.load(open(files["ladder"]))
    D["scores"] = pd.read_csv(files["scores"])
    D["nulls"] = json.load(open(files["nulls"]))
    D["sigs"] = json.load(open(files["signatures"]))
    D["comp"] = pd.read_csv(files["composition"], dtype={"cluster": str, "resolution": str})
    D["added"] = pd.read_csv(files["added_value"])
    D["strat"] = json.load(open(files["stratify"]))
    D["metabric"] = json.load(open(files["metabric"]))
    D["frozen"] = json.load(open(files["frozen"]))
    D["subtype"] = pd.read_csv(files["subtype"])
    D["donor_mixing"] = json.load(open(files["donor_mixing"]))
    return D


def picks(L: dict, analysis: str) -> pd.DataFrame:
    rows = L["ladder"] if analysis == "primary" else L["sensitivity"]["ladder"]
    df = pd.DataFrame(rows).set_index("model").loc[list(S.MODELS)]
    df["key"] = [f"{r:.1f}/{s}" for r, s in zip(df["resolution"], df["signature"])]
    return df


def programme(analysis: str, model: str, key: str, sigs: dict) -> str:
    p = PROGRAMMES[(analysis, model)]
    if p["pick"] != key:
        raise ValueError(f"{analysis}/{model}: pick is {key}, the programme label describes {p['pick']}")
    res, sig = key.split("/")
    missing = [g for g in p["genes"] if g not in sigs[model][res][sig]["genes"]]
    if missing:
        raise ValueError(f"{analysis}/{model} {key}: {missing} not in the signature")
    return p["label"]


def fmt_p(p: float) -> str:
    if p >= 0.1:
        return f"{p:.2f}"
    return f"{p:.2g}" if p >= 0.00095 else f"{p:.0e}".replace("e-0", "e-")


def signed(v, _=None) -> str:
    return "0" if abs(v) < 1e-9 else f"{v:+.2f}".replace("-", "\u2212")


def row_labels(ax, labels, ys, x=-0.02, size=S.FS):
    for y, t in zip(ys, labels):
        ax.text(x, y, t, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=size)


def group_header(ax, y, text, x=-0.02):
    ax.text(x, y, text, transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=S.FS_SMALL, color=S.MUTED, fontweight="bold")


# ---- 1. the ladder ----------------------------------------------------------------------------
def fig_ladder(D, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    L = D["L"]; ref = L["reference"]
    refs = [(ref["pam50"]["cindex_cv_mean"], "PAM50 subtype\nalone (CV)", S.COL["clinical"]),
            (ref["proliferation"]["cindex"], "Proliferation\nscore", S.COL["reference"]),
            (ref["clinical_full"]["cindex_cv_mean"], "Age + stage\n(CV)", S.COL["clinical"])]
    cols = [(1.02, "C of\npick"), (1.10, "pick's\np"), (1.205, "family\nmax ◇"), (1.285, "family-\nwise p")]
    mono = {"va": "center", "fontsize": S.FS, "family": "DejaVu Sans Mono"}
    fig, axes = plt.subplots(2, 1, figsize=(S.WIDE, 5.6), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.76, top=0.83, bottom=0.2, hspace=0.45)
    for ax, analysis in zip(axes, ("primary", "sensitivity")):
        P = picks(L, analysis); fw = L["family_wise_10k"][analysis]
        for i, m in enumerate(S.MODELS):
            r = P.loc[m]; y = -i
            ax.barh(y, r["cindex"] - 0.5, left=0.5, height=0.5, color=S.COL[m], zorder=2)
            ax.plot([r["ci_lo"], r["ci_hi"]], [y, y], color=S.INK, lw=1.1, zorder=3)
            for e in (r["ci_lo"], r["ci_hi"]):
                ax.plot([e, e], [y - 0.09, y + 0.09], color=S.INK, lw=1.1, zorder=3)
            ax.plot([r["floor_mean"]] * 2, [y - 0.36, y + 0.36], color=S.INK, lw=1.3, ls=(0, (1, 1.6)),
                    zorder=4)
            ax.plot([fw[m]["null_p95"]] * 2, [y - 0.36, y + 0.36], color=S.INK, lw=3.2, zorder=4,
                    solid_capstyle="butt")
            ax.plot(fw[m]["stat"], y, marker="D", ms=6.5, mfc=S.SURFACE, mec=S.INK, mew=1.2, zorder=5)
            vals = [f"{r['cindex']:.3f}", fmt_p(fw[m]["pick_fw_p"]), f"{fw[m]['stat']:.3f}",
                    fmt_p(fw[m]["fw_p"])]
            bold = [False, fw[m]["pick_fw_p"] < C.ALPHA, False, fw[m]["fw_p"] < C.ALPHA]
            for (x0, _), v, b in zip(cols, vals, bold):
                ax.text(x0, y, v, transform=ax.get_yaxis_transform(), fontweight="bold" if b else "normal",
                        **mono)
        ax.set_yticks([0, -1, -2]); ax.set_yticklabels([S.NAME[m] for m in S.MODELS])
        ax.set_ylim(-2.6, 0.6)
        for x, _, c in refs:
            ax.axvline(x, color=c, lw=1.0, ls=(0, (4, 3)), zorder=1)
        ax.axvline(0.5, color=S.MUTED, lw=0.8, zorder=1)
        S.clean(ax)
        title = ("Pre-specified: risk direction" if analysis == "primary"
                 else "Sensitivity (post hoc): either direction, C read in the direction the score acts")
        ax.set_title(title, fontsize=S.FS, pad=34 if analysis == "primary" else 8)
        for x0, t in cols:
            ax.text(x0, 0.62, t, transform=ax.get_yaxis_transform(), fontsize=S.FS_SMALL, color=S.MUTED,
                    va="bottom", linespacing=1.15)
    axes[1].set_xlim(0.44, 0.78)
    axes[1].set_xlabel("Concordance index (Harrell's C) of the signature score alone, TCGA-BRCA overall "
                       f"survival ({L['n_patients']:,} patients, {L['n_events']} deaths)")
    for x, t, c in refs:
        axes[0].text(x, 1.02, t + f" {x:.3f}", transform=axes[0].get_xaxis_transform(), ha="center",
                     va="bottom", fontsize=S.FS_SMALL, color=S.MUTED if c == S.COL["clinical"] else c,
                     linespacing=1.15)
    handles = [Patch(facecolor=S.SURFACE, edgecolor=S.MUTED, label="C of the pick (largest margin over its "
                     "floor; bar in the representation's colour)"),
               Line2D([], [], color=S.INK, lw=1.1, marker="|", ms=7,
                      label="95% bootstrap interval, pick held fixed (conditional)"),
               Line2D([], [], color=S.INK, lw=1.3, ls=(0, (1, 1.6)),
                      label="matched-random floor (mean of 200 random gene sets)"),
               Line2D([], [], color=S.INK, lw=3.2,
                      label=f"family-wise null, 95th percentile "
                            f"({L['family_wise_10k']['n_perm']:,} permutations)"),
               Line2D([], [], color=S.INK, lw=0, marker="D", ms=6.5, mfc=S.SURFACE, mew=1.2,
                      label="largest C in the family (family max): the statistic the family-wise p tests"),
               Line2D([], [], color=S.COL["clinical"], lw=1.0, ls=(0, (4, 3)),
                      label="clinical reference (cross-validated)"),
               Line2D([], [], color=S.COL["reference"], lw=1.0, ls=(0, (4, 3)),
                      label="published proliferation score"),
               Line2D([], [], lw=0, label="pick's p: the pick read against the same family-wise null")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=2, fontsize=S.FS_SMALL,
               handlelength=2.4, columnspacing=1.6, borderaxespad=0.4)
    return S.save(fig, out, "fig_ladder")


# ---- 2. the states ----------------------------------------------------------------------------
def fig_states(D, out):
    import matplotlib.pyplot as plt
    files = D["files"]
    if not all(files[f"cells_{m}"].exists() for m in S.MODELS):
        print("  fig_states skipped: results/states/cells_<model>.parquet not present (stage states)")
        return None
    L = D["L"]; comp = D["comp"].set_index(["model", "resolution", "cluster"])
    fig, axes = plt.subplots(1, 3, figsize=(S.WIDE, 4.3))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.01, wspace=0.04)
    corners = {"tl": (0.03, 0.97, "left", "top"), "tr": (0.97, 0.97, "right", "top"),
               "bl": (0.03, 0.03, "left", "bottom"), "br": (0.97, 0.03, "right", "bottom")}
    for ax, m in zip(axes, S.MODELS):
        cells = pd.read_parquet(files[f"cells_{m}"])
        xy = cells[["umap_1", "umap_2"]].to_numpy(float)
        ax.scatter(xy[:, 0], xy[:, 1], s=0.35, c=S.CELLS, lw=0, rasterized=True)
        lo, hi = xy.min(0), xy.max(0); span = hi - lo
        ax.set_xlim(lo[0] - 0.42 * span[0], hi[0] + 0.42 * span[0])
        ax.set_ylim(lo[1] - 0.30 * span[1], hi[1] + 0.30 * span[1])
        used = set()
        for analysis, col in (("primary", S.COL[m]), ("sensitivity", S.TINT[m])):
            r = picks(L, analysis).loc[m]
            res, sig = r["key"].split("/"); cl = D["sigs"][m][res][sig]["cluster"]
            sel = (cells[f"leiden_{res}"] == cl).to_numpy()
            if sel.sum() != D["sigs"][m][res][sig]["n_cells"]:
                raise ValueError(f"{m} {r['key']}: {sel.sum()} cells in the UMAP file, "
                                 f"{D['sigs'][m][res][sig]['n_cells']} in signatures.json")
            ax.scatter(xy[sel, 0], xy[sel, 1], s=1.6, c=col, lw=0, rasterized=True, zorder=3)
            c = comp.loc[(m, res, cl)]
            cen = np.median(xy[sel], axis=0)
            axc = ax.transAxes.inverted().transform(ax.transData.transform(cen))
            order = sorted(corners, key=lambda k: np.hypot(corners[k][0] - axc[0], corners[k][1] - axc[1]))
            k = next(k for k in order if k not in used and k[0] not in {u[0] for u in used}); used.add(k)
            x, y, ha, va = corners[k]
            nd, nds = int(c["n_donors"]), int(c["n_datasets"])
            donors = f"{nd} donor" if nd == 1 else f"{nd} donors ({c['top_donor_frac']:.0%} from one)"
            ds = "1 dataset" if nds == 1 else f"{nds} datasets ({c['top_dataset_frac']:.0%} from one)"
            text = (f"{ANALYSIS_NAME[analysis].split(' ')[0]} pick\n"
                    f"{programme(analysis, m, r['key'], D['sigs'])}\n"
                    f"{int(c['n_cells']):,} cells from {donors}\n{ds}")
            ax.annotate(text, xy=cen, xytext=(x, y), textcoords="axes fraction", ha=ha, va=va,
                        fontsize=S.FS_SMALL, linespacing=1.35, color=S.INK,
                        arrowprops={"arrowstyle": "-", "color": S.INK, "lw": 0.7, "shrinkA": 2,
                                    "shrinkB": 3})
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(f"{S.NAME[m]}   ·   {len(cells):,} cells", fontsize=S.FS, pad=4)
    return S.save(fig, out, "fig_states")


# ---- 3. added value ---------------------------------------------------------------------------
def _rows_added(D):
    A = D["added"]; L = D["L"]; rows = []
    for analysis in ("primary", "sensitivity"):
        P = picks(L, analysis)
        for m in S.MODELS:
            sub = A[(A["analysis"] == analysis) & (A["model"] == m)]
            if sub["signature"].nunique() != 1 or sub["signature"].iloc[0] != P.loc[m, "signature"]:
                raise ValueError(f"added_value.csv: {analysis}/{m} does not match the ladder pick")
            rows.append((analysis, m, f"{S.NAME[m]} · {programme(analysis, m, P.loc[m, 'key'], D['sigs'])}",
                         sub.set_index("base")))
    sub = A[A["analysis"] == "reference"]
    rows.append(("reference", "reference", "Proliferation score (11 genes)", sub.set_index("base")))
    return rows


def _layout(groups):
    """y positions top to bottom with a gap and a header slot before each group."""
    y, ys, heads = 0.0, [], []
    for name, n in groups:
        heads.append((y, name)); y -= 0.85
        for _ in range(n):
            ys.append(y); y -= 1.0
        y -= 0.25
    return ys, heads, y


def fig_added_value(D, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rows = _rows_added(D)
    ys, heads, ybot = _layout([("PRE-SPECIFIED PICKS", 3), ("SENSITIVITY PICKS (POST HOC)", 3),
                               ("PUBLISHED REFERENCE", 1)])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(S.WIDE, 5.0), sharey=True)
    fig.subplots_adjust(left=0.30, right=0.975, top=0.93, bottom=0.22, wspace=0.08)
    bases = [("age+stage", "o", 0.17), ("age+stage+PAM50", "o", -0.17)]
    for (_, m, _, sub), y in zip(rows, ys):
        col = S.COL[m]
        for base, mk, dy in bases:
            r = sub.loc[base]; filled = base == "age+stage"
            kw = {"marker": mk, "ms": 6.5, "mec": col, "mfc": col if filled else S.SURFACE, "mew": 1.4,
                  "zorder": 4}
            a1.plot([r["hr_ci_lo"], r["hr_ci_hi"]], [y + dy] * 2, color=col, lw=1.6, zorder=3)
            a1.plot(r["hr_per_sd"], y + dy, **kw)
            a2.plot([r["delta_cindex_cv_nested_lo"], r["delta_cindex_cv_nested_hi"]], [y + dy] * 2,
                    color=col, lw=1.6, zorder=3)
            a2.plot(r["delta_cindex_cv_nested_mean"], y + dy, **kw)
    a1.set_xscale("log")
    ticks = [0.5, 0.6, 0.8, 1.0, 1.25, 1.5]
    a1.set_xticks(ticks); a1.set_xticklabels([f"{t:g}" for t in ticks]); a1.minorticks_off()
    lo = min(r[3]["hr_ci_lo"].min() for r in rows); hi = max(r[3]["hr_ci_hi"].max() for r in rows)
    a1.set_xlim(lo / 1.08, hi * 1.08)
    a1.axvline(1.0, color=S.INK, lw=0.9, zorder=1); a2.axvline(0.0, color=S.INK, lw=0.9, zorder=1)
    lo = min(r[3]["delta_cindex_cv_nested_lo"].min() for r in rows)
    hi = max(r[3]["delta_cindex_cv_nested_hi"].max() for r in rows)
    a2.set_xlim(np.floor(lo * 100) / 100 - 0.002, np.ceil(hi * 100) / 100 + 0.002)
    a2.xaxis.set_major_locator(plt.MultipleLocator(0.01))
    a2.xaxis.set_major_formatter(plt.FuncFormatter(signed))
    for ax in (a1, a2):
        S.clean(ax); S.xgrid(ax); ax.set_yticks([])
    a1.set_ylim(ybot + 0.2, 0.5)
    row_labels(a1, [r[2] for r in rows], ys)
    for y, h in heads:
        group_header(a1, y, h)
    a1.set_title("HR per SD of the score", fontsize=S.FS)
    a2.set_title("Change in cross-validated C (nested)", fontsize=S.FS)
    a1.set_xlabel("Hazard ratio (HR) per standard deviation (SD),\nlog scale (below 1: protective)")
    a2.set_xlabel("ΔC over the same model without the score")
    n1 = rows[0][3].loc["age+stage"]; n2 = rows[0][3].loc["age+stage+PAM50"]
    handles = [Line2D([], [], color=S.MUTED, marker="o", ms=6.5, mfc=S.MUTED, lw=1.6,
                      label=f"added to age + stage ({int(n1['n']):,} patients, {int(n1['events'])} deaths)"),
               Line2D([], [], color=S.MUTED, marker="o", ms=6.5, mfc=S.SURFACE, mew=1.4, lw=1.6,
                      label=f"added to age + stage + PAM50 subtype ({int(n2['n']):,} patients, "
                            f"{int(n2['events'])} deaths)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=1, fontsize=S.FS_SMALL,
               borderaxespad=0.4)
    fig.text(0.715, 0.012, "HR: Cox model with the covariates, Wald 95%\ninterval, nominal (the picks were "
             "selected on\nthese outcomes). ΔC: the signature is re-selected\nin every training fold; "
             "whiskers span 5 repeats.", fontsize=S.FS_SMALL, color=S.MUTED, va="bottom", linespacing=1.3)
    return S.save(fig, out, "fig_added_value")


# ---- 4. stratification ------------------------------------------------------------------------
def fig_stratify(D, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    M = D["strat"]["models"]; F = D["strat"]["features"]
    spec = [("clinical", "Age + stage (Cox)", "clinical", "o"),
            ("clinical_prolif", "+ proliferation score (Cox)", "reference", "o")]
    spec += [(f"ridge_{m}", f"+ {F[m]} {S.NAME[m]} states (ridge Cox)", m, "o") for m in S.MODELS]
    spec += [(f"xgb_{m}", f"+ {F[m]} {S.NAME[m]} states (XGBoost)", m, "D") for m in S.MODELS]
    spec2 = [("clinical_on_pam50_set", "Age + stage (Cox)", "clinical", "o"),
             ("clinical_pam50", "+ PAM50 subtype (Cox)", "clinical", "o")]
    g1 = M["clinical"]; g2 = M["clinical_on_pam50_set"]
    ys, heads, ybot = _layout([(f"{g1['n']:,} PATIENTS WITH AGE, STAGE ({g1['events']} DEATHS)", len(spec)),
                               (f"{g2['n']:,} WITH A PAM50 CALL ({g2['events']} DEATHS)", len(spec2))])
    fig, ax = plt.subplots(figsize=(S.WIDE, 5.1))
    fig.subplots_adjust(left=0.36, right=0.9, top=0.95, bottom=0.2)
    for (key, _, colkey, mk), y in zip(spec + spec2, ys):
        r = M[key]; col = S.COL[colkey]
        v = np.asarray(r["cindex_per_repeat"], float)
        ax.plot([v.min(), v.max()], [y, y], color=col, lw=1.2, alpha=0.6, zorder=2)
        ax.scatter(v, np.full(len(v), y), s=16, color=col, alpha=0.55, lw=0, zorder=3, marker=mk)
        ax.plot([r["cindex_mean"]] * 2, [y - 0.32, y + 0.32], color=col, lw=2.6, zorder=4,
                solid_capstyle="butt")
        ax.text(1.01, y, f"{r['cindex_mean']:.3f}", transform=ax.get_yaxis_transform(), va="center",
                fontsize=S.FS, family="DejaVu Sans Mono")
    for (key, *_), y0 in ((spec[0], ys[0]), (spec2[0], ys[len(spec)])):
        n = len(spec) if key == "clinical" else len(spec2)
        yy = ys[ys.index(y0) + n - 1]
        ax.plot([M[key]["cindex_mean"]] * 2, [y0 + 0.45, yy - 0.45], color=S.MUTED, lw=0.8,
                ls=(0, (4, 3)), zorder=1)
    row_labels(ax, [s[1] for s in spec + spec2], ys)
    for y, h in heads:
        group_header(ax, y, h)
    ax.text(1.01, 0.3, "mean", transform=ax.get_yaxis_transform(), fontsize=S.FS_SMALL, color=S.MUTED)
    ax.set_yticks([]); ax.set_ylim(ybot + 0.2, 0.6)
    S.clean(ax); S.xgrid(ax)
    lo = min(min(M[k]["cindex_per_repeat"]) for k, *_ in spec + spec2)
    hi = max(max(M[k]["cindex_per_repeat"]) for k, *_ in spec + spec2)
    ax.set_xlim(np.floor(lo * 100) / 100 - 0.005, np.ceil(hi * 100) / 100 + 0.005)
    ax.set_xlabel("Cross-validated C, overall survival (5 × 5-fold; pooled out-of-fold risk)")
    handles = [Line2D([], [], color=S.MUTED, lw=0, marker="o", ms=4.5, alpha=0.6,
                      label="one cross-validation repeat"),
               Line2D([], [], color=S.MUTED, lw=2.6, label="mean of 5 repeats"),
               Line2D([], [], color=S.MUTED, lw=0.8, ls=(0, (4, 3)), label="age + stage, same patients"),
               Line2D([], [], color=S.MUTED, lw=0, marker="D", ms=4.5, label="XGBoost (survival:cox)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.36, 0.0), ncol=2, fontsize=S.FS_SMALL,
               borderaxespad=0.4, handlelength=2.0, columnspacing=1.5)
    fig.text(0.01, 0.02, "Models within a patient set share\none set of folds. Ridge penalty on\nthe "
             "state features only, chosen\nby inner cross-validation.", fontsize=S.FS_SMALL, color=S.MUTED,
             va="bottom", linespacing=1.3)
    return S.save(fig, out, "fig_stratify")


# ---- 5. replication ---------------------------------------------------------------------------
def _tcga_oriented(sig, D):
    """(oriented C, oriented floor p95) of a frozen signature in TCGA, read from translate."""
    if sig["model"] == "reference":
        r = D["nulls"]["reference"]["proliferation"]
        return r["cindex_or"], r["floor_or_p95"]
    sc = D["scores"]
    key = sig["tcga_keys"][0]; res, s = key.split("/")
    row = sc[(sc["model"] == sig["model"]) & (sc["resolution"].map(lambda v: f"{v:.1f}") == res)
             & (sc["signature"] == s)].iloc[0]
    c = row["cindex"] if sig["direction"] == 1 else 1 - row["cindex"]
    if int(row["direction"]) != sig["direction"]:
        raise ValueError(f"{sig['id']}: frozen direction differs from its TCGA direction")
    return float(c), float(row["floor_or_p95"])


def _replication_rows(D):
    frozen = {s["id"]: s for s in D["frozen"]["signatures"]}
    met = {s["id"]: s for s in D["metabric"]["signatures"]}
    L = D["L"]; rows = []
    groups = [("PRE-SPECIFIED PICKS", "primary"), ("SENSITIVITY PICKS (POST HOC)", "sensitivity"),
              ("FAMILY MAXIMA, PRE-SPECIFIED", "family_max_primary"),
              ("FAMILY MAXIMA, SENSITIVITY", "family_max_sensitivity"),
              ("PUBLISHED REFERENCE", "reference")]
    out = []
    for head, analysis in groups:
        ids = [i for i, s in frozen.items() if analysis in s["analyses"]]
        if analysis.startswith("family_max"):
            seen = {r[1] for r in rows}
            same = [S.NAME[frozen[i]["model"]] for i in ids if i in seen]
            ids = [i for i in ids if i not in seen]
            if same:
                head += f"\n({' and '.join(same)}: maximum = pick)"
        ids.sort(key=lambda i: (list(S.MODELS) + ["reference"]).index(frozen[i]["model"]))
        if not ids:
            continue
        g = []
        for i in ids:
            s = frozen[i]; m = s["model"]
            if analysis in ("primary", "sensitivity"):
                key = picks(L, analysis).loc[m, "key"]
                label = f"{S.NAME[m]} · {programme(analysis, m, key, D['sigs'])}"
            elif m == "reference":
                label = "Proliferation score (11 genes)"
            else:
                label = f"{S.NAME[m]} · max {s['resolution']}/{s['signature']}"
            rows.append((analysis, i)); g.append((i, label, s, met[i]))
        out.append((head, g))
    return out


def fig_replication(D, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    groups = _replication_rows(D)
    ys, heads, ybot = _layout([(h, len(g)) for h, g in groups])
    items = [it for _, g in groups for it in g]
    fig, axes = plt.subplots(1, 4, figsize=(S.WIDE, 6.9), sharey=True,
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.05]})
    fig.subplots_adjust(left=0.265, right=0.985, top=0.9, bottom=0.2, wspace=0.16)
    eps = D["metabric"]["endpoints"]
    titles = [f"TCGA-BRCA OS\ndiscovery, {D['L']['n_events']} deaths",
              f"METABRIC OS\n{eps['OS']['events']:,} deaths",
              f"METABRIC DSS\n{eps['DSS']['events']:,} deaths",
              "METABRIC OS, within the\nlead subtype (no floor)"]
    for (_, _, s, mr), y in zip(items, ys):
        col = S.COL[s["model"]]
        c, f95 = _tcga_oriented(s, D)
        cols = [(c, f95, True)]
        for ep in ("OS", "DSS"):
            e = mr["endpoints"][ep]
            cols.append((e["cindex_oriented"], e["floor_p95"], e["replicates"]))
        for ax, (cv, fl, ok) in zip(axes[:3], cols):
            ax.plot([fl, fl], [y - 0.34, y + 0.34], color=S.INK, lw=2.4, solid_capstyle="butt", zorder=2)
            ax.plot(cv, y, "o", ms=7, mec=col, mfc=col if ok else S.SURFACE, mew=1.6, zorder=4)
        for st, v in mr.get("p3", {}).items():
            e = v["OS"]; ok = e["perm_p"] < C.ALPHA
            axes[3].plot(e["cindex_oriented"], y, "o", ms=7, mec=col, mfc=col if ok else S.SURFACE, mew=1.6,
                         zorder=4)
            axes[3].text(e["cindex_oriented"] + 0.014, y, f"{st}, p {fmt_p(e['perm_p'])}", va="center",
                         fontsize=S.FS_SMALL)
    for ax, t in zip(axes, titles):
        S.clean(ax); S.xgrid(ax); ax.set_title(t, fontsize=S.FS_SMALL, fontweight="bold", linespacing=1.3)
        ax.axvline(0.5, color=S.MUTED, lw=0.8, zorder=1)
    axes[0].set_yticks([]); axes[0].set_ylim(ybot + 0.2, 0.5)
    for ax in axes:
        ax.set_xlim(0.47, 0.69); ax.set_xticks([0.5, 0.55, 0.6, 0.65])
    row_labels(axes[0], [it[1] for it in items], ys)
    for y, h in heads:
        group_header(axes[0], y, h)
    fig.supxlabel(f"C read in the direction fixed in TCGA (oriented C); METABRIC {eps['OS']['n']:,} patients",
                  fontsize=S.FS, x=0.625, y=0.13)
    handles = [Line2D([], [], color=S.MUTED, lw=0, marker="o", ms=7, mfc=S.MUTED,
                      label="replicates: C above the METABRIC floor p95, and age-adjusted, cohort-stratified "
                            "HR in the TCGA direction with p < 0.05 (TCGA column: discovery, always filled)"),
               Line2D([], [], color=S.MUTED, lw=0, marker="o", ms=7, mfc=S.SURFACE, mew=1.6,
                      label="does not replicate. Subtype panel: filled when the outcome-permutation "
                            "p < 0.05; no random-gene floor there, so not a replication test"),
               Line2D([], [], color=S.INK, lw=2.4, label="floor p95: 95th percentile of the matched-random "
                      "floor (200 random gene sets of the same size and expression, in that cohort)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 0.0), ncol=1, fontsize=S.FS_SMALL,
               borderaxespad=0.3)
    return S.save(fig, out, "fig_replication")


# ---- 6. Kaplan-Meier --------------------------------------------------------------------------
def km_choice(D) -> dict:
    """The lead replicated signature: a ladder pick that replicates on METABRIC OS, preferring one that
    passed P3 (within-subtype permutation test) in TCGA, then the larger METABRIC OS margin over its floor."""
    met = {s["id"]: s for s in D["metabric"]["signatures"]}
    p3 = D["L"]["p3_corrected"]["picks"]
    cand = []
    for s in D["frozen"]["signatures"]:
        for analysis in ("primary", "sensitivity"):
            if analysis not in s["analyses"] or s["model"] not in S.MODELS:
                continue
            e = met[s["id"]]["endpoints"]["OS"]
            if not e["replicates"]:
                continue
            q = p3[f"{analysis}/{s['model']}"]
            cand.append((bool(q["P3_pass"]), e["above_floor"], analysis, s, q))
    if not cand:
        raise ValueError("no ladder pick replicates on METABRIC OS")
    p3_pass, _, analysis, s, q = max(cand, key=lambda t: (t[0], t[1]))
    return {"analysis": analysis, "signature": s, "subtype": q["argmax_subtype"], "p3": q,
            "metabric": met[s["id"]]}


def _km_curve(t, e):
    from lifelines import KaplanMeierFitter
    k = KaplanMeierFitter().fit(t, e)
    sf = k.survival_function_.iloc[:, 0]
    return sf.index.to_numpy(float), sf.to_numpy(float)


def fig_km(D, out):
    import matplotlib.pyplot as plt
    ch = km_choice(D); s = ch["signature"]; m = s["model"]; analysis = ch["analysis"]
    L = D["L"]; P = picks(L, analysis)
    ev, tm = C.ENDPOINTS[L["endpoint"]]
    clin = pd.read_csv(D["files"]["clinical"], index_col=0)
    zexpr = zscore_genes(pd.read_parquet(D["files"]["expression"]))
    df = clin.assign(score=score_signature(zexpr, s["genes"])).dropna(subset=["score", tm, ev])
    df["years"] = df[tm] / 365.25
    sub = ch["subtype"]
    panels = [("All patients", df), (f"PAM50 {sub} tumours", df[df["subtype"] == sub])]
    shades = {"low": S.mix(S.COL[m], 0.55), "mid": S.mix(S.COL[m], 0.28), "high": S.COL[m]}
    fig, axes = plt.subplots(1, 2, figsize=(S.WIDE, 4.0), sharey=True)
    fig.subplots_adjust(left=0.085, right=0.89, top=0.72, bottom=0.27, wspace=0.3)
    horizon, ticks = 8, [0, 2, 4, 6, 8]
    fwm = L["family_wise"][analysis][m]
    stats = [f"C {P.loc[m, 'cindex']:.3f} ({s['direction_label']} direction); family-wise p "
             f"{fmt_p(fwm['pick_fw_p'])} for this pick,\n{fmt_p(fwm['fw_p'])} for the family maximum",
             f"within {sub}: C {ch['p3']['within_subtype_cindex'][sub]:.3f}, permutation p "
             f"{fmt_p(ch['p3']['p'])}"]
    for ax, (title, d), note in zip(axes, panels, stats):
        d = d.assign(tertile=pd.qcut(d["score"], 3, labels=["low", "mid", "high"]))
        ends = {}
        for t in ("low", "mid", "high"):
            g = d[d["tertile"] == t]
            x, yv = _km_curve(g["years"], g[ev])
            ax.step(x, yv, where="post", color=shades[t], lw=2.0)
            ends[t] = yv[x <= horizon][-1]
        prev = None
        for t in sorted(ends, key=ends.get, reverse=True):
            yl = ends[t] if prev is None else min(ends[t], prev - 0.06)
            ax.text(horizon + 0.25, yl, f"{t} score", va="center", fontsize=S.FS_SMALL)
            prev = yl
        ax.set_xlim(0, horizon); ax.set_ylim(0, 1.02); ax.set_xticks(ticks)
        ax.set_title(f"{title}: {len(d):,} patients, {int(d[ev].sum())} deaths", fontsize=S.FS, pad=34)
        ax.text(0, 1.03, note, transform=ax.transAxes, fontsize=S.FS_SMALL, color=S.MUTED, va="bottom")
        ax.set_xlabel("Years from diagnosis")
        S.clean(ax, left=True); ax.grid(axis="y", color=S.RULE, lw=0.6); ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=3, labelsize=S.FS_SMALL)
        for j, t in enumerate(("high", "mid", "low")):
            g = d[d["tertile"] == t]
            ax.text(-0.06, -0.30 - 0.075 * j, t, transform=ax.transAxes, ha="right", fontsize=S.FS_SMALL,
                    color=S.MUTED)
            for tk in ticks:
                ax.text(tk, -0.30 - 0.075 * j, f"{int((g['years'] >= tk).sum())}",
                        transform=ax.get_xaxis_transform(), ha="center", fontsize=S.FS_SMALL,
                        family="DejaVu Sans Mono")
        ax.text(-0.06, -0.225, "at risk", transform=ax.transAxes, ha="right", fontsize=S.FS_SMALL,
                color=S.MUTED, fontstyle="italic")
    axes[0].set_ylabel("Overall survival")
    lab = programme(analysis, m, P.loc[m, "key"], D["sigs"])
    fig.suptitle(f"TCGA-BRCA, the discovery cohort where it was selected: {S.NAME[m]} "
                 f"{ANALYSIS_NAME[analysis].lower()} pick ({lab}), score tertiles",
                 x=0.01, ha="left", fontsize=S.FS, fontweight="bold", y=0.975)
    me = ch["metabric"]["endpoints"]; mp = ch["metabric"].get("p3", {}).get(sub, {}).get("OS")
    line = (f"Frozen and re-scored in METABRIC: OS C {me['OS']['cindex_oriented']:.3f} (floor p95 "
            f"{me['OS']['floor_p95']:.3f}), DSS C {me['DSS']['cindex_oriented']:.3f} (floor p95 "
            f"{me['DSS']['floor_p95']:.3f})")
    if mp:
        line += f"; within {sub}, OS C {mp['cindex_oriented']:.3f}, permutation p {fmt_p(mp['perm_p'])}"
    fig.text(0.01, 0.905, line, fontsize=S.FS_SMALL, color=S.MUTED, ha="left", va="bottom")
    return S.save(fig, out, "fig_km"), ch


# ---- 7. the margin (P2) -----------------------------------------------------------------------
def fig_margin(D, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    L = D["L"]; comps = D["metabric"]["comparisons"]
    rows = []
    for analysis in ("primary", "sensitivity"):
        for m in ("scgpt", "geneformer"):
            t = L["p2_corrected"][analysis][m]
            fm_id = next(s["id"] for s in D["frozen"]["signatures"]
                         if s["model"] == m and analysis in s["analyses"])
            met = {c["endpoint"]: c for c in comps if c["analysis"] == analysis and c["fm"] == fm_id}
            rows.append((analysis, m, t, met))
    ys, heads, ybot = _layout([("PRE-SPECIFIED PICKS", 2), ("SENSITIVITY PICKS (POST HOC)", 2)])
    fig, ax = plt.subplots(figsize=(S.WIDE, 4.2))
    fig.subplots_adjust(left=0.25, right=0.80, top=0.95, bottom=0.27)
    marks = [(0.27, "o", True), (0.0, "s", False), (-0.27, "^", False)]    # TCGA OS, METABRIC OS, DSS
    for (_, m, t, met), y in zip(rows, ys):
        col = S.COL[m]
        vals = [(t["margin"], t["ci"]), (met["OS"]["margin"], met["OS"]["ci"]),
                (met["DSS"]["margin"], met["DSS"]["ci"])]
        for (dy, mk, filled), (v, ci) in zip(marks, vals):
            ax.plot(ci, [y + dy] * 2, color=col, lw=1.6, zorder=3)
            ax.plot(v, y + dy, mk, ms=6.5, mec=col, mfc=col if filled else S.SURFACE, mew=1.4, zorder=4)
            txt = f"{v:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]".replace("-", "\u2212")
            ax.text(1.01, y + dy, txt, transform=ax.get_yaxis_transform(),
                    va="center", fontsize=S.FS_SMALL, family="DejaVu Sans Mono")
    ax.axvline(0, color=S.INK, lw=0.9, zorder=1)
    labels = [f"{S.NAME[m]} − HVG-PCA" for _, m, _, _ in rows]
    row_labels(ax, labels, ys)
    for y, h in heads:
        group_header(ax, y, h)
    ax.set_yticks([]); ax.set_ylim(ybot + 0.25, 0.45)
    S.clean(ax); S.xgrid(ax)
    lo = min(min(t["ci"][0], met["OS"]["ci"][0], met["DSS"]["ci"][0]) for _, _, t, met in rows)
    hi = max(max(t["ci"][1], met["OS"]["ci"][1], met["DSS"]["ci"][1]) for _, _, t, met in rows)
    ax.set_xlim(np.floor(lo * 50) / 50 - 0.005, max(np.ceil(hi * 50) / 50, 0.02) + 0.005)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(signed))
    ax.set_xlabel("Margin over the matched-random floor, foundation-model pick minus HVG-PCA pick (ΔC)")
    nb = rows[0][2]["n_boot"]
    handles = [Line2D([], [], color=S.MUTED, marker="o", ms=6.5, mfc=S.MUTED, lw=1.6,
                      label=f"TCGA OS: picks re-selected in each of {nb:,} patient bootstraps (95%)"),
               Line2D([], [], color=S.MUTED, marker="s", ms=6.5, mfc=S.SURFACE, mew=1.4, lw=1.6,
                      label="METABRIC OS: frozen picks, paired patient bootstrap (95%)"),
               Line2D([], [], color=S.MUTED, marker="^", ms=6.5, mfc=S.SURFACE, mew=1.4, lw=1.6,
                      label="METABRIC DSS: frozen picks, paired patient bootstrap (95%)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.25, 0.0), ncol=1, fontsize=S.FS_SMALL,
               borderaxespad=0.4)
    fig.text(0.805, 0.955, "margin [95% interval]", fontsize=S.FS_SMALL, color=S.MUTED, va="bottom")
    fig.text(0.01, 0.04, "Right of 0: the foundation\nmodel's pick is further above\nits floor than the "
             "baseline's.", fontsize=S.FS_SMALL, color=S.MUTED, va="bottom", linespacing=1.3)
    return S.save(fig, out, "fig_margin")


# ---- 8. subtype association -------------------------------------------------------------------
SUBTYPE_ORDER = ["LumA", "LumB", "Her2", "Basal", "Normal"]
SUBTYPE_LABEL = {"LumA": "LumA", "LumB": "LumB", "Her2": "HER2", "Basal": "Basal", "Normal": "Normal"}


def _subtype_row(D, cohort: str, analysis: str, m: str, key: str) -> pd.Series:
    A = D["subtype"]
    hit = A[(A["cohort"] == cohort) & (A["model"] == m)
            & A["analyses"].str.split(";").map(lambda xs: analysis in xs)]
    if len(hit) != 1:
        raise ValueError(f"subtype_association.csv: {len(hit)} {cohort} rows for {analysis}/{m}")
    r = hit.iloc[0]
    if f"{float(r['resolution']):.1f}/{r['signature']}" != key:
        raise ValueError(f"subtype_association.csv: {analysis}/{m} is not the ladder pick {key}")
    return r


def fig_subtype(D, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    L = D["L"]
    fig, axes = plt.subplots(2, 3, figsize=(S.WIDE, 6.4), sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.82, bottom=0.14, hspace=0.78, wspace=0.07)
    x = np.arange(len(SUBTYPE_ORDER))
    for i, analysis in enumerate(("primary", "sensitivity")):
        P = picks(L, analysis)
        for j, m in enumerate(S.MODELS):
            ax = axes[i, j]; key = P.loc[m, "key"]; col = S.COL[m]
            t = _subtype_row(D, "TCGA", analysis, m, key)
            mb = _subtype_row(D, "METABRIC", analysis, m, key)
            for r, dx, filled in ((t, -0.14, True), (mb, 0.14, False)):
                med = [r[f"median_{g}"] for g in SUBTYPE_ORDER]
                lo = [r[f"q25_{g}"] for g in SUBTYPE_ORDER]
                hi = [r[f"q75_{g}"] for g in SUBTYPE_ORDER]
                for xx, a, b in zip(x + dx, lo, hi):
                    ax.plot([xx, xx], [a, b], color=col, lw=1.6, zorder=3)
                ax.plot(x + dx, med, "o" if filled else "s", ms=6, mec=col, mfc=col if filled else S.SURFACE,
                        mew=1.4, lw=0, zorder=4)
            ax.axhline(0, color=S.MUTED, lw=0.7, zorder=1)
            ax.set_xticks(x); ax.set_xticklabels([SUBTYPE_LABEL[g] for g in SUBTYPE_ORDER])
            ax.set_xlim(-0.5, len(x) - 0.5)
            S.clean(ax, left=(j == 0))
            if j == 0:
                ax.tick_params(axis="y", length=3, labelsize=S.FS_SMALL)
                ax.set_ylabel("signature score (mean z)", fontsize=S.FS_SMALL)
            ax.set_title(f"{S.NAME[m]} · {programme(analysis, m, key, D['sigs'])}", fontsize=S.FS, pad=30)
            g = t["tracked_group"]
            auc = t["auc_tracked"]
            word = "higher" if auc >= 0.5 else "lower"
            ax.text(0, 1.03, f"ε² {t['epsilon2']:.2f} TCGA, {mb['epsilon2']:.2f} METABRIC\n{word} in "
                    f"{SUBTYPE_LABEL.get(g, g)} tumours (TCGA AUC {auc:.2f})", transform=ax.transAxes,
                    fontsize=S.FS_SMALL, color=S.MUTED, va="bottom", linespacing=1.3)
        fig.text(0.075, 0.975 if i == 0 else 0.515,
                 "PRE-SPECIFIED PICKS" if i == 0 else "SENSITIVITY PICKS (POST HOC)",
                 fontsize=S.FS_SMALL, color=S.MUTED, fontweight="bold", va="top")
    n_t = int(_subtype_row(D, "TCGA", "primary", "scgpt", picks(L, "primary").loc["scgpt", "key"])["n"])
    n_m = int(_subtype_row(D, "METABRIC", "primary", "scgpt", picks(L, "primary").loc["scgpt", "key"])["n"])
    handles = [Line2D([], [], color=S.MUTED, marker="o", ms=6, mfc=S.MUTED, lw=1.6,
                      label=f"TCGA-BRCA, PAM50 call ({n_t:,} patients): median, interquartile range"),
               Line2D([], [], color=S.MUTED, marker="s", ms=6, mfc=S.SURFACE, mew=1.4, lw=1.6,
                      label=f"METABRIC ({n_m:,} patients; claudin-low not shown)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.07, 0.0), ncol=2, fontsize=S.FS_SMALL,
               borderaxespad=0.4)
    return S.save(fig, out, "fig_subtype")


# ---- 9. donor mixing --------------------------------------------------------------------------
def fig_donor(D, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scfm.stages.donors import states_with_margins
    d = states_with_margins(pd.read_csv(D["files"]["composition"]), D["scores"])
    summ = {r["model"]: r for r in D["donor_mixing"]["per_representation"]}
    marker = {"hvg_pca": "o", "scgpt": "s", "geneformer": "^"}
    for m in S.MODELS:
        g = d[d["model"] == m]
        if len(g) != summ[m]["n_states"] or int((g["top_donor_frac"] >= C.SINGLE_DONOR_FRAC).sum()) != \
                summ[m]["n_single_donor"]:
            raise ValueError(f"donor mixing of {m}: the states drawn differ from donor_mixing.json")
    fig = plt.figure(figsize=(S.WIDE, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1, 1], left=0.14, right=0.975, top=0.8, bottom=0.25,
                          wspace=0.36)
    a0 = fig.add_subplot(gs[0]); a1 = fig.add_subplot(gs[1]); a2 = fig.add_subplot(gs[2], sharey=a1)
    rng = np.random.default_rng(0)
    for i, m in enumerate(S.MODELS):
        g = d[d["model"] == m]; y = -i
        a0.scatter(g["top_donor_frac"], y + rng.uniform(-0.22, 0.22, len(g)), s=11, color=S.COL[m],
                   alpha=0.75, lw=0, marker=marker[m], zorder=3)
    a0.set_yticks([0, -1, -2])
    a0.set_yticklabels([f"{S.NAME[m]}\n{summ[m]['n_single_donor']} of {summ[m]['n_states']} single-donor\n"
                        f"({summ[m]['share_single_donor']:.0%})" for m in S.MODELS], fontsize=S.FS_SMALL,
                       linespacing=1.3)
    a0.set_ylim(-2.55, 0.75); a0.set_xlim(0, 1.0)
    a0.set_xlabel("share of the state's cells from its largest donor")
    a0.set_title("Donor mixing of every kept state", fontsize=S.FS, pad=44)
    for ax, col, title in ((a1, "above_floor", "Risk direction (pre-specified)"),
                           (a2, "above_floor_or", "Either direction (post hoc)")):
        for m in S.MODELS:
            g = d[d["model"] == m]
            ax.scatter(g["top_donor_frac"], g[col], s=13, color=S.COL[m], alpha=0.8, lw=0, marker=marker[m],
                       zorder=3)
            key = "best_margin_multi_donor" if col == "above_floor" else "best_margin_or_multi_donor"
            multi = g[g["top_donor_frac"] < C.MULTI_DONOR_FRAC]
            best = multi.loc[multi[col].idxmax()]
            if abs(best[col] - summ[m][key]) > 1e-12:
                raise ValueError(f"{m}: best multi-donor margin differs from donor_mixing.json")
            ax.plot(best["top_donor_frac"], best[col], marker=marker[m], ms=9, mfc="none", mec=S.INK, mew=1.1,
                    zorder=5)
            ax.plot([0, C.MULTI_DONOR_FRAC], [best[col]] * 2, color=S.COL[m], lw=1.1, ls=(0, (3, 2)),
                    zorder=2)
        key = "best_margin_multi_donor" if col == "above_floor" else "best_margin_or_multi_donor"
        v = {m: f"{summ[m][key]:+.3f}".replace("-", "\u2212") for m in S.MODELS}
        best_txt = f"HVG-PCA {v['hvg_pca']} · scGPT {v['scgpt']}\nGeneformer {v['geneformer']}"
        ax.text(0, 1.015, "best multi-donor state:\n" + best_txt, transform=ax.transAxes, fontsize=S.FS_SMALL,
                color=S.MUTED, va="bottom", linespacing=1.3)
        ax.axhline(0, color=S.MUTED, lw=0.7, zorder=1)
        ax.set_xlim(0, 1.0)
        ax.set_title(title, fontsize=S.FS, pad=44)
        ax.set_xlabel("share from the largest donor")
        S.clean(ax, left=True); ax.tick_params(axis="y", length=3, labelsize=S.FS_SMALL)
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.04))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(signed))
    a1.set_ylabel("C above the matched-random floor", fontsize=S.FS_SMALL)
    for ax in (a0, a1, a2):
        for v in (C.MULTI_DONOR_FRAC, C.SINGLE_DONOR_FRAC):
            ax.axvline(v, color=S.RULE if ax is not a0 else S.MUTED, lw=0.8, ls=(0, (2, 2)), zorder=1)
        ax.tick_params(axis="x", labelsize=S.FS_SMALL)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    S.clean(a0)
    handles = [Line2D([], [], color=S.COL[m], marker=marker[m], lw=0, ms=6, label=S.NAME[m])
               for m in S.MODELS]
    handles += [Line2D([], [], color=S.MUTED, lw=1.1, ls=(0, (3, 2)), marker="o", ms=9, mfc="none", mec=S.INK,
                       label=f"best multi-donor state (largest donor < {C.MULTI_DONOR_FRAC:.0%}) "
                             "and its margin"),
                Line2D([], [], color=S.MUTED, lw=0.8, ls=(0, (2, 2)),
                       label=f"{C.MULTI_DONOR_FRAC:.0%} and {C.SINGLE_DONOR_FRAC:.0%} (single-donor) "
                             "thresholds")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.1, 0.0), ncol=3, fontsize=S.FS_SMALL,
               borderaxespad=0.4, columnspacing=1.6)
    return S.save(fig, out, "fig_donor")


# ---- figures.json -----------------------------------------------------------------------------
def describe(D, paths: dict, km: dict | None) -> list[dict]:
    f = {k: relpath(str(v)) for k, v in D["files"].items()}
    L = D["L"]
    labels = {f"{a}/{m}": {"pick": picks(L, a).loc[m, "key"],
                           "programme": programme(a, m, picks(L, a).loc[m, "key"], D["sigs"])}
              for a in ("primary", "sensitivity") for m in S.MODELS}
    spec = {
        "fig_ladder": {
            "title": "The ladder: each representation's pick against its floor, its null and three "
                     "references",
            "question": "Does any representation's best cell state predict overall survival beyond chance "
                        "selection, and how far is it from routine clinical information?",
            "shows": "Two panels sharing the C axis. Top: pre-specified analysis (risk direction). Bottom: "
                     "post hoc sensitivity analysis (either direction, C read in the direction the score "
                     "acts). One row per representation.",
            "how_to_read": [
                "Bar: Harrell's C of the pick's score alone, unadjusted. Its whisker holds the pick fixed, "
                "so it is optimistic.",
                "A family clears its null when its diamond (the family maximum) lies right of the thick "
                "tick; the family-wise p belongs to that maximum, which can be a different signature from "
                "the pick.",
                "Right columns: the pick's C and its own p against the same family-wise null, then the "
                "family maximum's C and its family-wise p (bold when < 0.05); p-values and the thick tick "
                f"are from {L['family_wise_10k']['n_perm']:,} outcome permutations."],
            "data_files": [f["ladder"]]},
        "fig_states": {
            "title": "Where the picked states sit in each representation",
            "question": "What are the picked cell states, and are they shared programmes or single "
                        "patients' cells?",
            "shows": "Uniform manifold approximation and projection (UMAP) of all atlas cells per "
                     "representation (grey), the pre-specified pick in the representation's colour and the "
                     "sensitivity pick in a lighter tint, each labelled with its programme, cell count, "
                     "donors (top-donor share) and datasets (top-dataset share).",
            "how_to_read": ["Each panel is its own embedding; positions are not comparable across panels.",
                            "A top-donor share near 100% means the state is one patient's cells."],
            "programme_labels": labels,
            "data_files": [f[f"cells_{m}"] for m in S.MODELS] + [f["signatures"], f["composition"],
                                                                   f["ladder"]]},
        "fig_added_value": {
            "title": "Each pick added to a Cox model of age + stage",
            "question": "Does any pick add prognostic information to routine clinical variables?",
            "shows": "Forest plot, one row per pick (pre-specified, sensitivity) and the published "
                     "proliferation score. Left: Cox hazard ratio (HR) per standard deviation (SD) of the "
                     "score with its Wald 95% interval. Right: change in cross-validated C, with the "
                     "signature re-selected inside every training fold (nested).",
            "how_to_read": ["HRs are nominal: the picks were selected on these outcomes. Only the nested "
                            "change in C and the published score are free of that selection.",
                            "ΔC whiskers are the 2.5–97.5 percentile range over the 5 cross-validation "
                            "repeats on the same patients, not a confidence interval."],
            "data_files": [f["added_value"], f["ladder"], f["signatures"]]},
        "fig_stratify": {
            "title": "Multivariable models with and without the states",
            "question": "Do all of a representation's states together, in a penalised or tree model, "
                        "improve on age + stage?",
            "shows": "Cross-validated C (5 × 5-fold) of age + stage, + proliferation score, + all states of "
                     "each representation (ridge Cox and XGBoost), and, on the PAM50 subset, + PAM50 "
                     "subtype.",
            "how_to_read": ["The number before \"states\" in each label is the count of state features "
                            "added.",
                            "Differences of a few thousandths are within the spread across repeats."],
            "data_files": [f["stratify"], f["stratify_rows"]]},
        "fig_replication": {
            "title": "Frozen signatures in TCGA and in METABRIC",
            "question": "Do the signatures selected in TCGA replicate in an independent cohort?",
            "shows": "Every frozen signature (picks, family maxima, published proliferation score): oriented "
                     "C in TCGA OS, METABRIC OS and METABRIC disease-specific survival (DSS) against each "
                     "cohort's matched-random floor 95th percentile (p95), and, for the sensitivity picks, "
                     "METABRIC OS within the PAM50 subtype of the TCGA lead.",
            "how_to_read": ["The TCGA column is discovery: the C there is selected, and its dots are always "
                            "filled.",
                            "The subtype panel has no random-gene floor; its filled dots mark an "
                            "outcome-permutation p < 0.05, which is not the pre-specified replication "
                            "criterion."],
            "data_files": [f["metabric"], f["frozen"], f["scores"], f["nulls"], f["ladder"]]},
        "fig_km": {
            "title": "TCGA (discovery) survival curves for the signature that replicated in METABRIC",
            "question": "What does the lead replicated signature look like as survival curves, overall and "
                        "within its subtype?",
            "shows": "Descriptive: the pick and its direction were chosen on these TCGA-BRCA patients, so "
                     "the separation of the curves is optimistic; its replication is the METABRIC line in "
                     "the subtitle. Overall survival by tertile of the signature score (tertiles cut within "
                     "each panel), all patients and within the PAM50 subtype where its within-subtype C was "
                     "highest; numbers at risk below; curves stop at 8 years.",
            "how_to_read": ["Choice rule: a ladder pick that replicates on METABRIC OS, preferring one that "
                            "passed P3 in TCGA, then the larger METABRIC OS margin over its floor."],
            "data_files": [f["ladder"], f["metabric"], f["frozen"], f["clinical"], f["expression"]]},
        "fig_margin": {
            "title": "Foundation model minus baseline (P2, floor-adjusted)",
            "question": "Is either foundation model's pick further above its floor than the baseline's?",
            "shows": "Margin = (C − floor) of the foundation-model pick minus (C − floor) of the HVG-PCA "
                     "pick, per analysis and model: TCGA with the selection-aware bootstrap (picks "
                     "re-selected in each resample), METABRIC OS and DSS with frozen picks and a paired "
                     "patient bootstrap.",
            "how_to_read": ["P2 passes when the TCGA lower bound is above 0."],
            "data_files": [f["ladder"], f["metabric"], f["frozen"]]},
        "fig_subtype": {
            "title": "The picks' scores by PAM50 subtype, in TCGA and in METABRIC",
            "question": "Do the picked states' bulk scores track the intrinsic subtypes?",
            "shows": "One panel per pick: the median and interquartile range of the signature score in each "
                     "PAM50 subtype, TCGA-BRCA (filled) and METABRIC (open). No outcomes are used. Above "
                     "each panel: the rank-based effect size of subtype on the score (ε², Kruskal-Wallis "
                     "H / (n − 1)) "
                     "in both cohorts, and the one-vs-rest AUC of the TCGA subtype the score tracks most.",
            "how_to_read": ["ε² near 0: the score does not differ by subtype; the published proliferation "
                            "score, for comparison, has ε² above 0.5 in both cohorts.",
                            "AUC: the probability that a tumour of that subtype scores higher than one of "
                            "another subtype; below 0.5 means lower in that subtype."],
            "data_files": [f["subtype"], f["ladder"], f["signatures"]]},
        "fig_donor": {
            "title": "Donor mixing of the states, and their survival signal",
            "question": "Do the foundation models group cells across patients more than HVG-PCA, and are "
                        "the better-mixed states more prognostic?",
            "shows": "Left: the share of each kept state's cells that comes from its largest donor, per "
                     "representation. Middle and right: each state's above-floor C margin (its 50-gene "
                     "signature, TCGA overall survival) against that share, in the risk and in the either-"
                     "direction reading; dashed lines mark the best multi-donor state of each "
                     "representation.",
            "how_to_read": ["A state is single-donor when one donor supplies at least 80% of its cells, and "
                            "multi-donor when no donor supplies 50%.",
                            "The three resolutions are pooled, so nested states appear more than once; the "
                            "margins are selected on TCGA outcomes and carry no interval."],
            "data_files": [f["composition"], f["scores"], f["donor_mixing"]]},
    }
    out = []
    from PIL import Image
    for name in FIGURES:
        p = paths.get(name)
        e = {"name": name, "file": f"results/report/{name}.png", "produced": p is not None,
             "size": "wide", "display_width_px": 1140, **spec[name]}
        if p is not None:
            with Image.open(p) as im:
                e["pixels"] = list(im.size)
        if name == "fig_km" and km:
            a, m = km["analysis"], km["signature"]["model"]
            e["choice"] = {"analysis": a, "signature": km["signature"]["id"], "subtype": km["subtype"],
                           "programme": programme(a, m, picks(L, a).loc[m, "key"], D["sigs"])}
        out.append(e)
    return out


# ---- text outputs -----------------------------------------------------------------------------
def write_shortlist(D, out):
    L = D["L"]; scores = D["scores"]; sigs = D["sigs"]
    sl = scores[(scores["kind"] == "fm") & (scores["above_floor"] > 0)].copy()
    null95 = {m: v["null_p95"] for m, v in L["family_wise_10k"]["primary"].items()}
    sl = sl[sl.apply(lambda r: r["cindex"] > null95[r["model"]], axis=1)]
    sl = sl.sort_values("above_floor", ascending=False).head(10)
    md = ["# Candidate cell states", "",
          f"Endpoint {L['endpoint']}, {L['n_patients']} patients, {L['n_events']} events. "
          "Listed: foundation-model signatures above the matched-random floor and above the "
          "family-wise permutation null. Each is a hypothesis, not a result.", ""]
    for _, r in sl.iterrows():
        genes = sigs[r["model"]][str(r["resolution"])][r["signature"]]["genes"]
        md += [f"## {r['model']} · {r['signature']} — {r['top_cell_type']}",
               f"C = {r['cindex']:.3f} (floor {r['floor_mean']:.3f}, +{r['above_floor']:.3f}); "
               f"HR {r['hr']:.2f}, p = {r['p']:.2g}; {r['n_cells']:,} cells; "
               f"HVG fraction {r['hvg_frac']:.2f}",
               "", "`" + "`, `".join(genes[:25]) + ("`, …" if len(genes) > 25 else "`"), ""]
    (out / "shortlist.md").write_text("\n".join(md))
    (out / "validation.md").write_text("\n".join([
        "# Validation strategy — per shortlisted state", "",
        "1. **Independent bulk cohort.** Freeze the signature and its direction and re-score it in "
        "METABRIC under the criterion of DESIGN §13: oriented C above the METABRIC matched-random floor "
        "p95, and an age-adjusted, cohort-stratified hazard ratio in the TCGA direction with p < 0.05. A "
        "state that does not replicate there is dropped.",
        "2. **Back to the atlas.** Confirm the marker genes are expressed in the annotated cell type the "
        "state was assigned to, not in a contaminating population; check donor spread.",
        "3. **Orthogonal tissue assay.** Name the multiplexed IHC / spatial panel that would show the "
        "state's abundance in FFPE sections from a cohort with outcome, and the effect size it would "
        "need to reach to matter clinically.", "",
        "Nothing on the shortlist is a target or biomarker until step 1 has passed.", ""]))
    return len(sl)


def atlas_duplicates(D, data_dir: Path, out: Path):
    """Atlas cells whose raw count vector is identical to another atlas cell's (descriptive QC)."""
    path = data_dir / "atlas.h5ad"
    if not path.exists():
        print("  atlas_duplicates.json skipped: data/atlas.h5ad not present (stage data)")
        return None
    import hashlib
    from collections import Counter, defaultdict

    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(path)
    X = sp.csr_matrix(a.layers["counts"])
    groups = defaultdict(list)
    for i in range(X.shape[0]):
        lo, hi = X.indptr[i], X.indptr[i + 1]
        key = X.indices[lo:hi].tobytes() + np.asarray(X.data[lo:hi], dtype=np.float64).tobytes()
        groups[hashlib.sha1(key).hexdigest()].append(i)
    dup = [g for g in groups.values() if len(g) > 1]
    obs = a.obs
    pairs = Counter(" + ".join(sorted({str(d)[:8] for d in obs["dataset_id"].iloc[g]})) for g in dup)
    in_dup = {a.obs_names[i] for g in dup for i in g}
    res = {"definition": "cells whose raw count vector (layers['counts']) is identical to that of another "
                         "atlas cell; not removed from the analysis",
           "n_cells": int(X.shape[0]), "n_groups": len(dup),
           "group_sizes": dict(Counter(len(g) for g in dup)),
           "n_redundant_cells": int(sum(len(g) - 1 for g in dup)),
           "n_cells_in_groups": int(sum(len(g) for g in dup)),
           "n_unique_cells": int(X.shape[0] - sum(len(g) - 1 for g in dup)),
           "all_groups_same_donor": all(obs["donor_id"].iloc[g].nunique() == 1 for g in dup),
           "all_groups_across_datasets": all(obs["dataset_id"].iloc[g].nunique() == len(g) for g in dup),
           "dataset_pairs": dict(pairs.most_common()), "picks": {}}
    for analysis in ("primary", "sensitivity"):
        P = picks(D["L"], analysis)
        for m in S.MODELS:
            f = D["files"][f"cells_{m}"]
            if not f.exists():
                continue
            res_, sig = P.loc[m, "key"].split("/")
            cells = pd.read_parquet(f, columns=["obs_name", f"leiden_{res_}"])
            sel = cells[cells[f"leiden_{res_}"] == D["sigs"][m][res_][sig]["cluster"]]["obs_name"].astype(str)
            res["picks"][f"{analysis}/{m}"] = {"pick": P.loc[m, "key"], "n_cells": int(len(sel)),
                                               "n_cells_with_duplicate": int(sel.isin(in_dup).sum())}
    (out / "atlas_duplicates.json").write_text(json.dumps(res, indent=2) + "\n")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--translate-dir", default=str(C.RESULTS / "translate"))
    ap.add_argument("--ladder-dir", default=str(C.RESULTS / "ladder"))
    ap.add_argument("--stratify-dir", default=str(C.RESULTS / "stratify"))
    ap.add_argument("--replicate-dir", default=str(C.RESULTS / "replicate"))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "report"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    out = Path(args.out_dir)
    if args.dry_run:
        print(json.dumps({"figures": FIGURES, "docs": ["figures.json", "shortlist.md", "validation.md",
                                                                   "atlas_duplicates.json"]},
                         indent=2)); return
    out.mkdir(parents=True, exist_ok=True)
    S.setup()
    D = load(args)
    paths = {"fig_ladder": fig_ladder(D, out), "fig_states": fig_states(D, out),
             "fig_added_value": fig_added_value(D, out), "fig_stratify": fig_stratify(D, out),
             "fig_replication": fig_replication(D, out)}
    paths["fig_km"], km = fig_km(D, out)
    paths["fig_margin"] = fig_margin(D, out)
    paths["fig_subtype"] = fig_subtype(D, out)
    paths["fig_donor"] = fig_donor(D, out)
    figs = describe(D, paths, km)
    (out / "figures.json").write_text(json.dumps({"figures": figs, "style": {
        "palette": {S.NAME.get(k, k): v for k, v in S.COL.items()}, "ink": S.INK, "muted": S.MUTED,
        "rule": S.RULE, "dpi": S.DPI, "inches_wide": S.WIDE, "display_px_per_inch": 114}}, indent=2) + "\n")
    n = write_shortlist(D, out)
    atlas_duplicates(D, Path(args.data_dir), out)
    made = [f["name"] for f in figs if f["produced"]]
    record_params(out, args, extra={"n_shortlisted": n, "figures": made})
    print(f"  {len(made)} figures + figures.json + shortlist ({n}) -> {out}")


if __name__ == "__main__":
    main()
