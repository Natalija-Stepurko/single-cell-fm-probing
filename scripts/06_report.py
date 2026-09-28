"""Stage 06 — figures, the candidate shortlist, and the validation-strategy template.

Figures (matplotlib, one colour per arm, values printed at bar ends):
  fig_ladder.png     C-index per representation with patient-bootstrap intervals; floor, null
                     and PAM50 reference drawn as rungs
  fig_margin.png     FM − baseline margin distribution (P2) — the figure the study is about
  fig_km.png         Kaplan–Meier, best FM signature, tertiles, overall and within subtype (P3)
  fig_hvg.png        HVG fraction of every signature against its above-floor margin (P4)

shortlist.md        FM signatures passing P1–P2, with their marker genes and cell-type annotation
validation.md       the three-step validation each shortlisted state would need
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import config as C
import qc_common as qc
from surv_common import score_signature, zscore_genes

COL = {"hvg_pca": "#4A4F55", "scgpt": "#1F6F6B", "geneformer": "#C2681A"}


def fig_ladder(ladder, ref, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    y = np.arange(len(ladder))
    for i, (m, r) in enumerate(ladder.iterrows()):
        ax.barh(i, r["cindex"] - 0.5, left=0.5, color=COL.get(m, "#888"), height=0.55)
        ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], color="#16191D", lw=1.2)
        ax.text(r["ci_hi"] + 0.005, i, f"{r['cindex']:.3f}", va="center", fontsize=9,
                family="monospace")
        ax.plot([r["floor_mean"]] * 2, [i - .3, i + .3], color="#888", lw=1, ls=":")
    ax.axvline(0.5, color="#B9C0C6", lw=1)
    if np.isfinite(ref):
        ax.axvline(ref, color="#16191D", lw=1, ls="--"); ax.text(ref, len(y) - .4, "PAM50",
                                                                 fontsize=8, ha="center")
    ax.set_yticks(y); ax.set_yticklabels(ladder.index); ax.set_xlabel("concordance index")
    ax.set_title("Best signature per representation — patient bootstrap 95%; dotted = matched-random floor",
                 fontsize=9, loc="left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out / "fig_ladder.png", dpi=150); plt.close(fig)


def fig_margin(P, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    for i, (m, p) in enumerate(P.items()):
        lo, hi = p["P2_margin_ci"]; mu = p["P2_margin_over_baseline"]
        ax.plot([lo, hi], [i, i], color=COL.get(m, "#888"), lw=3)
        ax.plot(mu, i, "o", color=COL.get(m, "#888"))
        ax.text(hi + 0.003, i, f"{mu:+.3f}", va="center", fontsize=9, family="monospace")
    ax.axvline(0, color="#16191D", lw=1)
    ax.set_yticks(range(len(P))); ax.set_yticklabels(list(P))
    ax.set_xlabel("C-index margin over the HVG-PCA baseline (paired patient bootstrap)")
    ax.set_title("P2 — does the foundation model add prognostic value over the linear baseline?",
                 fontsize=9, loc="left")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out / "fig_margin.png", dpi=150); plt.close(fig)


def fig_km(clin, zexpr, genes, ev, tm, model, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    df = clin.assign(score=score_signature(zexpr, genes)).dropna(subset=["score", tm, ev])
    df["tertile"] = pd.qcut(df["score"], 3, labels=["low", "mid", "high"])
    subtypes = [s for s, g in df.dropna(subset=["subtype"]).groupby("subtype") if g[ev].sum() >= 10]
    fig, axes = plt.subplots(1, 1 + len(subtypes), figsize=(3.4 * (1 + len(subtypes)), 3.2),
                             squeeze=False)
    panels = [("all patients", df)] + [(s, df[df["subtype"] == s]) for s in subtypes]
    for ax, (title, d) in zip(axes[0], panels):
        for t, c in zip(["low", "mid", "high"], ["#4A4F55", "#B9C0C6", COL.get(model, "#1F6F6B")]):
            g = d[d["tertile"] == t]
            if len(g) > 5:
                KaplanMeierFitter().fit(g[tm], g[ev], label=t).plot_survival_function(
                    ax=ax, color=c, ci_show=False, lw=1.6)
        ax.set_title(f"{title} (n={len(d)})", fontsize=9); ax.set_xlabel("days"); ax.set_ylim(0, 1)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.suptitle(f"Best {model} signature — score tertiles", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(); fig.savefig(out / "fig_km.png", dpi=150); plt.close(fig)


def fig_hvg(scores, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    for m, g in scores.groupby("model"):
        ax.scatter(g["hvg_frac"], g["above_floor"], s=14, alpha=.7, color=COL.get(m, "#888"), label=m)
    ax.axhline(0, color="#16191D", lw=1); ax.axvline(0.5, color="#B9C0C6", lw=1, ls=":")
    ax.set_xlabel("fraction of signature genes that are HVGs")
    ax.set_ylabel("C-index above matched-random floor")
    ax.set_title("P4 — do the informative signatures live outside the HVG set?", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out / "fig_hvg.png", dpi=150); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default=str(C.DATA))
    ap.add_argument("--states-dir", default=str(C.RESULTS / "states"))
    ap.add_argument("--translate-dir", default=str(C.RESULTS / "translate"))
    ap.add_argument("--ladder-dir", default=str(C.RESULTS / "ladder"))
    ap.add_argument("--out-dir", default=str(C.RESULTS / "report"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    if args.dry_run:
        print(json.dumps({"figures": ["fig_ladder", "fig_margin", "fig_km", "fig_hvg"],
                          "docs": ["shortlist.md", "validation.md"]}, indent=2)); return

    L = json.load(open(Path(args.ladder_dir) / "ladder.json"))
    ladder = pd.DataFrame(L["ladder"]).set_index("model")
    scores = pd.read_csv(Path(args.translate_dir) / "scores.csv")
    sigs = json.load(open(Path(args.states_dir) / "signatures.json"))
    clin = pd.read_csv(Path(args.data_dir) / "bulk_clinical.csv", index_col=0)
    zexpr = zscore_genes(pd.read_parquet(Path(args.data_dir) / "bulk_expr.parquet"))
    ev, tm = C.ENDPOINTS[L["endpoint"]]

    fig_ladder(ladder, L["reference_pam50"], out)
    fig_margin(L["predictions"], out)
    fig_hvg(scores, out)
    fm = [m for m in ladder.index if ladder.loc[m, "kind"] == "fm"]
    if fm:
        top = max(fm, key=lambda m: ladder.loc[m, "above_floor"])
        r = ladder.loc[top]
        fig_km(clin, zexpr, sigs[top][str(scores.set_index("signature").loc[r["signature"], "resolution"])]
               [r["signature"]]["genes"], ev, tm, top, out)

    # shortlist: FM signatures above floor AND above the family-wise null
    null95 = {m: L["predictions"].get(m, {}) for m in fm}
    sl = scores[(scores["kind"] == "fm") & (scores["above_floor"] > 0)].copy()
    sl = sl[sl.apply(lambda r: r["cindex"] > ladder.loc[r["model"], "null_p95"], axis=1)]
    sl = sl.sort_values("above_floor", ascending=False).head(10)
    md = ["# Candidate cell states", "",
          f"Endpoint {L['endpoint']}, {L['n_patients']} patients, {L['n_events']} events. "
          "Listed: foundation-model signatures above the matched-random floor and above the "
          "family-wise permutation null. Each is a hypothesis, not a result.", ""]
    for _, r in sl.iterrows():
        genes = sigs[r["model"]][str(r["resolution"])][r["signature"]]["genes"]
        md += [f"## {r['model']} · {r['signature']} — {r['top_cell_type']}",
               f"C = {r['cindex']:.3f} (floor {r['floor_mean']:.3f}, +{r['above_floor']:.3f}); "
               f"HR {r['hr']:.2f}, p = {r['p']:.2g}; {r['n_cells']:,} cells; HVG fraction {r['hvg_frac']:.2f}",
               "", "`" + "`, `".join(genes[:25]) + ("`, …" if len(genes) > 25 else "`"), ""]
    (out / "shortlist.md").write_text("\n".join(md))

    (out / "validation.md").write_text("\n".join([
        "# Validation strategy — per shortlisted state", "",
        "1. **Independent bulk cohort.** Re-score the signature in METABRIC (n≈2,000, OS) with the same "
        "Cox model. A state that does not replicate there is dropped.",
        "2. **Back to the atlas.** Confirm the marker genes are expressed in the annotated cell type the "
        "state was assigned to, not in a contaminating population; check donor spread.",
        "3. **Orthogonal tissue assay.** Name the multiplexed IHC / spatial panel that would show the "
        "state's abundance in FFPE sections from a cohort with outcome, and the effect size it would "
        "need to reach to matter clinically.", "",
        "Nothing on the shortlist is a target or biomarker until step 1 has passed.", ""]))
    qc.record_params(out, args, extra={"n_shortlisted": len(sl)})
    print(f"  figures + shortlist ({len(sl)}) + validation template -> {out}")


if __name__ == "__main__":
    main()
