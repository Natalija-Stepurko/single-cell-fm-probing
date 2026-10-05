"""Build the project page as one self-contained HTML file: docs/index.html (served by GitHub Pages).

Every number on the page is read from results/ when the page is built, and every pass/fail verb in the
prose is chosen from the booleans in results/ladder/ladder.json and results/replicate/metabric.json and
checked with an assert, so a change in the results either changes the text or stops the build.
Figures are the PNGs that `scfm run report` writes, embedded as data URIs.

    python docs/site/build.py                -> docs/index.html
"""
import base64
import csv
import html
import json
import re
import sys
from datetime import date
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from scfm import config as C  # noqa: E402
from scfm.cli import STAGES  # noqa: E402

OUT = ROOT / "docs" / "index.html"
RES = C.RESULTS
REPORT = RES / "report"
REPO = "https://github.com/Natalija-Stepurko/single-cell-fm-probing"
SITE = "https://natalija-stepurko.com"
DESIGN_URL = f"{REPO}/blob/main/docs/DESIGN.md"

INK, MUTED, RULE, AXIS, PANEL = "#16191D", "#5B646E", "#DDE1E4", "#B9C0C6", "#FFFFFF"
MONO = "ui-monospace,SFMono-Regular,Menlo,Consolas,monospace"
SANS = "ui-sans-serif,system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif"
MODELS = ("hvg_pca", "scgpt", "geneformer")
FMS = ("scgpt", "geneformer")
NAME = {"hvg_pca": "HVG-PCA", "scgpt": "scGPT", "geneformer": "Geneformer"}
SUBTYPE = {"LumA": "luminal A", "LumB": "luminal B", "Her2": "HER2-enriched", "Basal": "basal-like"}
MINUS = "−"
WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}


# ----------------------------------------------------------------------------- results
def jload(rel):
    return json.loads((RES / rel).read_text())


def rows(rel):
    with open(RES / rel, newline="") as fh:
        return list(csv.DictReader(fh))


L = jload("ladder/ladder.json")
LP = jload("ladder_pfi/ladder.json")
M = jload("replicate/metabric.json")
FROZ = jload("replicate/frozen_signatures.json")
TPAR = jload("translate/params.json")
SPAR = jload("states/params.json")
SSUM = jload("states/summary.json")
STR = jload("stratify/cv_summary.json")
FIGS = jload("report/figures.json")
RUNLOG = jload("run_log.json")
AV = rows("ladder/added_value.csv")
COMP = rows("states/state_composition.csv")
SCORES = rows("translate/scores.csv")
ATLAS = rows("atlas_cells.csv")
SIGS = jload("states/signatures.json")
DUPS = jload("report/atlas_duplicates.json")

FIG = {f["name"]: f for f in FIGS["figures"]}
PAL = FIGS["style"]["palette"]
COL = {"hvg_pca": PAL["HVG-PCA"], "scgpt": PAL["scGPT"], "geneformer": PAL["Geneformer"],
       "prolif": PAL["Proliferation score"]}


# ----------------------------------------------------------------------------- formatting
def f3(x, d=3):
    return f"{x:.{d}f}".replace("-", MINUS)


def sgn(x, d=3):
    return f"{x:+.{d}f}".replace("-", MINUS)


def n_(x):
    return f"{int(x):,}"


def pv(p):
    """A p-value as printed on the page."""
    if p >= 0.1:
        return f"{p:.2f}"
    if round(p, 3) >= 0.01:
        return f"{p:.3f}"
    if p >= 0.00095:
        return f"{p:.2g}"
    m, e = f"{p:.1e}".split("e")
    return f"{m}&nbsp;×&nbsp;10<sup>{MINUS}{abs(int(e))}</sup>"


def ci(lo, hi, d=3):
    return f"{sgn(lo, d)} to {sgn(hi, d)}"


def pct(x):
    return f"{100 * x:.0f}%"


def esc(s):
    return html.escape(str(s), quote=False)


def check(cond, msg):
    """A claim the prose makes about the results; the build stops if the results disagree."""
    if not cond:
        raise AssertionError(f"page text no longer matches results/: {msg}")


def cite(*keys):
    return "(" + "; ".join(REFS[k][0] for k in keys) + ")"


# ----------------------------------------------------------------------------- derived facts
N_PAT, N_EV = L["n_patients"], L["n_events"]
check(L["endpoint"] == "OS", "ladder.json is the OS ladder")
REF = L["reference"]
PROLIF = REF["proliferation"]
FW = L["family_wise"]
P2C = L["p2_corrected"]
P3C = L["p3_corrected"]
VER = L["verdicts"]

# atlas
N_CELLS = len(ATLAS)
N_DONORS = len({r["donor_id"] for r in ATLAS})
N_DATASETS = len({r["dataset_id"] for r in ATLAS})
N_CELLTYPES = len({r["cell_type"] for r in ATLAS})
# the four datasets of the Chen et al. 2026 collection (CELLxGENE collection records)
CHEN = ("de5416ef", "ed880090", "5a9cfb44", "75011e96")
CHEN_FRAC = sum(r["dataset_id"].startswith(CHEN) for r in ATLAS) / N_CELLS
N_SIGS = {m: SSUM[m]["n_signatures"] for m in MODELS}
N_STATES = {m: sum(SSUM[m]["n_states"].values()) for m in MODELS}
UNIVERSE, UNIVERSE_HVG = SSUM["universe_n"], SSUM["universe_n_hvg"]
RESOLUTIONS = SPAR["args"]["resolutions"]
TOPK = SPAR["args"]["topk"]
MIN_CELLS = SPAR["args"]["min_cells"]
N_PERM = TPAR["args"]["n_perm"]
N_FLOOR = TPAR["args"]["n_floor"]
MIN_C = min(float(r["cindex"]) for r in SCORES)
check(len(SCORES) == sum(N_SIGS.values()), "scores.csv holds every signature")


def comp(model, res, cluster):
    for r in COMP:
        if r["model"] == model and float(r["resolution"]) == float(res) and r["cluster"] == str(cluster):
            return r
    raise KeyError((model, res, cluster))


def pick_rows(analysis):
    src = L["ladder"] if analysis == "primary" else L["sensitivity"]["ladder"]
    return {r["model"]: r for r in src}


PICKS = {a: pick_rows(a) for a in ("primary", "sensitivity")}


def pick_key(r):
    return f"{r['resolution']:.1f}/{r['signature']}"


def pick_id(r):
    """The id of a pick in results/replicate (direction suffix for protective picks)."""
    return f"{r['model']}:{pick_key(r)}" + (":protective" if r.get("direction") == "protective" else "")


def programme(analysis, model):
    lab = FIG["fig_states"]["programme_labels"][f"{analysis}/{model}"]
    check(lab["pick"] == pick_key(PICKS[analysis][model]), f"programme label of {analysis}/{model} is for its pick")
    return lab["programme"]


def genes(r):
    return SIGS[r["model"]][f"{r['resolution']:.1f}"][r["signature"]]["genes"]


def state_comp(r):
    cl = int(r["signature"].split("_")[0][1:])
    return comp(r["model"], r["resolution"], cl)


# the genes that name each programme, checked against the pick's signature at build time
PROGRAMME_GENES = {
    "primary/hvg_pca": ["CEACAM7", "TKTL1", "FAM9C", "INSL3"],
    "primary/scgpt": ["ERBB2", "EPCAM", "KRT7", "KRT19", "SPDEF", "TFAP2A"],
    "primary/geneformer": ["KRT5", "KRT14", "KRT17", "KRT6A", "KRT6B", "TACSTD2", "SFRP1"],
    "sensitivity/hvg_pca": ["GATA3", "XBP1", "PGR", "CA12", "STC2", "ANKRD30A", "AZGP1", "TRPS1"],
    "sensitivity/scgpt": ["GATA3", "XBP1", "TFF3", "AGR3", "SCGB2A2", "AZGP1", "TRPS1"],
    "sensitivity/geneformer": ["MALAT1", "NEAT1", "XIST", "JUND", "FOSB", "PPP1R15A"],
}
for key, gl in PROGRAMME_GENES.items():
    a, m = key.split("/")
    missing = [g for g in gl if g not in genes(PICKS[a][m])]
    check(not missing, f"{key} signature contains {missing}")
T_CELL_GENES = ["CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "PTPRC"]
check(not set(T_CELL_GENES) & set(genes(PICKS["primary"]["hvg_pca"])), "baseline primary pick has no canonical T-cell marker")
BASE_EXTRA_GENES = {"treg": ["IKZF4"], "cycle": ["CDK2", "ERCC6L"]}
check(all(g in genes(PICKS["primary"]["hvg_pca"]) for gl in BASE_EXTRA_GENES.values() for g in gl),
      "baseline primary pick contains IKZF4, CDK2 and ERCC6L")


def is_ribo(g):
    return bool(re.match(r"^(RPL|RPS|MRPL|MRPS)\d|^EEF1|^MT-", g))


def top_markers(r, k=9):
    return [g for g in genes(r) if not is_ribo(g)][:k]


# added value (results/ladder/added_value.csv)
def av_row(analysis, model, base):
    out = [r for r in AV if r["analysis"] == analysis and r["model"] == model and r["base"] == base]
    check(len(out) == 1, f"one added-value row for {analysis}/{model}/{base}")
    return {k: (float(v) if re.fullmatch(r"-?[\d.]+(e-?\d+)?", v or "x") else v) for k, v in out[0].items()}


AV_PICKS = {(a, m, b): av_row(a, m, b) for a in ("primary", "sensitivity") for m in MODELS
            for b in ("age+stage", "age+stage+PAM50")}
AV_PROLIF = {b: av_row("reference", "proliferation", b) for b in ("age+stage", "age+stage+PAM50")}
AV_N = {b: (int(AV_PROLIF[b]["n"]), int(AV_PROLIF[b]["events"])) for b in AV_PROLIF}


def adds(r):
    """Nested-CV change in C above zero in every repeat."""
    return r["delta_cindex_cv_nested_lo"] > 0


# METABRIC
MSIG = {s["id"]: s for s in M["signatures"]}
MEND = M["endpoints"]
MCMP = M["comparisons"]


def mrep(sid, ep):
    return MSIG[sid]["endpoints"][ep]


# run log: the latest successful run of each stage
def stage_seconds(name, pfi=False):
    hits = [r for r in RUNLOG if r.get("rc") == 0 and not r.get("dry_run")
            and r["stage"].split("_")[-1] == name and (("PFI" in " ".join(r["cmd"])) == pfi)]
    hits = [r for r in hits if not any(a in " ".join(r["cmd"]) for a in ("--freeze", "--shuffle"))]
    check(hits, f"run_log.json has a run of {name}{' PFI' if pfi else ''}")
    return hits[-1]["seconds"]


def hm(sec):
    h, m = divmod(round(sec / 60), 60)
    return f"{h} h {m} min" if h else f"{m} min"


T_TIER_A = sum(stage_seconds(s, pfi) for s, pfi in [("translate", False), ("translate", True), ("ladder", False),
                                                    ("ladder", True), ("stratify", False), ("report", False)])
T_EMBED = stage_seconds("embed")
T_STATES = stage_seconds("states")
T_DATA = stage_seconds("data")
T_REPL = stage_seconds("replicate")
T_TIER_C = T_DATA + T_EMBED + T_STATES + T_TIER_A + T_REPL


def about(sec, step=5):
    m = max(step, step * round(sec / 60 / step))
    return f"about {m // 60} h {m % 60:02d} min" if m >= 60 else f"about {m} min"


# ----------------------------------------------------------------------------- references
REFS = {
    "ahlmann": ("Ahlmann-Eltze et al. 2025", "Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. <em>Nature Methods</em> 22:1657–1661 (2025).", "10.1038/s41592-025-02772-6"),
    "bendidi": ("Bendidi et al. 2024", "Bendidi I, Whitfield S, Kenyon-Dean K, Ben Yedder H, El Mesbahi Y, Noutahi E, Denton AK. Benchmarking transcriptomics foundation models for perturbation analysis: one PCA still rules them all. arXiv:2410.13956 (2024).", "arXiv:2410.13956"),
    "chenA": ("Chen A et al. 2026", "Chen A, Kroehling L, Ennis CS, Denis GV, Monti S. A highly resolved integrated single-cell atlas of human breast cancers. <em>NAR Genomics and Bioinformatics</em> 8:lqaf217 (2026).", "10.1093/nargab/lqaf217"),
    "cui": ("Cui et al. 2024", "Cui H, Wang C, Maan H, Pang K, Luo F, Duan N, Wang B. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. <em>Nature Methods</em> 21:1470–1480 (2024).", "10.1038/s41592-024-02201-0"),
    "curtis": ("Curtis et al. 2012", "Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. <em>Nature</em> 486:346–352 (2012).", "10.1038/nature10983"),
    "davidson": ("Davidson-Pilon 2019", "Davidson-Pilon C. lifelines: survival analysis in Python. <em>Journal of Open Source Software</em> 4:1317 (2019).", "10.21105/joss.01317"),
    "goldman": ("Goldman et al. 2020", "Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. <em>Nature Biotechnology</em> 38:675–678 (2020).", "10.1038/s41587-020-0546-8"),
    "harrell": ("Harrell et al. 1982", "Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA. Evaluating the yield of medical tests. <em>JAMA</em> 247:2543–2546 (1982).", "10.1001/jama.1982.03320430047030"),
    "kedzierska": ("Kedzierska et al. 2025", "Kedzierska KZ, Crawford L, Amini AP, Lu AX. Zero-shot evaluation reveals limitations of single-cell foundation models. <em>Genome Biology</em> 26:101 (2025).", "10.1186/s13059-025-03574-x"),
    "liuJ": ("Liu J et al. 2018", "Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. <em>Cell</em> 173:400–416.e11 (2018).", "10.1016/j.cell.2018.02.052"),
    "liuW": ("Liu W et al. 2026", "Liu W, Wang Q, Long L, Wang W. Leveraging single-cell foundation models for accurate survival outcome prediction. <em>Bioinformatics Advances</em> 6:vbag076 (2026).", "10.1093/bioadv/vbag076"),
    "luca": ("Luca et al. 2021", "Luca BA, Steen CB, Matusiak M, et al. Atlas of clinically distinct cell states and ecosystems across human solid tumors. <em>Cell</em> 184:5482–5496.e28 (2021).", "10.1016/j.cell.2021.09.014"),
    "nielsen": ("Nielsen et al. 2010", "Nielsen TO, Parker JS, Leung S, et al. A comparison of PAM50 intrinsic subtyping with immunohistochemistry and clinical prognostic factors in tamoxifen-treated estrogen receptor-positive breast cancer. <em>Clinical Cancer Research</em> 16:5222–5232 (2010).", "10.1158/1078-0432.CCR-10-1282"),
    "parker": ("Parker et al. 2009", "Parker JS, Mullins M, Cheang MCU, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. <em>Journal of Clinical Oncology</em> 27:1160–1167 (2009).", "10.1200/JCO.2008.18.1370"),
    "pereira": ("Pereira et al. 2016", "Pereira B, Chin SF, Rueda OM, et al. The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. <em>Nature Communications</em> 7:11479 (2016).", "10.1038/ncomms11479"),
    "roman": ("Roman et al. 2025", "Roman A, Johri S, Conci R, Van Allen EM, Elmarakeby H. Empirical evaluation of single-cell foundation models for predicting cancer outcomes. bioRxiv (2025; v2 2026).", "10.1101/2025.10.31.685892"),
    "sun": ("Sun et al. 2022", "Sun D, Guan X, Moran AE, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. <em>Nature Biotechnology</em> 40:527–538 (2022).", "10.1038/s41587-021-01091-3"),
    "tcga": ("Cancer Genome Atlas Network 2012", "Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. <em>Nature</em> 490:61–70 (2012).", "10.1038/nature11412"),
    "theodoris": ("Theodoris et al. 2023", "Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. <em>Nature</em> 618:616–624 (2023).", "10.1038/s41586-023-06139-9"),
    "tirosh": ("Tirosh et al. 2016", "Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. <em>Science</em> 352:189–196 (2016).", "10.1126/science.aad0501"),
    "traag": ("Traag et al. 2019", "Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. <em>Scientific Reports</em> 9:5233 (2019).", "10.1038/s41598-019-41695-z"),
    "venet": ("Venet et al. 2011", "Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. <em>PLoS Computational Biology</em> 7:e1002240 (2011).", "10.1371/journal.pcbi.1002240"),
    "wirapati": ("Wirapati et al. 2008", "Wirapati P, Sotiriou C, Kunkel S, et al. Meta-analysis of gene expression profiles in breast cancer: toward a unified understanding of breast cancer subtyping and prognosis signatures. <em>Breast Cancer Research</em> 10:R65 (2008).", "10.1186/bcr2124"),
    "wolf": ("Wolf et al. 2018", "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. <em>Genome Biology</em> 19:15 (2018).", "10.1186/s13059-017-1382-0"),
    "wuSZ": ("Wu et al. 2021", "Wu SZ, Al-Eryani G, Roden DL, et al. A single-cell and spatially resolved atlas of human breast cancers. <em>Nature Genetics</em> 53:1334–1347 (2021).", "10.1038/s41588-021-00911-1"),
}


# ----------------------------------------------------------------------------- svg helpers
def T(x, y, s, size=12, fill=INK, anchor="start", weight=400, mono=False, dy=0, italic=False, halo=False):
    fam = MONO if mono else SANS
    s = s.replace("& ", "&amp; ")
    st = f"font-family:{fam};font-size:{size}px;fill:{fill};font-weight:{weight}"
    if italic:
        st += ";font-style:italic"
    if halo:
        st += f";paint-order:stroke;stroke:{PANEL};stroke-width:5px;stroke-linejoin:round"
    return (f'<text x="{x}" y="{y + dy}" text-anchor="{anchor}" dominant-baseline="middle" '
            f'style="{st}">{s}</text>')


def box(x, y, w, h, fill=PANEL, stroke=RULE, r=3, sw=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>'


def arrow(x1, y1, x2, y2, stroke=MUTED, sw=1.2, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{sw}" '
            f'marker-end="url(#arr)"{d}/>')


def line(x1, y1, x2, y2, stroke=RULE, sw=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{sw}"{d}/>'


def svg(w, h, body, label):
    return (f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" xmlns="http://www.w3.org/2000/svg">'
            f'<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{MUTED}"/></marker></defs>{body}</svg>')


def cells_glyph(x, y, n=18, seed=3, r=2.2, fill=MUTED):
    """A scatter of small dots standing in for cells."""
    import random
    rnd = random.Random(seed)
    return "".join(f'<circle cx="{x + rnd.uniform(6, 54):.1f}" cy="{y + rnd.uniform(6, 34):.1f}" r="{r}" fill="{fill}" opacity=".8"/>'
                   for _ in range(n))


def grid_glyph(x, y, cols, rows_, step=14, size=9):
    return "".join(f'<rect x="{x + i * step}" y="{y + j * step}" width="{size}" height="{size}" rx="1" fill="{MUTED}" opacity=".55"/>'
                   for i in range(cols) for j in range(rows_))


# ----------------------------------------------------------------------------- schematics
def fig_cohorts():
    """Three cohort cards, each its own small SVG so they stack on a narrow screen."""
    ms_os = MEND["OS"]
    cols = [
        ("single-cell atlas", [f"{n_(N_CELLS)} cells, {N_DONORS} donors,", f"{N_DATASETS} datasets; cell states", "are defined here"],
         f"CELLxGENE Census {C.CENSUS_VERSION}", "atlas"),
        ("TCGA-BRCA", [f"{n_(N_PAT)} primary tumours,", f"{N_EV} deaths; picks are", "selected and tested here"],
         "UCSC Xena · OS, PFI · age, stage, PAM50", "bulk"),
        ("METABRIC", [f"{n_(ms_os['n'])} patients, {n_(ms_os['events'])}", "deaths; frozen signatures", "are re-scored here"],
         "cBioPortal datahub · OS, DSS · age, PAM50", "bulk"),
    ]
    cards = []
    for title, lines, foot, kind in cols:
        b = [box(1, 1, 298, 118)]
        if kind == "atlas":
            b.append(cells_glyph(7, 15, n=34, seed=2, fill=MUTED))
            b.append(cells_glyph(31, 41, n=24, seed=7, fill="#8A9298"))
        else:
            b.append(grid_glyph(17, 17, 5, 6, step=13, size=9))
        b.append(T(113, 25, title, 13, INK, weight=600))
        for i, s in enumerate(lines):
            b.append(T(113, 49 + i * 18, s, 11.5, MUTED))
        b.append(T(150, 136, foot, 10.5, MUTED, anchor="middle", mono=True))
        cards.append(svg(300, 148, "".join(b), f"{title}: {' '.join(lines)}"))
    return ('<div class="cardrow">' + "".join(f'<div class="card">{c}</div>' for c in cards) + '</div>'
            '<p class="cardnote">The cohorts are independent collections; the data stage also asserts that no atlas donor '
            'carries a TCGA barcode.</p>')


def fig_flow():
    """Atlas -> three representations -> states -> signatures -> bulk cohort -> ladder."""
    import random
    W, H = 980, 262
    b = []
    X = {"atlas": 20, "rep": 165, "states": 350, "sig": 535, "bulk": 680, "out": 850}
    rows_ = {"hvg_pca": 44, "scgpt": 114, "geneformer": 184}   # top of each 56-high row
    RH = 56

    def cy(y):
        return y + RH / 2
    mid = 44 + 70 + 28   # vertical centre of the middle row

    def head(x, w, s):
        return T(x + w / 2, 18, s, 10, MUTED, anchor="middle", mono=True, weight=600)
    b.append(head(X["atlas"], 110, "1 · ONE INPUT"))
    b.append(box(X["atlas"], 44, 110, 196))
    b.append(cells_glyph(X["atlas"] + 2, 56, n=36, seed=5, r=2.2, fill=MUTED))
    b.append(cells_glyph(X["atlas"] + 2, 106, n=30, seed=9, r=2.2, fill="#8A9298"))
    b.append(cells_glyph(X["atlas"] + 2, 156, n=24, seed=13, r=2.2, fill=MUTED))
    b.append(T(X["atlas"] + 55, 205, "breast tumour", 11.5, INK, anchor="middle", weight=600))
    b.append(T(X["atlas"] + 55, 221, "single-cell atlas", 11.5, INK, anchor="middle", weight=600))
    b.append(head(X["rep"], 150, "2 · THREE VIEWS"))
    subs = {"hvg_pca": "linear · unsupervised", "scgpt": "foundation model", "geneformer": "foundation model"}
    for m, y in rows_.items():
        b.append(box(X["rep"], y, 150, RH, stroke=COL[m], sw=1.5))
        b.append(f'<rect x="{X["rep"]}" y="{y}" width="5" height="{RH}" rx="2" fill="{COL[m]}"/>')
        b.append(T(X["rep"] + 16, y + 20, NAME[m], 12.5, INK, weight=600))
        b.append(T(X["rep"] + 16, y + 39, subs[m], 10, MUTED, mono=True))
        b.append(arrow(X["atlas"] + 110, mid, X["rep"] - 3, cy(y)))
    b.append(head(X["states"], 150, "3 · CELL STATES"))
    seeds = {"hvg_pca": 1, "scgpt": 2, "geneformer": 3}
    for m, y in rows_.items():
        b.append(box(X["states"], y, 150, RH))
        rnd = random.Random(seeds[m])
        for k in range(3):
            cx = X["states"] + 20 + k * 21
            for _ in range(7):
                b.append(f'<circle cx="{cx + rnd.uniform(-6, 6):.1f}" cy="{cy(y) + rnd.uniform(-11, 11):.1f}" r="2" fill="{COL[m]}" opacity=".8"/>')
        b.append(T(X["states"] + 86, cy(y) - 8, "Leiden", 10.5, MUTED))
        b.append(T(X["states"] + 86, cy(y) + 8, "clusters", 10.5, MUTED))
        b.append(arrow(X["rep"] + 150, cy(y), X["states"] - 3, cy(y)))
    b.append(head(X["sig"], 100, "4 · MARKER GENES"))
    for m, y in rows_.items():
        b.append(box(X["sig"], y, 100, RH))
        for i, wl in enumerate([68, 50, 60, 38]):
            b.append(f'<rect x="{X["sig"] + 16}" y="{y + 11 + i * 10}" width="{wl}" height="4" rx="2" fill="{COL[m]}" opacity="{.95 - i * .18:.2f}"/>')
        b.append(arrow(X["states"] + 150, cy(y), X["sig"] - 3, cy(y)))
    b.append(head(X["bulk"], 120, "5 · BULK COHORT"))
    b.append(box(X["bulk"], 44, 120, 196))
    b.append(grid_glyph(X["bulk"] + 17, 58, 6, 6, step=15, size=10))
    b.append(T(X["bulk"] + 60, 168, "TCGA breast", 11.5, INK, anchor="middle", weight=600))
    b.append(T(X["bulk"] + 60, 184, "tumours", 11.5, INK, anchor="middle", weight=600))
    b.append(T(X["bulk"] + 60, 208, f"{n_(N_PAT)} patients", 10, MUTED, anchor="middle", mono=True))
    b.append(T(X["bulk"] + 60, 223, f"{N_EV} deaths", 10, MUTED, anchor="middle", mono=True))
    for y in rows_.values():
        b.append(arrow(X["sig"] + 100, cy(y), X["bulk"] - 3, mid))
    b.append(head(X["out"], 110, "6 · READ"))
    b.append(box(X["out"], 44, 110, 196))
    for lab, yy, c, dash in [("reference", 80, INK, "6 3"), ("test", 120, COL["scgpt"], ""),
                             ("baseline", 160, COL["hvg_pca"], ""), ("floor", 200, MUTED, "2 3")]:
        b.append(T(X["out"] + 55, yy - 12, lab, 10, MUTED, anchor="middle", mono=True))
        b.append(line(X["out"] + 18, yy, X["out"] + 92, yy, stroke=c, sw=2.2 if not dash else 1.5, dash=dash))
    b.append(arrow(X["bulk"] + 120, mid, X["out"] - 3, mid))
    return svg(W, H, "".join(b), "Design flow from one atlas through three representations to survival in the bulk cohort")


def fig_ladder_schematic():
    W, H = 980, 262
    x0, x1, ay = 60, 920, 196

    def X(c):
        return x0 + (c - 0.5) / 0.3 * (x1 - x0)
    b = [line(x0, ay, x1, ay, stroke=INK, sw=1.2)]
    for c in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8):
        b.append(line(X(c), ay, X(c), ay + 6, stroke=INK))
        b.append(T(X(c), ay + 20, f"{c:.2f}", 11, MUTED, anchor="middle", mono=True))
    b.append(T((x0 + x1) / 2, ay + 46, "concordance index: of two patients, how often the higher-scored one dies first (0.50 = chance)", 12, MUTED, anchor="middle"))
    rungs = [("null", 0.53, AXIS, "6 4", "outcomes permuted"),
             ("floor", 0.57, MUTED, "3 3", "matched random genes"),
             ("baseline", 0.63, COL["hvg_pca"], "", "HVG-PCA states"),
             ("test", 0.69, COL["scgpt"], "", "scGPT / Geneformer"),
             ("reference", 0.78, INK, "8 4", "clinical information")]
    for lab, c, col, dash, desc in rungs:
        x = X(c)
        b.append(T(x, 16, lab.upper(), 10.5, MUTED if col == AXIS else col, anchor="middle", mono=True, weight=600))
        b.append(T(x, 34, desc, 10.5, MUTED, anchor="middle"))
        b.append(line(x, 48, x, ay - 1, stroke=col, sw=2.4 if not dash else 1.6, dash=dash))
    p1y, p2y = 92, 150
    b.append(f'<line x1="{X(0.53) + 3}" y1="{p1y}" x2="{X(0.69) - 3}" y2="{p1y}" stroke="{INK}" stroke-width="1.6" marker-end="url(#arr)" marker-start="url(#arr)"/>')
    b.append(T(X(0.69) + 10, p1y, "P1 · clears null and floor", 11, INK, mono=True, weight=500, halo=True))
    b.append(f'<line x1="{X(0.63) + 3}" y1="{p2y}" x2="{X(0.69) - 3}" y2="{p2y}" stroke="{INK}" stroke-width="1.6" marker-end="url(#arr)" marker-start="url(#arr)"/>')
    b.append(T(X(0.69) + 10, p2y, "P2 · margin over baseline", 11, INK, mono=True, weight=600, halo=True))
    return svg(W, H, "".join(b), "Schematic ladder of rungs: null, floor, baseline, test, reference")


def fig_shared_input():
    import random
    W, H = 980, 206
    b = [box(30, 14, 150, 132, fill="#F2F4F5")]
    rnd = random.Random(11)
    for i in range(10):
        for j in range(9):
            b.append(f'<rect x="{40 + i * 13}" y="{22 + j * 13}" width="11" height="11" fill="{INK}" opacity="{rnd.choice([.1, .2, .35, .55, .8])}"/>')
    b.append(T(105, 165, "the same cells × genes", 11.5, INK, anchor="middle", weight=600))
    b.append(T(105, 182, "expression matrix", 11.5, INK, anchor="middle", weight=600))
    for m, y in {"hvg_pca": 20, "scgpt": 70, "geneformer": 120}.items():
        b.append(arrow(180, 80, 300, y + 15))
        b.append(box(304, y, 160, 30, stroke=COL[m], sw=1.5))
        b.append(T(384, y + 15, NAME[m], 12, INK, anchor="middle", weight=600))
        b.append(arrow(464, y + 15, 570, y + 15))
        b.append(box(574, y, 120, 30))
        b.append(T(634, y + 15, "signature", 11, MUTED, anchor="middle", mono=True))
        b.append(arrow(694, y + 15, 760, 80))
    b.append(box(764, 50, 190, 60, fill="#F2F4F5"))
    b.append(T(859, 71, "scored in the same", 11.5, INK, anchor="middle"))
    b.append(T(859, 89, "bulk tumour genes", 11.5, INK, anchor="middle"))
    return svg(W, H, "".join(b), "All three representations are built from one expression matrix and scored on the same genes")


def fig_outcomes():
    """The three outcome shapes of DESIGN §7, each its own small SVG; the TCGA result lies between shapes 2 and 3."""
    scen = [("1 · Foundation model ahead", {"floor": .56, "base": .62, "test": .70, "ref": .74}, "test above baseline, both above floor"),
            ("2 · Equal, below the reference", {"floor": .56, "base": .64, "test": .65, "ref": .74}, "test ≈ baseline, both far below the reference"),
            ("3 · Both fail", {"floor": .56, "base": .575, "test": .585, "ref": .74}, "test ≈ baseline ≈ floor")]
    cards = []
    for i, (title, v, desc) in enumerate(scen):
        hit = i in (1, 2)
        b = [box(1, 1, 288, 128, stroke=INK if hit else RULE, sw=1.6 if hit else 1, dash="" if hit else "")]
        b.append(T(16, 21, title, 13, INK, weight=600))
        ay = 91
        b.append(line(20, ay, 270, ay, stroke=INK))

        def Xc(c):
            return 20 + (c - .5) / .3 * 250
        for c, col, dash, sw_ in [(v["floor"], MUTED, "3 3", 1.6), (v["base"], COL["hvg_pca"], "", 2.4),
                                  (v["test"], COL["scgpt"], "", 2.4), (v["ref"], INK, "8 4", 1.6)]:
            b.append(line(Xc(c), 37, Xc(c), ay, stroke=col, sw=sw_, dash=dash))
        b.append(T(16, 112, desc, 11, MUTED))
        cards.append(svg(290, 130, "".join(b), f"Outcome {title}: {desc}"))
    legend = "".join(
        f'<span class="lg"><svg viewBox="0 0 28 10" width="28" height="10" aria-hidden="true">{line(0, 5, 28, 5, stroke=col, sw=sw_, dash=dash)}</svg>{lab}</span>'
        for lab, col, dash, sw_ in [("floor", MUTED, "3 3", 1.6), ("baseline", COL["hvg_pca"], "", 2.4),
                                    ("test", COL["scgpt"], "", 2.4), ("reference", INK, "8 4", 1.6)])
    return ('<div class="cardrow">' + "".join(f'<div class="card">{c}</div>' for c in cards) + '</div>'
            f'<p class="cardnote legend">{legend}</p>'
            '<p class="cardnote">Outlined: the TCGA result lies between shapes 2 and 3. Test ≈ baseline, both above their '
            'random-gene floors but not beyond the family-wise null, both far below age + stage.</p>')


# ----------------------------------------------------------------------------- page pieces
def data_uri(p: Path):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


def figure(name, label):
    """A result figure from results/report with its caption from figures.json."""
    meta = FIG[name]
    p = ROOT / meta["file"]
    reads = meta["how_to_read"]
    reads = reads if isinstance(reads, list) else [reads]
    cap = (f'<div class="rtext"><p><strong>{esc(meta["question"])}</strong> {esc(meta["shows"])}</p>'
           f'<ul>{"".join(f"<li>{esc(r)}</li>" for r in reads)}</ul></div>')
    if p.exists() and meta.get("produced", True):
        w, h = meta["pixels"]
        img = (f'<div class="imgscroll"><img src="{data_uri(p)}" width="{w}" height="{h}" '
               f'alt="{esc(meta["title"])}"></div>')
    else:
        img = (f'<p class="missing">{esc(meta["file"])} is not in this build; <code>scfm run report</code> writes it '
               f'when its inputs are present.</p>')
    return (f'<figure class="result" id="{name}"><figcaption><span class="fk">{label}</span> {esc(meta["title"])}</figcaption>'
            f'{img}<p class="swipe">Wide figure: scroll sideways to see all of it.</p>{cap}</figure>')


def panel(fig_html, caption, kind="scroll"):
    if kind == "cards":
        return f'<figure class="panel fig"><figcaption>{caption}</figcaption>{fig_html}</figure>'
    return (f'<figure class="panel fig"><figcaption>{caption}</figcaption><div class="svgscroll">{fig_html}</div>'
            f'<p class="swipe">Scroll sideways to see all of it.</p></figure>')


def sw(m):
    return f'<span class="swatch" style="background:{COL[m]}"></span>'


def yesno(b, yes="yes", no="no"):
    return f'<span class="pill {"on" if b else "off"}">{yes if b else no}</span>'


# ----------------------------------------------------------------------------- claims (asserted)
# P1, pre-specified direction
fw_prim = {m: FW["primary"][m] for m in MODELS}
fw_sens = {m: FW["sensitivity"][m] for m in MODELS}
for fwd in (fw_prim, fw_sens):
    for m in MODELS:
        check(not fwd[m]["P1_corrected"] or fwd[m]["fw_p"] < C.ALPHA, f"P1_corrected of {m} matches its family-wise p")
NONE_CLEARS_PRIM = not any(fw_prim[m]["P1_corrected"] for m in MODELS)
check(NONE_CLEARS_PRIM, "no representation clears its family-wise null in the pre-specified direction")
for m in FMS:
    check(VER["primary"][m]["P1_corrected"] == fw_prim[m]["P1_corrected"], f"verdict P1 of {m}")
SENS_CLEAR = [m for m in MODELS if fw_sens[m]["P1_corrected"]]
SENS_NOT = [m for m in MODELS if not fw_sens[m]["P1_corrected"]]
check(set(SENS_CLEAR) == {"hvg_pca", "scgpt"}, "either direction: HVG-PCA and scGPT clear their nulls")
check(SENS_NOT == ["geneformer"], "either direction: Geneformer does not clear its null")

# P2
P2_ALL_FAIL = not any(P2C[a][m]["P2_pass"] for a in ("primary", "sensitivity") for m in FMS)
check(P2_ALL_FAIL, "P2 fails for both foundation models in both readings")
check(all(P2C["primary"][m]["margin"] < 0 for m in FMS), "pre-specified P2 margins are negative")
for a in ("primary", "sensitivity"):
    for m in FMS:
        check(VER[a][m]["P2_corrected"] == P2C[a][m]["P2_pass"], f"verdict P2 of {a}/{m}")
P2_BELOW = [m for m in FMS if P2C["primary"][m]["ci"][1] < 0]

# added value
for a in ("primary", "sensitivity"):
    for m in MODELS:
        check(not adds(AV_PICKS[(a, m, "age+stage")]), f"{a}/{m} pick does not add to age + stage (nested)")
check(all(AV_PICKS[(a, m, "age+stage")]["delta_cindex_cv_nested_lo"] <= 0 <= AV_PICKS[(a, m, "age+stage")]["delta_cindex_cv_nested_hi"]
          for a in ("primary", "sensitivity") for m in MODELS), "every pick's nested change in C spans zero")
check(adds(AV_PROLIF["age+stage"]), "the proliferation score adds to age + stage (nested)")
check(not adds(AV_PROLIF["age+stage+PAM50"]), "with PAM50 the proliferation score does not add")
NESTED_PICK = [AV_PICKS[(a, m, "age+stage")]["delta_cindex_cv_nested_mean"] for a in ("primary", "sensitivity") for m in MODELS]
FAMMAX_ADDS = [r for r in AV if r["analysis"].startswith("family_max") and r["base"] == "age+stage"
               and float(r["delta_cindex_cv_nested_lo"]) > 0]
FAMMAX_ADDS_MODELS = sorted({(r["analysis"], r["model"]) for r in FAMMAX_ADDS})
check(FAMMAX_ADDS_MODELS == [("family_max_sensitivity", "hvg_pca")],
      "among family maxima only HVG-PCA's either-direction maximum adds")
FAMMAX_HVG = [r for r in AV if r["analysis"] == "family_max_sensitivity" and r["model"] == "hvg_pca"]
FAMMAX_HVG_AS = float([r for r in FAMMAX_HVG if r["base"] == "age+stage"][0]["delta_cindex_cv_nested_mean"])
FAMMAX_HVG_PAM = [r for r in FAMMAX_HVG if r["base"] == "age+stage+PAM50"][0]
check(float(FAMMAX_HVG_PAM["delta_cindex_cv_nested_lo"]) <= 0, "HVG-PCA family maximum does not add with PAM50")

# stratify
SM = STR["models"]
CLIN_CV = SM["clinical"]["cindex_mean"]
check(all(SM[f"xgb_{m}"]["cindex_mean"] < CLIN_CV for m in MODELS),
      "gradient-boosted models with any representation's states fall below age + stage alone")
XGB_RANGE = (min(SM[f"xgb_{m}"]["cindex_mean"] for m in MODELS), max(SM[f"xgb_{m}"]["cindex_mean"] for m in MODELS))
check(SM["clinical_pam50"]["cindex_mean"] <= SM["clinical_on_pam50_set"]["cindex_mean"],
      "PAM50 subtype adds nothing to age + stage on the PAM50 set")
# the same models under the fold seeds of different sections
AGE_STAGE_CVS = [REF["clinical_full"]["cindex_cv_mean"], CLIN_CV, float(AV_PROLIF["age+stage"]["cindex_base_cv"])]
PAM_CLIN_CVS = [REF["pam50_clinical"]["cindex_cv_mean"], SM["clinical_pam50"]["cindex_mean"]]
check(all(abs(SM[f"ridge_{m}"]["cindex_mean"] - CLIN_CV) < 0.01 and SM[f"ridge_{m}"]["vs_clinical"]["min"] <= 0
          for m in MODELS), "no ridge model improves on age + stage in every repeat")

# METABRIC
check(M["outcomes"] == "real", "metabric.json is the real-outcome run")
MARGINS = [c["margin"] for c in MCMP]
check(len(MCMP) == 8 and all(x < 0 for x in MARGINS), "every METABRIC foundation-model-minus-baseline margin is negative")
N_BELOW = sum(c["ci"][1] < 0 for c in MCMP)
CMP_PRIM = [c for c in MCMP if c["analysis"] == "primary"]
CMP_SENS = [c for c in MCMP if c["analysis"] == "sensitivity"]
BELOW_PRIM = [c for c in CMP_PRIM if c["ci"][1] < 0]
check(len(CMP_PRIM) == 4 and len(CMP_SENS) == 4, "four pre-specified and four post hoc METABRIC comparisons")
check(all(c["ci"][1] < 0 for c in CMP_SENS), "every post hoc METABRIC interval lies below zero")
BELOW_PRIM_TXT = ", ".join(f"{NAME[c['fm'].split(':')[0]]}, {c['endpoint']}" for c in BELOW_PRIM)
PROLIF_ID = "reference:proliferation"


def gap(sid, ep):
    e = mrep(sid, ep)
    return e["cindex_oriented"] - e["floor_p95"]


def floor_se(sid, ep):
    """Approximate Monte-Carlo SE of a floor p95 from n random sets (normal approximation from mean and p95)."""
    e = mrep(sid, ep)
    sd = (e["floor_p95"] - e["floor_mean"]) / 1.645
    return (0.95 * 0.05 / FROZ["n_floor_sets"]) ** 0.5 / 0.10314 * sd


def floor_borderline(sid, ep):
    return abs(gap(sid, ep)) < 2 * floor_se(sid, ep)


for ep in ("OS", "DSS"):
    check(mrep(PROLIF_ID, ep)["replicates"], f"proliferation replicates on {ep}")
    check(max(MSIG, key=lambda s: gap(s, ep) if MSIG[s]["scored"] else -1) == PROLIF_ID,
          f"proliferation has the largest margin over its METABRIC floor on {ep}")
PID = {a: {m: pick_id(PICKS[a][m]) for m in MODELS} for a in PICKS}
for a in PICKS:
    for m in MODELS:
        check(PID[a][m] in MSIG, f"{PID[a][m]} is frozen and scored")
HVG_REPL = all(mrep(PID[a]["hvg_pca"], ep)["replicates"] for a in PICKS for ep in ("OS", "DSS"))
check(HVG_REPL, "both HVG-PCA picks replicate on OS and DSS")
SC_SENS = PID["sensitivity"]["scgpt"]
GF_SENS = PID["sensitivity"]["geneformer"]
check(not any(mrep(SC_SENS, ep)["replicates"] for ep in ("OS", "DSS")), "scGPT protective pick does not replicate overall")
SC_P3 = MSIG[SC_SENS]["p3"]
check(list(SC_P3) == ["LumA"] and SC_P3["LumA"]["OS"]["perm_p"] < C.ALPHA,
      "scGPT protective pick: within luminal A the outcome-permutation p is below 0.05")
check(all("floor" not in k for k in SC_P3["LumA"]["OS"]), "no random-gene floor was computed within subtypes")
check(not any(mrep(GF_SENS, ep)["replicates"] for ep in ("OS", "DSS")), "Geneformer protective pick does not replicate")
GF_P3 = MSIG[GF_SENS]["p3"]
check(all(v["OS"]["perm_p"] >= C.ALPHA for v in GF_P3.values()), "Geneformer protective pick fails within its subtype")
HIGH_HR_NO = [(sid, ep) for sid, s in MSIG.items() if s["scored"] for ep in ("OS", "DSS")
              if s["endpoints"][ep]["p"] < 1e-5 and not s["endpoints"][ep]["replicates"]]
FLOOR_RANGE = (min(s["endpoints"][ep]["floor_p95"] for s in MSIG.values() if s["scored"] for ep in ("OS", "DSS")),
               max(s["endpoints"][ep]["floor_p95"] for s in MSIG.values() if s["scored"] for ep in ("OS", "DSS")))

# P3
P3_PASS = {k: v["P3_pass"] for k, v in P3C["picks"].items()}
check(not any(P3_PASS[f"primary/{m}"] for m in MODELS), "no pre-specified pick passes P3")
check([k for k, v in P3_PASS.items() if v] == ["sensitivity/hvg_pca"], "only the baseline's sensitivity pick passes P3")
for m in FMS:
    check(VER["primary"][m]["P3_corrected"] == P3_PASS[f"primary/{m}"], f"verdict P3 of {m}")
    check(VER["sensitivity"][m]["P3_corrected"] == P3_PASS[f"sensitivity/{m}"], f"sensitivity verdict P3 of {m}")
check(P3C["picks"]["sensitivity/scgpt"]["borderline"], "scGPT's luminal A P3 p is borderline")

# outcome shapes of DESIGN §7, read against the TCGA ladder
PICKS_ABOVE_FLOOR = all(r["cindex"] > r["floor_mean"] for a in PICKS for r in PICKS[a].values())
check(PICKS_ABOVE_FLOOR and NONE_CLEARS_PRIM, "every pick is above its floor and no pre-specified family clears its null")

# composition flags
SINGLE_DONOR_PICKS = [(a, m) for a in PICKS for m in MODELS if state_comp(PICKS[a][m])["single_donor"] == "True"]
check(SINGLE_DONOR_PICKS == [("primary", "hvg_pca"), ("primary", "geneformer"), ("sensitivity", "hvg_pca")],
      "the single-donor picks are the baseline's two and Geneformer's pre-specified pick")
BASE_PRIM_COMP = state_comp(PICKS["primary"]["hvg_pca"])
check(BASE_PRIM_COMP["n_datasets"] == "1", "the baseline's pre-specified pick comes from one dataset")
GF_SENS_COMP = state_comp(PICKS["sensitivity"]["geneformer"])
check(float(GF_SENS_COMP["top_cell_type_frac"]) < 0.5, "Geneformer's protective pick is mostly not one cell type")
KEPT = [r for r in COMP if r["kept"] == "True"]
SINGLE_DONOR_STATES = {m: (sum(r["single_donor"] == "True" for r in KEPT if r["model"] == m),
                           sum(r["model"] == m for r in KEPT)) for m in MODELS}

P3_NPERM = P3C["n_perm"]
check(P3_NPERM == L["monte_carlo"]["p3_permutations"], "P3 permutation count agrees within ladder.json")


DEV_TERMS = [("FM ", "foundation-model "), ("FM-derived", "foundation-model-derived"),
             ("Monte-Carlo SE", "Monte-Carlo standard error (SE)"), ("pre-registered", "pre-specified"),
             ("pre-specified proliferation reference", "published proliferation reference"),
             ("clinical_full (age + stage, all patients with both)", "age + stage on all patients with both"),
             ("P1_corrected = fw_p < 0.05", "P1 passes if the family-wise p < 0.05"),
             ("P2_corrected = lower", "P2 passes if the lower"), ("P3_corrected = p <", "P3 passes if p <"),
             ("unless P2_corrected passes", "unless corrected P2 passes"), ("hvg_frac < 0.5", "HVG fraction < 0.5"),
             (" (baseline_rank_p)", ""),
             (" (clinical, pam50, pam50_clinical)", ""), ("reference_pam50", "a separate field"),
             (" (selected_on_outcome)", ""), ("lrt_p_nominal", "the likelihood-ratio p"),
             ("P2_point_as_coded", "a separate field"), ("P2_margin_ci_as_coded", "the as-coded P2 interval"),
             ("P2_margin_as_coded", "the as-coded P2 margin")]


def fix_dev_text(s):
    """Deviation text from ladder.json, with field names glossed for readers."""
    m = re.search(r"against ([\d,]+) within-subtype outcome permutations", s)
    check(not m or int(m.group(1).replace(",", "")) == P3_NPERM, "deviation text gives the P3 permutation count")
    for a, b in DEV_TERMS:
        s = s.replace(a, b)
    return s


# ----------------------------------------------------------------------------- sections
def header_html(today):
    pp = {m: pv(fw_prim[m]["fw_p"]) for m in MODELS}
    m2 = {m: sgn(P2C["primary"][m]["margin"]) for m in FMS}
    prol = AV_PROLIF["age+stage"]["delta_cindex_cv_nested_mean"]
    h1 = "Foundation-model cell states predict breast-cancer survival no better than a linear baseline, in TCGA and in METABRIC"
    check(P2_ALL_FAIL and all(x < 0 for x in MARGINS), "headline")
    commit = f'<a href="{REPO}/commit/8d2c2d0b7391f48f22ff92be15fff72b8a76cd50">commit&nbsp;8d2c2d0</a>'
    abstract = (
        f"Cell states were defined in {n_(N_CELLS)} cells from {N_DONORS} donors of a breast-cancer single-cell atlas with two "
        f"single-cell foundation models, scGPT and Geneformer, and with a linear baseline, HVG-PCA: principal component "
        f"analysis (PCA) of the highly variable genes (HVGs). Each state's marker genes were scored as a signature in {n_(N_PAT)} tumours of The "
        f"Cancer Genome Atlas breast cohort (TCGA-BRCA; {N_EV} deaths) and tested for overall survival (OS) against outcome "
        f"permutations, random gene sets matched for size and expression, and clinical references: age, stage and the PAM50 "
        f"(Prediction Analysis of Microarray, 50 genes) intrinsic subtype. In the pre-specified risk "
        f"direction no representation's best signature cleared its family-wise permutation null (p = {pp['hvg_pca']} HVG-PCA, "
        f"{pp['scgpt']} scGPT, {pp['geneformer']} Geneformer), and each foundation model's margin over the baseline in "
        f"concordance index (C-index) above the random-gene floor was negative (scGPT {m2['scgpt']}, Geneformer {m2['geneformer']}), "
        f"measured against a baseline pick that is a single-dataset state with a low floor. "
        f"Added to age and stage, no pick changed the cross-validated C-index by more than {sgn(min(NESTED_PICK))} to "
        f"{sgn(max(NESTED_PICK))} under nested selection (only HVG-PCA's either-direction family maximum added, "
        f"{sgn(FAMMAX_HVG_AS)}); the published 11-gene proliferation score of PAM50 added {sgn(prol)}. With the signatures "
        f"frozen before outcomes were read, all {len(MCMP)} foundation-model-minus-baseline margins were negative in "
        f"{n_(MEND['OS']['n'])} patients of the Molecular Taxonomy of Breast Cancer International Consortium (METABRIC); the "
        f"interval lies below zero for {WORD[len(BELOW_PRIM)]} of the {WORD[len(CMP_PRIM)]} pre-specified comparisons "
        f"({BELOW_PRIM_TXT}) and for all {WORD[len(CMP_SENS)]} post hoc ones.")
    return f"""<header>
<p class="eyebrow"><span>Single-cell foundation models · breast cancer · design committed before data ({commit})</span></p>
<h1>{h1}</h1>
<p class="byline"><span class="nw">Natalija Stepurko</span> · <span class="nw">{today:%B %Y}</span> · <a class="nw" href="{SITE}">natalija-stepurko.com</a> · <a href="{REPO}">code</a></p>
<p class="abstract">{abstract}</p>
</header>"""


def result_html():
    pp = {m: pv(fw_prim[m]["fw_p"]) for m in MODELS}
    ps = {m: pv(fw_sens[m]["fw_p"]) for m in MODELS}
    kp = {m: pv(fw_prim[m]["pick_fw_p"]) for m in MODELS}
    ks = {m: pv(fw_sens[m]["pick_fw_p"]) for m in MODELS}
    check(not any(fw[m]["pick_fw_p"] < C.ALPHA for fw in (fw_prim, fw_sens) for m in MODELS),
          "no margin-selected pick clears the family-wise null")
    m2 = {m: sgn(P2C["primary"][m]["margin"]) for m in FMS}
    ci2 = {m: ci(*P2C["primary"][m]["ci"]) for m in FMS}
    near0 = [m for m in P2_BELOW if P2C["primary"][m]["ci"][1] > -0.005]
    below_txt = "".join(f"; {NAME[m]}'s interval lies below zero, with its upper bound {sgn(P2C['primary'][m]['ci'][1])} "
                        f"near zero" if m in near0 else f"; {NAME[m]}'s interval lies below zero" for m in P2_BELOW)
    borderline = [m for m in MODELS if fw_sens[m]["fw_p_borderline"]]
    bl = (f" ({' and '.join(NAME[m] for m in borderline)} within two Monte-Carlo standard errors of 0.05)" if borderline else "")
    bp, sp_ = PICKS["primary"]["hvg_pca"], PICKS["primary"]["scgpt"]
    check(sp_["cindex"] > bp["cindex"] and sp_["floor_mean"] > bp["floor_mean"], "scGPT's pick: higher C and higher floor")
    ms = {m: sgn(P2C["sensitivity"][m]["margin"]) for m in FMS}
    check(all(P2C["sensitivity"][m]["ci"][0] < 0 < P2C["sensitivity"][m]["ci"][1] for m in FMS),
          "sensitivity margins span zero")
    clin = REF["clinical_full"]
    pam = REF["pam50"]
    prol = AV_PROLIF["age+stage"]
    os_p = mrep(PROLIF_ID, "OS")
    ds_p = mrep(PROLIF_ID, "DSS")
    lumA = SC_P3["LumA"]["OS"]
    pam_set, pam_base = SM["clinical_pam50"], SM["clinical_on_pam50_set"]
    gf_dss = PID["primary"]["geneformer"]
    check([(a, m, ep) for a in PICKS for m in FMS for ep in ("OS", "DSS") if mrep(PID[a][m], ep)["replicates"]]
          == [("primary", "geneformer", "DSS")], "the only foundation-model pick that replicates is Geneformer's risk pick on DSS")
    check(floor_borderline(gf_dss, "DSS"), "Geneformer's DSS replication is within Monte-Carlo error of its floor p95")
    f1 = (f"<p><strong>No representation beats chance selection in the pre-specified direction.</strong> The largest C-index in "
          f"each representation's family of {min(N_SIGS.values())}–{max(N_SIGS.values())} signatures stays inside its "
          f"family-wise null (p = {pp['hvg_pca']} HVG-PCA, {pp['scgpt']} scGPT, {pp['geneformer']} Geneformer). Read in the "
          f"direction each signature acts (post hoc), the family maxima of HVG-PCA (p = {ps['hvg_pca']}) and scGPT "
          f"(p = {ps['scgpt']}) clear their nulls and Geneformer's (p = {ps['geneformer']}) does not{bl}. The margin-selected "
          f"picks shown in the figures do not clear the null in either reading (p = {ks['hvg_pca']}, {ks['scgpt']} and "
          f"{ks['geneformer']} either direction; {kp['hvg_pca']}, {kp['scgpt']} and {kp['geneformer']} pre-specified). "
          f"No foundation-model pick is further above its random-gene floor than the baseline's pick; in the pre-specified "
          f"reading the margins are {m2['scgpt']} for scGPT (95% interval {ci2['scgpt']}) and {m2['geneformer']} for "
          f"Geneformer ({ci2['geneformer']}){below_txt}. scGPT's pick has the higher C-index ({f3(sp_['cindex'])} against "
          f"{f3(bp['cindex'])}) and also the higher floor ({f3(sp_['floor_mean'])} against {f3(bp['floor_mean'])}); the "
          f"baseline pick it is measured against is the single-dataset state described below and does not itself beat "
          f"chance (p = {kp['hvg_pca']}). Against the baseline's luminal sensitivity pick the margins are {ms['scgpt']} and "
          f"{ms['geneformer']}, with intervals spanning zero.</p>")
    f2 = (f"<p><strong>No pick adds to age and stage.</strong> Age and stage alone reach a "
          f"cross-validated C-index of {f3(clin['cindex_cv_mean'])} ({n_(clin['n'])} patients); the PAM50 subtype call alone "
          f"reaches {f3(pam['cindex_cv_mean'])} ({n_(pam['n'])} patients). With the signature re-selected inside every training "
          f"fold, adding a pick changes the cross-validated C-index by {sgn(min(NESTED_PICK), 4)} to {sgn(max(NESTED_PICK), 4)}, "
          f"and every pick's change spans zero across the five repeats; the published 11-gene proliferation score adds "
          f"{sgn(prol['delta_cindex_cv_nested_mean'])}. In this cohort PAM50 subtype adds nothing to age and stage either "
          f"({f3(pam_set['cindex_mean'])} against {f3(pam_base['cindex_mean'])} on the {n_(pam_set['n'])} patients with a "
          f"call), so the test has little room to show added value.</p>")
    f3_ = (f"<p><strong>The replication in METABRIC points the same way.</strong> With signatures and directions frozen before "
           f"its outcomes were read, every foundation-model-minus-baseline margin is negative (intervals below zero: "
           f"{len(BELOW_PRIM)} of {len(CMP_PRIM)} pre-specified, {len(CMP_SENS)} of {len(CMP_SENS)} post hoc). The published "
           f"proliferation score replicates most strongly (OS C = {f3(os_p['cindex_oriented'])} against a floor 95th "
           f"percentile of {f3(os_p['floor_p95'])}; disease-specific survival (DSS) C = {f3(ds_p['cindex_oriented'])}). Both "
           f"HVG-PCA picks replicate. Of the foundation-model picks only Geneformer's risk pick meets the criterion, on DSS "
           f"and by {sgn(gap(gf_dss, 'DSS'), 4)} (borderline). Within luminal A tumours scGPT's protective pick has "
           f"C = {f3(lumA['cindex_oriented'])} (outcome-permutation p = {pv(lumA['perm_p'])}; no random-gene floor, and not a "
           f"pre-specified replication criterion); Geneformer's protective pick does not replicate.</p>")
    return f"""<h2 id="result">Key result</h2>
{figure("fig_ladder", "Figure 1")}
<div class="findings">{f1}{f2}{f3_}</div>"""


def states_html():
    rows_ = []
    for a in ("primary", "sensitivity"):
        for m in MODELS:
            r = PICKS[a][m]
            c = state_comp(r)
            key = f"{a}/{m}"
            direction = r.get("direction", "risk")
            prog_genes = ", ".join(PROGRAMME_GENES[key])
            k = r["signature"].split("_k")[1]
            rows_.append(
                f'<tr><td>{sw(m)}<strong>{NAME[m]}</strong><br><span class="dim2">{"pre-specified" if a == "primary" else "sensitivity"}'
                f' · <span class="mono">{r["resolution"]:.1f}/{r["signature"]}</span></span></td>'
                f'<td><strong>{esc(programme(a, m))}</strong><br><span class="dim2">{esc(prog_genes)}</span></td>'
                f'<td>{n_(c["n_cells"])} cells<br><span class="dim2">{c["n_donors"]} donors ({pct(float(c["top_donor_frac"]))}), '
                f'{c["n_datasets"]} dataset{"" if c["n_datasets"] == "1" else "s"} ({pct(float(c["top_dataset_frac"]))})<br>'
                f'{esc(c["top_cell_type"])} {pct(float(c["top_cell_type_frac"]))}<br>'
                f'ribosomal genes: {pct(float(c[f"ribo_frac_k{k}"]))} of the signature</span></td>'
                f'<td class="genes">{", ".join(top_markers(r))}</td>'
                f'<td>{direction}</td><td class="num">{f3(r["cindex"])}</td></tr>')
    table = ('<div class="scroll"><table class="states"><thead><tr><th>pick</th>'
             '<th>programme, and the genes that name it</th><th>composition in the atlas</th>'
             '<th>top-ranked markers (ribosomal genes omitted)</th>'
             '<th>direction</th><th>TCGA C</th></tr></thead><tbody>' + "".join(rows_) + "</tbody></table></div>")
    b = PICKS["primary"]["hvg_pca"]
    sds = SINGLE_DONOR_STATES
    share = {k: pct(float(state_comp(PICKS[a][m])["top_donor_frac"])) for k, (a, m) in
             {"gf": ("primary", "geneformer"), "hs": ("sensitivity", "hvg_pca"), "hp": ("primary", "hvg_pca")}.items()}
    frac = {m: sds[m][0] / sds[m][1] for m in MODELS}
    check(all(frac[m] < frac["hvg_pca"] for m in FMS), "the foundation models have fewer single-donor states than HVG-PCA")
    return f"""<h2 id="states">What the cell states are</h2>
<p class="sub">Each representation's pick is the state whose signature lies furthest above its own random-gene floor: in the risk direction for the pre-specified analysis, in either direction for the sensitivity analysis. Composition is from the {n_(N_CELLS)}-cell atlas; markers are ranked by a Wilcoxon test of the state against all other cells.</p>
{table}
<p class="note">Composition: donors and datasets with, in brackets, the share of the state's cells from the largest one, then the most frequent annotated cell type and its share. A state is flagged single-donor when one donor supplies at least {pct(C.SINGLE_DONOR_FRAC)} of its cells. ER, PR: oestrogen and progesterone receptor; HER2: the ERBB2 gene product.</p>
<p><strong>Most picks are tumour-cell programmes that follow the intrinsic subtypes.</strong> In the risk direction, scGPT's pick is an ERBB2-containing luminal epithelial programme and Geneformer's a basal keratin programme. In the protective direction, HVG-PCA and scGPT both recover a luminal, oestrogen-responsive programme (GATA3, XBP1, AZGP1, TRPS1 in both lists). Their marker lists contain genes of the subtype axes that carry most prognostic information in breast cancer {cite("wirapati")}, so their bulk scores are expected to track PAM50 subtype; that association was not measured here. Malignant cells cluster by patient in single-cell data {cite("tirosh")}: Geneformer's pre-specified pick ({share['gf']} of its cells from one donor) and HVG-PCA's sensitivity pick ({share['hs']}) are essentially one patient's tumour cells, and HVG-PCA's pre-specified pick ({share['hp']}) is the single-dataset state below. Across all kept states, {sds['hvg_pca'][0]} of {sds['hvg_pca'][1]} HVG-PCA, {sds['scgpt'][0]} of {sds['scgpt'][1]} scGPT and {sds['geneformer'][0]} of {sds['geneformer'][1]} Geneformer states are single-donor: the foundation models group cells across patients more often than HVG-PCA without batch correction does, and this did not make their signatures more prognostic in bulk.</p>
<p><strong>The baseline's pre-specified pick is a single-dataset artefact.</strong> Its cells are annotated as exhausted T cells, yet its markers are testis and colon genes (CEACAM7, TKTL1, FAM9C, INSL3), and none of the canonical T-cell markers ({", ".join(T_CELL_GENES)}) is among its 50 genes; it does include {", ".join(BASE_EXTRA_GENES["treg"])}, a regulatory-T-cell transcription factor, and the cell-cycle genes {" and ".join(BASE_EXTRA_GENES["cycle"])}. It was picked because its random-gene floor is low ({f3(b["floor_mean"])}), which makes its margin large although its C-index ({f3(b["cindex"])}) is not. The pre-specified P2 margins are measured against it.</p>
<p><strong>Geneformer's protective pick is a mixed-lineage state.</strong> Its markers are nuclear-retained transcripts and immediate-early genes, the profile of low-quality or dissociation-stressed cells.</p>
<p>Figure 2 places the picks in each representation's two-dimensional uniform manifold approximation and projection (UMAP) of all {n_(N_CELLS)} cells.</p>
{figure("fig_states", "Figure 2")}"""


def clinical_html():
    nas, eas = AV_N["age+stage"]
    nap, eap = AV_N["age+stage+PAM50"]
    prol = AV_PROLIF["age+stage"]
    prolp = AV_PROLIF["age+stage+PAM50"]
    lum = [AV_PICKS[("sensitivity", m, b)] for m in ("hvg_pca", "scgpt") for b in ("age+stage", "age+stage+PAM50")]
    check(all(r["hr_ci_hi"] < 1 for r in lum), "the luminal protective picks keep nominal HR < 1 with and without PAM50")
    nfeat = STR["features"]
    spread = max(max(AGE_STAGE_CVS) - min(AGE_STAGE_CVS), max(PAM_CLIN_CVS) - min(PAM_CLIN_CVS))
    return f"""<h2 id="clinical">Does anything add to clinical information?</h2>
<p class="sub">Each pick is added to a Cox proportional-hazards model of age and stage ({n_(nas)} patients, {eas} deaths) and to age, stage and PAM50 subtype ({n_(nap)} patients, {eap} deaths). The ladder's C-indices are of the score alone; this section asks whether the score adds to what a clinician already has.</p>
{figure("fig_added_value", "Figure 3")}
<p>The hazard ratios are nominal: every pick was chosen on the outcomes of these patients. The change in cross-validated C-index repeats the selection, and the choice of direction, inside every training fold, so it carries no selection. On that measure every pick's change in C spans zero across the five repeats. The luminal protective picks of HVG-PCA and scGPT keep nominal hazard ratios below 1 with PAM50 in the model ({f3(lum[1]["hr_per_sd"], 2)} and {f3(lum[3]["hr_per_sd"], 2)} per standard deviation), but under nested selection they do not raise the cross-validated C-index. The published proliferation score (a fixed gene list, chosen after the first run to implement the design's published-signature rung, and pre-specified for METABRIC) is free of selection; it adds {sgn(prol["delta_cindex_cv_nested_mean"])} to age and stage and nothing once PAM50 is in the model ({sgn(prolp["delta_cindex_cv_nested_mean"])}). Of the family maxima, only HVG-PCA's either-direction maximum adds under nested selection ({sgn(FAMMAX_HVG_AS)}), and it too stops adding once PAM50 is in the model.</p>
<h3>All of a representation's states at once</h3>
<p>A ridge-penalised Cox model with every state of a representation (50-gene signatures: {nfeat['hvg_pca']} HVG-PCA, {nfeat['scgpt']} scGPT, {nfeat['geneformer']} Geneformer) plus age and stage does not improve on age and stage alone (Figure 4); differences of a few thousandths are within the spread across fold seeds. Gradient-boosted survival models (XGBoost, Cox objective, fixed hyperparameters) fall below age and stage alone with every representation's states ({f3(XGB_RANGE[0])}–{f3(XGB_RANGE[1])}), so they are not read as a comparison between representations. Adding PAM50 subtype to age and stage does not raise the cross-validated C either.</p>
<p class="note">Models within one section share their cross-validation folds; sections use different fold seeds, so the same model can differ by up to {f3(spread)} between sections (age + stage {f3(min(AGE_STAGE_CVS))}–{f3(max(AGE_STAGE_CVS))}; PAM50 + age + stage {f3(min(PAM_CLIN_CVS))} and {f3(max(PAM_CLIN_CVS))}).</p>
{figure("fig_stratify", "Figure 4")}"""


def replication_html():
    os_, ds_ = MEND["OS"], MEND["DSS"]
    check(all(k in M["criterion"] for k in ("floor p95", "age-adjusted, cohort-stratified HR", "TCGA direction", f"p < {C.ALPHA}")),
          "the replication criterion stated on the page is the one in metabric.json")
    BL = ' <span class="dim2">(borderline)</span>'

    def ep_cells(sid, ep):
        e = mrep(sid, ep)
        return (f'<td class="num">{f3(e["cindex_oriented"])}</td><td class="num">{f3(e["floor_p95"])}</td>'
                f'<td class="num">{f3(e["hr_per_sd"], 2)} ({pv(e["p"])})</td>'
                f'<td>{yesno(e["replicates"])}{BL if floor_borderline(sid, ep) else ""}</td>')

    sigrows = [(PID[a][m], f'{sw(m)}{NAME[m]}', f'{"pre-specified" if a == "primary" else "sensitivity"} · {esc(programme(a, m))}')
               for a in ("primary", "sensitivity") for m in MODELS]
    sigrows.append((PROLIF_ID, f'{sw("prolif")}Proliferation score', "published reference"))
    os_rows, ds_rows = [], []
    for sid, rep, lab in sigrows:
        s_ = MSIG[sid]
        cls = ' class="ref"' if sid == PROLIF_ID else ""
        os_rows.append(f'<tr{cls}><td>{rep}</td><td>{lab}</td><td class="num">{pct(s_["coverage"])}</td>{ep_cells(sid, "OS")}</tr>')
        sub = ""
        if s_.get("p3"):
            st, v = next(iter(s_["p3"].items()))
            sub = f'{SUBTYPE.get(st, st)}: C {f3(v["OS"]["cindex_oriented"])}, p {pv(v["OS"]["perm_p"])}'
        ds_rows.append(f'<tr{cls}><td>{rep}</td><td>{lab}</td>{ep_cells(sid, "DSS")}<td>{sub}</td></tr>')
    head = '<th>C</th><th>floor p95</th><th>HR/SD (p)</th><th>replicates</th>'
    table = ('<details class="every"><summary>Numbers behind figure 5 for the six picks and the reference</summary><div class="inner">'
             '<p class="tcap">Overall survival</p>'
             '<div class="scroll"><table><thead><tr><th>representation</th><th>pick</th><th>genes measured</th>'
             + head + '</tr></thead><tbody>' + "".join(os_rows) + '</tbody></table></div>'
             '<p class="tcap">Disease-specific survival, and overall survival within the TCGA lead subtype</p>'
             '<div class="scroll"><table><thead><tr><th>representation</th><th>pick</th>'
             + head + '<th>within subtype, OS (no floor)</th></tr></thead><tbody>' + "".join(ds_rows) + '</tbody></table></div>'
             '<p class="note">C is read in the direction fixed in TCGA. Floor p95: 95th percentile of '
             f'{FROZ["n_floor_sets"]} matched random gene sets in METABRIC. HR/SD: age-adjusted hazard ratio per standard '
             'deviation of the score in that direction, stratified by METABRIC cohort. Borderline: C within two approximate '
             'Monte-Carlo standard errors of the floor p95 (normal approximation from the floor mean and p95). Within-subtype '
             f'p: {n_(FROZ["n_perm_p3"])} outcome permutations; {pv(1 / (FROZ["n_perm_p3"] + 1))} is the smallest value '
             'attainable.</p></div></details>')
    gf = PID["primary"]["geneformer"]
    lumA = SC_P3["LumA"]["OS"]
    hi = sorted({sid for sid, _ in HIGH_HR_NO})
    check(hi and all(sid.split(":")[0] in NAME for sid in hi), "high-HR signatures that fail the floor are named")
    pick_ids = {PID[a][m] for a in PID for m in MODELS}
    hi_txt = ", ".join(f"{NAME[sid.split(':')[0]]}'s {'protective' if sid.endswith('protective') else 'risk'} "
                       f"{'pick' if sid in pick_ids else 'family maximum'} "
                       f"({' and '.join(ep for s2, ep in HIGH_HR_NO if s2 == sid)})" for sid in hi)
    base = PID["primary"]["hvg_pca"]
    check(all(mrep(base, ep)["replicates"] for ep in ("OS", "DSS")), "the baseline's artefact pick replicates")
    return f"""<h2 id="replication">Independent replication in METABRIC</h2>
<p class="sub">Specified in <a href="{DESIGN_URL}#13-independent-replication-in-metabric-specified-before-the-run">DESIGN §13</a> and committed (<a href="{REPO}/commit/2ebb73e3f3a51326894186d317c47e963723686a">2ebb73e</a>) before any frozen signature was scored against METABRIC outcomes; the results followed in <a href="{REPO}/commit/35c9233ed8e95f4fe3512472c795dcf2bbafbfb7">35c9233</a>. The stage refuses to run on real outcomes unless the frozen file is committed and unmodified, and it was first run on permuted outcomes.</p>
<p>METABRIC {cite("curtis", "pereira")} is an expression-microarray cohort with long follow-up: {n_(os_['n'])} patients, {n_(os_['events'])} deaths (OS). Disease-specific survival (DSS) is the secondary endpoint; it censors deaths from other causes, leaving {n_(ds_['events'])} events. Twelve signatures were frozen with their genes and TCGA directions: the six picks, the family maxima and the proliferation reference. Nothing was re-selected or re-oriented. A signature replicates on an endpoint if its C, read in the TCGA direction, exceeds the 95th percentile of its METABRIC floor and its age-adjusted hazard ratio, stratified by METABRIC cohort, points in the TCGA direction with p &lt; {C.ALPHA}.</p>
{figure("fig_replication", "Figure 5")}
{table}
<p><strong>Verdict.</strong> The replication shows no foundation-model advantage (figure 6). Geneformer's risk pick meets the criterion on DSS by {sgn(gap(gf, 'DSS'), 4)}, within Monte-Carlo error of the floor p95. Within luminal A tumours scGPT's protective pick has C = {f3(lumA['cindex_oriented'])} with outcome-permutation p = {pv(lumA['perm_p'])}; no random-gene floor was computed within subtypes, so this is a result, not a pre-specified replication.</p>
<p><strong>Random gene sets predict survival in METABRIC.</strong> With {n_(os_['events'])} deaths, the matched-random floors have 95th percentiles of {f3(FLOOR_RANGE[0], 2)}–{f3(FLOOR_RANGE[1], 2)}. Some signatures have hazard-ratio p-values below 10<sup>{MINUS}5</sup> and still do not beat random genes: {hi_txt}. This is the effect described by {REFS["venet"][0]}.</p>
<p><strong>The baseline's artefact pick replicates on both endpoints.</strong> Why it replicates was not examined in the tracked analysis; its 50 genes include the cell-cycle genes {" and ".join(BASE_EXTRA_GENES["cycle"])}.</p>
{figure("fig_margin", "Figure 6")}"""


def approach_html():
    hv = UNIVERSE_HVG
    steps = [
        ("Assemble the atlas",
         f"CELLxGENE Census release {C.CENSUS_VERSION}: every breast-carcinoma disease label, primary data only, droplet-based "
         f"3′ and 5′ chemistries. After quality control, {n_(N_CELLS)} cells from {N_DONORS} donors in {N_DATASETS} datasets "
         f"were sampled in proportion to their {N_CELLTYPES} annotated cell types. No batch correction was applied; "
         f"{pct(CHEN_FRAC)} of the cells come from one collection {cite('chenA')}. {n_(DUPS['n_redundant_cells'])} cells "
         f"({100 * DUPS['n_redundant_cells'] / N_CELLS:.1f}%) duplicate another atlas cell exactly (identical raw counts, same "
         f"donor), between that collection's Global Atlas and its compartment datasets; they were not removed "
         f"({n_(DUPS['n_unique_cells'])} unique cells)."),
        ("Represent every cell three ways",
         f"HVG-PCA: the {n_(C.N_HVG)} most variable genes, log-normalised, reduced to {C.MODELS['hvg_pca']['n_comps']} principal "
         f"components fitted on this atlas. scGPT {cite('cui')}: the whole-human checkpoint, pretrained on non-malignant human "
         f"cells, cell embedding from its classification token. Geneformer {cite('theodoris')}: the 6-layer V1 model, "
         f"pretrained on a corpus that excluded malignant cells, mean-pooled cell embedding. Both models run in bfloat16 on "
         f"CPU; the model revisions in the configuration were pinned after the embedding run and match the only "
         f"downloaded snapshot of each model."),
        ("Find cell states",
         f"Leiden clustering {cite('traag')} of each representation's 15-nearest-neighbour graph at resolutions "
         f"{', '.join(str(r) for r in RESOLUTIONS)}; states with at least {MIN_CELLS} cells are kept: "
         f"{N_STATES['hvg_pca']} HVG-PCA, {N_STATES['scgpt']} scGPT and {N_STATES['geneformer']} Geneformer states."),
        ("Turn states into signatures",
         f"The top {', '.join(str(k) for k in TOPK)} marker genes of each state (Wilcoxon test against all other cells), drawn "
         f"only from the {n_(UNIVERSE)} genes also measured in the bulk cohort ({n_(hv)} of them highly variable): "
         f"{N_SIGS['hvg_pca']}, {N_SIGS['scgpt']} and {N_SIGS['geneformer']} signatures."),
        ("Score the signatures in TCGA-BRCA",
         f"Each of {n_(N_PAT)} primary tumours {cite('tcga')} gets one score per signature: the mean of the signature genes' "
         f"expression, each z-scored across patients. The ladder uses Harrell's C {cite('harrell')} of that score alone, without "
         f"covariates. A Cox model with age and stage gives each pick's hazard ratio, reported separately (Figure 3)."),
        ("Pick and read",
         f"Each representation's pick is the signature furthest above its own random-gene floor across all nine settings "
         f"(three resolutions × three marker counts), chosen on TCGA outcomes. The family-wise null tests the family "
         f"maximum, the largest C among a representation's signatures; read against the same null, the picks have "
         f"p = {', '.join(pv(fw_prim[m]['pick_fw_p']) for m in MODELS)} (pre-specified; HVG-PCA, scGPT, Geneformer) and "
         f"{', '.join(pv(fw_sens[m]['pick_fw_p']) for m in MODELS)} (either direction). A patient bootstrap that re-selects "
         f"the pick in every resample gives the P2 margins. P1–P4 are then read off the ladder."),
    ]
    steps_html = "".join(f'<div class="step"><h3>{t}</h3><p>{d}</p></div>' for t, d in steps)
    return f"""<h2 id="approach">Approach</h2>
<p class="sub">The three representations share everything except the representation itself: the cells, the clustering method, the marker ranking, the gene universe, the scoring and the survival statistics.</p>
{panel(fig_cohorts(), "Three cohorts (OS: overall survival; PFI: progression-free interval; DSS: disease-specific survival)", "cards")}
{panel(fig_flow(), "One atlas, three views, one patient cohort, one ladder")}
<div class="steps">{steps_html}</div>"""


def controls_html():
    clin, pam, pc, cf = REF["clinical"], REF["pam50"], REF["pam50_clinical"], REF["clinical_full"]
    rungs = [
        ("null", f"Survival outcomes permuted across patients {N_PERM} times. In each permutation the largest C-index among all of a "
                 f"representation's {min(N_SIGS.values())}–{max(N_SIGS.values())} signatures is kept, so the null includes the "
                 f"selection of a winner. The family-wise p of the observed maximum is read against it."),
        ("floor", f"{N_FLOOR} random gene sets per signature, matched for size and for mean-expression decile. In breast cancer "
                  f"most random gene sets are associated with outcome, largely through proliferation {cite('venet')}; the floor "
                  f"is what a gene list of that composition scores."),
        ("baseline", "The HVG-PCA states, derived, scored and selected identically."),
        ("test", "The scGPT and Geneformer states."),
        ("references", f"Age + stage: cross-validated C {f3(cf['cindex_cv_mean'])} ({n_(cf['n'])} patients; "
                       f"{f3(clin['cindex_cv_mean'])} on the {n_(clin['n'])} with a PAM50 call). PAM50 subtype alone "
                       f"{cite('parker')}: {f3(pam['cindex_cv_mean'])} cross-validated, {f3(pam['cindex_in_sample'])} in sample. "
                       f"PAM50 + age + stage: {f3(pc['cindex_cv_mean'])}. The published 11-gene proliferation score of PAM50 "
                       f"{cite('nielsen')}, chosen after the first run to fill the design's published-signature rung: C {f3(PROLIF['cindex'])} against a floor mean of {f3(PROLIF['floor_mean'])}, "
                       f"permutation p = {pv(PROLIF['perm_p'])} ({n_(PROLIF['n_perm'])} permutations). Cross-validation is "
                       f"5 × 5-fold, with this section's fold seed."),
    ]
    rungs_html = "".join(f"<tr><td>{r}</td><td>{d}</td></tr>" for r, d in rungs)
    return f"""<h2 id="controls">Controls</h2>
<p class="sub">Every C-index sits on a ladder of controls, and findings are stated as distances between rungs.</p>
{panel(fig_ladder_schematic(), "The ladder (rung positions illustrative)")}
<table class="rungs"><tbody>{rungs_html}</tbody></table>
<h3>What the baseline controls for</h3>
{panel(fig_shared_input(), "One input, three views, one scoring set")}
<p>All three representations are computed from the same expression matrix, and their signatures are scored on the same bulk genes. Prognostic signal present in all three can come from that shared input: cell-type composition, tumour purity or proliferation. The baseline measures what the shared input yields through a linear, unsupervised representation, and every claim about the foundation models is stated as the margin over it.</p>"""


def predictions_html():
    v = VER["primary"]
    pf = {True: "pass", False: "fail", None: "n/a"}

    def two(key):
        return "<br>".join(f'{NAME[m]} <span class="pill {"on" if v[m][key] else "off"}">{pf[v[m][key]]}</span>' for m in FMS)
    check(all(r["cindex"] >= REF["pam50"]["cindex_cv_mean"] and r["cindex"] < REF["clinical_full"]["cindex_cv_mean"] - 0.1
              for a in PICKS for r in PICKS[a].values()), "every pick is below age + stage and not below PAM50 alone")
    MIN_PICK_C = min(r["cindex"] for a in PICKS for r in PICKS[a].values())
    MAX_PICK_C = max(r["cindex"] for a in PICKS for r in PICKS[a].values())
    p4c = L["p4_corrected"]["primary"]
    check(all(p4c[m]["P4_corrected"] is None for m in FMS), "P4 corrected is not applicable")
    rows_ = [
        ("P1", "Foundation-model signatures reach a C-index above the floor (§3), each read against the family-wise "
               "permutation null (§6).",
         two("P1_as_coded"), two("P1_corrected"),
         "The first run compared the margin-selected pick with a null built for the family maximum. Corrected: the family "
         "maximum's family-wise p, with the maximising signature above its floor."),
        ("P2", "Foundation-model signatures exceed HVG-PCA signatures by a margin whose bootstrap interval excludes zero.",
         two("P2_as_coded"), two("P2_corrected"),
         "The design reads both picks above their floors; the first run used raw C-indices with the picks fixed. Corrected: "
         "floor-adjusted margin, picks re-selected in each of "
         f"{n_(P2C['primary']['scgpt']['n_boot'])} patient bootstraps."),
        ("P3", "At least one foundation-model state is prognostic within a PAM50 subtype.",
         two("P3_as_coded"), two("P3_corrected"),
         f"The first run applied a C &gt; {C.P3_WITHIN_SUBTYPE_CINDEX} threshold with no null and did not test the baseline. "
         f"Corrected: best within-subtype C against {n_(P3_NPERM)} within-subtype outcome permutations, every representation "
         f"tested."),
        ("P4", "Foundation-model signatures that pass P2 are enriched for genes outside the HVG set.",
         two("P4_as_coded"), "<br>".join(f'{NAME[m]} <span class="pill off">n/a</span>' for m in FMS),
         f"The design makes P4 conditional on P2, which fails. Only {n_(UNIVERSE_HVG)} of the {n_(UNIVERSE)} marker "
         f"candidates are HVGs, so the as-coded cut (HVG fraction below {C.P4_MAX_HVG_FRAC}) is met by most signatures of "
         f"every representation, the baseline's included."),
    ]
    body = "".join(f'<tr><td class="mono">{a}</td><td data-label="as specified">{b}</td><td data-label="as coded, first run">{c}</td>'
                   f'<td data-label="corrected">{d}</td><td data-label="why they differ">{e}</td></tr>' for a, b, c, d, e in rows_)
    table = ('<div class="scroll"><table class="wrap preds cards"><thead><tr><th></th><th>as specified (DESIGN §3)</th>'
             '<th>as coded, first run</th><th>corrected</th><th>why they differ</th></tr></thead><tbody>'
             + body + "</tbody></table></div>")
    return f"""<h2 id="predictions">Predictions and how they came out</h2>
<p class="sub">Four pass/fail predictions were committed to the repository with the analysis code on 2026-09-28 (<a href="{REPO}/commit/8d2c2d0b7391f48f22ff92be15fff72b8a76cd50">8d2c2d0</a>), before any data were downloaded; this is what "pre-specified" means on this page (there is no external registry). A review of the code against the design found that the first run computed some predictions differently from the text; both versions are kept in <code>results/ladder/ladder.json</code>, and the corrections draw from their own random streams, so no first-run number changed.</p>
{table}
<p>The design also wrote down, before the run, what each outcome would look like (<a href="{DESIGN_URL}#7-what-a-positive-and-a-negative-result-each-look-like">DESIGN §7</a>). The TCGA result matches none of the three shapes exactly. Test ≈ baseline, as in the second and third shapes. Every pick lies above its own random-gene floor, which the third shape ("both fail", test ≈ baseline ≈ floor) excludes, but no family clears its family-wise null, a rung the shapes do not name. On the reference, the data fit the second shape: every pick's C-index ({f3(MIN_PICK_C)}–{f3(MAX_PICK_C)}, in sample and selected on these outcomes) is far below age and stage ({f3(REF['clinical_full']['cindex_cv_mean'])}, cross-validated). METABRIC, where every foundation-model margin is negative, points to the second.</p>
{panel(fig_outcomes(), "The three outcomes specified before the run (DESIGN §7)", "cards")}"""


def flag(b, on, off, borderline):
    return yesno(b, on, off) + (' <span class="dim2">(borderline)</span>' if borderline else "")


def sensitivity_html():
    rows_ = []
    for m in MODELS:
        r = PICKS["sensitivity"][m]
        fw = fw_sens[m]
        p3 = P3C["picks"][f"sensitivity/{m}"]
        rows_.append(f'<tr><td>{sw(m)}{NAME[m]}</td><td>{esc(programme("sensitivity", m))} · {r["direction"]}</td>'
                     f'<td class="num">{f3(r["cindex"])}</td><td class="num">{f3(r["floor_mean"])}</td>'
                     f'<td class="num">{pv(fw["pick_fw_p"])}</td>'
                     f'<td class="num">{f3(fw["stat"])}</td><td class="num">{pv(fw["fw_p"])}</td>'
                     f'<td>{flag(fw["P1_corrected"], "clears", "does not clear", fw["fw_p_borderline"])}</td>'
                     f'<td>{SUBTYPE[p3["argmax_subtype"]]}</td><td class="num">{f3(p3["stat"])}</td><td class="num">{pv(p3["p"])}</td>'
                     f'<td>{flag(p3["P3_pass"], "passes", "fails", p3["borderline"])}</td></tr>')
    table = ('<div class="scroll"><table class="sens"><thead><tr><th>representation</th><th>pick</th><th>C</th><th>floor</th>'
             '<th>pick\'s p</th><th>family max C</th><th>family-wise p</th><th>family null</th>'
             '<th>best subtype</th><th>C within</th><th>p</th><th>P3</th>'
             '</tr></thead><tbody>' + "".join(rows_) + "</tbody></table></div>")
    prim_p3 = [P3C["picks"][f"primary/{m}"]["p"] for m in MODELS]
    km = FIG["fig_km"]["choice"]
    check(km["signature"] == PID["sensitivity"]["hvg_pca"], "the Kaplan-Meier figure shows the baseline's sensitivity pick")
    check(P3C["picks"]["sensitivity/scgpt"]["argmax_subtype"] == "LumA", "scGPT's P3 lead is in luminal A")
    argmax_differs = [m for m in ("hvg_pca", "scgpt") if fw_sens[m]["stat"] != fw_sens[m]["pick_stat"]]
    check(argmax_differs == ["hvg_pca", "scgpt"], "for HVG-PCA and scGPT the family maximum is a different state from the pick")
    return f"""<h2 id="sensitivity">Sensitivity analysis: direction of effect (post hoc)</h2>
<p class="sub">Added on 2026-09-30, after the primary run (<a href="{DESIGN_URL}#11-sensitivity-analysis-direction-of-effect-added-after-the-primary-run">DESIGN §11</a>). The pre-specified ladder is unchanged.</p>
<p>The design reads the C-index in one direction: a signature scores high only when a higher score means worse survival. Protective signatures turned out to be common and include the largest departures from 0.5 (C down to {f3(MIN_C)}). The sensitivity analysis reads every signature in the direction it acts on the full cohort; its floor is read the same way, and the family-wise null takes the largest max(C, 1 − C) over the same {N_PERM} permutations. P3 is tested for every representation, the baseline included, with {n_(P3_NPERM)} within-subtype permutations; pick and direction were fixed on the full cohort, which contains these patients.</p>
{table}
<p>The signatures that clear the null are the family maxima; for HVG-PCA and scGPT these are different states from the picks carried into added value, replication and Figure 7. P2 still fails for both foundation models in this reading, and none of the pre-specified picks passes P3 (p = {pv(min(prim_p3))}–{pv(max(prim_p3))}). Figure 7 shows the lead replicated signature as survival curves in TCGA, where it was selected: the {NAME[km['signature'].split(':')[0]]} sensitivity pick, chosen by the rule in its caption.</p>
{figure("fig_km", "Figure 7")}"""


def pfi_html():
    rows_ = []
    for m in MODELS:
        fp_, fs_ = LP["family_wise"]["primary"][m], LP["family_wise"]["sensitivity"][m]
        p2 = LP["p2_corrected"]["primary"].get(m)
        p2s = (f'{sgn(p2["margin"])} ({ci(*p2["ci"])})' if p2 else "baseline")
        rows_.append(f'<tr><td>{sw(m)}{NAME[m]}</td><td class="num">{pv(fp_["fw_p"])}</td>'
                     f'<td>{flag(fp_["P1_corrected"], "clears", "does not clear", fp_["fw_p_borderline"])}</td>'
                     f'<td class="num">{pv(fs_["fw_p"])}</td>'
                     f'<td>{flag(fs_["P1_corrected"], "clears", "does not clear", fs_["fw_p_borderline"])}</td>'
                     f'<td class="num">{p2s}</td></tr>')
    check(not any(LP["family_wise"]["primary"][m]["P1_corrected"] for m in MODELS), "PFI: no pre-specified family clears")
    check(not any(LP["p2_corrected"][a][m]["P2_pass"] for a in ("primary", "sensitivity") for m in FMS), "PFI: P2 fails")
    diff = [m for m in MODELS if LP["family_wise"]["sensitivity"][m]["P1_corrected"] != fw_sens[m]["P1_corrected"]]
    check(diff == ["geneformer"] and LP["family_wise"]["sensitivity"]["geneformer"]["P1_corrected"]
          and LP["family_wise"]["sensitivity"]["geneformer"]["fw_p_borderline"],
          "PFI differs from OS only for Geneformer either-direction, which clears (borderline)")
    table = ('<div class="scroll"><table><thead><tr><th>representation</th><th>family-wise p, pre-specified</th>'
             '<th>null (pre-specified)</th><th>family-wise p, either direction</th><th>null (either direction)</th>'
             '<th>P2 margin, pre-specified (95% interval)</th></tr></thead>'
             '<tbody>' + "".join(rows_) + "</tbody></table></div>")
    return f"""<h2 id="pfi">Secondary endpoint: progression-free interval</h2>
<p class="sub">The progression-free interval (PFI) was specified as the secondary endpoint for power, but in this cohort it has fewer events than OS ({LP['n_events']} against {N_EV}). It runs through the same stages with the same statistics.</p>
{table}
<p>The pattern is the same as for OS, with one difference: on PFI, Geneformer's either-direction family also clears its null (borderline), which it does not on OS.</p>"""


def pipeline_html():
    def tidy(d):
        return d.replace("hvg_pca, scgpt, geneformer", "HVG-PCA, scGPT, Geneformer").replace("->", "→")
    stages = "".join(f"<tr><td class=\"mono\">{n}</td><td>{esc(tidy(d))}</td></tr>" for n, d in STAGES)
    xena = "".join(f"<li><code>{k}</code> <span class=\"hash\">{v[:12]}…</span></li>" for k, v in C.XENA_SHA256.items())
    mlab = {"data_clinical_patient.txt": "patient clinical table", "data_clinical_sample.txt": "sample clinical table"}
    check(len(FROZ["files_sha256"]) == 3 and set(mlab) < set(FROZ["files_sha256"]), "three METABRIC files")
    mb = "".join(f"<li>{mlab.get(k, 'expression matrix')} <span class=\"hash\">{v[:12]}…</span></li>"
                 for k, v in FROZ["files_sha256"].items())
    rev = TPAR["model_revisions"]
    n_runs = sum(1 for r in RUNLOG if not r.get("dry_run"))
    return f"""<h2 id="pipeline">How the analysis runs</h2>
<p class="sub">The pipeline is a Python package, <code>scfm</code>, whose stages are named commands. Each prints its plan with <code>--dry-run</code> before it touches data, so a person or a coding agent can inspect a run before starting it.</p>
<div class="scroll"><table class="wrap"><thead><tr><th>stage</th><th>does</th></tr></thead><tbody>{stages}</tbody></table></div>
<pre><code>uv run scfm list                          # stages, interpreter, root and log path
uv run scfm run all --dry-run             # every stage's plan; dry runs are not logged
uv run scfm run translate -- --endpoint PFI
make reproduce &amp;&amp; make verify             # tier A, then compare with results/MANIFEST.sha256</code></pre>
<div class="conv">
<div><h3>Run log and provenance</h3><p>Every real run appends its command, start time, duration and exit code to <code>results/run_log.json</code> ({n_runs} runs so far). Every stage writes <code>params.json</code> beside its outputs: arguments, git commit and dirty flag, library versions and the pinned model revisions (scGPT <span class="hash">{rev['scgpt'][:8]}</span>, Geneformer <span class="hash">{rev['geneformer'][:8]}</span>).</p></div>
<div><h3>Checked downloads</h3><p>Every file download (Xena, METABRIC) is checked against a sha256 before use; the atlas is pinned to Census release {C.CENSUS_VERSION} and the models to repository revisions. TCGA-BRCA from UCSC Xena {cite('goldman')}, survival from the TCGA Pan-Cancer Clinical Data Resource {cite('liuJ')}:</p><ul class="hashes">{xena}</ul><p>METABRIC, cBioPortal datahub commit <span class="hash">{C.METABRIC_COMMIT[:8]}</span>:</p><ul class="hashes">{mb}</ul></div>
<div><h3>Three reproduction tiers</h3><p><strong>A, from the tracked results</strong> ({about(T_TIER_A)}, core environment): rebuilds translate, ladder, stratify and report for OS and PFI from the tracked signatures. <strong>B, from the embeddings</strong>: the per-cell embeddings are not distributed; with them, <code>states</code> onward reruns from the tracked cell list. <strong>C, from scratch</strong> ({about(T_TIER_C, 15)} on CPU, both environments): embedding took {hm(T_EMBED)} on 4 cores.</p></div>
<div><h3>Tests, continuous integration and numerics</h3><p>pytest covers the statistics, the stages' helpers and the Geneformer tokeniser, checked against the official tokeniser on a fixture; ruff lints. Continuous integration runs both on every pull request. A smoke run executes every stage except <code>replicate</code>, both foundation models included, on 2,000 cells. Embeddings are computed in bfloat16; <code>scfm run embed -- --precision fp32 --limit 300</code> recomputes a subset in fp32, and no such comparison is tracked in <code>results/</code>.</p></div>
</div>
<p class="note">Times are from <code>results/run_log.json</code> on a shared CPU machine. The atlas UMAP (figure 2) needs per-cell tables from <code>states</code> that are not tracked, so tier A does not rebuild it.</p>"""


DEV_LABEL = {
    "reference_rung": "Reference rung",
    "P1_statistic": "P1 statistic",
    "P2_margin": "P2 margin",
    "P3_null": "P3 null",
    "monte_carlo_error": "Monte-Carlo error",
    "P4_gate": "P4 gate",
    "per_setting_ladder": "Ladder per setting",
    "covariate_adjustment": "Covariate adjustment",
    "selection_conditional_intervals": "Bootstrap intervals",
    "clustering_stability": "Clustering stability",
    "reference_penalty": "Reference penalty",
    "secondary_endpoint": "Secondary endpoint",
}
VERDICT_KEYS = [(a, m, p) for a in ("primary", "sensitivity") for m in FMS for p in ("P1", "P2", "P3", "P4")]


def verdict_changes():
    pf = {True: "pass", False: "fail", None: "not applicable"}
    ch = [(a, m, p, VER[a][m][f"{p}_as_coded"], VER[a][m][f"{p}_corrected"]) for a, m, p in VERDICT_KEYS
          if VER[a][m][f"{p}_as_coded"] != VER[a][m][f"{p}_corrected"]]
    check(not [c for c in ch if c[0] == "primary" and c[2] != "P4"], "no pre-specified P1-P3 verdict changed")
    out = {}
    for a, m, p, x, y in ch:
        out.setdefault((a, p, pf[x], pf[y]), []).append(NAME[m])
    parts = [f"{'post hoc ' if a == 'sensitivity' else 'pre-specified '}{p} for {' and '.join(ms)} from {x} to {y}"
             for (a, p, x, y), ms in out.items()]
    return "No pre-specified P1–P3 verdict changed. The corrections moved " + "; ".join(parts) + "."


def deviations_html():
    devs = L["deviations"]
    unknown = [d["id"] for d in devs if d["id"] not in DEV_LABEL]
    check(not unknown, f"every deviation has a label: {unknown}")
    body = "".join(f'<tr><td>{DEV_LABEL[d["id"]]}</td><td data-label="design">{esc(fix_dev_text(d["design"]))}</td>'
                   f'<td data-label="first run">{esc(fix_dev_text(d["as_coded"]))}</td>'
                   f'<td data-label="reported now">{esc(fix_dev_text(d["now_reported"]))}</td></tr>' for d in devs)
    sd = {d["id"]: d for d in STR["deviations"]}
    rp = sd["stratify_ridge_penalty"]
    body += (f'<tr><td>Stratification penalty</td><td data-label="design">{esc(rp["design"])}</td>'
             f'<td data-label="first run">{esc(rp["why"])}</td><td data-label="reported now">{esc(rp["as_run"])}</td></tr>')
    gf_cov = MSIG[GF_SENS]
    ex_os, ex_ds = MEND["OS"], MEND["DSS"]
    check(ex_os["excluded_time_nonpositive"] == 1 and ex_ds["excluded_status_unknown"] == 1, "one exclusion each")
    metabric_notes = (
        f"METABRIC ran as specified in DESIGN §13. {WORD[ex_os['excluded_time_nonpositive']].capitalize()} patient with "
        f"non-positive OS time is excluded, and {WORD[ex_ds['excluded_status_unknown']]} more from DSS with unknown status. Geneformer's protective pick is "
        f"scored on {gf_cov['n_genes_measured']} of its {gf_cov['n_genes']} genes ({', '.join(gf_cov['genes_missing'])} have "
        f"no probe), above the {pct(C.REPLICATE_MIN_COVERAGE)} minimum. METABRIC's subtype call is PAM50 plus claudin-low, "
        f"taken as it is; stage is sparse there, so the Cox model adjusts for age and is stratified by METABRIC cohort.")
    pairs = ", ".join(f"{n} between {k.replace(' + ', ' and ')}" for k, n in DUPS["dataset_pairs"].items())
    check(DUPS["all_groups_same_donor"] and DUPS["all_groups_across_datasets"] and set(DUPS["group_sizes"]) == {"2"},
          "duplicates are pairs from one donor in two datasets")
    dup_txt = (f"{n_(DUPS['n_redundant_cells'])} atlas cells ({100 * DUPS['n_redundant_cells'] / N_CELLS:.1f}%) have a raw "
               f"count vector identical to another atlas cell from the same donor in a different dataset of the Chen et al. "
               f"collection ({pairs}): its Global Atlas repeats cells of its compartment datasets. The design did not "
               f"foresee this, and the duplicates were not removed; the atlas holds {n_(DUPS['n_unique_cells'])} unique cells "
               f"(<code>results/report/atlas_duplicates.json</code>, which also counts them per pick).")
    labels_txt = ("The first run labelled each state by its most frequent annotated cell type; the baseline's pre-specified "
                  "pick was shown as an exhausted T-cell state. Donor, dataset and cell-type composition per state are now in "
                  "<code>results/states/state_composition.csv</code>, and the page names programmes from marker genes.")
    return f"""<h2 id="deviations">Deviations from the design and corrections</h2>
<p class="sub">Where the first run computed something other than the design describes, or did not compute it. The full account, with every number, is <a href="{DESIGN_URL}#12-deviations-from-the-design-and-corrections">DESIGN §12</a>. The rows are read from <code>results/ladder/ladder.json</code> and <code>results/stratify/cv_summary.json</code>.</p>
<p>{verdict_changes()}</p>
<div class="scroll"><table class="wrap devs cards"><thead><tr><th>item</th><th>design</th><th>first run</th><th>reported now</th></tr></thead><tbody>{body}</tbody></table></div>
<p><strong>State labels.</strong> {labels_txt}</p>
<p><strong>Replication.</strong> {esc(metabric_notes)}</p>
<p><strong>Duplicate cells.</strong> {dup_txt}</p>
<p><strong>Model revisions.</strong> The tracked embeddings were computed before the model revisions were pinned in the configuration, and their <code>params.json</code> records no revision. The pinned revisions (scGPT <span class="hash">{TPAR['model_revisions']['scgpt'][:8]}</span>, Geneformer <span class="hash">{TPAR['model_revisions']['geneformer'][:8]}</span>) are the only snapshots of each model in the download cache used for that run.</p>"""


def limits_html():
    sds = SINGLE_DONOR_STATES
    return f"""<h2 id="limits">Limits</h2>
<ul class="limits">
<li><strong>Neither model saw malignant cells in pretraining.</strong> The scGPT whole-human checkpoint was pretrained on non-diseased cells and Geneformer's corpus excluded malignant cells, while most of the picks here are malignant-cell states. Cancer-adapted checkpoints were not tested.</li>
<li><strong>Light checkpoints.</strong> scGPT whole-human and the 6-layer Geneformer V1 model, run on CPU. A negative result is a result about these checkpoints.</li>
<li><strong>One atlas, no batch correction.</strong> {N_DATASETS} datasets, 3′ and 5′ chemistries, pooled without integration; one baseline pick is confined to a single dataset.</li>
<li><strong>Single-donor states.</strong> {WORD[len(SINGLE_DONOR_PICKS)].capitalize()} of the six picks, and {sds['hvg_pca'][0]}, {sds['scgpt'][0]} and {sds['geneformer'][0]} of the HVG-PCA, scGPT and Geneformer states, draw at least {pct(C.SINGLE_DONOR_FRAC)} of their cells from one donor.</li>
<li><strong>Duplicate cells.</strong> {n_(DUPS['n_redundant_cells'])} cells ({100 * DUPS['n_redundant_cells'] / N_CELLS:.1f}%) duplicate other cells across Chen et al.'s Global Atlas and its compartment datasets; they were not removed.</li>
<li><strong>Power.</strong> {N_EV} deaths among {n_(N_PAT)} TCGA patients, with short follow-up; the TCGA Pan-Cancer Clinical Data Resource flags breast-cancer OS as needing longer follow-up {cite('liuJ')}.</li>
<li><strong>Selection in one cohort.</strong> The picks were chosen on TCGA outcomes. The family-wise null, the selection-aware bootstrap and the nested cross-validation pay for that choice within TCGA; METABRIC is the only out-of-sample test.</li>
<li><strong>A bulk score is a proxy.</strong> The mean z-score of a state's markers mixes the state's abundance with tumour purity, proliferation and subtype.</li>
<li><strong>One indication.</strong> Breast cancer was chosen for cohort size and established references; nothing here extends to other tumours until it is run there.</li>
</ul>"""


def related_html():
    return f"""<h2 id="related">Related work</h2>
<p>Carrying single-cell states into bulk cohorts to test prognosis is established: EcoTyper defined prognostic cell states and ecosystems across carcinomas {cite('luca')}; in breast cancer, {REFS['wuSZ'][0]} deconvolved bulk cohorts into ecotypes with distinct outcomes, and {REFS['chenA'][0]}, whose collection supplies {pct(CHEN_FRAC)} of this atlas's cells, found that the states whose associations held across cohorts were protective. Scissor {cite('sun')} uses the bulk phenotype to choose cells, so its associations are not independent tests of outcome; here state discovery never sees outcomes. Random gene sets predict breast-cancer outcome, mostly through proliferation {cite('venet', 'wirapati')}, which is the reason for the matched-random floor. Foundation models have been evaluated on cancer outcomes within single-cell cohorts {cite('roman')} and on TCGA survival from bulk embeddings {cite('liuW')}, and simple baselines match them on zero-shot cell-level tasks {cite('kedzierska')} and perturbation prediction {cite('ahlmann', 'bendidi')}. We found no prior study that combines foundation-model-defined states in a tumour atlas, a linear arm through the identical pipeline, matched-random and family-wise controls, predictions committed before the data, and a replication with frozen signatures. The related-work note in the repository (<a href="{REPO}/blob/main/research/literature.md">research/literature.md</a>) covers the field and how each reference was verified.</p>"""


def refs_html():
    items = sorted(REFS.values(), key=lambda t: t[1])
    out = []
    for _, text, doi in items:
        url = f"https://arxiv.org/abs/{doi.split(':')[1]}" if doi.startswith("arXiv") else f"https://doi.org/{doi}"
        lab = doi if doi.startswith("arXiv") else f"doi:{doi}"
        out.append(f'<li>{text} <a href="{url}">{lab.replace("/", "/<wbr>")}</a></li>')
    return f"""<h2 id="refs">References</h2>
<ol class="refs">{"".join(out)}</ol>"""


# ----------------------------------------------------------------------------- page
# the ladder chart in miniature: three bars in the view colours and the null tick
FAVICON = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32" aria-hidden="true">'
           '<rect width="32" height="32" rx="6" fill="#16191D"/>'
           '<rect x="6" y="7" width="15" height="4.5" rx="1" fill="#C9CED3"/>'
           '<rect x="6" y="14" width="18" height="4.5" rx="1" fill="#2FA39B"/>'
           '<rect x="6" y="21" width="13" height="4.5" rx="1" fill="#E08A3C"/>'
           '<rect x="24.5" y="5" width="2" height="22" rx="1" fill="#FFFFFF"/></svg>')

NAV_JS = """
(()=>{
  const nav=document.querySelector('.topnav'), btn=nav.querySelector('.navtoggle');
  const links=[...nav.querySelectorAll('ul a')];
  const byId=new Map(links.map(a=>[a.getAttribute('href').slice(1),a]));
  const setOpen=o=>{nav.classList.toggle('open',o);btn.setAttribute('aria-expanded',String(o));};
  btn.addEventListener('click',()=>setOpen(!nav.classList.contains('open')));
  links.forEach(a=>a.addEventListener('click',()=>setOpen(false)));
  document.addEventListener('keydown',e=>{if(e.key==='Escape')setOpen(false);});
  document.addEventListener('click',e=>{if(!nav.contains(e.target))setOpen(false);});
  const heads=[...byId.keys()].map(id=>document.getElementById(id)).filter(Boolean);
  const mark=()=>{
    const y=nav.offsetHeight+24; let cur=null;
    for(const h of heads){ if(h.getBoundingClientRect().top<=y) cur=h; else break; }
    links.forEach(a=>{a.classList.remove('active');a.removeAttribute('aria-current');});
    if(cur){const a=byId.get(cur.id);a.classList.add('active');a.setAttribute('aria-current','true');
      btn.textContent=a.textContent;} else btn.textContent='Sections';
  };
  addEventListener('scroll',mark,{passive:true}); addEventListener('resize',mark); mark();
})();
(()=>{
  const boxes=[...document.querySelectorAll('.scroll,.imgscroll,.svgscroll')];
  boxes.forEach(b=>{let n=b.nextElementSibling;
    if(!(n&&n.classList.contains('swipe'))){n=document.createElement('p');n.className='swipe';
      n.textContent='Wide: scroll sideways to see all of it.';b.after(n);}});
  const hint=()=>boxes.forEach(b=>b.nextElementSibling.classList.toggle('on',b.scrollWidth>b.clientWidth+2));
  addEventListener('resize',hint); addEventListener('load',hint); hint();
})();
"""

CSS = f"""
:root{{
  --paper:#F7F8F9; --panel:#FFFFFF; --ink:{INK}; --muted:{MUTED}; --rule:{RULE}; --axis:{AXIS};
  --base:{COL['hvg_pca']}; --scgpt:{COL['scgpt']}; --gf:{COL['geneformer']}; --band:#EEF1F2;
  --sans:{SANS}; --mono:{MONO};}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);line-height:1.6;
  -webkit-font-smoothing:antialiased;padding-inline:16px}}
.wrap{{max-width:1140px;margin:0 auto;padding:34px 0 96px}}
html{{scroll-behavior:smooth}}
h2[id]{{scroll-margin-top:74px}}
a{{color:var(--ink)}}
.topnav{{position:sticky;top:0;z-index:10;background:color-mix(in srgb,var(--paper) 90%,transparent);
  backdrop-filter:blur(8px);border-bottom:1px solid var(--rule);margin-inline:-16px;padding-inline:16px}}
.topnav .inner{{max-width:1140px;margin:0 auto;display:flex;align-items:center;justify-content:space-between;gap:16px;height:54px;position:relative}}
.topnav .brand{{display:flex;align-items:center;gap:8px;font-family:var(--mono);font-size:12px;color:var(--ink);white-space:nowrap;text-decoration:none}}
.topnav .brand svg{{width:18px;height:18px;flex:none}}
.topnav ul{{list-style:none;margin:0;padding:0;display:flex;gap:10px;font-family:var(--mono);font-size:10px;
  letter-spacing:.02em;text-transform:uppercase}}
.topnav ul a{{display:block;color:var(--muted);text-decoration:none;white-space:nowrap;padding:17px 0 15px;border-bottom:2px solid transparent}}
.topnav ul a:hover{{color:var(--ink)}}
.topnav ul a.active{{color:var(--ink);border-bottom-color:var(--ink)}}
.navtoggle{{display:none;font:inherit;font-family:var(--mono);font-size:11px;letter-spacing:.06em;text-transform:uppercase;
  color:var(--ink);background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:6px 10px;cursor:pointer;
  max-width:60vw;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.navtoggle::before{{content:"\\2261\\00a0\\00a0";font-size:13px}}
.navtoggle:focus-visible{{outline:2px solid var(--scgpt);outline-offset:2px}}
@media (max-width:1200px){{
  .navtoggle{{display:block}}
  .topnav ul{{display:none;position:absolute;top:54px;right:0;flex-direction:column;gap:0;min-width:230px;
    max-height:calc(100vh - 70px);overflow-y:auto;
    background:var(--panel);border:1px solid var(--rule);border-radius:3px;box-shadow:0 8px 24px rgba(22,25,29,.12);padding:6px 0}}
  .topnav.open ul{{display:flex}}
  .topnav ul a{{padding:8px 16px;border-bottom:0;border-left:2px solid transparent;font-size:10.5px}}
  .topnav ul a.active{{border-left-color:var(--ink);background:var(--band)}}
}}
header{{border-bottom:2px solid var(--ink);padding-bottom:22px;margin-bottom:34px}}
.eyebrow{{font-family:var(--mono);font-size:11px;letter-spacing:.12em;text-transform:uppercase;color:var(--muted);margin:0 0 12px}}
.eyebrow a{{color:var(--muted);white-space:nowrap}}
.nw{{white-space:nowrap}}
h1{{font-size:clamp(25px,3.4vw,36px);line-height:1.15;letter-spacing:-.02em;margin:0 0 12px;text-wrap:balance;font-weight:640}}
.byline{{font-family:var(--mono);font-size:12.5px;color:var(--muted);margin:0 0 18px}}
.abstract{{font-size:16.5px;line-height:1.62;margin:0;max-width:980px}}
h2{{font-size:21px;letter-spacing:-.01em;margin:60px 0 8px;font-weight:640;text-wrap:balance}}
h3{{font-size:14.5px;font-weight:640;margin:26px 0 8px}}
.sub{{color:var(--muted);font-size:15px;margin:0 0 18px}}
p{{margin:0 0 14px}}
.note{{color:var(--muted);font-size:12.5px;margin:6px 0 16px}}
.findings{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px;margin:18px 0 0}}
@media (max-width:900px){{.findings{{grid-template-columns:minmax(0,1fr)}}}}
.findings p{{background:var(--panel);border:1px solid var(--rule);border-top:3px solid var(--ink);border-radius:3px;padding:14px 16px;font-size:14px;margin:0}}
.findings strong{{display:block;margin-bottom:6px;font-size:14.5px}}
.panel{{margin:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:14px 12px 8px}}
.panel figcaption{{font-family:var(--mono);font-size:11.5px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);margin:0 0 6px;padding-left:4px}}
.svgscroll{{overflow-x:auto}}
@media (max-width:760px){{.svgscroll svg{{min-width:980px}}}}
.cardrow{{display:flex;flex-wrap:wrap;gap:12px;justify-content:space-between}}
.cardrow .card{{flex:1 1 260px;max-width:340px}}
.cardrow svg{{width:100%;height:auto;display:block}}
@media (max-width:760px){{.cardrow .card{{max-width:none}}}}
.cardnote{{font-size:13px;color:var(--ink);margin:10px 4px 6px;font-weight:600}}
.cardnote.legend{{font-weight:400;color:var(--muted);font-family:var(--mono);font-size:11.5px;display:flex;flex-wrap:wrap;gap:6px 22px}}
.lg{{display:inline-flex;align-items:center;gap:8px}}
.tcap{{font-family:var(--mono);font-size:11px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);margin:12px 0 6px}}
.panel svg{{width:100%;height:auto;display:block}}
.fig{{margin:18px 0 14px}}
.result{{margin:18px 0 16px;background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:12px 14px 10px}}
.result figcaption{{font-size:14px;font-weight:600;margin:0 0 10px;padding-left:2px}}
.result .fk{{font-family:var(--mono);font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);font-weight:500;margin-right:6px}}
.result .imgscroll{{overflow-x:auto}}
.result img{{width:100%;height:auto;display:block}}
@media (max-width:760px){{.result img{{min-width:860px}}}}
.swipe{{display:none;font-family:var(--mono);font-size:10.5px;letter-spacing:.04em;color:var(--muted);margin:6px 0 0}}
.swipe.on{{display:block}}
.result .rtext{{border-top:1px solid var(--rule);margin-top:10px;padding-top:8px}}
.result .rtext p{{font-size:13px;color:var(--muted);margin:0 0 6px}}
.result .rtext strong{{color:var(--ink);font-weight:600}}
.result .rtext ul{{margin:0;padding-left:18px;font-size:12.5px;color:var(--muted)}}
.result .rtext li{{margin:0 0 3px}}
.missing{{color:var(--muted);font-size:13px;border:1.5px dashed var(--axis);padding:24px;text-align:center}}
.steps{{counter-reset:step;margin:26px 0 0}}
.step{{counter-increment:step;position:relative;padding:0 0 18px 46px;border-left:1px solid var(--rule);margin-left:12px}}
.step:last-child{{border-left-color:transparent;padding-bottom:4px}}
.step::before{{content:counter(step);position:absolute;left:-13px;top:-2px;width:26px;height:26px;border-radius:50%;
  background:var(--panel);border:1px solid var(--rule);color:var(--muted);font-family:var(--mono);font-size:11px;
  display:flex;align-items:center;justify-content:center}}
.step h3{{margin:0 0 4px;font-size:15px}}
.step p{{font-size:14.5px}}
.scroll{{overflow-x:auto;border:1px solid var(--rule);border-radius:3px;margin:0 0 16px;
  background:linear-gradient(to right,var(--panel) 30%,rgba(255,255,255,0)) left/40px 100% no-repeat local,
  linear-gradient(to left,var(--panel) 30%,rgba(255,255,255,0)) right/40px 100% no-repeat local,
  radial-gradient(farthest-side at 0 50%,rgba(22,25,29,.16),transparent) left/12px 100% no-repeat scroll,
  radial-gradient(farthest-side at 100% 50%,rgba(22,25,29,.16),transparent) right/12px 100% no-repeat scroll,var(--panel)}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{padding:7px 12px;text-align:left;border-bottom:1px solid var(--rule);white-space:nowrap;vertical-align:top}}
th{{font-family:var(--mono);font-size:10.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:500;background:#F2F4F5;white-space:normal;min-width:52px}}
@media (max-width:760px){{
  table.cards,table.cards tbody,table.cards tr,table.cards td{{display:block;width:100%}}
  table.cards thead{{display:none}}
  table.cards tr{{border-bottom:1px solid var(--axis);padding:6px 0}}
  table.cards td{{border-bottom:none;min-width:0!important;padding:4px 12px}}
  table.cards td[data-label]::before{{content:attr(data-label);display:block;font-family:var(--mono);font-size:10px;
    letter-spacing:.06em;text-transform:uppercase;color:var(--muted)}}
  table.cards td:first-child{{font-weight:600}}
}}
tbody tr:last-child td{{border-bottom:none}}
tr.ref td{{border-top:1px solid var(--axis)}}
table.wrap td{{white-space:normal;min-width:150px}}
table.wrap td:first-child{{min-width:0;white-space:nowrap}}
table.preds td:nth-child(2){{min-width:220px}}
table.preds td:nth-child(3),table.preds td:nth-child(4){{min-width:150px;line-height:1.9}}
table.preds td:nth-child(5){{min-width:280px}}
table.sens td:nth-child(2){{white-space:normal;min-width:150px}}
table.sens th,table.sens td{{padding:7px 8px}}
table.sens td .pill{{white-space:nowrap}}
table.devs td{{font-size:12.5px;min-width:230px}}
table.states td{{white-space:normal}}
table.states td:nth-child(1){{min-width:150px}}
table.states td:nth-child(2){{min-width:200px}}
table.states td:nth-child(3){{min-width:190px}}
table.states td.genes{{min-width:220px;font-family:var(--mono);font-size:11.5px}}
.num{{font-family:var(--mono);font-variant-numeric:tabular-nums;text-align:right}}
.mono{{font-family:var(--mono);font-size:12px}}
.dim2{{color:var(--muted);font-size:12px}}
.swatch{{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:8px;vertical-align:-1px}}
.pill{{font-family:var(--mono);font-size:10.5px;letter-spacing:.04em;border-radius:2px;padding:1px 6px;border:1px solid var(--rule);white-space:nowrap}}
.pill.on{{background:#EEF3F3;color:var(--scgpt);border-color:var(--scgpt)}}
.pill.off{{color:var(--muted)}}
.rungs{{width:100%;border-collapse:collapse;font-size:13.5px;margin:14px 0 6px}}
.rungs td{{padding:9px 12px;vertical-align:top;border-bottom:1px solid var(--rule);white-space:normal}}
.rungs td:first-child{{font-family:var(--mono);font-size:11px;letter-spacing:.1em;text-transform:uppercase;white-space:nowrap;width:110px;padding-top:12px}}
.rungs tr:last-child td{{border-bottom:none}}
.conv{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin:18px 0}}
@media (max-width:720px){{.conv{{grid-template-columns:minmax(0,1fr)}}}}
.conv>div{{background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:14px 16px 8px;min-width:0}}
.conv h3{{margin:0 0 6px;font-size:14px}}
.conv p{{font-size:13.5px;color:var(--muted);margin:0 0 8px}}
.hashes{{margin:0 0 8px;padding-left:18px;font-size:12.5px;color:var(--muted)}}
.hash{{font-family:var(--mono);font-size:11.5px;color:var(--muted)}}
pre{{background:#EDEFF1;border-radius:3px;padding:12px 14px;overflow-x:auto;font-size:12.5px;line-height:1.5;margin:14px 0}}
pre code{{background:none;padding:0}}
.limits li,.refs li{{margin:0 0 8px;font-size:14px}}
.refs{{padding-left:22px;font-size:13.5px}}
.refs li{{color:var(--muted);overflow-wrap:break-word}}
details.every{{border:1px solid var(--rule);border-radius:3px;background:var(--panel);margin:8px 0 16px}}
details.every>summary{{cursor:pointer;padding:9px 13px;font-size:12.5px;color:var(--muted);list-style:none;display:flex;justify-content:space-between;align-items:center;gap:12px}}
details.every>summary::-webkit-details-marker{{display:none}}
details.every>summary::after{{content:"show";font-family:var(--mono);font-size:10px;letter-spacing:.06em;color:var(--muted);border:1px solid var(--rule);border-radius:2px;padding:1px 6px}}
details.every[open]>summary::after{{content:"hide"}}
details.every>summary:focus-visible,a:focus-visible{{outline:2px solid var(--scgpt);outline-offset:2px}}
details.every .inner{{padding:0 13px 13px}}
code{{font-family:var(--mono);font-size:.9em;background:#EDEFF1;padding:1px 5px;border-radius:3px;overflow-wrap:anywhere}}
footer{{margin-top:64px;padding-top:18px;border-top:1px solid var(--rule);font-family:var(--mono);font-size:11.5px;color:var(--muted);line-height:1.8}}
footer a{{color:var(--ink)}}
@media (prefers-reduced-motion:reduce){{*{{scroll-behavior:auto!important}}}}
"""

NAV = [("result", "Result"), ("states", "States"), ("clinical", "Clinical"), ("replication", "Replication"),
       ("approach", "Approach"), ("controls", "Controls"), ("predictions", "Predictions"),
       ("sensitivity", "Sensitivity"), ("pfi", "PFI"), ("pipeline", "Pipeline"), ("deviations", "Deviations"),
       ("limits", "Limits"), ("related", "Related"), ("refs", "Refs")]


def build():
    today = date.today()
    body = "\n".join([header_html(today), result_html(), states_html(), clinical_html(), replication_html(),
                      approach_html(), controls_html(), predictions_html(), sensitivity_html(), pfi_html(),
                      pipeline_html(), deviations_html(), limits_html(), related_html(), refs_html()])
    ids = re.findall(r'<h2 id="([^"]+)"', body)
    check(ids == [a for a, _ in NAV], f"navbar sections match the page sections: {ids}")
    nav = "".join(f'<li><a href="#{a}">{t}</a></li>' for a, t in NAV)
    xena = ", ".join(C.XENA_SHA256)
    html_ = f"""<title>Cell states and survival</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="icon" type="image/svg+xml" href="data:image/svg+xml,{quote(FAVICON)}">
<style>{CSS}</style>
<nav class="topnav" aria-label="Sections"><div class="inner">
<a class="brand" href="#top">{FAVICON}<span>cell states → survival</span></a>
<button class="navtoggle" type="button" aria-label="Jump to section" aria-expanded="false" aria-controls="navlist">Sections</button>
<ul id="navlist">{nav}</ul></div></nav>
<span id="top"></span>
<div class="wrap">
{body}
<footer>
<div>Data: CELLxGENE Census release {C.CENSUS_VERSION} (breast-carcinoma donors, primary data) · TCGA-BRCA from UCSC Xena ({xena}) · METABRIC from the cBioPortal datahub, commit {C.METABRIC_COMMIT[:8]}. Data are downloaded at run time and not redistributed.</div>
<div>Code: <a href="{REPO}">single-cell-fm-probing</a> · design: <a href="{DESIGN_URL}">docs/DESIGN.md</a> · this page is built by <code>docs/site/build.py</code> from <code>results/</code>.</div>
<div>Page built {today.isoformat()}.</div>
</footer>
</div>
<script>{NAV_JS}</script>
"""
    OUT.write_text('<!doctype html>\n<html lang="en">\n<meta charset="utf-8">\n'
                   '<meta name="description" content="Cell states from scGPT, Geneformer and an HVG-PCA baseline, turned into '
                   'signatures and tested for breast-cancer survival in TCGA-BRCA and METABRIC.">\n'
                   + html_ + "\n</html>\n")
    print(f"wrote {OUT}  {OUT.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    build()
