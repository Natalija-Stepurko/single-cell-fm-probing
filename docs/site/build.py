"""Build the portfolio page for this study as a single self-contained HTML file.

The page is the pre-registered design, written before any result exists. Result panels are
slots: when `results/report/fig_*.png` and `results/ladder/ladder.json` exist they are embedded
(PNGs as data URIs, the ladder as a table); until then each slot shows what the panel will
contain and how to read it. Re-run after stage 06 to fill the page; nothing else changes.

    python3 docs/site/build.py                -> docs/index.html (served by GitHub Pages)
"""
import base64
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import config as C  # noqa: E402

OUT = ROOT / "docs" / "index.html"
REPORT = C.RESULTS / "report"
LADDER = C.RESULTS / "ladder" / "ladder.json"
REPO = "https://github.com/Natalija-Stepurko/single-cell-fm-probing"

COL = {"hvg_pca": "#4A4F55", "scgpt": "#1F6F6B", "geneformer": "#C2681A"}
NAME = {"hvg_pca": "HVG-PCA (baseline)", "scgpt": "scGPT", "geneformer": "Geneformer"}
INK, MUTED, RULE, AXIS, PANEL = "#16191D", "#5B646E", "#DDE1E4", "#B9C0C6", "#FFFFFF"
MONO = "ui-monospace,SFMono-Regular,Menlo,Consolas,monospace"
SANS = "ui-sans-serif,system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif"


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
    out = []
    for _ in range(n):
        out.append(f'<circle cx="{x + rnd.uniform(6, 54):.1f}" cy="{y + rnd.uniform(6, 34):.1f}" r="{r}" fill="{fill}" opacity=".8"/>')
    return "".join(out)


# ----------------------------------------------------------------------------- figures
def fig_dispute():
    W, H = 980, 196
    b = []
    b.append(box(20, 10, 390, 122, fill="#F2F4F5"))
    b.append(T(40, 32, "SCEPTICS", 10.5, MUTED, mono=True, weight=600))
    b.append(T(40, 58, "A linear summary of ~2,000 variable genes", 14))
    b.append(T(40, 78, "matches or beats the foundation models", 14))
    b.append(T(40, 108, "Kedzierska 2025 · Boiarsky 2023 · Souza &amp; Mehta 2026", 10, MUTED, mono=True))
    b.append(box(570, 10, 390, 122, fill="#F2F4F5"))
    b.append(T(590, 32, "OPTIMISTS", 10.5, MUTED, mono=True, weight=600))
    b.append(T(590, 58, "The models learn cell biology that", 14))
    b.append(T(590, 78, "transfers zero-shot to unseen data", 14))
    b.append(T(590, 108, "Cui 2024 · Theodoris 2023 · Rosen 2026", 10, MUTED, mono=True))
    b.append(T(490, 52, "so far argued on", 11, MUTED, anchor="middle", italic=True))
    b.append(T(490, 74, "cluster maps and", 12, INK, anchor="middle", weight=600))
    b.append(T(490, 92, "cell-type labels", 12, INK, anchor="middle", weight=600))
    b.append(arrow(490, 140, 490, 158))
    b.append(T(490, 176, "This study asks the same question on a clinical endpoint: patient survival.", 13.5, INK, anchor="middle", weight=600))
    return svg(W, H, "".join(b), "The two poles of the dispute and where this study moves it")


def fig_cohorts():
    W, H = 980, 186
    b = []
    b.append(box(20, 10, 440, 112))
    b.append(cells_glyph(30, 26, n=40, seed=2, fill=MUTED))
    b.append(cells_glyph(58, 44, n=28, seed=7, fill="#8A9298"))
    b.append(T(170, 34, "single-cell atlas", 13, INK, weight=600))
    b.append(T(170, 58, "tens of thousands of cells", 12, MUTED))
    b.append(T(170, 76, "from a few dozen tumours,", 12, MUTED))
    b.append(T(170, 94, "used to find cell states", 12, MUTED))
    b.append(T(240, 138, "CELLxGENE Census · breast cancer · primary data", 10, MUTED, anchor="middle", mono=True))
    b.append(box(520, 10, 440, 112))
    for i in range(6):
        for j in range(5):
            b.append(f'<rect x="{540 + i * 17}" y="{26 + j * 17}" width="11" height="12" rx="1" fill="{MUTED}" opacity=".55"/>')
    b.append(T(660, 34, "bulk tumour cohort", 13, INK, weight=600))
    b.append(T(660, 58, "~1,100 patients, one RNA profile", 12, MUTED))
    b.append(T(660, 76, "each, with follow-up, used only", 12, MUTED))
    b.append(T(660, 94, "to test signatures on survival", 12, MUTED))
    b.append(T(740, 138, "TCGA-BRCA via UCSC Xena · OS, PFI · age, stage, PAM50", 10, MUTED, anchor="middle", mono=True))
    b.append(T(490, 60, "∅", 24, INK, anchor="middle"))
    b.append(T(490, 168, "No patient appears in both cohorts.", 12.5, INK, anchor="middle", weight=600))
    return svg(W, H, "".join(b), "Two cohorts with no shared patients")


def fig_flow():
    """Atlas -> three representations -> states -> signatures -> bulk cohort -> ladder."""
    import random
    W, H = 980, 262
    b = []
    X = {"atlas": 20, "rep": 165, "states": 350, "sig": 535, "bulk": 680, "out": 850}
    rows = {"hvg_pca": 44, "scgpt": 114, "geneformer": 184}   # top of each 56-high row
    RH = 56
    cy = lambda y: y + RH / 2
    mid = 44 + 70 + 28   # 142: vertical centre of the middle row
    head = lambda x, w, s: T(x + w / 2, 18, s, 10, MUTED, anchor="middle", mono=True, weight=600)
    # 1 atlas
    b.append(head(X["atlas"], 110, "1 · ONE INPUT"))
    b.append(box(X["atlas"], 44, 110, 196))
    b.append(cells_glyph(X["atlas"] + 2, 56, n=36, seed=5, r=2.2, fill=MUTED))
    b.append(cells_glyph(X["atlas"] + 2, 106, n=30, seed=9, r=2.2, fill="#8A9298"))
    b.append(cells_glyph(X["atlas"] + 2, 156, n=24, seed=13, r=2.2, fill=MUTED))
    b.append(T(X["atlas"] + 55, 205, "breast tumour", 11.5, INK, anchor="middle", weight=600))
    b.append(T(X["atlas"] + 55, 221, "single-cell atlas", 11.5, INK, anchor="middle", weight=600))
    # 2 representations
    b.append(head(X["rep"], 150, "2 · THREE VIEWS"))
    subs = {"hvg_pca": "linear · no learning", "scgpt": "foundation model", "geneformer": "foundation model"}
    for m, y in rows.items():
        b.append(box(X["rep"], y, 150, RH, stroke=COL[m], sw=1.5))
        b.append(f'<rect x="{X["rep"]}" y="{y}" width="5" height="{RH}" rx="2" fill="{COL[m]}"/>')
        b.append(T(X["rep"] + 16, y + 20, NAME[m].replace(" (baseline)", ""), 12.5, INK, weight=600))
        b.append(T(X["rep"] + 16, y + 39, subs[m], 10, MUTED, mono=True))
        b.append(arrow(X["atlas"] + 110, mid, X["rep"] - 3, cy(y)))
    # 3 states
    b.append(head(X["states"], 150, "3 · CELL STATES"))
    seeds = {"hvg_pca": 1, "scgpt": 2, "geneformer": 3}
    for m, y in rows.items():
        b.append(box(X["states"], y, 150, RH))
        rnd = random.Random(seeds[m])
        for k in range(3):
            cx = X["states"] + 20 + k * 21
            for _ in range(7):
                b.append(f'<circle cx="{cx + rnd.uniform(-6, 6):.1f}" cy="{cy(y) + rnd.uniform(-11, 11):.1f}" r="2" fill="{COL[m]}" opacity=".8"/>')
        b.append(T(X["states"] + 86, cy(y) - 8, "groups of", 10.5, MUTED))
        b.append(T(X["states"] + 86, cy(y) + 8, "similar cells", 10.5, MUTED))
        b.append(arrow(X["rep"] + 150, cy(y), X["states"] - 3, cy(y)))
    # 4 signatures
    b.append(head(X["sig"], 100, "4 · MARKER GENES"))
    for m, y in rows.items():
        b.append(box(X["sig"], y, 100, RH))
        for i, wl in enumerate([68, 50, 60, 38]):
            b.append(f'<rect x="{X["sig"] + 16}" y="{y + 11 + i * 10}" width="{wl}" height="4" rx="2" fill="{COL[m]}" opacity="{.95 - i * .18:.2f}"/>')
        b.append(arrow(X["states"] + 150, cy(y), X["sig"] - 3, cy(y)))
    # 5 bulk
    b.append(head(X["bulk"], 120, "5 · SECOND COHORT"))
    b.append(box(X["bulk"], 44, 120, 196))
    for i in range(6):
        for j in range(6):
            b.append(f'<rect x="{X["bulk"] + 17 + i * 15}" y="{58 + j * 15}" width="10" height="10" rx="1" fill="{MUTED}" opacity=".55"/>')
    b.append(T(X["bulk"] + 60, 168, "TCGA breast", 11.5, INK, anchor="middle", weight=600))
    b.append(T(X["bulk"] + 60, 184, "tumours", 11.5, INK, anchor="middle", weight=600))
    b.append(T(X["bulk"] + 60, 208, "~1,100 patients", 10, MUTED, anchor="middle", mono=True))
    b.append(T(X["bulk"] + 60, 223, "with follow-up", 10, MUTED, anchor="middle", mono=True))
    for m, y in rows.items():
        b.append(arrow(X["sig"] + 100, cy(y), X["bulk"] - 3, mid))
    # 6 ladder
    b.append(head(X["out"], 110, "6 · READ"))
    b.append(box(X["out"], 44, 110, 196))
    for lab, yy, c, dash in [("reference", 80, INK, "6 3"), ("test", 120, COL["scgpt"], ""),
                             ("baseline", 160, COL["hvg_pca"], ""), ("floor", 200, MUTED, "2 3")]:
        b.append(T(X["out"] + 55, yy - 12, lab, 10, MUTED, anchor="middle", mono=True))
        b.append(line(X["out"] + 18, yy, X["out"] + 92, yy, stroke=c, sw=2.2 if not dash else 1.5, dash=dash))
    b.append(arrow(X["bulk"] + 120, mid, X["out"] - 3, mid))
    return svg(W, H, "".join(b), "Design flow from one atlas through three representations to survival in a second cohort")


def fig_ladder_schematic():
    W, H = 980, 262
    x0, x1, ay = 60, 920, 196

    def X(c):
        return x0 + (c - 0.5) / 0.3 * (x1 - x0)
    b = [line(x0, ay, x1, ay, stroke=INK, sw=1.2)]
    for c in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8):
        b.append(line(X(c), ay, X(c), ay + 6, stroke=INK))
        b.append(T(X(c), ay + 20, f"{c:.2f}", 11, MUTED, anchor="middle", mono=True))
    b.append(T((x0 + x1) / 2, ay + 46, "concordance index: how often, of two patients, the higher-scored one dies first (0.50 = coin flip)", 12, MUTED, anchor="middle"))
    rungs = [("null", 0.53, AXIS, "6 4", "outcomes shuffled"),
             ("floor", 0.57, MUTED, "3 3", "matched random genes"),
             ("baseline", 0.63, COL["hvg_pca"], "", "HVG-PCA states"),
             ("test", 0.69, COL["scgpt"], "", "scGPT / Geneformer"),
             ("reference", 0.74, INK, "8 4", "PAM50 subtype")]
    for lab, c, col, dash, desc in rungs:
        x = X(c)
        b.append(T(x, 16, lab.upper(), 10.5, MUTED if col == AXIS else col, anchor="middle", mono=True, weight=600))
        b.append(T(x, 34, desc, 10.5, MUTED, anchor="middle"))
        b.append(line(x, 48, x, ay - 1, stroke=col, sw=2.4 if not dash else 1.6, dash=dash))
    for (a, c2), yy, lab, w in [((0.57, 0.69), 100, "P1 · test above floor", 500), ((0.63, 0.69), 150, "P2 · margin over baseline", 600)]:
        b.append(f'<line x1="{X(a) + 3}" y1="{yy}" x2="{X(c2) - 3}" y2="{yy}" stroke="{INK}" stroke-width="1.6" marker-end="url(#arr)" marker-start="url(#arr)"/>')
        b.append(T((X(a) + X(c2)) / 2, yy - 13, lab, 11, INK, anchor="middle", mono=True, weight=w, halo=True))
    return svg(W, H, "".join(b), "Schematic ladder of concordance rungs: null, floor, baseline, test, reference")


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
        b.append(T(384, y + 15, NAME[m].replace(" (baseline)", ""), 12, INK, anchor="middle", weight=600))
        b.append(arrow(464, y + 15, 570, y + 15))
        b.append(box(574, y, 120, 30))
        b.append(T(634, y + 15, "signature", 11, MUTED, anchor="middle", mono=True))
        b.append(arrow(694, y + 15, 760, 80))
    b.append(box(764, 50, 190, 60, fill="#F2F4F5"))
    b.append(T(859, 71, "scored in the same", 11.5, INK, anchor="middle"))
    b.append(T(859, 89, "bulk tumour genes", 11.5, INK, anchor="middle"))
    return svg(W, H, "".join(b), "All three representations are built from one expression matrix and scored on the same genes")


def fig_prediction(kind):
    W, H = 220, 108
    b = []
    if kind == "P1":
        b.append(line(10, 68, 210, 68, stroke=INK))
        b.append(T(70, 20, "floor", 10, MUTED, anchor="middle", mono=True))
        b.append(line(70, 30, 70, 68, stroke=MUTED, sw=1.4, dash="3 3"))
        b.append(T(150, 20, "FM state", 10, COL["scgpt"], anchor="middle", mono=True))
        b.append(line(150, 30, 150, 68, stroke=COL["scgpt"], sw=2.4))
        b.append(arrow(76, 50, 144, 50, stroke=INK, sw=1.4))
        b.append(T(110, 90, "concordance →", 10, MUTED, anchor="middle", mono=True))
    elif kind == "P2":
        b.append(T(110, 10, "FM − baseline", 10, MUTED, anchor="middle", mono=True))
        b.append(line(110, 22, 110, 72, stroke=INK))
        b.append(T(110, 88, "0 = no gain", 10, MUTED, anchor="middle", mono=True))
        b.append(line(124, 40, 200, 40, stroke=COL["scgpt"], sw=3)); b.append(f'<circle cx="162" cy="40" r="4" fill="{COL["scgpt"]}"/>')
        b.append(T(8, 40, "pass", 10, COL["scgpt"], mono=True))
        b.append(line(30, 62, 130, 62, stroke=AXIS, sw=3)); b.append(f'<circle cx="80" cy="62" r="4" fill="{AXIS}"/>')
        b.append(T(150, 62, "fail", 10, MUTED, mono=True))
    elif kind == "P3":
        for i, lab in enumerate(["LumA", "LumB", "Basal"]):
            x = 8 + i * 72
            b.append(T(x + 32, 16, lab, 10, MUTED, anchor="middle", mono=True))
            b.append(box(x, 26, 64, 46, fill="#F2F4F5"))
            gap = 15 if i == 2 else 3
            b.append(f'<path d="M{x + 7},34 L{x + 28},44 L{x + 46},52 L{x + 58},62" stroke="{INK}" fill="none" stroke-width="1.2"/>')
            b.append(f'<path d="M{x + 7},34 L{x + 28},{44 - gap} L{x + 46},{52 - gap} L{x + 58},{62 - gap}" stroke="{COL["scgpt"] if i == 2 else AXIS}" fill="none" stroke-width="1.8"/>')
        b.append(T(110, 92, "survival split within a subtype", 10, MUTED, anchor="middle", mono=True))
    else:
        import random
        b.append(box(10, 14, 200, 44, fill="#F2F4F5"))
        b.append(f'<rect x="10" y="14" width="90" height="44" rx="3" fill="{AXIS}" opacity=".5"/>')
        rnd = random.Random(4)
        for _ in range(10):
            b.append(f'<circle cx="{rnd.uniform(112, 200):.0f}" cy="{rnd.uniform(20, 52):.0f}" r="3" fill="{COL["geneformer"]}"/>')
        for _ in range(3):
            b.append(f'<circle cx="{rnd.uniform(18, 92):.0f}" cy="{rnd.uniform(20, 52):.0f}" r="3" fill="{COL["geneformer"]}"/>')
        b.append(T(55, 72, "HVG set", 10, MUTED, anchor="middle", mono=True))
        b.append(T(155, 72, "other genes", 10, MUTED, anchor="middle", mono=True))
        b.append(T(110, 96, "where the signature's genes fall", 10, MUTED, anchor="middle", mono=True))
    return svg(W, H, "".join(b), f"Prediction {kind} glyph")


def fig_outcomes():
    W, H = 980, 196
    b = []
    scen = [("FM wins", {"floor": .56, "base": .62, "test": .70, "ref": .74}, "test clearly above baseline, both above floor"),
            ("Sceptics win", {"floor": .56, "base": .64, "test": .65, "ref": .74}, "test ≈ baseline, both below the reference"),
            ("Both fail", {"floor": .56, "base": .575, "test": .585, "ref": .74}, "test ≈ baseline ≈ floor: signal lost")]
    for i, (title, v, desc) in enumerate(scen):
        x0 = 30 + i * 320
        b.append(box(x0, 10, 290, 128))
        b.append(T(x0 + 16, 30, title, 13, INK, weight=600))
        ay = 100
        b.append(line(x0 + 20, ay, x0 + 270, ay, stroke=INK))
        Xc = lambda c: x0 + 20 + (c - .5) / .3 * 250
        for c, col, dash, sw in [(v["floor"], MUTED, "3 3", 1.6), (v["base"], COL["hvg_pca"], "", 2.4),
                                 (v["test"], COL["scgpt"], "", 2.4), (v["ref"], INK, "8 4", 1.6)]:
            b.append(line(Xc(c), 46, Xc(c), ay, stroke=col, sw=sw, dash=dash))
        b.append(T(x0 + 16, 120, desc, 10.5, MUTED))
    lx = 190
    for lab, col, dash, sw in [("floor", MUTED, "3 3", 1.6), ("baseline", COL["hvg_pca"], "", 2.4),
                               ("test", COL["scgpt"], "", 2.4), ("reference", INK, "8 4", 1.6)]:
        b.append(line(lx, 170, lx + 26, 170, stroke=col, sw=sw, dash=dash))
        b.append(T(lx + 34, 170, lab, 10.5, MUTED, mono=True))
        lx += 170
    return svg(W, H, "".join(b), "Three possible outcomes shown as miniature ladders")


# ----------------------------------------------------------------------------- results
def data_uri(p: Path):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


def result_slot(name, title, what, read, wide=False):
    p = REPORT / f"{name}.png"
    if p.exists():
        img = f'<img src="{data_uri(p)}" alt="{title}">'
        if wide:
            img = f'<div class="imgscroll">{img}</div>'
        state = ""
    else:
        img = (f'<div class="pending"><span class="pending-k">{name}.png</span>'
               f'<span>Result panel — fills in when stage 06 runs.</span></div>')
        state = ' data-pending="true"'
    return (f'<figure class="result{" wide" if wide else ""}" data-slot="{name}"{state}>'
            f'<figcaption><span class="rk">{title}</span></figcaption>{img}'
            f'<div class="rtext"><p><strong>What it shows.</strong> {what}</p><p><strong>How to read it.</strong> {read}</p></div>'
            f'</figure>')


def ladder_table(sensitivity=False):
    rows = []
    have = LADDER.exists()
    L = json.load(open(LADDER)) if have else None
    order = ["hvg_pca", "scgpt", "geneformer"]
    src = (L["sensitivity"]["ladder"] if sensitivity else L["ladder"]) if have else []
    by = {r["model"]: r for r in src}
    for m in order:
        r = by.get(m)
        sw = f'<span class="swatch" style="background:{COL[m]}"></span>'
        if r:
            ct = r["top_cell_type"] + (f' · {r["direction"]}' if sensitivity else "")
            rows.append(f'<tr><td>{sw}{NAME[m]}</td><td>{ct}</td>'
                        f'<td class="num">{r["cindex"]:.3f}</td><td class="num">{r["ci_lo"]:.3f}–{r["ci_hi"]:.3f}</td>'
                        f'<td class="num">{r["floor_mean"]:.3f}</td><td class="num">{r["above_floor"]:+.3f}</td>'
                        f'<td class="num">{r["null_p95"]:.3f}</td><td class="num">{r["hvg_frac"]:.2f}</td></tr>')
        else:
            rows.append(f'<tr><td>{sw}{NAME[m]}</td><td class="dim">—</td>' + '<td class="num dim">—</td>' * 6 + '</tr>')
    ref = f'{L["reference_pam50"]:.3f}' if have and L.get("reference_pam50") is not None else "—"
    head = ("<thead><tr><th>representation</th><th>best state (cell type)</th><th>C-index</th><th>95% patient bootstrap</th>"
            "<th>floor</th><th>above floor</th><th>null 95th pct</th><th>HVG fraction</th></tr></thead>")
    foot = f'<tr class="ref"><td>PAM50 reference</td><td class="dim">intrinsic subtype call</td><td class="num">{ref}</td>' + '<td class="num dim">—</td>' * 5 + '</tr>'
    meta = (f'{L["endpoint"]} · {L["n_patients"]} patients · {L["n_events"]} events' if have
            else f'{C.PRIMARY_ENDPOINT} · patients and events fill in with the data')
    return (f'<div class="scroll"><table>{head}<tbody>{"".join(rows)}{foot}</tbody></table></div>'
            f'<p class="sub" style="margin-top:8px">{meta}. One row per representation: its best cell state, chosen by margin above its own floor, not by raw score.</p>')


def predictions_status(sensitivity=False):
    if not LADDER.exists():
        return ""
    L = json.load(open(LADDER))
    P = L["sensitivity"]["predictions"] if sensitivity else L["predictions"]
    out = ['<div class="pstatus">']
    for m, p in P.items():
        cells = "".join(f'<span class="pill {"on" if p[k] else "off"}">{k.split("_")[0]} {"pass" if p[k] else "fail"}</span>'
                        for k in ("P1_above_floor", "P2_pass", "P3_pass", "P4_pass"))
        out.append(f'<div><span class="swatch" style="background:{COL[m]}"></span><strong>{NAME[m]}</strong> {cells}</div>')
    out.append("</div>")
    return "".join(out)


def findings_html():
    if not LADDER.exists():
        return ""
    L = json.load(open(LADDER))
    P = L["predictions"]
    fm = [m for m in ("scgpt", "geneformer") if m in P]
    margins = "; ".join(f'{NAME[m]} {P[m]["P2_margin_over_baseline"]:+.3f} '
                        f'(95% interval {P[m]["P2_margin_ci"][0]:+.3f} to {P[m]["P2_margin_ci"][1]:+.3f})'
                        for m in fm)
    return (f'<div class="call"><p><strong>Neither foundation model beats the linear baseline.</strong> '
            f'In {L["n_patients"]:,} TCGA-BRCA patients with {L["n_events"]} deaths, the best '
            f'signature of every representation stays below its permutation null, so P1 fails for both '
            f'models. The margin over HVG-PCA is {margins}; both intervals include zero, so P2 fails. '
            f'PAM50 subtype alone reaches C = {L["reference_pam50"]:.3f}. The sensitivity analysis below, '
            f'which also credits protective signatures, gives the same answer.</p></div>')


def sensitivity_html():
    if not LADDER.exists() or not json.load(open(LADDER)).get("sensitivity"):
        return ""
    L = json.load(open(LADDER))
    P = L["sensitivity"]["predictions"]
    p3 = []
    for m in ("scgpt", "geneformer"):
        if m in P and P[m]["P3_pass"]:
            st, v = max(P[m]["P3_within_subtype_cindex"].items(), key=lambda kv: kv[1])
            p3.append(f"{NAME[m]} in {st} (C = {v:.3f})")
    p3_txt = (" P3 passes for " + " and ".join(p3) + ". P3 is a single threshold of 0.6 applied "
              "across four subtypes, with no null of its own, and these signatures do not pass P1; "
              "the within-subtype values are leads for replication, not findings.") if p3 else ""
    fig = REPORT / "fig_ladder_sensitivity.png"
    img = f'<figure class="result narrow"><img src="{data_uri(fig)}" alt="Sensitivity ladder"></figure>' if fig.exists() else ""
    return (f'<h2 id="sensitivity">Sensitivity analysis: protective signatures</h2>'
            f'<p class="sub">Added after the primary run. The primary ladder is unchanged.</p>'
            f'<p>The pre-registered concordance index credits a signature only when a higher score means '
            f'worse survival. Many signatures work the other way: a higher score means longer survival. '
            f'They include the largest departures from 0.5 in the run, down to C = 0.40. Here every '
            f'signature is read in the direction it acts. Its floor is oriented the same way, and the '
            f'permutation null takes the best of max(C, 1 − C) over the same 500 permutations.</p>'
            f'{img}{ladder_table(sensitivity=True)}{predictions_status(sensitivity=True)}'
            f'<p>The best state of every representation is a protective malignant-cell state, and every one '
            f'sits just below its null. HVG-PCA is still the strongest.{p3_txt}</p>')


# ----------------------------------------------------------------------------- page
CSS = f"""
:root{{
  --paper:#F7F8F9; --panel:#FFFFFF; --ink:{INK}; --muted:{MUTED}; --rule:{RULE}; --axis:{AXIS};
  --base:{COL['hvg_pca']}; --scgpt:{COL['scgpt']}; --gf:{COL['geneformer']}; --band:#EEF1F2;
  --sans:{SANS}; --mono:{MONO};}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);line-height:1.6;
  -webkit-font-smoothing:antialiased;padding-inline:16px}}
.wrap{{max-width:1140px;margin:0 auto;padding:34px 0 96px}}
.topnav{{position:sticky;top:0;z-index:10;background:color-mix(in srgb,var(--paper) 88%,transparent);
  backdrop-filter:blur(8px);border-bottom:1px solid var(--rule);margin-inline:-16px;padding-inline:16px}}
.topnav .inner{{max-width:1140px;margin:0 auto;display:flex;align-items:center;justify-content:space-between;gap:16px;height:54px}}
.topnav .brand{{font-family:var(--mono);font-size:12px;color:var(--ink);white-space:nowrap}}
.topnav ul{{list-style:none;margin:0;padding:0;display:flex;gap:18px;font-family:var(--mono);font-size:11px;
  letter-spacing:.06em;text-transform:uppercase;overflow-x:auto;scrollbar-width:none}}
.topnav ul::-webkit-scrollbar{{display:none}}
.topnav a{{color:var(--muted);text-decoration:none;white-space:nowrap}}
.topnav a:hover,.topnav a:focus-visible{{color:var(--ink)}}
header{{border-bottom:2px solid var(--ink);padding-bottom:22px;margin-bottom:34px}}
.eyebrow{{font-family:var(--mono);font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);margin:0 0 10px;
  display:flex;gap:14px;flex-wrap:wrap;align-items:center}}
.status{{display:inline-block;border:1px solid var(--gf);color:var(--gf);border-radius:2px;padding:1px 7px;letter-spacing:.1em}}
h1{{font-size:clamp(26px,3.6vw,38px);line-height:1.12;letter-spacing:-.02em;margin:0 0 14px;text-wrap:balance;font-weight:640}}
.lede{{font-size:17px;color:var(--muted);margin:0 0 16px}}
.repo{{font-family:var(--mono);font-size:12px;color:var(--muted)}}
.repo a{{color:var(--ink)}}
h2{{font-size:21px;letter-spacing:-.01em;margin:56px 0 8px;font-weight:640;text-wrap:balance}}
h3{{font-size:14.5px;font-weight:640;margin:24px 0 8px;display:flex;align-items:center;gap:8px}}
.sub{{color:var(--muted);font-size:15px;margin:0 0 20px}}
p{{margin:0 0 14px}}
.panel{{margin:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:14px 12px 8px}}
.panel figcaption{{font-family:var(--mono);font-size:11.5px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);margin:0 0 6px;padding-left:4px}}
.svgscroll{{overflow-x:auto}}
.svgscroll svg{{min-width:760px}}
.panel svg,.result img{{width:100%;height:auto;display:block}}
.fig{{margin:18px 0 10px}}
.arms{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px;margin:18px 0}}
@media (max-width:760px){{.arms{{grid-template-columns:minmax(0,1fr)}}}}
.arm{{background:var(--panel);border:1px solid var(--rule);border-top:4px solid var(--c);border-radius:3px;padding:14px 16px 10px}}
.arm h3{{margin:0 0 6px}}
.arm .k{{font-family:var(--mono);font-size:10.5px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:0 0 8px}}
.arm p{{font-size:13.5px;margin:0 0 8px}}
.arm .meta{{font-family:var(--mono);font-size:11px;color:var(--muted);line-height:1.7}}
.call{{border-left:3px solid var(--scgpt);background:#EEF3F3;padding:16px 22px;margin:22px 0;border-radius:0 3px 3px 0}}
.call.warn{{border-left-color:var(--gf);background:#FAF2EA}}
.call p:last-child{{margin:0}}
.steps{{counter-reset:step;margin:26px 0 0}}
.step{{counter-increment:step;position:relative;padding:0 0 22px 46px;border-left:1px solid var(--rule);margin-left:12px}}
.step:last-child{{border-left-color:transparent;padding-bottom:4px}}
.step::before{{content:counter(step);position:absolute;left:-13px;top:-2px;width:26px;height:26px;border-radius:50%;
  background:var(--panel);border:1px solid var(--rule);color:var(--muted);font-family:var(--mono);font-size:11px;
  display:flex;align-items:center;justify-content:center}}
.step h3{{margin:0 0 6px;font-size:15px}}
.why{{color:var(--muted);font-size:14px;border-left:2px solid var(--rule);padding-left:12px;margin:10px 0 0}}
.preds{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin:18px 0}}
@media (max-width:720px){{.preds{{grid-template-columns:minmax(0,1fr)}}}}
.pred{{background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:14px 16px 10px;display:grid;
  grid-template-columns:minmax(0,1fr) 220px;gap:6px 16px;align-items:start}}
@media (max-width:520px){{.pred{{grid-template-columns:minmax(0,1fr)}}}}
.pred .id{{font-family:var(--mono);font-size:11px;letter-spacing:.1em;color:var(--muted);grid-column:1/-1}}
.pred h3{{margin:0 0 4px;font-size:14.5px}}
.pred p{{font-size:13.5px;margin:0 0 6px}}
.pred .ifnot{{font-size:12.5px;color:var(--muted)}}
.pred svg{{width:220px;max-width:100%;height:auto;display:block;min-width:0}}
.results{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin:18px 0}}
@media (max-width:820px){{.results{{grid-template-columns:minmax(0,1fr)}}}}
.result{{margin:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:12px 12px 10px}}
.result figcaption{{font-family:var(--mono);font-size:11.5px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);margin:0 0 8px;padding-left:2px}}
.result .rtext p{{font-size:12.5px;color:var(--muted);margin:8px 0 0}}
.result .rtext strong{{color:var(--ink);font-weight:600}}
.result.wide{{grid-column:1/-1}}
.result.narrow{{max-width:640px;margin-bottom:14px}}
.result.wide .imgscroll{{overflow-x:auto}}
.result.wide .imgscroll img{{min-width:900px}}
.pending{{aspect-ratio:2/1;border:1.5px dashed var(--axis);border-radius:3px;background:repeating-linear-gradient(135deg,transparent 0 10px,#F2F4F5 10px 20px);
  display:flex;flex-direction:column;align-items:center;justify-content:center;gap:6px;color:var(--muted);font-size:12.5px;text-align:center;padding:12px}}
.pending-k{{font-family:var(--mono);font-size:11px;letter-spacing:.06em;color:var(--ink)}}
.scroll{{overflow-x:auto;border:1px solid var(--rule);border-radius:3px;background:var(--panel);margin:0 0 8px}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{padding:7px 14px;text-align:left;border-bottom:1px solid var(--rule);white-space:nowrap}}
th{{font-family:var(--mono);font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;color:var(--muted);font-weight:500;background:#F2F4F5}}
tbody tr:last-child td{{border-bottom:none}}
tr.ref td{{border-top:1px solid var(--axis);color:var(--muted)}}
.num{{font-family:var(--mono);font-variant-numeric:tabular-nums;text-align:right}}
.dim{{color:var(--axis)}}
.swatch{{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:8px;vertical-align:-1px}}
.pstatus{{display:flex;flex-direction:column;gap:6px;margin:12px 0;font-size:13.5px}}
.pill{{font-family:var(--mono);font-size:10.5px;letter-spacing:.06em;border-radius:2px;padding:1px 7px;margin-left:8px;border:1px solid var(--rule)}}
.pill.on{{background:#EEF3F3;color:var(--scgpt);border-color:var(--scgpt)}}
.pill.off{{color:var(--muted)}}
.rungs{{width:100%;border-collapse:collapse;font-size:13.5px;margin:14px 0 6px}}
.rungs td{{padding:9px 12px;vertical-align:top;border-bottom:1px solid var(--rule);white-space:normal}}
.rungs td:first-child{{font-family:var(--mono);font-size:11px;letter-spacing:.1em;text-transform:uppercase;white-space:nowrap;width:110px;padding-top:12px}}
.rungs tr:last-child td{{border-bottom:none}}
.conv{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin:18px 0}}
@media (max-width:720px){{.conv{{grid-template-columns:minmax(0,1fr)}}}}
.conv div{{background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:14px 16px 8px}}
.conv h3{{margin:0 0 6px;font-size:14px}}
.conv p{{font-size:13.5px;color:var(--muted);margin:0 0 8px}}
.limits li,.refs li{{margin:0 0 8px;font-size:14px}}
.refs{{padding-left:20px;font-size:13.5px}}
.refs li{{color:var(--muted)}}
.refs a{{color:var(--ink)}}
details.every{{border:1px solid var(--rule);border-radius:3px;background:var(--panel);margin:14px 0}}
details.every>summary{{cursor:pointer;padding:9px 13px;font-size:12px;color:var(--muted);list-style:none;display:flex;justify-content:space-between;align-items:center;gap:12px}}
details.every>summary::-webkit-details-marker{{display:none}}
details.every>summary::after{{content:"show";font-family:var(--mono);font-size:10px;letter-spacing:.06em;color:var(--muted);border:1px solid var(--rule);border-radius:2px;padding:1px 6px}}
details.every[open]>summary::after{{content:"hide"}}
details.every>summary:focus-visible,a:focus-visible{{outline:2px solid var(--scgpt);outline-offset:2px}}
details.every .inner{{padding:0 13px 13px}}
code{{font-family:var(--mono);font-size:.9em;background:#EDEFF1;padding:1px 5px;border-radius:3px}}
footer{{margin-top:64px;padding-top:18px;border-top:1px solid var(--rule);font-family:var(--mono);font-size:11.5px;color:var(--muted);line-height:1.8}}
footer a{{color:var(--ink)}}
@media (prefers-reduced-motion:reduce){{*{{scroll-behavior:auto!important}}}}
"""


def panel(fig_html, caption):
    return f'<figure class="panel fig"><figcaption>{caption}</figcaption><div class="svgscroll">{fig_html}</div></figure>'


def build():
    today = date.today().isoformat()
    have_results = LADDER.exists()
    status = "results in · negative" if have_results else "design complete · results pending"

    arms = [
        ("hvg_pca", "The baseline: no learning",
         "Keep the ~2,000 genes that vary most between cells and summarise them with principal components. This is what a careful analyst did before foundation models existed, and what the sceptics say is still enough.",
         "linear · fitted on this atlas only · 50 components"),
        ("scgpt", "A foundation model",
         "A transformer trained on 33 million human cells to fill in masked gene expression. It was never told what a cell type, a tumour, or a patient is. We ask what it picked up anyway.",
         "Cui et al. 2024 · whole-human checkpoint · 512-dim cell embedding"),
        ("geneformer", "A foundation model",
         "A transformer trained on ~30 million cells to predict which genes rank highest in each cell. Same self-supervised idea as scGPT, different training signal and gene encoding — so agreement between the two is informative.",
         "Theodoris et al. 2023 · 6-layer checkpoint · mean-pooled cell embedding"),
    ]
    arms_html = "".join(
        f'<div class="arm" style="--c:{COL[m]}"><p class="k">{k}</p><h3>{NAME[m]}</h3><p>{d}</p><p class="meta">{meta}</p></div>'
        for m, k, d, meta in arms)

    preds = [
        ("P1", "Foundation-model cell states carry prognostic biology at all",
         "The best cell state from scGPT or Geneformer predicts survival better than random gene sets of the same size and expression level, and better than the shuffled-outcome null.",
         "If not: the cell states these models find in a tumour say nothing about how the patient fares. Everything downstream is moot."),
        ("P2", "They add something over the linear baseline",
         "The gap between the best foundation-model state and the best HVG-PCA state, recomputed on 200 resamples of the patients, stays above zero.",
         "If not: on a clinical endpoint the linear method is sufficient — the sceptics' claim, stated on a harder test than cluster maps."),
        ("P3", "At least one state is new biology, not PAM50 rediscovered",
         "Within a single PAM50 subtype — luminal A, luminal B, HER2-enriched, basal — the state still separates patients who do well from those who do not.",
         "If not: the model has found the subtypes the field has used since 2009 and nothing more. Useful as a sanity check, not as a discovery."),
        ("P4", "The informative genes lie outside the highly variable set",
         "Signatures that pass P2 draw most of their marker genes from outside the ~2,000 HVGs the baseline was built on.",
         "If so: the claim that these models learn from genes a linear method discards has a concrete, checkable instance."),
    ]
    preds_html = "".join(
        f'<div class="pred"><span class="id">{pid}</span><div><h3>{t}</h3><p>{d}</p><p class="ifnot">{n}</p></div><div>{fig_prediction(pid)}</div></div>'
        for pid, t, d, n in preds)

    steps = [
        ("Assemble one tumour atlas and one patient cohort",
         "Breast-cancer single-cell data from CELLxGENE Census: quality-filter the cells, flag the highly variable genes, subsample to at most 50,000 cells keeping every cell type represented. Separately, TCGA-BRCA: ~1,100 primary tumours with bulk RNA, age, stage, PAM50 call and follow-up.",
         "The two cohorts share no patients, so every survival result is an out-of-cohort test. Breast cancer is the largest TCGA cohort with survival and the one with a reference the field already accepts."),
        ("Summarise every cell three ways",
         "Each cell gets three coordinates: an HVG-PCA position, a scGPT embedding and a Geneformer embedding. Same cells, same counts in — three views out.",
         "Everything that follows is done identically for the three views. The comparison is only fair because nothing else differs."),
        ("Find cell states in each view",
         "Group similar cells in each view (Leiden clustering at three resolutions), keep groups of at least 200 cells, and record which annotated cell type each group mostly contains.",
         "A cell state is a group of cells the view considers alike. Whether the groups mean anything is exactly what the later steps test."),
        ("Turn each state into a gene signature",
         "For each group, the 25, 50 or 100 genes most specifically expressed in it, compared with all other cells in the same view.",
         "A signature is the portable form of a cell state: a list of genes that can be scored in any tumour, including a bulk one where cells were never separated."),
        ("Score the signatures in the patient cohort and ask about survival",
         "Each patient's tumour gets a score per signature — the mean expression of its genes, standardised. A Cox model adjusted for age and stage relates score to overall survival; its concordance index says how well the score orders patients by outcome.",
         "Bulk expression of a signature is a proxy for how much of that cell state a tumour contains. It is stated as a proxy throughout."),
        ("Read every number against the ladder",
         "For each view take its best signature — best by margin over its own random-gene floor, not by raw score — and place it on the ladder with a patient-bootstrap interval. Test P1–P4.",
         "The study reports one number: the margin of the foundation models over the baseline, both read above their floors. No absolute C-index is a finding on its own."),
    ]
    steps_html = "".join(f'<div class="step"><h3>{t}</h3><p>{d}</p><p class="why">{w}</p></div>' for t, d, w in steps)

    rungs_rows = [
        ("null", "Survival outcomes shuffled across patients, 500 times. What any signature scores by luck. Per view, the null is the distribution of the <em>best</em> of its signatures, so picking the winner of many is already corrected for."),
        ("floor", "200 random gene sets matched to each signature for size and mean expression level. What any gene list of that shape scores. This is the honest zero: highly expressed genes predict survival for reasons that have nothing to do with cell states."),
        ("baseline", "Cell states from HVG-PCA. What the linear, learning-free method achieves. The sceptics' position, made concrete."),
        ("test", "Cell states from scGPT and Geneformer. The claim under test."),
        ("reference", "PAM50 intrinsic subtype, the classifier clinicians already use. Where the field is. Test against reference says whether either method is clinically interesting; it does not decide the question."),
    ]
    rungs_html = "".join(f'<tr><td>{r}</td><td>{d}</td></tr>' for r, d in rungs_rows)

    results_html = "".join([
        result_slot("fig_ladder", "A · the ladder",
                    "One bar per representation: the concordance index of its best cell-state signature, with a 95% interval from resampling patients. Dotted tick: that signature's matched-random floor. Thick tick: the 95th percentile of the permutation null for the best of that representation's signatures. Dashed line: PAM50.",
                    "A bar that clears both ticks passes P1. The distance between the scGPT / Geneformer bars and the HVG-PCA bar is the study's result."),
        result_slot("fig_margin", "B · the margin over the baseline",
                    "For each foundation model, the difference in concordance from the HVG-PCA baseline, computed on the same 200 patient resamples so the comparison is paired.",
                    "An interval wholly to the right of zero passes P2. An interval straddling zero is the sceptics' result, and is reported as such."),
        result_slot("fig_hvg", "C · where the informative genes come from",
                    "Every signature from every representation: the fraction of its genes that are HVGs (x) against its concordance above the matched-random floor (y).",
                    "Foundation-model points high on the y-axis and left of 0.5 on the x-axis are P4: prognostic signal built from genes the baseline never saw."),
        result_slot("fig_km", "D · survival curves, overall and within subtype",
                    "Kaplan–Meier curves for patients split into thirds by the best foundation-model signature: all patients, then each PAM50 subtype with enough events.",
                    "Separated curves inside a single subtype panel is P3: the state orders patients that PAM50 puts in one box. Curves that separate only in the all-patients panel mean the state is tracking subtype. Thirds are cut on all patients, so a small subtype can lack one of them.", wide=True),
    ])

    refs = [
        ('Kedzierska KZ, Crawford L, Amini AP, Lu AX. Assessing the limits of zero-shot foundation models in single-cell biology. <em>Genome Biology</em> 2025.', "https://doi.org/10.1186/s13059-025-03574-x"),
        ('Boiarsky R, Singh N, Buendia A, Getz G, Uhler C. A deep dive into single-cell RNA sequencing foundation models. bioRxiv 2023.', "https://doi.org/10.1101/2023.10.19.563100"),
        ('Souza &amp; Mehta. Parameter-free representations outperform single-cell foundation models on downstream benchmarks. arXiv 2026.', "https://arxiv.org/abs/2602.16696"),
        ('Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. <em>Nature Methods</em> 2024.', "https://doi.org/10.1038/s41592-024-02201-0"),
        ('Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. <em>Nature</em> 2023.', "https://doi.org/10.1038/s41586-023-06139-9"),
        ('Rosen Y, Roohani Y, Agarwal A, et al. Universal cell embedding provides a foundation model for cell biology. <em>Nature</em> 2026.', "https://doi.org/10.1038/s41586-026-10689-z"),
        ('Parker JS, Mullins M, Cheang MCU, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. <em>J Clin Oncol</em> 2009.', "https://doi.org/10.1200/JCO.2008.18.1370"),
        ('The Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. <em>Nature</em> 2012.', "https://doi.org/10.1038/nature11412"),
        ('Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. <em>Cell</em> 2018. (TCGA-CDR)', "https://doi.org/10.1016/j.cell.2018.02.052"),
        ('Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. <em>Nature Biotechnology</em> 2020.', "https://doi.org/10.1038/s41587-020-0546-8"),
        ('CZI Single-Cell Biology Program. CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. <em>Nucleic Acids Research</em> 2025.', "https://doi.org/10.1093/nar/gkae1142"),
        ('Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA. Evaluating the yield of medical tests. <em>JAMA</em> 1982. (concordance index)', "https://doi.org/10.1001/jama.1982.03320430047030"),
        ('Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. <em>Scientific Reports</em> 2019.', "https://doi.org/10.1038/s41598-019-41695-z"),
    ]
    refs_html = "".join(f'<li>{t} <a href="{u}">{u.replace("https://", "")}</a></li>' for t, u in refs)

    nav = "".join(f'<li><a href="#{a}">{t}</a></li>' for a, t in [
        ("question", "Question"), ("compared", "What is compared"), ("design", "Design"), ("ladder", "Ladder"),
        ("predictions", "Predictions"), ("results", "Results"), ("sensitivity", "Sensitivity"),
        ("conventions", "Conventions"),
        ("outcomes", "Outcomes"), ("validation", "Validation"), ("limits", "Limits"), ("refs", "References")])

    html = f"""<title>Cell states versus survival</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>{CSS}</style>
<nav class="topnav"><div class="inner"><span class="brand">cell states → survival</span><ul>{nav}</ul></div></nav>
<div class="wrap">
<header>
<p class="eyebrow"><span>Single-cell foundation models · breast cancer · pre-registered design</span><span class="status">{status}</span></p>
<h1>Do the cell states a foundation model finds in a tumour predict how patients fare — better than a plain linear method?</h1>
<p class="lede">Two transformer models trained on tens of millions of cells, one linear baseline, one breast-cancer atlas, and ~1,100 patients they have never seen. Every result is read against a ladder of controls, and the study reports a single margin.</p>
<p class="repo">Code, design document and pipeline: <a href="{REPO}">github.com/Natalija-Stepurko/single-cell-fm-probing</a>. Every stage writes its parameters, command and git commit beside its outputs.</p>
</header>

<h2 id="question">The question, and why it is still open</h2>
<p class="sub">Single-cell foundation models learn to reconstruct hidden gene expression across enormous cell collections. Whether they learn any biology beyond what a linear summary of the most variable genes already captures is disputed — and so far the dispute has been fought on how well each method's cell map matches annotated cell types.</p>
{panel(fig_dispute(), "The dispute and where this study moves it")}
<div class="call"><p><strong>The move.</strong> Cell-type labels are a soft target: the annotations were made from the same kind of data the models see, and a method can match them by capturing what is obvious. Whether a patient lives or dies is not soft, is recorded independently of any expression data, and is what a tumour cell state would need to inform to matter. If the foundation models have learned something a linear method has not, it should show here. If it does not show here, that is a finding.</p></div>

<h2 id="compared">What is compared</h2>
<p class="sub">Three ways of summarising a cell. All three are handed the same cells; the only difference is what each does with them.</p>
<div class="arms">{arms_html}</div>
{panel(fig_cohorts(), "Two cohorts, no shared patients")}
<p>The atlas is where cell states are found. The bulk cohort is where they are tested. A bulk tumour sample is a mixture of many cells, so a cell state is tested there as a gene signature: the more of a state a tumour contains, the higher its signature scores.</p>

<h2 id="design">What is done, step by step</h2>
<p class="sub">Six steps, run identically for the three representations. The figure reads left to right; the numbered notes below say what each step is for.</p>
{panel(fig_flow(), "Design — one atlas, three views, one patient cohort, one ladder")}
<p><strong>Every view sees the same cells and is scored in the same patients.</strong> The only thing that differs between the three rows is how the cells were summarised, so any difference in how well their cell states predict survival is a difference in what each method learned. The baseline row is what no learning at all achieves.</p>
<div class="steps">{steps_html}</div>

<h2 id="ladder">The ladder: what every number is read against</h2>
<p class="sub">A concordance index, on its own, means nothing. 0.65 might be a strong result or what any random handful of highly expressed genes achieves. So every score sits on a ladder, and the findings are gaps between rungs, not rungs.</p>
{panel(fig_ladder_schematic(), "The ladder — rung positions are illustrative until the data land")}
<table class="rungs"><tbody>{rungs_html}</tbody></table>
<h3>The trap the baseline rung guards against</h3>
{panel(fig_shared_input(), "One input, three views, one scoring set")}
<p>The three representations are not independent witnesses. They were built from the same expression matrix and their signatures are scored in the same bulk genes. Prognostic signal that appears in all three can come from that shared input — from cell-type composition, from tumour purity, from proliferation, and none of that reflects what a model learned. The baseline rung is the control: it is what the shared input yields when no model learns anything. Every claim about the foundation models is therefore the margin over that rung.</p>

<h2 id="predictions">Four predictions, written down before running</h2>
<p class="sub">Each is a yes-or-no test with the failing outcome spelled out. The study can lose.</p>
<div class="preds">{preds_html}</div>
{predictions_status()}

<h2 id="results">Results</h2>
{findings_html()}
<p class="sub">{"The four panels stage 06 writes, and the ladder table. Everything in this block is the pre-registered analysis." if have_results else "The full run is in progress. The four panels below are the ones stage 06 will write, each with what it will show and how to read it; the table has the shape of the final ladder. When the data land, this section fills in and nothing else on the page changes."}</p>
<div class="results">{results_html}</div>
<h3>The ladder as a table</h3>
{ladder_table()}
{sensitivity_html()}

<h2 id="conventions">Statistical conventions</h2>
<p class="sub">Each of these was learned the hard way on an earlier study and is fixed here before any data are seen.</p>
<div class="conv">
<div><h3>The patient is the unit</h3><p>All intervals come from resampling patients. Resampling genes or cells would count the same patient many times and shrink every interval below its true width.</p></div>
<div><h3>Above the null and above the floor</h3><p>A signature must beat both the shuffled-outcome null and random gene sets matched for size and expression. The floor is matched because prognostic signal correlates with both.</p></div>
<div><h3>Best-of-many is corrected at the family level</h3><p>Each representation yields dozens of signatures. The null for each is the distribution of its <em>best</em> signature under permuted outcomes, so the winner is not judged against a null for a single, pre-chosen one.</p></div>
<div><h3>Adjust for age and stage; stratify by subtype separately</h3><p>The primary model adjusts for age and stage. Subtype-specific models are P3 and are reported on their own, so a state that merely tracks subtype cannot pass as new biology.</p></div>
<div><h3>Best by margin, not by raw score</h3><p>The signature that represents each view is the one furthest above its own floor. A high raw score built from highly expressed genes does not qualify.</p></div>
<div><h3>Signature choices are swept, not tuned</h3><p>Three clustering resolutions × three marker-gene counts. The ladder is reported at every setting; no setting is chosen after seeing the survival result.</p></div>
</div>

<h2 id="outcomes">What each outcome would mean</h2>
<p class="sub">Three shapes the ladder can take. Each is a result; none is a failure of the study.</p>
{panel(fig_outcomes(), "Three possible ladders")}
<table class="rungs"><tbody>
<tr><td>FM wins</td><td>Test above baseline with an interval excluding zero, both above floor; at least one state prognostic within a PAM50 subtype; its genes drawn from outside the HVG set. The shortlisted states go to validation.</td></tr>
<tr><td>Sceptics win</td><td>Test ≈ baseline, both below the PAM50 reference. On a clinical endpoint the linear method is sufficient. Reported as the finding, with the margin and its interval.</td></tr>
<tr><td>Both fail</td><td>Test ≈ baseline ≈ floor. The route from cell state to bulk signature loses the signal; this design cannot separate the two poles, and the report says so and does not look for a setting where it can.</td></tr>
</tbody></table>

<h2 id="validation">What a positive result would still need</h2>
<p class="sub">A cell state that passes P1–P3 in TCGA is a hypothesis, not a biomarker. Stage 06 writes a validation plan for each shortlisted state with three fixed steps.</p>
<div class="steps">
<div class="step"><h3>Replicate in an independent patient cohort</h3><p>Re-score the signature in METABRIC (~2,000 patients, overall survival) with the same age- and stage-adjusted Cox model. A state that does not replicate is dropped.</p><p class="why">TCGA is one cohort with one sequencing platform and ~15% events. A single-cohort survival result is where most published signatures stop, and where most of them fail.</p></div>
<div class="step"><h3>Go back to the atlas</h3><p>Check that the signature's genes are expressed in the annotated cell type the state was assigned to, not in a contaminating population, and that the state is present in several donors.</p><p class="why">A signature can score well in bulk for a reason that has nothing to do with the cells it was named after.</p></div>
<div class="step"><h3>Name the tissue assay</h3><p>Specify the multiplexed immunohistochemistry or spatial panel that would show the state's abundance in sections from a cohort with outcome, and the effect size it would need to reach to matter clinically.</p><p class="why">A gene list becomes a finding when it can be seen in tissue and changes a decision. Writing the assay down is the test of whether the state is concrete enough to be one.</p></div>
</div>

<h2 id="limits">Limits, stated in advance</h2>
<ul class="limits">
<li><strong>A bulk signature score is a proxy.</strong> It stands in for the abundance of a cell state and is contaminated by everything else that changes bulk expression — purity, proliferation, stroma.</li>
<li><strong>Signature derivation carries choices.</strong> Cluster resolution and marker-gene count set the signatures. Both are swept and reported; neither is chosen after seeing survival.</li>
<li><strong>Power.</strong> TCGA-BRCA has ~15% deaths in follow-up. Progression-free interval is the secondary endpoint with more events; both are reported.</li>
<li><strong>The atlas and the cohort differ in subtype mix.</strong> P3 addresses this directly by testing within subtype.</li>
<li><strong>Light checkpoints on CPU.</strong> The scGPT whole-human and Geneformer 6-layer models, not the largest released. A negative result is a result about these checkpoints.</li>
<li><strong>One indication.</strong> Breast cancer was chosen for cohort size and a fixed reference. Nothing here generalises to other tumours until it is run there.</li>
</ul>

<h2 id="refs">References</h2>
<ol class="refs">{refs_html}</ol>

<footer>
<div>Data: CELLxGENE Census release 2025-11-08 (breast-carcinoma donors, primary data, 10x assays) · TCGA-BRCA via UCSC Xena (HiSeqV2, clinical matrix, TCGA-CDR survival). Both public.</div>
<div>Pipeline: <a href="{REPO}">single-cell-fm-probing</a> · stages 01–06 · Python, scanpy, lifelines · result panels are written by <code>scripts/06_report.py</code> and embedded by <code>docs/site/build.py</code>.</div>
<div>Page built {today}. {"Results embedded." if have_results else "No results embedded."}</div>
</footer>
</div>
"""
    # head and body tags are optional in HTML; the page stays one self-contained file
    OUT.write_text('<!doctype html>\n<html lang="en">\n<meta charset="utf-8">\n'
                   '<meta name="description" content="Do single-cell foundation models find cell '
                   'states that predict breast-cancer survival better than a linear baseline?">\n'
                   + html + "\n</html>\n")
    print(f"wrote {OUT}  {OUT.stat().st_size / 1024:.0f} KB  results={'yes' if have_results else 'no'}")


if __name__ == "__main__":
    build()
