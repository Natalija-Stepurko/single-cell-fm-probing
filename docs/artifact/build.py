"""Build the portfolio page for this study as a single self-contained HTML file.

The page is the pre-registered design, written before any result exists. Result panels are
slots: when `results/report/fig_*.png` and `results/ladder/ladder.json` exist they are embedded
(PNGs as data URIs, the ladder as a table); until then each slot shows what the panel will
contain and how to read it. Re-run after stage 06 to fill the page; nothing else changes.

    python3 docs/artifact/build.py            -> docs/artifact/cell_states_survival.html
"""
import base64
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import config as C  # noqa: E402

OUT = Path(__file__).parent / "cell_states_survival.html"
REPORT = C.RESULTS / "report"
LADDER = C.RESULTS / "ladder" / "ladder.json"
REPO = "https://github.com/Natalija-Stepurko/single-cell-fm-probing"

COL = {"hvg_pca": "#4A4F55", "scgpt": "#1F6F6B", "geneformer": "#C2681A"}
NAME = {"hvg_pca": "HVG-PCA (baseline)", "scgpt": "scGPT", "geneformer": "Geneformer"}
INK, MUTED, RULE, AXIS, PANEL = "#16191D", "#5B646E", "#DDE1E4", "#B9C0C6", "#FFFFFF"
MONO = "ui-monospace,SFMono-Regular,Menlo,Consolas,monospace"
SANS = "ui-sans-serif,system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif"


# ----------------------------------------------------------------------------- svg helpers
def T(x, y, s, size=12, fill=INK, anchor="start", weight=400, mono=False, dy=0, italic=False):
    fam = MONO if mono else SANS
    s = s.replace("& ", "&amp; ")
    st = f"font-family:{fam};font-size:{size}px;fill:{fill};font-weight:{weight}"
    if italic:
        st += ";font-style:italic"
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
    W, H = 980, 210
    b = []
    # two poles
    b.append(box(20, 24, 430, 118, fill="#F2F4F5"))
    b.append(T(38, 46, "SCEPTICS", 10.5, MUTED, mono=True, weight=600))
    b.append(T(38, 72, "A linear summary of ~2,000 variable genes", 14))
    b.append(T(38, 92, "matches or beats the foundation models", 14))
    b.append(T(38, 118, "Kedzierska 2025 · Boiarsky 2023 · Souza & Mehta 2026", 11, MUTED, mono=True))
    b.append(box(530, 24, 430, 118, fill="#F2F4F5"))
    b.append(T(548, 46, "OPTIMISTS", 10.5, MUTED, mono=True, weight=600))
    b.append(T(548, 72, "The models learn cell biology that", 14))
    b.append(T(548, 92, "transfers zero-shot to unseen data", 14))
    b.append(T(548, 118, "Cui 2024 · Theodoris 2023 · Rosen 2026", 11, MUTED, mono=True))
    # the arena so far vs here
    b.append(T(490, 60, "argued on", 11, MUTED, anchor="middle", italic=True))
    b.append(T(490, 80, "cluster maps", 12.5, INK, anchor="middle", weight=600))
    b.append(T(490, 98, "cell-type labels", 12.5, INK, anchor="middle", weight=600))
    b.append(arrow(490, 150, 490, 172))
    b.append(T(490, 190, "This study asks the same question on a clinical endpoint: patient survival.", 13.5, INK, anchor="middle", weight=600))
    return svg(W, H, "".join(b), "The two poles of the dispute and where this study moves it")


def fig_flow():
    """The whole design, left to right: atlas -> three representations -> states -> signatures -> bulk -> survival -> ladder."""
    W, H = 980, 400
    b = []
    # column x positions
    xs = {"atlas": 20, "rep": 190, "states": 380, "sig": 560, "bulk": 720, "out": 880}
    # atlas
    b.append(box(xs["atlas"], 120, 120, 150))
    b.append(cells_glyph(xs["atlas"] + 30, 140, n=34, seed=5, fill=MUTED))
    b.append(cells_glyph(xs["atlas"] + 30, 175, n=22, seed=9, fill="#8A9298"))
    b.append(T(xs["atlas"] + 60, 232, "breast tumour", 12, INK, anchor="middle", weight=600))
    b.append(T(xs["atlas"] + 60, 248, "single-cell atlas", 12, INK, anchor="middle", weight=600))
    b.append(T(xs["atlas"] + 60, 100, "1  ONE INPUT", 10.5, MUTED, anchor="middle", mono=True, weight=600))
    b.append(T(xs["atlas"] + 60, 290, "CELLxGENE Census", 10.5, MUTED, anchor="middle", mono=True))
    b.append(T(xs["atlas"] + 60, 304, "≤50,000 cells", 10.5, MUTED, anchor="middle", mono=True))
    # three representations
    b.append(T(xs["rep"] + 70, 100, "2  THREE WAYS TO SEE IT", 10.5, MUTED, anchor="middle", mono=True, weight=600))
    ys = {"hvg_pca": 120, "scgpt": 190, "geneformer": 260}
    for m, y in ys.items():
        b.append(box(xs["rep"], y, 140, 52, stroke=COL[m], sw=1.5))
        b.append(f'<rect x="{xs["rep"]}" y="{y}" width="5" height="52" rx="2" fill="{COL[m]}"/>')
        b.append(T(xs["rep"] + 16, y + 18, NAME[m].replace(" (baseline)", ""), 12.5, INK, weight=600))
        sub = {"hvg_pca": "linear · no learning", "scgpt": "foundation model", "geneformer": "foundation model"}[m]
        b.append(T(xs["rep"] + 16, y + 36, sub, 10.5, MUTED, mono=True))
        b.append(arrow(xs["atlas"] + 120, 195, xs["rep"] - 4, y + 26))
    # states
    b.append(T(xs["states"] + 70, 100, "3  CELL STATES", 10.5, MUTED, anchor="middle", mono=True, weight=600))
    for m, y in ys.items():
        b.append(box(xs["states"], y, 140, 52))
        # clusters glyph
        import random
        rnd = random.Random({"hvg_pca": 1, "scgpt": 2, "geneformer": 3}[m])
        for k in range(3):
            cx, cy = xs["states"] + 22 + k * 24, y + 26
            for _ in range(7):
                b.append(f'<circle cx="{cx + rnd.uniform(-7, 7):.1f}" cy="{cy + rnd.uniform(-9, 9):.1f}" r="2" fill="{COL[m]}" opacity=".75"/>')
        b.append(T(xs["states"] + 96, y + 20, "groups of", 10.5, MUTED))
        b.append(T(xs["states"] + 96, y + 34, "similar cells", 10.5, MUTED))
        b.append(arrow(xs["rep"] + 140, y + 26, xs["states"] - 4, y + 26))
    # signatures
    b.append(T(xs["sig"] + 55, 100, "4  MARKER GENES", 10.5, MUTED, anchor="middle", mono=True, weight=600))
    for m, y in ys.items():
        b.append(box(xs["sig"], y, 110, 52))
        for i in range(4):
            wl = [70, 52, 62, 40][i]
            b.append(f'<rect x="{xs["sig"] + 14}" y="{y + 10 + i * 9}" width="{wl}" height="4" rx="2" fill="{COL[m]}" opacity="{.9 - i * .18}"/>')
        b.append(arrow(xs["states"] + 140, y + 26, xs["sig"] - 4, y + 26))
    # bulk cohort
    b.append(T(xs["bulk"] + 60, 100, "5  A SECOND COHORT", 10.5, MUTED, anchor="middle", mono=True, weight=600))
    b.append(box(xs["bulk"], 120, 120, 150))
    # patients glyph: grid of small squares
    for i in range(6):
        for j in range(5):
            b.append(f'<rect x="{xs["bulk"] + 18 + i * 15}" y="{138 + j * 15}" width="9" height="9" rx="1" fill="{MUTED}" opacity=".55"/>')
    b.append(T(xs["bulk"] + 60, 232, "TCGA breast", 12, INK, anchor="middle", weight=600))
    b.append(T(xs["bulk"] + 60, 248, "tumours · survival", 12, INK, anchor="middle", weight=600))
    b.append(T(xs["bulk"] + 60, 290, "~1,100 patients", 10.5, MUTED, anchor="middle", mono=True))
    b.append(T(xs["bulk"] + 60, 304, "no overlap with atlas", 10.5, MUTED, anchor="middle", mono=True))
    for m, y in ys.items():
        b.append(arrow(xs["sig"] + 110, y + 26, xs["bulk"] - 4, 195))
    # outcome
    b.append(T(xs["out"] + 40, 100, "6  READ", 10.5, MUTED, anchor="middle", mono=True, weight=600))
    b.append(box(xs["out"], 120, 80, 150))
    # tiny ladder
    for i, (lab, yy, c) in enumerate([("ref", 150, INK), ("test", 180, COL["scgpt"]), ("base", 210, COL["hvg_pca"]), ("floor", 240, AXIS)]):
        b.append(line(xs["out"] + 14, yy, xs["out"] + 66, yy, stroke=c, sw=2 if lab != "floor" else 1.2, dash="" if lab != "floor" else "3 3"))
        b.append(T(xs["out"] + 40, yy - 9, lab, 9.5, MUTED, anchor="middle", mono=True))
    b.append(arrow(xs["bulk"] + 120, 195, xs["out"] - 4, 195))
    b.append(T(xs["out"] + 40, 290, "against a", 10.5, MUTED, anchor="middle", mono=True))
    b.append(T(xs["out"] + 40, 304, "ladder", 10.5, MUTED, anchor="middle", mono=True))
    # bottom caption
    b.append(line(20, 340, 960, 340))
    b.append(T(20, 362, "Every arm sees the same cells and is scored in the same patients. The only thing that differs between the three rows is how the cells were summarised.", 12.5, INK))
    b.append(T(20, 382, "So any difference in how well their cell states predict survival is a difference in what each method learned — and the baseline row is what no learning at all achieves.", 12.5, MUTED))
    return svg(W, H, "".join(b), "Design flow from one atlas through three representations to survival in a second cohort")


def fig_ladder_schematic():
    """The ladder as a horizontal axis of concordance with rungs and the reported gaps."""
    W, H = 980, 300
    x0, x1 = 90, 900
    axis_y = 200

    def X(c):  # map concordance 0.5..0.8
        return x0 + (c - 0.5) / 0.3 * (x1 - x0)
    b = []
    b.append(line(x0, axis_y, x1, axis_y, stroke=INK, sw=1.2))
    for c in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8):
        b.append(line(X(c), axis_y, X(c), axis_y + 6, stroke=INK))
        b.append(T(X(c), axis_y + 20, f"{c:.2f}", 11, MUTED, anchor="middle", mono=True))
    b.append(T((x0 + x1) / 2, axis_y + 44, "concordance index — chance that, of two patients, the one the signature calls higher-risk is the one who dies first", 12, MUTED, anchor="middle"))
    b.append(T(X(0.5), axis_y + 62, "0.50 = coin flip", 11, MUTED, anchor="middle", mono=True))
    rungs = [
        ("null", 0.53, AXIS, "6 4", "outcomes shuffled across patients"),
        ("floor", 0.57, MUTED, "3 3", "random gene sets of the same size and expression"),
        ("baseline", 0.63, COL["hvg_pca"], "", "HVG-PCA cell states"),
        ("test", 0.69, COL["scgpt"], "", "scGPT / Geneformer cell states"),
        ("reference", 0.74, INK, "8 4", "PAM50 subtype — where the field already is"),
    ]
    for i, (lab, c, col, dash, desc) in enumerate(rungs):
        x = X(c)
        top = 40 + i * 0  # all rungs same height
        b.append(line(x, 60, x, axis_y - 2, stroke=col, sw=2.2 if not dash else 1.6, dash=dash))
        b.append(T(x, 44, lab.upper(), 10.5, col if col != AXIS else MUTED, anchor="middle", mono=True, weight=600))
        b.append(T(x, 30, desc if lab in ("baseline", "test") else "", 10.5, MUTED, anchor="middle"))
    # descriptions along a second row for the dashed rungs
    b.append(T(X(0.53), 80, "outcomes shuffled", 10, MUTED, anchor="middle"))
    b.append(T(X(0.57), 80, "matched random genes", 10, MUTED, anchor="middle"))
    b.append(T(X(0.74), 80, "PAM50 subtype", 10, MUTED, anchor="middle"))
    # gaps
    yg = 130
    b.append(f'<path d="M{X(0.57)},{yg} L{X(0.69)},{yg}" stroke="{COL["scgpt"]}" stroke-width="1.4" marker-end="url(#arr)" marker-start="url(#arr)"/>')
    b.append(T((X(0.57) + X(0.69)) / 2, yg - 12, "P1 · test above floor", 11, INK, anchor="middle", mono=True))
    yg2 = 165
    b.append(f'<path d="M{X(0.63)},{yg2} L{X(0.69)},{yg2}" stroke="{INK}" stroke-width="1.8" marker-end="url(#arr)" marker-start="url(#arr)"/>')
    b.append(T((X(0.63) + X(0.69)) / 2, yg2 - 12, "P2 · the margin the study reports", 11, INK, anchor="middle", mono=True, weight=600))
    b.append(T(x1 + 4, 60, "positions are", 9.5, MUTED, mono=True))
    b.append(T(x1 + 4, 72, "illustrative", 9.5, MUTED, mono=True))
    return svg(W, H, "".join(b), "Schematic ladder of concordance rungs: null, floor, baseline, test, reference")


def fig_shared_input():
    W, H = 980, 230
    b = []
    # matrix
    b.append(box(40, 40, 150, 140, fill="#F2F4F5"))
    import random
    rnd = random.Random(11)
    for i in range(10):
        for j in range(9):
            o = rnd.choice([.1, .2, .35, .55, .8])
            b.append(f'<rect x="{50 + i * 13}" y="{50 + j * 13}" width="11" height="11" fill="{INK}" opacity="{o}"/>')
    b.append(T(115, 196, "the same cells × genes", 11.5, INK, anchor="middle", weight=600))
    b.append(T(115, 212, "expression matrix", 11.5, INK, anchor="middle", weight=600))
    ys = {"hvg_pca": 50, "scgpt": 100, "geneformer": 150}
    for m, y in ys.items():
        b.append(arrow(190, 110, 300, y + 15))
        b.append(box(304, y, 150, 30, stroke=COL[m], sw=1.5))
        b.append(T(379, y + 15, NAME[m].replace(" (baseline)", ""), 12, INK, anchor="middle", weight=600))
        b.append(arrow(454, y + 15, 560, y + 15))
        b.append(box(564, y, 120, 30))
        b.append(T(624, y + 15, "signature", 11.5, MUTED, anchor="middle", mono=True))
        b.append(arrow(684, y + 15, 760, 110))
    b.append(box(764, 80, 190, 60, fill="#F2F4F5"))
    b.append(T(859, 100, "scored in the same", 11.5, INK, anchor="middle"))
    b.append(T(859, 118, "bulk tumour genes", 11.5, INK, anchor="middle"))
    b.append(T(500, 200, "Signal that appears in all three rows can come from the shared input, not from anything a model learned.", 12.5, INK, anchor="middle"))
    b.append(T(500, 218, "The baseline row is what the shared input yields with no learned model — so every claim is the margin over it.", 12.5, MUTED, anchor="middle"))
    return svg(W, H, "".join(b), "All three representations are built from one expression matrix and scored on the same genes")


def fig_cohorts():
    W, H = 980, 150
    b = []
    b.append(box(40, 20, 400, 100))
    b.append(cells_glyph(60, 30, n=40, seed=2, fill=MUTED))
    b.append(cells_glyph(110, 36, n=30, seed=7, fill="#8A9298"))
    b.append(T(200, 50, "single-cell atlas", 13, INK, weight=600))
    b.append(T(200, 70, "tens of thousands of cells from a few dozen tumours", 11.5, MUTED))
    b.append(T(200, 90, "used to find cell states and their marker genes", 11.5, MUTED))
    b.append(T(240, 135, "CELLxGENE Census · disease = breast cancer · primary data", 10.5, MUTED, anchor="middle", mono=True))
    b.append(box(540, 20, 400, 100))
    for i in range(10):
        for j in range(4):
            b.append(f'<rect x="{560 + i * 13}" y="{34 + j * 17}" width="9" height="12" rx="1" fill="{MUTED}" opacity=".55"/>')
    b.append(T(700, 50, "bulk tumour cohort", 13, INK, weight=600))
    b.append(T(700, 70, "~1,100 patients · one RNA profile each · follow-up", 11.5, MUTED))
    b.append(T(700, 90, "used only to test the signatures against survival", 11.5, MUTED))
    b.append(T(740, 135, "TCGA-BRCA via UCSC Xena · OS and PFI · age, stage, PAM50", 10.5, MUTED, anchor="middle", mono=True))
    b.append(T(490, 60, "∅", 22, INK, anchor="middle"))
    b.append(T(490, 84, "no shared", 10, MUTED, anchor="middle", mono=True))
    b.append(T(490, 96, "patients", 10, MUTED, anchor="middle", mono=True))
    return svg(W, H, "".join(b), "Two cohorts with no shared patients")


def fig_prediction(kind):
    """Small pass/fail glyph for each prediction."""
    W, H = 300, 96
    b = []
    if kind == "P1":
        b.append(line(20, 70, 280, 70, stroke=INK))
        b.append(line(90, 30, 90, 70, stroke=MUTED, sw=1.4, dash="3 3")); b.append(T(90, 20, "floor", 10, MUTED, anchor="middle", mono=True))
        b.append(line(190, 30, 190, 70, stroke=COL["scgpt"], sw=2.2)); b.append(T(190, 20, "FM state", 10, COL["scgpt"], anchor="middle", mono=True))
        b.append(f'<path d="M96,50 L184,50" stroke="{INK}" stroke-width="1.4" marker-end="url(#arr)"/>')
        b.append(T(150, 86, "concordance →", 10, MUTED, anchor="middle", mono=True))
    elif kind == "P2":
        b.append(line(150, 20, 150, 70, stroke=INK))
        b.append(T(150, 86, "0", 10, MUTED, anchor="middle", mono=True))
        b.append(line(170, 40, 250, 40, stroke=COL["scgpt"], sw=3)); b.append(f'<circle cx="210" cy="40" r="4" fill="{COL["scgpt"]}"/>')
        b.append(T(262, 40, "pass", 10, COL["scgpt"], mono=True))
        b.append(line(110, 58, 200, 58, stroke=AXIS, sw=3)); b.append(f'<circle cx="155" cy="58" r="4" fill="{AXIS}"/>')
        b.append(T(60, 58, "fail", 10, MUTED, anchor="middle", mono=True))
        b.append(T(150, 10, "FM − baseline, patient bootstrap", 10, MUTED, anchor="middle", mono=True))
    elif kind == "P3":
        for i, lab in enumerate(["LumA", "LumB", "Basal"]):
            x = 40 + i * 90
            b.append(box(x, 26, 70, 44, fill="#F2F4F5"))
            b.append(T(x + 35, 18, lab, 10, MUTED, anchor="middle", mono=True))
            # two tiny survival curves
            hi = COL["scgpt"] if i == 2 else AXIS
            b.append(f'<path d="M{x + 8},34 L{x + 30},44 L{x + 50},52 L{x + 62},60" stroke="{INK}" fill="none" stroke-width="1.2"/>')
            gap = 14 if i == 2 else 3
            b.append(f'<path d="M{x + 8},34 L{x + 30},{44 - gap} L{x + 50},{52 - gap} L{x + 62},{60 - gap}" stroke="{hi}" fill="none" stroke-width="1.6"/>')
        b.append(T(150, 86, "survival split within a PAM50 subtype", 10, MUTED, anchor="middle", mono=True))
    elif kind == "P4":
        b.append(box(30, 26, 240, 44, fill="#F2F4F5"))
        b.append(f'<rect x="30" y="26" width="110" height="44" rx="3" fill="{AXIS}" opacity=".5"/>')
        b.append(T(85, 48, "HVG set", 10.5, INK, anchor="middle", mono=True))
        b.append(T(205, 48, "other genes", 10.5, MUTED, anchor="middle", mono=True))
        import random
        rnd = random.Random(4)
        for _ in range(10):
            b.append(f'<circle cx="{150 + rnd.uniform(6, 112):.0f}" cy="{rnd.uniform(30, 66):.0f}" r="3" fill="{COL["geneformer"]}"/>')
        for _ in range(3):
            b.append(f'<circle cx="{rnd.uniform(36, 130):.0f}" cy="{rnd.uniform(30, 66):.0f}" r="3" fill="{COL["geneformer"]}"/>')
        b.append(T(150, 86, "where the signature's genes fall", 10, MUTED, anchor="middle", mono=True))
    return svg(W, H, "".join(b), f"Prediction {kind} glyph")


def fig_outcomes():
    """Three miniature ladders: FM wins / sceptics win / both fail."""
    W, H = 980, 190
    b = []
    scen = [
        ("FM wins", {"floor": .56, "base": .62, "test": .70, "ref": .74}, "test clearly above baseline, both above floor"),
        ("Sceptics win", {"floor": .56, "base": .64, "test": .645, "ref": .74}, "test ≈ baseline, both below the reference"),
        ("Both fail", {"floor": .56, "base": .565, "test": .57, "ref": .74}, "test ≈ baseline ≈ floor — the translation lost the signal"),
    ]
    for i, (title, v, desc) in enumerate(scen):
        x0 = 30 + i * 320
        b.append(box(x0, 16, 290, 150))
        b.append(T(x0 + 14, 36, title, 13, INK, weight=600))
        ax = 120
        b.append(line(x0 + 20, ax, x0 + 270, ax, stroke=INK))

        def X(c):
            return x0 + 20 + (c - .5) / .3 * 250
        for lab, c, col, dash in [("floor", v["floor"], MUTED, "3 3"), ("base", v["base"], COL["hvg_pca"], ""), ("test", v["test"], COL["scgpt"], ""), ("ref", v["ref"], INK, "8 4")]:
            b.append(line(X(c), 56, X(c), ax, stroke=col, sw=2 if not dash else 1.4, dash=dash))
            yl = 48 if lab != "test" else 66
            if lab == "test" and abs(v["test"] - v["base"]) < .02:
                b.append(T(X(c) + 4, 66, "test", 9.5, COL["scgpt"], mono=True))
                continue
            b.append(T(X(c), 48, lab, 9.5, col if col != MUTED else MUTED, anchor="middle", mono=True))
        b.append(T(x0 + 14, 146, desc, 10.5, MUTED))
    return svg(W, H, "".join(b), "Three possible outcomes shown as miniature ladders")


# ----------------------------------------------------------------------------- results
def data_uri(p: Path):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


def result_slot(name, title, what, read):
    p = REPORT / f"{name}.png"
    if p.exists():
        img = f'<img src="{data_uri(p)}" alt="{title}">'
        state = ""
    else:
        img = (f'<div class="pending"><span class="pending-k">{name}.png</span>'
               f'<span>Result panel — fills in when stage 06 runs.</span></div>')
        state = ' data-pending="true"'
    return (f'<figure class="result" data-slot="{name}"{state}>'
            f'<figcaption><span class="rk">{title}</span></figcaption>{img}'
            f'<div class="rtext"><p><strong>What it shows.</strong> {what}</p><p><strong>How to read it.</strong> {read}</p></div>'
            f'</figure>')


def ladder_table():
    rows = []
    have = LADDER.exists()
    L = json.load(open(LADDER)) if have else None
    order = ["hvg_pca", "scgpt", "geneformer"]
    by = {r["model"]: r for r in L["ladder"]} if have else {}
    for m in order:
        r = by.get(m)
        sw = f'<span class="swatch" style="background:{COL[m]}"></span>'
        if r:
            rows.append(f'<tr><td>{sw}{NAME[m]}</td><td>{r["top_cell_type"]}</td>'
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


def predictions_status():
    if not LADDER.exists():
        return ""
    L = json.load(open(LADDER))
    out = ['<div class="pstatus">']
    for m, p in L["predictions"].items():
        cells = "".join(f'<span class="pill {"on" if p[k] else "off"}">{k.split("_")[0]} {"pass" if p[k] else "fail"}</span>'
                        for k in ("P1_above_floor", "P2_pass", "P3_pass", "P4_pass"))
        out.append(f'<div><span class="swatch" style="background:{COL[m]}"></span><strong>{NAME[m]}</strong> {cells}</div>')
    out.append("</div>")
    return "".join(out)


# ----------------------------------------------------------------------------- page
CSS = f"""
:root{{
  --paper:#F7F8F9; --panel:#FFFFFF; --ink:{INK}; --muted:{MUTED}; --rule:{RULE}; --axis:{AXIS};
  --base:{COL['hvg_pca']}; --scgpt:{COL['scgpt']}; --gf:{COL['geneformer']}; --band:#EEF1F2;
  --sans:{SANS}; --mono:{MONO};}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);line-height:1.6;
  -webkit-font-smoothing:antialiased;padding-inline:16px}}
.wrap{{max-width:1040px;margin:0 auto;padding:34px 0 96px}}
.topnav{{position:sticky;top:0;z-index:10;background:color-mix(in srgb,var(--paper) 88%,transparent);
  backdrop-filter:blur(8px);border-bottom:1px solid var(--rule);margin-inline:-16px;padding-inline:16px}}
.topnav .inner{{max-width:1040px;margin:0 auto;display:flex;align-items:center;justify-content:space-between;gap:16px;height:54px}}
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
.lede{{font-size:17px;color:var(--muted);margin:0 0 16px;max-width:64ch}}
.repo{{font-family:var(--mono);font-size:12px;color:var(--muted)}}
.repo a{{color:var(--ink)}}
h2{{font-size:21px;letter-spacing:-.01em;margin:56px 0 8px;font-weight:640;text-wrap:balance}}
h3{{font-size:14.5px;font-weight:640;margin:24px 0 8px;display:flex;align-items:center;gap:8px}}
.sub{{color:var(--muted);font-size:15px;margin:0 0 20px;max-width:70ch}}
p{{margin:0 0 14px;max-width:72ch}}
.panel{{margin:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:14px 12px 8px}}
.panel figcaption{{font-family:var(--mono);font-size:11.5px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);margin:0 0 6px;padding-left:4px}}
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
.why{{color:var(--muted);font-size:14px;border-left:2px solid var(--rule);padding-left:12px;margin:10px 0 0;max-width:68ch}}
.preds{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin:18px 0}}
@media (max-width:720px){{.preds{{grid-template-columns:minmax(0,1fr)}}}}
.pred{{background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:14px 16px 10px;display:grid;
  grid-template-columns:minmax(0,1fr) 150px;gap:6px 14px;align-items:start}}
@media (max-width:520px){{.pred{{grid-template-columns:minmax(0,1fr)}}}}
.pred .id{{font-family:var(--mono);font-size:11px;letter-spacing:.1em;color:var(--muted);grid-column:1/-1}}
.pred h3{{margin:0 0 4px;font-size:14.5px}}
.pred p{{font-size:13.5px;margin:0 0 6px}}
.pred .ifnot{{font-size:12.5px;color:var(--muted)}}
.pred svg{{width:100%;height:auto;display:block}}
.results{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin:18px 0}}
@media (max-width:820px){{.results{{grid-template-columns:minmax(0,1fr)}}}}
.result{{margin:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:12px 12px 10px}}
.result figcaption{{font-family:var(--mono);font-size:11.5px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);margin:0 0 8px;padding-left:2px}}
.result .rtext p{{font-size:12.5px;color:var(--muted);margin:8px 0 0;max-width:none}}
.result .rtext strong{{color:var(--ink);font-weight:600}}
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
.conv p{{font-size:13.5px;color:var(--muted);margin:0 0 8px;max-width:none}}
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
    return f'<figure class="panel fig"><figcaption>{caption}</figcaption>{fig_html}</figure>'


def build():
    today = date.today().isoformat()
    have_results = LADDER.exists()
    status = "results in" if have_results else "design complete · results pending"

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
        f'<div class="pred"><span class="id">{pid}</span><div><h3>{t}</h3><p>{d}</p><p class="ifnot">{n}</p></div>{fig_prediction(pid)}</div>'
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
                    "One bar per representation: the concordance index of its best cell-state signature, with a 95% interval from resampling patients. Dotted tick: that signature's matched-random floor. Dashed line: PAM50.",
                    "Bars that clear their dotted tick pass P1. The distance between the scGPT / Geneformer bars and the HVG-PCA bar is the study's result."),
        result_slot("fig_margin", "B · the margin over the baseline",
                    "For each foundation model, the difference in concordance from the HVG-PCA baseline, computed on the same 200 patient resamples so the comparison is paired.",
                    "An interval wholly to the right of zero passes P2. An interval straddling zero is the sceptics' result, and is reported as such."),
        result_slot("fig_km", "C · survival curves, overall and within subtype",
                    "Kaplan–Meier curves for patients split into thirds by the best foundation-model signature: all patients, then each PAM50 subtype with enough events.",
                    "Separated curves inside a single subtype panel is P3: the state orders patients that PAM50 puts in one box. Curves that separate only in the all-patients panel mean the state is tracking subtype."),
        result_slot("fig_hvg", "D · where the informative genes come from",
                    "Every signature from every representation: the fraction of its genes that are HVGs (x) against its concordance above the matched-random floor (y).",
                    "Foundation-model points high on the y-axis and left of 0.5 on the x-axis are P4: prognostic signal built from genes the baseline never saw."),
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
        ("predictions", "Predictions"), ("results", "Results"), ("conventions", "Conventions"),
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
<div class="steps">{steps_html}</div>

<h2 id="ladder">The ladder: what every number is read against</h2>
<p class="sub">A concordance index, on its own, means nothing. 0.65 might be a strong result or what any random handful of highly expressed genes achieves. So every score sits on a ladder, and the findings are gaps between rungs, not rungs.</p>
{panel(fig_ladder_schematic(), "The ladder — rung positions are illustrative until the data land")}
<table class="rungs"><tbody>{rungs_html}</tbody></table>
<h3>The trap the baseline rung guards against</h3>
{panel(fig_shared_input(), "One input, three views, one scoring set")}
<p>The three representations are not independent witnesses. They were built from the same expression matrix and their signatures are scored in the same bulk genes. Prognostic signal that appears in all three can come from that shared input — from cell-type composition, from tumour purity, from proliferation — rather than from anything a model learned. The baseline rung is the control: it is what the shared input yields when no model learns anything. Every claim about the foundation models is therefore the margin over that rung.</p>

<h2 id="predictions">Four predictions, written down before running</h2>
<p class="sub">Each is a yes-or-no test with the failing outcome spelled out. The study can lose.</p>
<div class="preds">{preds_html}</div>
{predictions_status()}

<h2 id="results">Results</h2>
<p class="sub">{"The four panels stage 06 writes, and the ladder table." if have_results else "Nothing has run. The four panels below are the ones stage 06 will write, each with what it will show and how to read it; the table has the shape of the final ladder. When the data land, this section fills in and nothing else on the page changes."}</p>
<div class="results">{results_html}</div>
<h3>The ladder as a table</h3>
{ladder_table()}

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
<tr><td>Both fail</td><td>Test ≈ baseline ≈ floor. The route from cell state to bulk signature loses the signal; this design cannot separate the two poles, and the report says so rather than looking for a setting where it can.</td></tr>
</tbody></table>

<h2 id="validation">What a positive result would still need</h2>
<p class="sub">A cell state that passes P1–P3 in TCGA is a hypothesis, not a biomarker. Stage 06 writes a validation plan for each shortlisted state with three fixed steps.</p>
<div class="steps">
<div class="step"><h3>Replicate in an independent patient cohort</h3><p>Re-score the signature in METABRIC (~2,000 patients, overall survival) with the same age- and stage-adjusted Cox model. A state that does not replicate is dropped.</p><p class="why">TCGA is one cohort with one sequencing platform and ~15% events. A single-cohort survival result is where most published signatures stop, and where most of them fail.</p></div>
<div class="step"><h3>Go back to the atlas</h3><p>Check that the signature's genes are expressed in the annotated cell type the state was assigned to, not in a contaminating population, and that the state is present across donors rather than in one.</p><p class="why">A signature can score well in bulk for a reason that has nothing to do with the cells it was named after.</p></div>
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
<div>Data: CELLxGENE Census (breast cancer, primary data) · TCGA-BRCA via UCSC Xena (HiSeqV2, clinical matrix, TCGA-CDR survival). Both public.</div>
<div>Pipeline: <a href="{REPO}">single-cell-fm-probing</a> · stages 01–06 · Python, scanpy, lifelines · result panels are written by <code>scripts/06_report.py</code> and embedded by <code>docs/artifact/build.py</code>.</div>
<div>Page built {today}. {"Results embedded." if have_results else "No results embedded."}</div>
</footer>
</div>
"""
    OUT.write_text(html)
    print(f"wrote {OUT}  {OUT.stat().st_size / 1024:.0f} KB  results={'yes' if have_results else 'no'}")


if __name__ == "__main__":
    build()
