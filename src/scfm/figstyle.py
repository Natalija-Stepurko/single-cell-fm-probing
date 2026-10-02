"""Shared figure style: the page palette, display names and matplotlib settings.

Figures are drawn at 114 display pixels per inch: a wide figure is 10 in (the 1,140 px page column), a
half figure 5 in. They are saved at 240 dpi (2,400 px wide), enough for a 2x screen at the column width.
Fonts are the DejaVu Sans that ships with matplotlib, so the PNGs do not depend on system fonts.
"""
from pathlib import Path

INK = "#16191D"
MUTED = "#5B646E"
RULE = "#DDE1E4"
CELLS = "#E6E9EB"
SURFACE = "#FFFFFF"

COL = {"hvg_pca": "#4A4F55", "scgpt": "#1F6F6B", "geneformer": "#C2681A", "reference": "#8A6FA8",
       "clinical": "#9AA1A8"}
TINT = {"hvg_pca": "#848B92", "scgpt": "#7FB5B0", "geneformer": "#E5AE7C"}
NAME = {"hvg_pca": "HVG-PCA", "scgpt": "scGPT", "geneformer": "Geneformer",
        "reference": "Proliferation score"}
MODELS = ("hvg_pca", "scgpt", "geneformer")

WIDE = 10.0
HALF = 5.0
DPI = 240
FS = 9.0          # body text, 14 px on the page
FS_SMALL = 8.0    # annotations, 12.7 px on the page


def mix(hex_color: str, white: float) -> str:
    """`hex_color` blended toward white by the fraction `white`."""
    rgb = [int(hex_color[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(c + (255 - c) * white):02X}" for c in rgb)


def setup():
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
        "xtick.labelsize": FS_SMALL, "ytick.labelsize": FS, "legend.fontsize": FS_SMALL,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED, "axes.linewidth": 0.7,
        "xtick.color": MUTED, "ytick.color": INK, "xtick.major.width": 0.7, "ytick.major.width": 0,
        "xtick.major.size": 3, "ytick.major.size": 0, "axes.titlelocation": "left",
        "axes.titleweight": "bold", "axes.titlepad": 8, "axes.grid": False,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "legend.frameon": False, "svg.hashsalt": "scfm", "path.simplify": True,
    })


def clean(ax, left=False):
    for s in ("top", "right") + (() if left else ("left",)):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)


def xgrid(ax):
    ax.grid(axis="x", color=RULE, lw=0.6)
    ax.set_axisbelow(True)


def save(fig, out: Path, name: str) -> Path:
    path = Path(out) / f"{name}.png"
    fig.savefig(path, dpi=DPI, metadata={"Software": None})
    import matplotlib.pyplot as plt
    plt.close(fig)
    return path
