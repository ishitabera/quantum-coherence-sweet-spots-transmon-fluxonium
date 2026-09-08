"""Consistent figure styling and the standard plots for this study."""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap

# --- Publication style -----------------------------------------------------
# Applied at import so every notebook gets identical typography and spine
# treatment without repeating rcParams cells. fonttype 42 embeds TrueType
# rather than Type-3 fonts, which many journals (IEEE, APS) reject outright.
mpl.rcParams.update({
    "savefig.dpi": 300,
    "figure.dpi": 110,
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 10.5,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.autolayout": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# --- Categorical channel colours -------------------------------------------
# CVD audit (Vienot/Brettel & Mollon 1999 LMS dichromat projection, cross-checked
# against Machado et al. 2009 severity-1.0; pairwise distances as CIEDE2000):
#
#                        normal   deuteranopia   protanopia
#   blue   vs orange      48.5        59.7          52.7
#   blue   vs green       45.2        48.6          49.7
#   orange vs green       56.1        11.0  <--     14.4  <--
#
# The orange/green pair COLLAPSES under both dichromacies: #d95f0e -> ~#919100
# and #31a354 -> ~#8d8d57 for a deuteranope, i.e. two olives differing by
# dL* = 0.8. dE00 ~ 11 is below the ~15-20 needed for confident categorical
# discrimination in adjacent map patches, so "charge" and "crit-current"
# regions are effectively indistinguishable for ~6% of male readers.
#
# FIX APPLIED (2026-09, alongside the notebook-03 regeneration pass): the green
# is lightened to ColorBrewer Set2 green, which raises the worst-case pair to
# dE00 = 26.0 (deuteranopia) / 28.5 (protanopia) while still reading as
# blue/orange/green to normal vision (worst pair 35.4). #7fcdbb is even safer
# (30.9 / 33.0) but visibly paler.
#
# NOTE a categorical map still has no cue other than hue, so `dominant_map`
# should additionally encode the channel redundantly (hatching or direct
# labels) if it is ever used for more than three categories.
CHANNEL_COLORS = ListedColormap(["#2c7fb8", "#d95f0e", "#66c2a5"])

DPI = 300


def save(fig_or_plt, path):
    """Save a figure at publication DPI as both raster and vector.

    Writes the requested file (usually a PNG) at ``DPI`` and, unless the
    requested file is already a PDF, a ``.pdf`` sibling alongside it, so every
    figure has a scalable version for the LaTeX build.

    Returns the list of paths written.
    """
    written = []
    fig_or_plt.savefig(path, dpi=DPI, bbox_inches="tight")
    written.append(path)

    base, ext = os.path.splitext(path)
    pdf_path = base + ".pdf"
    if ext.lower() != ".pdf":
        fig_or_plt.savefig(pdf_path, bbox_inches="tight")
        written.append(pdf_path)
    return written


def heatmap(df, value="T2_us", index="EC", columns="EJ", log=True,
            title="", xlabel="", ylabel="", cbar_label="T2 (us)", path=None):
    """2D map of a swept quantity. log=True is usually required because
    coherence spans orders of magnitude and a linear scale hides structure."""
    pivot = df.pivot(index=index, columns=columns, values=value)
    plt.figure(figsize=(8, 5))
    norm = LogNorm() if log else None
    plt.pcolormesh(pivot.columns, pivot.index, pivot.values,
                   shading="auto", norm=norm)
    plt.colorbar(label=cbar_label + (" , log" if log else ""))
    plt.xlabel(xlabel or columns); plt.ylabel(ylabel or index)
    plt.title(title)
    written = save(plt, path) if path else []
    plt.show()
    return written


def dominant_map(df, code_col="dom_code", labels=("relaxation", "charge", "crit-current"),
                 index="EC", columns="EJ", title="", path=None):
    """Categorical map showing which noise channel dominates where."""
    grid = df.pivot(index=index, columns=columns, values=code_col)
    plt.figure(figsize=(8, 5))
    plt.pcolormesh(grid.columns, grid.index, grid.values,
                   cmap=CHANNEL_COLORS, vmin=0, vmax=2, shading="auto")
    cbar = plt.colorbar(ticks=[0.33, 1.0, 1.67])
    cbar.ax.set_yticklabels(list(labels))
    plt.xlabel(columns); plt.ylabel(index); plt.title(title)
    written = save(plt, path) if path else []
    plt.show()
    return written
