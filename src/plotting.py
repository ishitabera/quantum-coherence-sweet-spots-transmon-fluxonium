"""Consistent figure styling and the standard plots for this study."""
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap

CHANNEL_COLORS = ListedColormap(["#2c7fb8", "#d95f0e", "#31a354"])
DPI = 150


def save(fig_or_plt, path):
    """Save a figure at publication DPI with tight bounding box."""
    fig_or_plt.savefig(path, dpi=DPI, bbox_inches="tight")


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
    if path:
        plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.show()


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
    if path:
        plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.show()
