---
name: paper-figure
description: Produce or revise a publication figure for this study so it matches the existing house style (log color scales, categorical channel maps, dpi=150 PNGs into figures/). Use when asked to plot, re-plot, restyle, or add a panel to any figure destined for figures/, docs/, or the paper.
---

# Making a figure for this study

## Use the shared helpers

`src/plotting.py` already encodes the house style — `heatmap()`, `dominant_map()`, `save()`,
`CHANNEL_COLORS`, `DPI = 150`. Reach for these before writing raw matplotlib. If a new figure
type is genuinely needed, add it to `plotting.py` rather than leaving one-off styling inside a
notebook, so the next figure inherits it.

## Style rules that are not negotiable

**Log color scale for coherence.** T1/T2/Tφ span orders of magnitude across these parameter
spaces; a linear scale flattens the structure into a single bright corner and hides the sweet
spot entirely. `heatmap(..., log=True)` is the default for a reason. Linear is acceptable only
for a deliberately narrow range, and then say why in the caption.

**Categorical maps use the fixed channel colors.** `CHANNEL_COLORS` is
`["#2c7fb8" relaxation, "#d95f0e" charge/flux, "#31a354" crit-current]`, and the same channel
keeps the same color across every figure so panels can be read against each other. If a fourth
channel exists, extend the colormap deliberately — do not let matplotlib recycle.

**Axis labels carry units.** `EJ (GHz)`, `EC (GHz)`, `T2 (µs)`, `Φ/Φ₀`. A bare `EJ` is
incomplete. Colorbars state the quantity and that the scale is log.

**Save at `dpi=150` with `bbox_inches="tight"`, as PNG, into `figures/`.** Reuse the existing
filename when revising a figure so the paper's includes keep resolving; pick a new descriptive
name only for a genuinely new figure.

## Plot from the CSV, not from a live sweep

Read `data/*.csv`. The committed CSVs are the results of record, and re-running a sweep just to
draw a figure risks producing a plot that disagrees with the data everything else cites. If the
CSV you need does not exist, that is a sweep task first (see the `new-sweep` skill).

## Multi-panel figures

`summary_figure.png` is a 4-panel layout. When adding panels, keep a shared color meaning
across panels, label panels `(a) (b) (c) (d)`, and let each panel keep its own colorbar when
the quantities differ — a shared colorbar across different quantities is misleading.

## Before finishing

Display the figure and actually look at it. Check that the sweet spot is visible rather than
saturated, that no label is clipped, that the colorbar range spans the real data, and that
nothing renders as a solid block of one color — that last one usually means a linear scale
where a log scale was needed, or an `inf` poisoning the range.

State which CSV the figure came from and what filter was applied, so the caption can be written
truthfully.
