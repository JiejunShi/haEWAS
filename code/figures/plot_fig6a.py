"""Plot Fig6a using only its prepared CSV."""
from _fig6_plot_common import *
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator, NullLocator

CELL_ORDER = ("Hep", "Chol", "EC", "Kup", "Lym")
PLOT_METHODS = ("haEWAS-specific", "Common", "EWAS-specific")
METHOD_COLORS = {
    "haEWAS-specific": "#B64342",
    "Common": "#8F8F8F",
    "EWAS-specific": "#0F4D92",
}
METHOD_LABELS = {
    "haEWAS-specific": "haEWAS-\nspecific",
    "Common": "Common",
    "EWAS-specific": "EWAS-\nspecific",
}

MM = 1 / 25.4
PANEL_W_MM, PANEL_H_MM = 82, 44
FS_PANEL, FS_TITLE, FS_AXIS = 7.0, 6.0, 6.0
FS_TICK, FS_LEGEND, FS_LEGEND_TITLE, FS_ANNOT = 5.0, 5.0, 5.0, 5.0
LW_AXIS, LW_GRID, LW_BAR_EDGE, LW_ZERO = 0.55, 0.30, 0.35, 0.50
COLOR_AXIS, COLOR_GRID = "#272727", "#C7C7C7"




def plot_dataset(dataset, subset):
    expected = {(cell, category) for cell in CELL_ORDER for category in PLOT_METHODS}
    actual = set(subset[["cell_type", "category"]].itertuples(index=False, name=None))
    if actual != expected:
        raise ValueError("Incomplete cell/category grid for {}".format(dataset))
    if subset.duplicated(["cell_type", "category"]).any():
        raise ValueError("Duplicate cell/category counts for {}".format(dataset))

    indexed = subset.set_index(["cell_type", "category"])
    x_indices = np.arange(len(CELL_ORDER), dtype=float)
    max_label_len = max(map(len, CELL_ORDER))
    label_rotation = 22 if max_label_len <= 6 else 35
    bottom_margin = 0.185 if label_rotation <= 22 else 0.235

    legend_left = 0.810
    axes_left = 0.105
    axes_top = 0.855
    axes_width = max(0.560, legend_left - axes_left - 0.030)
    axes_height = axes_top - bottom_margin

    fig = plt.figure(figsize=(PANEL_W_MM * MM, PANEL_H_MM * MM),
                     facecolor="white", constrained_layout=False)
    ax = fig.add_axes([axes_left, bottom_margin, axes_width, axes_height])
    ax.minorticks_off()
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.grid(True, linestyle=(0, (2.0, 2.0)), linewidth=LW_GRID,
                  alpha=0.50, color=COLOR_GRID, zorder=0)
    ax.set_axisbelow(True)

    total_width_per_cell = min(0.78, max(0.52, 0.88 - 0.035 * max(len(CELL_ORDER) - 2, 0)))
    bar_width = total_width_per_cell / len(PLOT_METHODS)
    max_pos = pd.to_numeric(subset["positive"], errors="coerce").max()
    max_neg = pd.to_numeric(subset["negative"], errors="coerce").max()
    if not np.isfinite(max_pos) or not np.isfinite(max_neg) or min(max_pos, max_neg) < 0:
        raise ValueError("Invalid directional counts for {}".format(dataset))
    y_lim = max(max_pos, max_neg, 1.0) * 1.38
    ax.set_ylim(-y_lim, y_lim)
    ax.set_xlim(-0.50, len(CELL_ORDER) - 0.50)

    y_span = 2 * y_lim
    text_offset = max(0.35, y_span * 0.012)
    top_inner_limit = y_lim - y_span * 0.060
    bottom_inner_limit = -y_lim + y_span * 0.060

    for j, category in enumerate(PLOT_METHODS):
        color = METHOD_COLORS[category]
        offset = (j - (len(PLOT_METHODS) - 1) / 2) * bar_width
        x_pos = x_indices + offset
        rows = indexed.loc[[(cell, category) for cell in CELL_ORDER]]
        pos_vals = pd.to_numeric(rows["positive"], errors="raise").to_numpy(float)
        neg_vals = -pd.to_numeric(rows["negative"], errors="raise").to_numpy(float)

        ax.bar(x_pos, pos_vals, width=bar_width, facecolor=color,
               edgecolor=COLOR_AXIS, linewidth=LW_BAR_EDGE, zorder=3)
        ax.bar(x_pos, neg_vals, width=bar_width, facecolor="white",
               edgecolor=color, linewidth=LW_BAR_EDGE, hatch="///", zorder=3)

        for x, positive, negative in zip(x_pos, pos_vals, neg_vals):
            if positive > 0:
                y_text = min(positive + text_offset, top_inner_limit)
                va = "bottom" if positive + text_offset <= top_inner_limit else "top"
                ax.text(x, y_text, str(int(positive)), ha="center", va=va,
                        fontsize=FS_ANNOT, color=COLOR_AXIS, rotation=90,
                        zorder=5, clip_on=True)
            if negative < 0:
                y_text = max(negative - text_offset, bottom_inner_limit)
                va = "top" if negative - text_offset >= bottom_inner_limit else "bottom"
                ax.text(x, y_text, str(int(abs(negative))), ha="center", va=va,
                        fontsize=FS_ANNOT, color=COLOR_AXIS, rotation=90,
                        zorder=5, clip_on=True)

    ax.axhline(0, color=COLOR_AXIS, linewidth=LW_ZERO, zorder=4)
    ax.set_xticks(x_indices)
    ax.set_xticklabels(CELL_ORDER, rotation=label_rotation, ha="right",
                       fontsize=FS_TICK, color=COLOR_AXIS)
    ax.set_ylabel("DMC count", fontsize=FS_AXIS, labelpad=0.8, color=COLOR_AXIS)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.yaxis.set_major_formatter(FuncFormatter(
        lambda y, _: "0" if abs(y) < 1e-9 else str(abs(int(y)))))
    ax.tick_params(axis="both", which="major", labelsize=FS_TICK,
                   width=LW_AXIS, length=2.0, direction="out", pad=1.0,
                   colors=COLOR_AXIS)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(LW_AXIS)
        spine.set_color(COLOR_AXIS)

    fig.text(0.018, 0.982, "a", ha="left", va="top", fontsize=FS_PANEL,
             fontweight="bold", color=COLOR_AXIS)
    fig.text(0.455, 0.975, "{} DMC counts".format(dataset), ha="center", va="top",
             fontsize=FS_TITLE, color=COLOR_AXIS)
    method_handles = [Patch(facecolor=METHOD_COLORS[m], edgecolor=COLOR_AXIS,
        linewidth=LW_BAR_EDGE, label=METHOD_LABELS[m]) for m in PLOT_METHODS]
    leg1 = fig.legend(handles=method_handles, title="Loci category", loc="upper left",
        bbox_to_anchor=(legend_left, 0.790), frameon=False, fontsize=FS_LEGEND,
        title_fontsize=FS_LEGEND_TITLE, handlelength=0.82, handleheight=0.54,
        handletextpad=0.26, labelspacing=0.15, borderaxespad=0.0)
    leg1._legend_box.align = "left"
    direction_handles = [
        Patch(facecolor="#BDBDBD", edgecolor=COLOR_AXIS, linewidth=LW_BAR_EDGE,
              label="Hyper (+)"),
        Patch(facecolor="white", edgecolor=COLOR_AXIS, linewidth=LW_BAR_EDGE,
              hatch="///", label="Hypo (-)"),
    ]
    leg2 = fig.legend(handles=direction_handles, title="Direction", loc="upper left",
        bbox_to_anchor=(legend_left, 0.335), frameon=False, fontsize=FS_LEGEND,
        title_fontsize=FS_LEGEND_TITLE, handlelength=0.82, handleheight=0.54,
        handletextpad=0.26, labelspacing=0.16, borderaxespad=0.0)
    leg2._legend_box.align = "left"
    return fig

def run():
    frame = pd.read_csv(ROOT / 'data/Fig6a.csv').rename(columns={
        'CellType': 'cell_type', 'Category': 'category', 'Positive': 'positive', 'Negative': 'negative'})
    style()
    finish(plot_dataset('GSE157341', frame), 'a', tight=True)

if __name__ == '__main__':
    run()
