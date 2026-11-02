"""Plot Fig6b using only its prepared CSV."""
from _fig6_plot_common import *
from matplotlib.patches import Patch

CELL_ORDER = ("Hep", "Chol", "EC", "Kup", "Lym")
MAX_INTERSECTIONS = None
INTERSECTION_Y_MAX = 1000
INTERSECTION_Y_TICKS = (0, 200, 400, 600, 800, 1000)
INTERSECTION_HEADROOM_FRAC = 0.14
TOTAL_X_PAD_FRAC = 0.40

MM = 1 / 25.4
PANEL_W_MM, PANEL_H_MM = 82, 44
FS_PANEL, FS_TITLE, FS_AXIS = 7.0, 6.0, 6.0
FS_TICK, FS_COUNT, FS_CELL = 5.0, 5.0, 5.2
LW_AXIS, LW_BAR_EDGE, LW_MATRIX_LINE, LW_DOT_EDGE = 0.55, 0.35, 0.55, 0.35
COLOR_AXIS, COLOR_INACTIVE, COLOR_BAR = "#272727", "#D8D8D8", "#3A3A3A"
CELL_COLOR_MAP = {
    "Hep": "#5A9B4B", "Chol": "#7A68A6", "EC": "#2F9BA3",
    "Kup": "#9A6A42", "Lym": "#D28A32",
}
COLOR_FALLBACK = "#7F7F7F"




def plot_dataset(dataset, ranked, totals):
    if not ranked:
        raise ValueError("Empty Fig6b intersection table")
    n_cells, n_cols = len(CELL_ORDER), len(ranked)
    x = np.arange(n_cols)
    fig = plt.figure(figsize=(PANEL_W_MM * MM, PANEL_H_MM * MM),
                     facecolor="white", constrained_layout=False)
    ax_int = fig.add_axes([0.135, 0.525, 0.585, 0.350])
    ax_mat = fig.add_axes([0.135, 0.145, 0.585, 0.330])
    ax_tot = fig.add_axes([0.735, 0.145, 0.200, 0.330])

    x = np.arange(n_cols)
    raw_counts = [count for _, count in ranked]
    display_counts = [min(int(value), INTERSECTION_Y_MAX) for value in raw_counts]
    bar_width = min(0.62, max(0.34, 5.6 / max(n_cols, 1)))
    count_rotation = 90 if n_cols > 9 else 0
    top_ylim = INTERSECTION_Y_MAX * (1 + INTERSECTION_HEADROOM_FRAC)
    ax_int.bar(x, display_counts, width=bar_width, color=COLOR_BAR,
               edgecolor=COLOR_AXIS, linewidth=LW_BAR_EDGE, zorder=3)
    ax_int.set_ylim(0, top_ylim)
    ax_int.set_yticks(INTERSECTION_Y_TICKS)
    ax_int.spines["left"].set_bounds(0, INTERSECTION_Y_MAX)
    label_offset = INTERSECTION_Y_MAX * 0.028
    top_safe_limit = INTERSECTION_Y_MAX * 1.095
    for xi, raw_value, display_value in zip(x, raw_counts, display_counts):
        y_text = min(display_value + label_offset, top_safe_limit)
        ax_int.text(xi, y_text, str(int(raw_value)), ha="center", va="bottom",
                    fontsize=FS_COUNT, rotation=count_rotation, color=COLOR_AXIS,
                    clip_on=False, zorder=5)
    ax_int.set_xlim(-0.6, n_cols - 0.4)
    ax_int.set_ylabel("Intersection\nsize", fontsize=FS_AXIS, labelpad=0.8,
                      color=COLOR_AXIS)
    ax_int.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_int.tick_params(axis="y", labelsize=FS_TICK, width=LW_AXIS, length=2.0,
                       pad=0.8, colors=COLOR_AXIS)
    for spine in ("top", "right", "bottom"):
        ax_int.spines[spine].set_visible(False)
    ax_int.spines["left"].set_linewidth(LW_AXIS)
    ax_int.spines["left"].set_color(COLOR_AXIS)
    ax_int.grid(False)

    ax_mat.set_xlim(-0.6, n_cols - 0.4)
    ax_mat.set_ylim(n_cells - 0.5, -0.5)
    dot_size = 17 if n_cols <= 14 else 14
    for col_idx, (pattern, _) in enumerate(ranked):
        active_y = []
        pattern_set = set(pattern)
        for y_idx, cell in enumerate(CELL_ORDER):
            active = cell in pattern_set
            color = CELL_COLOR_MAP.get(cell, COLOR_FALLBACK) if active else COLOR_INACTIVE
            if active:
                active_y.append(y_idx)
            ax_mat.scatter(col_idx, y_idx, s=dot_size, color=color,
                edgecolor=COLOR_AXIS if active else "none",
                linewidth=LW_DOT_EDGE if active else 0,
                zorder=3 if active else 2)
        if len(active_y) >= 2:
            ax_mat.plot([col_idx, col_idx], [min(active_y), max(active_y)],
                        color=COLOR_AXIS, linewidth=LW_MATRIX_LINE, zorder=1)
    ax_mat.set_yticks(range(n_cells))
    ax_mat.set_yticklabels(CELL_ORDER, fontsize=FS_CELL, fontweight="normal")
    for label, cell in zip(ax_mat.get_yticklabels(), CELL_ORDER):
        label.set_color(CELL_COLOR_MAP.get(cell, COLOR_FALLBACK))
    ax_mat.tick_params(axis="y", length=0, pad=1.4)
    ax_mat.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)
    ax_mat.grid(False)

    total_values = [totals[cell] for cell in CELL_ORDER]
    max_total = max(total_values) if total_values else 1
    x_right = max(max_total * (1 + TOTAL_X_PAD_FRAC), max_total + 5, 1)
    ax_tot.set_ylim(n_cells - 0.5, -0.5)
    ax_tot.set_xlim(0, x_right)
    for y_idx, cell in enumerate(CELL_ORDER):
        value = totals[cell]
        ax_tot.barh(y_idx, value, height=0.56,
            color=CELL_COLOR_MAP.get(cell, COLOR_FALLBACK), edgecolor=COLOR_AXIS,
            linewidth=LW_BAR_EDGE, zorder=3)
        label_x = value + x_right * 0.035
        if label_x > x_right * 0.92:
            label_x, ha = x_right * 0.92, "right"
        else:
            ha = "left"
        ax_tot.text(label_x, y_idx, str(int(value)), ha=ha, va="center",
                    fontsize=FS_COUNT, color=COLOR_AXIS, clip_on=False, zorder=5)
    ax_tot.set_yticks([])
    ax_tot.set_xticks([])
    ax_tot.set_xlabel("")
    ax_tot.tick_params(axis="both", which="both", left=False, right=False,
        bottom=False, top=False, labelleft=False, labelbottom=False)
    for spine in ax_tot.spines.values():
        spine.set_visible(False)
    ax_tot.grid(False)

    fig.text(0.018, 0.982, "b", ha="left", va="top", fontsize=FS_PANEL,
             fontweight="bold", color=COLOR_AXIS)
    fig.text(0.53, 0.965, "{} haEWAS-specific DMC events".format(dataset),
             ha="center", va="top", fontsize=FS_TITLE, color=COLOR_AXIS)
    return fig

def run():
    frame = pd.read_csv(ROOT / 'data/Fig6b.csv')
    ranked = [(tuple(row.Pattern.split(';')), int(row.Count)) for row in frame.itertuples(index=False)]
    if any(not set(pattern).issubset(CELL_ORDER) or count <= 0 for pattern, count in ranked):
        raise ValueError('Invalid Fig6b intersection table')
    totals = {cell: sum(count for pattern, count in ranked if cell in pattern) for cell in CELL_ORDER}
    style()
    finish(plot_dataset('GSE157341', ranked, totals), 'b', tight=True)

if __name__ == '__main__':
    run()
