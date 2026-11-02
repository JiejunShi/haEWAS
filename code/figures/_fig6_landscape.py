"""Shared Fig6c/d plotting and label/rug selection."""
from _fig6_plot_common import *
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

X_MIN, X_MAX = -12, 12
X_TICKS = (-12, -6, 0, 6, 12)
FDR_CUTOFF = 0.05
TRACK_TOP_K, LABEL_TOP_N = 50, 5
LABEL_CATEGORY_PRIORITY = ("haEWAS-specific", "Common", "EWAS-specific")

COLOR_TRAD = "#0F4D92"
COLOR_HA = "#B64342"
COLOR_COMMON = "#8F8F8F"
COLOR_AXIS = "#272727"
COLOR_FDR = "#4F7F3A"
COLOR_REF = "#7A7A7A"
CATEGORY_COLORS = {
    "EWAS-specific": COLOR_TRAD,
    "haEWAS-specific": COLOR_HA,
    "Common": COLOR_COMMON,
}
PANEL_W_MM, PANEL_H_MM = 82, 68
MM = 1 / 25.4
FS_PANEL, FS_TITLE, FS_AXIS = 7.0, 6.2, 6.0
FS_TICK, FS_TRACK, FS_GENE = 5.0, 5.0, 5.0
LW_AXIS, LW_REF, LW_FDR, LW_TRACK = 0.55, 0.50, 0.55, 0.35
LABEL_MIN_DY_FRAC, LABEL_SAFE_X_FRAC = 0.055, 0.030
LABEL_SAFE_Y_FRAC, LABEL_TEXT_GAP_FRAC = 0.065, 0.010
LABEL_FONTWEIGHT, LABEL_OUTLINE_LW = "normal", 0.90
USE_LABEL_OUTLINE = True




def _plot_rug(ax, frame):
    rows = (("haEWAS-specific", COLOR_HA, 2), ("Common", COLOR_COMMON, 1),
            ("EWAS-specific", COLOR_TRAD, 0))
    for category, color, y in rows:
        subset = frame.loc[frame["Category"].eq(category) &
                             frame["is_track_topk"] &
                             frame["Tstat"].between(X_MIN, X_MAX)]
        ax.add_patch(Rectangle((X_MIN, y - .38), X_MAX - X_MIN, .76,
            facecolor="#F3F3F3", edgecolor="none", alpha=.82, zorder=0))
        if len(subset):
            ax.eventplot(subset["Tstat"].to_numpy(float), lineoffsets=y,
                linelengths=.56, linewidths=LW_TRACK, colors=color,
                alpha=.88, zorder=5)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(-.5, 2.5)
    ax.set_yticks((0, 1, 2))
    ax.set_yticklabels(("EWAS-specific", "Common", "haEWAS-specific"), fontsize=FS_TRACK)
    for label, color in zip(ax.get_yticklabels(), (COLOR_TRAD, COLOR_COMMON, COLOR_HA)):
        label.set_color(color)
    ax.tick_params(axis="y", length=0, pad=1.2)
    for side in ("top", "left", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(LW_AXIS)
    ax.spines["bottom"].set_color(COLOR_AXIS)
    ax.spines["bottom"].set_bounds(X_MIN, X_MAX)
    ax.axvline(0, color=COLOR_REF, linestyle=(0, (2, 2)), linewidth=LW_REF, alpha=.75)
    ax.set_xlabel("t-statistic", fontsize=FS_AXIS, labelpad=1.0, color=COLOR_AXIS)
    ax.set_xticks(X_TICKS)
    ax.tick_params(axis="x", labelsize=FS_TICK, width=LW_AXIS,
                   length=2.0, pad=.8, colors=COLOR_AXIS)

def _spread_label_y(values, min_gap, y_min, y_max):
    values = np.asarray(values, dtype=float)
    if len(values) <= 1:
        return np.clip(values, y_min, y_max)
    order = np.argsort(values)
    y = values[order].copy()
    for i in range(1, len(y)):
        y[i] = max(y[i], y[i - 1] + min_gap)
    available = y_max - y_min
    span = y[-1] - y[0]
    if span > available and span > 0:
        center = (y[-1] + y[0]) / 2
        y = (y - center) * (available / span) + center
    if y[0] < y_min:
        y += y_min - y[0]
    if y[-1] > y_max:
        y -= y[-1] - y_max
    result = np.empty_like(y)
    result[order] = y
    return result

def _italic_gene_font():
    import matplotlib.font_manager as fm
    candidates = (
        ARIAL_TTF.with_name("ariali.ttf"),
        ARIAL_TTF.with_name("Arial Italic.ttf"),
        ARIAL_TTF.with_name("ArialItalic.ttf"),
        ARIAL_TTF.with_name("Arial_Italic.ttf"),
    )
    for path in candidates:
        if path.is_file():
            fm.fontManager.addfont(str(path))
            return fm.FontProperties(fname=str(path), weight=LABEL_FONTWEIGHT)
    return fm.FontProperties(
        family="DejaVu Sans", style="italic", weight=LABEL_FONTWEIGHT)

def _add_global_repel_labels(ax, labels):
    """Original one-dimensional label repulsion and straight connectors."""
    if labels.empty:
        return
    labels = labels.loc[
        labels["gene"].notna() & labels["gene"].astype(str).str.strip().ne("") &
        labels["Tstat"].between(X_MIN, X_MAX)
    ].copy()
    if labels.empty:
        return
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x_range, y_range = x1 - x0, y1 - y0
    safe_x_left = x0 + x_range * LABEL_SAFE_X_FRAC
    safe_x_right = x1 - x_range * LABEL_SAFE_X_FRAC
    safe_y_min = y0 + y_range * LABEL_SAFE_Y_FRAC
    safe_y_max = y1 - y_range * LABEL_SAFE_Y_FRAC
    min_dy = y_range * LABEL_MIN_DY_FRAC
    text_gap = x_range * LABEL_TEXT_GAP_FRAC
    gene_font = _italic_gene_font()
    for direction, positive_t in (("+", True), ("-", False)):
        subset = labels.loc[labels["Direction"].eq(direction)].sort_values(
            ["mlog10FDR", "CpG", "Method"], ascending=[False, True, True], kind="stable")
        subset = subset.drop_duplicates("gene", keep="first")
        if subset.empty:
            continue
        y_positions = _spread_label_y(
            subset["mlog10FDR"].to_numpy(float), min_dy, safe_y_min, safe_y_max)
        for (_, row), y_text in zip(subset.iterrows(), y_positions):
            xi = float(row["Tstat"])
            yi = float(row["mlog10FDR"])
            if positive_t:
                x_text, ha = safe_x_right, "right"
                x_label_edge = x_text - text_gap
            else:
                x_text, ha = safe_x_left, "left"
                x_label_edge = x_text + text_gap
            ax.plot([xi, x_label_edge],
                    [yi, y_text],
                    color=COLOR_AXIS, linewidth=0.34, alpha=0.62,
                    zorder=3, solid_capstyle="round", clip_on=True)
            text = ax.text(x_text, y_text, str(row["gene"]),
                ha=ha, va="center", fontsize=FS_GENE,
                fontproperties=gene_font, color=row["label_color"],
                zorder=21, clip_on=True)
            # Explicit patch clipping prevents long symbols from painting into
            # the figure margin even if a backend ignores the axes clip flag.
            text.set_clip_path(ax.patch)
            if USE_LABEL_OUTLINE:
                text.set_path_effects([
                    pe.withStroke(linewidth=LABEL_OUTLINE_LW, foreground="white"),
                    pe.Normal(),
                ])

def _draw(cell, panel, frame):
    finite_y = frame["mlog10FDR"].replace([np.inf, -np.inf], np.nan).dropna()
    y_upper = max(5, int(np.ceil(float(finite_y.max()) * 1.06))) if len(finite_y) else 5
    fig = plt.figure(figsize=(PANEL_W_MM * MM, PANEL_H_MM * MM), facecolor="white")
    ax_vol = fig.add_axes([.265, .305, .500, .620])
    ax_rug = fig.add_axes([.265, .130, .500, .145], sharex=ax_vol)

    for category, color, size, zorder in (
        ("EWAS-specific", COLOR_TRAD, 8.2, 6),
        ("haEWAS-specific", COLOR_HA, 8.2, 7),
        ("Common", COLOR_COMMON, 8.2, 8),
    ):
        subset = frame.loc[frame["Category"].eq(category) &
                           frame["Tstat"].between(X_MIN, X_MAX)]
        if len(subset):
            ax_vol.scatter(subset["Tstat"], subset["mlog10FDR"], s=size, c=color,
                alpha=.94, edgecolors="white", linewidths=.16, zorder=zorder)

    ax_vol.axvline(0, color=COLOR_REF, linewidth=LW_REF,
                   linestyle=(0, (2, 2)), alpha=.72)
    y_threshold = -np.log10(FDR_CUTOFF)
    ax_vol.axhline(y_threshold, color=COLOR_FDR, linewidth=LW_FDR,
                   linestyle=(0, (2, 2)), alpha=.85)
    ax_vol.set_xlim(X_MIN, X_MAX)
    ax_vol.set_ylim(0, y_upper)
    ax_vol.yaxis.set_major_locator(MultipleLocator(5))
    ax_vol.set_ylabel(r"$-\log_{10}(\mathrm{FDR})$", fontsize=FS_AXIS,
                      labelpad=1.0, color=COLOR_AXIS)
    ax_vol.tick_params(axis="y", labelsize=FS_TICK, width=LW_AXIS,
                       length=2.0, pad=.8, colors=COLOR_AXIS)
    ax_vol.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ax_vol.spines.values():
        spine.set_linewidth(LW_AXIS)
        spine.set_color(COLOR_AXIS)
        spine.set_visible(True)

    labels = frame.loc[frame["is_label"] & frame["gene"].ne("")].copy()
    label_colors = {"EWAS-specific": COLOR_TRAD, "haEWAS-specific": COLOR_HA, "Common": COLOR_COMMON}
    labels["label_color"] = labels["Category"].map(label_colors).fillna(COLOR_AXIS)
    _add_global_repel_labels(ax_vol, labels)

    _plot_rug(ax_rug, frame)
    fig.text(.515, .975, "{} t-statistic landscape".format(cell),
             ha="center", va="top", fontsize=FS_TITLE, color=COLOR_AXIS)
    fig.text(.018, .985, panel, ha="left", va="top", fontsize=FS_PANEL,
             fontweight="bold", color=COLOR_AXIS)
    return fig, frame

def _style_original():
    import matplotlib as mpl
    import matplotlib.font_manager as fm
    if ARIAL_TTF.is_file():
        fm.fontManager.addfont(str(ARIAL_TTF))
        family = fm.FontProperties(fname=str(ARIAL_TTF)).get_name()
    else:
        family = "DejaVu Sans"
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [family, "Helvetica", "DejaVu Sans", "Liberation Sans"],
        "font.size": 5.5,
        "axes.labelsize": 6.0,
        "axes.titlesize": 6.0,
        "xtick.labelsize": 5.0,
        "ytick.labelsize": 5.0,
        "axes.linewidth": 0.55,
        "lines.linewidth": 0.70,
        "patch.linewidth": 0.45,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.unicode_minus": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
    })

def _mark_annotations(frame):
    frame = frame.copy()
    frame["is_track_topk"] = False
    frame["is_label"] = False
    fdr_for_rank = np.power(10.0, -frame["mlog10FDR"])
    for category in ("EWAS-specific", "haEWAS-specific", "Common"):
        for direction in ("+", "-"):
            subset = frame.loc[
                frame["Category"].eq(category) & frame["Direction"].eq(direction) &
                frame["Tstat"].between(X_MIN, X_MAX)
            ].sort_values("mlog10FDR", ascending=False, kind="stable")
            frame.loc[subset.head(TRACK_TOP_K).index, "is_track_topk"] = True

    for direction in ("+", "-"):
        used = set()
        for category in LABEL_CATEGORY_PRIORITY:
            subset = frame.loc[
                frame["Category"].eq(category) & frame["Direction"].eq(direction) &
                frame["gene"].notna() & frame["gene"].astype(str).str.strip().ne("") &
                np.isfinite(frame["mlog10FDR"]) & fdr_for_rank.between(0, 1) &
                np.isfinite(frame["Tstat"]) & frame["Tstat"].between(X_MIN, X_MAX) &
                ~frame["gene"].isin(used)
            ].sort_values(
                ["mlog10FDR", "CpG", "Method"], ascending=[False, True, True],
                kind="stable")
            chosen = subset.drop_duplicates("gene", keep="first").head(LABEL_TOP_N)
            frame.loc[chosen.index, "is_label"] = True
            used.update(chosen["gene"])

    red = frame.loc[
        frame["Category"].eq("haEWAS-specific") & frame["gene"].ne("") &
        frame["Tstat"].between(X_MIN, X_MAX)
    ].sort_values("mlog10FDR", ascending=False, kind="stable")
    if len(red):
        frame.loc[red.index[0], "is_track_topk"] = True
    return frame


def run(cell, panel):
    frame = pd.read_csv(ROOT / "data" / ("Fig6" + panel + ".csv")).rename(columns={
        "Mlog10FDR": "mlog10FDR", "Gene": "gene"})
    frame["gene"] = frame["gene"].fillna("")
    required = {"gene", "CpG", "Method", "Tstat", "mlog10FDR", "Category", "Direction"}
    if not required.issubset(frame.columns):
        raise ValueError("Prepared landscape table lacks required columns")
    if (frame.empty or not np.isfinite(frame[["Tstat", "mlog10FDR"]]).all().all()
            or not frame["Category"].isin(
                ("EWAS-specific", "haEWAS-specific", "Common")).all()
            or not frame["Direction"].isin(("+", "-", "0")).all()):
        raise ValueError("Invalid Fig6 landscape data")
    frame = _mark_annotations(frame)
    _style_original()
    fig, _ = _draw(cell, panel, frame)
    finish(fig, panel)