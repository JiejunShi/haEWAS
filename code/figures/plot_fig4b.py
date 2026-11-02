"""Draw Fig4b using only its prepared CSV; export PDF and display inline."""
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(os.environ["REPRO_ROOT"])
PDF = ROOT / "figures" / "Fig4b.pdf"

def style():
    import matplotlib as mpl
    import matplotlib.font_manager as fm
    font_file = ROOT.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
    if font_file.is_file():
        fm.fontManager.addfont(str(font_file))
        family = fm.FontProperties(fname=str(font_file)).get_name()
    else:
        family = "Arial"
    mpl.rcParams.update({
        "font.family": family, "font.size": 6, "axes.labelsize": 6,
        "axes.titlesize": 7, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
        "legend.fontsize": 5.5, "axes.linewidth": 0.5, "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })

def finish(fig, bbox_inches=None):
    import matplotlib.pyplot as plt
    PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF, format="pdf", facecolor="white", bbox_inches=bbox_inches)
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None:
            display(fig)
    except ImportError:
        pass
    plt.close(fig)

from matplotlib.patches import Patch

GROUPS = ("ewas_specific", "common", "haewas_specific")

TYPE_LABELS = {
    "ewas_specific": "EWAS-specific",
    "common": "Common",
    "haewas_specific": "haEWAS-specific",
}

TYPE_COLORS = {
    "ewas_specific": "#0F4D92",
    "common": "#8F8F8F",
    "haewas_specific": "#B64342",
}

LEGEND_ORDER = ("haewas_specific", "common", "ewas_specific")

Y_ORDER = (
    "OpenSea", "Shelf", "Shore", "Island", "3'UTR", "Gene body",
    "1stExon", "5'UTR", "TSS1500", "TSS200",
)

COLOR_AXIS = "#272727"

COLOR_GRID = "#BDBDBD"

def get_star(p_value):
    if p_value < 1e-10:
        return "***"
    if p_value < 1e-4:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""

def axis_parameters(data):
    values = data["Log2_odds_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
    max_abs = values.abs().max() if not values.empty else 1.0
    if max_abs <= 2.0:
        step = 0.5
        limit = np.ceil(max_abs * 2) / 2
    else:
        step = 1.0
        limit = np.ceil(max_abs)
    if limit < 1.0:
        limit = 1.0
    margin = step * 0.18
    return (-limit - margin, limit + margin), np.arange(-limit, limit + 0.0001, step)

def draw_original_style(data):
    """Reproduce the original Fig. 4b single-axis grouped horizontal bars."""
    import matplotlib.pyplot as plt

    style()  # Use the current shared Arial file, then retain original figure typography.
    plt.rcParams.update({
        "font.size": 5.5,
        "axes.labelsize": 6.0,
        "axes.titlesize": 6.0,
        "xtick.labelsize": 5.2,
        "ytick.labelsize": 5.2,
        "legend.fontsize": 5.0,
        "axes.linewidth": 0.55,
        "lines.linewidth": 0.7,
        "patch.linewidth": 0.45,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    panel_width_mm, panel_height_mm = 79.0, 64.0
    fs_title, fs_axis, fs_tick, fs_legend = 6.0, 6.0, 5.2, 5.0
    fs_legend_title, fs_star, fs_panel = 5.4, 5.4, 7.0
    lw_axis, lw_grid, lw_zero, lw_bar_edge = 0.55, 0.30, 0.50, 0.20
    xlim, xticks = axis_parameters(data)

    frame = data.copy()
    frame["Region"] = pd.Categorical(frame["Region"], categories=Y_ORDER, ordered=True)
    frame["Group"] = pd.Categorical(frame["Group"], categories=GROUPS, ordered=True)
    frame = frame.sort_values(["Region", "Group"])
    template = pd.DataFrame({
        "Region": pd.Categorical(Y_ORDER, categories=Y_ORDER, ordered=True)
    })

    fig, ax = plt.subplots(
        figsize=(panel_width_mm / 25.4, panel_height_mm / 25.4),
        constrained_layout=False,
        facecolor="white",
    )
    fig.subplots_adjust(left=0.25, right=0.69, bottom=0.18, top=0.86)

    y_base = np.arange(len(Y_ORDER))
    n_groups = len(GROUPS)
    bar_height = 0.72 / n_groups
    for index, y in enumerate(y_base):
        ax.axhspan(y - 0.5, y + 0.5,
                   facecolor="#F3F5F6" if index % 2 == 0 else "#FFFFFF",
                   edgecolor="none", zorder=0)

    for group_index, group in enumerate(GROUPS):
        subset = frame.loc[frame["Group"].eq(group),
                           ["Region", "Log2_odds_ratio", "Fisher_p"]].copy()
        subset = pd.merge(template, subset, on="Region", how="left")
        subset["Log2_odds_ratio"] = pd.to_numeric(
            subset["Log2_odds_ratio"], errors="coerce").fillna(0.0)
        subset["Fisher_p"] = pd.to_numeric(
            subset["Fisher_p"], errors="coerce").fillna(1.0)
        y_locs = y_base + (group_index - n_groups / 2 + 0.5) * bar_height
        bars = ax.barh(
            y_locs, subset["Log2_odds_ratio"].to_numpy(float),
            height=bar_height, color=TYPE_COLORS[group],
            edgecolor="white", linewidth=lw_bar_edge, zorder=3,
        )
        for bar, p_value in zip(bars, subset["Fisher_p"].to_numpy(float)):
            star = get_star(p_value)
            if not star:
                continue
            width = bar.get_width()
            y = bar.get_y() + bar.get_height() / 2
            ax.annotate(
                star, xy=(width, y),
                xytext=(1.5 if width >= 0 else -1.5, 0),
                textcoords="offset points",
                ha="left" if width >= 0 else "right", va="center",
                fontsize=fs_star, fontweight="bold", color=COLOR_AXIS,
                zorder=5, clip_on=False,
            )

    ax.set_title("Genomic locations", fontsize=fs_title,
                 fontweight="normal", color=COLOR_AXIS, pad=2.0)
    ax.axvline(0, color=COLOR_AXIS, linewidth=lw_zero, zorder=2)
    ax.grid(axis="x", linestyle=(0, (2.0, 2.0)), linewidth=lw_grid,
            color=COLOR_GRID, alpha=0.55, zorder=1)
    ax.set_xlabel(r"$\log_2$(Odds ratio)", fontsize=fs_axis,
                  fontweight="normal", labelpad=1.5)
    ax.set_yticks(y_base)
    ax.set_yticklabels(Y_ORDER, fontsize=5.2, fontweight="normal")
    ax.set_xlim(xlim)
    ax.set_xticks(xticks)
    ax.tick_params(axis="x", labelsize=fs_tick, width=lw_axis, length=2.0,
                   pad=1.0, colors=COLOR_AXIS)
    ax.tick_params(axis="y", width=lw_axis, length=2.0,
                   pad=1.2, colors=COLOR_AXIS)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(lw_axis)
        spine.set_color(COLOR_AXIS)

    legend_handles = [
        Patch(facecolor=TYPE_COLORS[group], edgecolor="none",
              label=TYPE_LABELS[group])
        for group in LEGEND_ORDER
    ]
    legend = ax.legend(
        handles=legend_handles, title="Loci category",
        loc="upper left", bbox_to_anchor=(0.705, 0.78),
        bbox_transform=fig.transFigure, frameon=False,
        fontsize=fs_legend, title_fontsize=fs_legend_title,
        handlelength=1.0, handleheight=0.55, handletextpad=0.35,
        labelspacing=0.32, borderaxespad=0,
    )
    legend._legend_box.align = "left"
    fig.text(0.018, 0.972, "b", ha="left", va="top",
             fontsize=fs_panel, fontweight="bold", color=COLOR_AXIS)
    return fig
data = pd.read_csv(ROOT / "data" / "Fig4b.csv")
finish(draw_original_style(data))
