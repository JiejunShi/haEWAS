"""Draw Fig4a using only its prepared CSV; export PDF and display inline."""
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(os.environ["REPRO_ROOT"])
PDF = ROOT / "figures" / "Fig4a.pdf"

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

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

TYPE_COLORS = {
    "EWAS-specific": "#0F4D92",
    "Common": "#8F8F8F",
    "haEWAS-specific": "#B64342",
}

FIG_W_MM, FIG_H_MM = 84.0, 64.0

MIN_BUBBLE, MAX_BUBBLE = 10.0, 38.0

def _bubble_size(value, low, high):
    if high == low:
        return (MIN_BUBBLE + MAX_BUBBLE) / 2
    return MIN_BUBBLE + (float(value) - low) / (high - low) * (MAX_BUBBLE - MIN_BUBBLE)

def _draw_figure(data):
    import matplotlib as mpl
    style()
    mpl.rcParams.update({
        "font.size": 5.5, "axes.labelsize": 6.0, "axes.titlesize": 6.0,
        "xtick.labelsize": 5.2, "ytick.labelsize": 5.0, "legend.fontsize": 5.0,
        "axes.linewidth": 0.55, "lines.linewidth": 0.7, "patch.linewidth": 0.45,
    })
    fig, ax = plt.subplots(figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4),
                           constrained_layout=False, facecolor="white")
    fig.subplots_adjust(left=0.42, right=0.72, bottom=0.17, top=0.86)

    low, high = float(data["Overlap_genes"].min()), float(data["Overlap_genes"].max())
    data = data.copy()
    data["bubble_size"] = data["Overlap_genes"].map(lambda value: _bubble_size(value, low, high))
    y = np.arange(len(data))
    for index, row in data.iterrows():
        ax.scatter(row["Log_fdr"], y[index], s=row["bubble_size"],
                   color=TYPE_COLORS[row["Category"]], alpha=0.88,
                   edgecolors="white", linewidth=0.25, zorder=3)

    ax.axvline(-np.log10(0.05), color="#272727", linestyle=(0, (2.2, 2.2)),
               linewidth=0.50, alpha=0.85, zorder=1)
    ax.grid(axis="x", linestyle=(0, (2.0, 2.0)), linewidth=0.30,
            color="#BDBDBD", alpha=0.55, zorder=0)
    counts = [int(data["Category"].eq(category).sum())
              for category in ("EWAS-specific", "Common", "haEWAS-specific")]
    cumulative = 0
    for count in counts[:-1]:
        cumulative += count
        ax.axhline(cumulative - 0.5, color="#BDBDBD", linewidth=0.30,
                   alpha=0.75, zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels(data["Term"])
    for label, category in zip(ax.get_yticklabels(), data["Category"]):
        label.set_color(TYPE_COLORS[category])
        label.set_fontsize(5.0)
    ax.set_xlabel(r"$-\log_{10}(\mathrm{FDR})$", fontsize=6.0, labelpad=1.5)
    ax.tick_params(axis="x", labelsize=5.2, width=0.55, length=2.0,
                   pad=1.0, colors="#272727")
    ax.tick_params(axis="y", width=0.55, length=2.0, pad=1.5)
    xmax = max(6, int(np.ceil((float(data["Log_fdr"].max()) + 0.5) / 5) * 5))
    ax.set_xlim(0, xmax)
    ax.xaxis.set_major_locator(MultipleLocator(5 if xmax > 20 else 2))
    ax.set_title("GO enrichment", fontsize=6.0, fontweight="normal", pad=2.0,
                 color="#272727")
    fig.text(0.018, 0.972, "a", ha="left", va="top", fontsize=7.0,
             fontweight="bold", color="#272727")

    color_handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=3,
                            markerfacecolor=TYPE_COLORS[category], markeredgecolor="none",
                            label=category)
                     for category in ("EWAS-specific", "Common", "haEWAS-specific")]
    leg1 = ax.legend(handles=color_handles, title="Loci category", loc="upper left",
                     bbox_to_anchor=(1.02, 1.0), frameon=False,
                     fontsize=5.0, title_fontsize=5.4, borderaxespad=0.0)
    ax.add_artist(leg1)
    sizes = sorted(set([int(low), int(round((low + high) / 2)), int(high)]))
    size_handles = [plt.scatter([], [], s=_bubble_size(value, low, high),
                                color="#A0A0A0", edgecolors="none", label=str(value))
                    for value in sizes]
    ax.legend(handles=size_handles, title="Overlap genes", loc="lower left",
              bbox_to_anchor=(1.02, 0.0), frameon=False,
              fontsize=5.0, title_fontsize=5.4, borderaxespad=0.0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.55)
        spine.set_color("#272727")
    return fig
data = pd.read_csv(ROOT / "data" / "Fig4a.csv")
fig = _draw_figure(data)
finish(fig, bbox_inches="tight")
