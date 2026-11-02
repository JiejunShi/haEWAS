"""Draw the TCGA discovery attribution panels in FigS20."""

import gzip
import importlib
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA_FILE = ROOT / "data" / "FigS20.csv"
OUTPUT = ROOT / "figures" / "FigS20.pdf"
FIG_W_MM, FIG_H_MM = 100.6666666667, 96.0
DATA_LEFT_MM = 28.0
DATA_WIDTH_MM = 88.0 * (2.0 / 3.0)
FRACTION_COLUMNS = ("CAMDA only", "CHALM only", "CHALM-CAMDA shared")
COUNT_COLUMNS = ("haEWAS-specific", "CAMDA", "CHALM", "EWAS")
COUNT_LABELS = COUNT_COLUMNS
FRACTION_LABELS = ("CAMDA only", "CHALM only", "CHALM–CAMDA shared")
COLORS = ("#3B9DAA", "#E58C36", "#A16B99")


def read_data():
    path = DATA_FILE if DATA_FILE.is_file() else Path(str(DATA_FILE) + ".gz")
    if not path.is_file():
        raise FileNotFoundError("Prepare FigS20 data first")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        data = pd.read_csv(handle, low_memory=False)
    expected = ["Cancer"] + list(FRACTION_COLUMNS) + list(COUNT_COLUMNS)
    if list(data.columns) != expected or len(data) != 15:
        raise ValueError("FigS20 data has an unexpected format")
    values = data[list(FRACTION_COLUMNS + COUNT_COLUMNS)].to_numpy(float)
    if (not np.isfinite(values).all() or (values < 0).any()
            or (values != np.floor(values)).any()):
        raise ValueError("FigS20 counts must be nonnegative integers")
    if not np.array_equal(data[list(FRACTION_COLUMNS)].sum(axis=1),
                          data["haEWAS-specific"]):
        raise ValueError("haEWAS-specific components do not sum to their total")
    return data


def run():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    importlib.invalidate_caches()
    from figure_workflow_common import style

    style()
    mpl.rcParams.update({
        "font.size": 6.5,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 6.0,
        "ytick.labelsize": 6.0,
        "legend.fontsize": 6.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    })
    data = read_data()
    fractions = np.divide(
        data[list(FRACTION_COLUMNS)].to_numpy(float),
        data["haEWAS-specific"].to_numpy(float)[:, None],
        out=np.full((len(data), len(FRACTION_COLUMNS)), np.nan),
        where=data["haEWAS-specific"].to_numpy(float)[:, None] > 0,
    )
    counts = data[list(COUNT_COLUMNS)].to_numpy(float).T

    figure = plt.figure(figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4),
                        facecolor="white")

    def rectangle(left, bottom, width, height):
        return [left / FIG_W_MM, bottom / FIG_H_MM,
                width / FIG_W_MM, height / FIG_H_MM]

    left, width = DATA_LEFT_MM, DATA_WIDTH_MM
    top = figure.add_axes(rectangle(left, 46, width, 32))
    x = np.arange(len(data))
    bottom = np.zeros(len(data))
    for i, color in enumerate(COLORS):
        height = np.nan_to_num(fractions[:, i], nan=0.0) * 100
        top.bar(x, height, bottom=bottom, width=.91, color=color,
                edgecolor="white", linewidth=.35, zorder=3)
        bottom += height
    for i, total in enumerate(data["haEWAS-specific"]):
        if total == 0:
            top.bar(i, 100, width=.91, color="#EEEEEE", edgecolor="white",
                    linewidth=.35, zorder=2)
            top.text(i, 50, "NA", ha="center", va="center", rotation=90,
                     fontsize=6)
    top.set_xlim(-.5, len(data) - .5)
    top.set_ylim(0, 100)
    top.set_xticks([])
    top.set_yticks([0, 25, 50, 75, 100])
    top.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    top.set_ylabel("haEWAS-specific loci (%)", labelpad=2)
    top.set_axisbelow(True)
    top.grid(axis="y", color="#E3E3E3", linewidth=.25)
    for side in ("top", "right"):
        top.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        top.spines[side].set_linewidth(.4)
    top.tick_params(axis="y", width=.4, length=1.5, pad=1)

    heat = figure.add_axes(rectangle(left, 12, width, 32))
    cmap = LinearSegmentedColormap.from_list(
        "figs20_count_blues", ("#EEF2F7", "#C8D6EA", "#8FAFD3",
                               "#4F8FC0", "#0F5B8C"), N=256,
    )
    cmap.set_bad("#EEEEEE")
    transformed = np.log1p(counts)
    norm = Normalize(vmin=0, vmax=max(float(transformed.max()), 1.0))
    image = heat.pcolormesh(
        np.arange(len(data) + 1) - .5, np.arange(5) - .5,
        np.ma.masked_where(counts == 0, transformed), cmap=cmap, norm=norm,
        edgecolors="white", linewidth=.35, shading="flat",
    )
    for row in range(4):
        for col in range(len(data)):
            heat.text(col, row, "{:,}".format(int(counts[row, col])),
                      ha="center", va="center", rotation=90, fontsize=6,
                      color="white" if norm(transformed[row, col]) > .55
                      else "#272727")
    heat.set_xlim(-.5, len(data) - .5)
    heat.set_ylim(3.5, -.5)
    heat.set_yticks(np.arange(4))
    heat.set_yticklabels(COUNT_LABELS)
    heat.set_xticks(x)
    heat.set_xticklabels(data["Cancer"], rotation=90, ha="center", va="top")
    heat.tick_params(axis="both", length=0, pad=2, colors="#272727")
    for spine in heat.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.35)

    cax = figure.add_axes(rectangle(left + width + 4, 15, 1.8, 24))
    colorbar = figure.colorbar(image, cax=cax)
    colorbar.set_ticks([])
    colorbar.outline.set_linewidth(.35)
    cax.text(1.7, 1, "More", transform=cax.transAxes,
             ha="left", va="center", fontsize=6)
    cax.text(1.7, 0, "Less", transform=cax.transAxes,
             ha="left", va="center", fontsize=6)
    cax.text(.5, 1.12, "Loci", transform=cax.transAxes,
             ha="center", va="bottom", fontsize=6)
    handles = [Patch(facecolor=color, edgecolor="none", label=label)
               for color, label in zip(COLORS, FRACTION_LABELS)]
    figure.legend(handles=handles, loc="upper center",
                  bbox_to_anchor=((left + width / 2) / FIG_W_MM, .95),
                  ncol=3, frameon=False, handlelength=.9,
                  handletextpad=.35, columnspacing=1.0, borderaxespad=0)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", facecolor="white")
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None and hasattr(get_ipython(), "kernel"):
            display(figure)
    except ImportError:
        pass
    plt.close(figure)


if __name__ == "__main__":
    run()
