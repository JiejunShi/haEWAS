"""Draw the positive and negative liver contribution panels in FigS22."""

import gzip
import importlib
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA_FILE = ROOT / "data" / "FigS22.csv"
OUTPUT = ROOT / "figures" / "FigS22.pdf"
CELLS = ("Hep", "Chol", "EC", "Kup", "Lym")
DIRECTIONS = ("positive", "negative")
ROWS = ("haEWAS-specific", "CAMDA", "CHALM", "EWAS")
PARTS = ("CAMDA-only", "CHALM-only", "Shared")
COLORS = ("#3B9DAA", "#E58C36", "#A16B99")
FIG_W_MM, FIG_H_MM = 88.0, 75.0


def read_data():
    path = DATA_FILE if DATA_FILE.is_file() else Path(str(DATA_FILE) + ".gz")
    if not path.is_file():
        raise FileNotFoundError("Prepare FigS22 data first")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        data = pd.read_csv(handle, low_memory=False)
    expected = [
        "Dataset", "CellType", "Direction", "haEWAS-specific",
        "CAMDA-only", "CHALM-only", "Shared", "CAMDA", "CHALM", "EWAS",
    ]
    if list(data.columns) != expected or data.empty:
        raise ValueError("FigS22 data has an unexpected format")
    counts = data[list(ROWS)].to_numpy(float)
    if (not np.isfinite(counts).all() or (counts < 0).any()
            or (counts != np.floor(counts)).any()):
        raise ValueError("FigS22 counts must be nonnegative integers")
    if not data["haEWAS-specific"].eq(
            data[["CAMDA-only", "CHALM-only", "Shared"]].sum(axis=1)).all():
        raise ValueError("haEWAS-specific counts do not match their components")
    return data


def style():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    importlib.invalidate_caches()
    from figure_workflow_common import style as apply_style
    apply_style()
    mpl.rcParams.update({
        "font.size": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.titlesize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def plot_dataset(data, dataset):
    subset = data.loc[data["Dataset"].eq(dataset)]
    expected = {(cell, direction) for cell in CELLS for direction in DIRECTIONS}
    actual = set(subset[["CellType", "Direction"]].itertuples(
        index=False, name=None
    ))
    if actual != expected or subset.duplicated(["CellType", "Direction"]).any():
        raise ValueError("Incomplete cell/direction grid: {}".format(dataset))

    values = subset[list(ROWS)].to_numpy(float)
    norm = Normalize(0, max(1, np.log1p(values).max()))
    cmap = plt.get_cmap("Blues")
    figure = plt.figure(figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4))

    for panel, direction in enumerate(DIRECTIONS):
        left, width = (29 + 29 * panel) / 88, 20 / 88
        top = figure.add_axes([left, 39 / 75, width, 24 / 75])
        heat = figure.add_axes([left, 15 / 75, width, 20 / 75])
        block = (subset.loc[subset["Direction"].eq(direction)]
                 .set_index("CellType").loc[list(CELLS)])
        denominator = block["haEWAS-specific"].to_numpy(float)
        bottom = np.zeros(len(CELLS))
        for part, color in zip(PARTS, COLORS):
            heights = np.divide(
                block[part].to_numpy(float), denominator,
                out=np.zeros(len(CELLS)), where=denominator > 0,
            )
            top.bar(np.arange(5), heights, bottom=bottom, color=color, width=.82,
                    edgecolor="white", linewidth=.35)
            bottom += heights
        for i, total in enumerate(denominator):
            if not total:
                top.text(i, .03, "NA", ha="center", va="bottom",
                         rotation=90, fontsize=6)
        top.set(xlim=(-.5, 4.5), ylim=(0, 1), xticks=[],
                yticks=[0, .25, .5, .75, 1])
        top.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        top.set_title(direction.capitalize() + " association", pad=4)
        if panel == 0:
            top.set_ylabel("haEWAS-specific DMCTs", fontsize=6, labelpad=2)

        matrix = block[list(ROWS)].to_numpy(float).T
        heat.imshow(np.log1p(matrix), cmap=cmap, norm=norm, aspect="auto",
                    interpolation="nearest")
        for row in range(4):
            for col in range(5):
                heat.text(col, row, str(int(matrix[row, col])), ha="center",
                          va="center", rotation=90, fontsize=6,
                          color="white" if norm(np.log1p(matrix[row, col])) > .58
                          else "#222222")
        heat.set_xticks(range(5))
        heat.set_xticklabels(CELLS, rotation=90, ha="center", va="top")
        heat.set_yticks(range(4))
        heat.set_yticklabels(ROWS if panel == 0 else [""] * 4)
        heat.set_xticks(np.arange(-.5, 5, 1), minor=True)
        heat.set_yticks(np.arange(-.5, 4, 1), minor=True)
        heat.grid(which="minor", color="white", linewidth=.5)
        heat.tick_params(which="both", length=0)
        for ax in (top, heat):
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(.5)
            ax.tick_params(width=.5, length=2, pad=1.5)

    figure.legend(
        handles=[Patch(facecolor=color, label=label)
                 for color, label in zip(
                     COLORS, ("CAMDA-only", "CHALM-only", "CHALM–CAMDA shared")
                 )],
        loc="upper center", bbox_to_anchor=(.54, .997), ncol=3,
        frameon=False, handlelength=1, handletextpad=.3, columnspacing=.8,
    )
    bar = figure.add_axes([80 / 88, 18 / 75, 1.5 / 88, 13 / 75])
    figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=bar, ticks=[])
    bar.text(.5, 1.06, "More", transform=bar.transAxes,
             ha="center", fontsize=6)
    bar.text(.5, -.06, "Less", transform=bar.transAxes,
             ha="center", va="top", fontsize=6)

    return figure


def run():
    style()
    data = read_data()
    datasets = sorted(data["Dataset"].astype(str).unique())
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    for dataset in datasets:
        figure = plot_dataset(data, dataset)
        figure.savefig(OUTPUT, format="pdf")
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
