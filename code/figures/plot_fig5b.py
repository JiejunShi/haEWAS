"""Draw Fig5b from the prepared C-index summary."""

import gzip
import importlib.util
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get(
    "REPRO_ROOT", str(Path(__file__).resolve().parents[2])
))
REVISION_ROOT = Path(os.environ.get("REVISION_ROOT", str(ROOT.parent)))
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA_FILE = ROOT / "data" / "Fig5b.csv"
OUTPUT = ROOT / "figures" / "Fig5b.pdf"
METHODS = ("m", "hm", "m+hm")
TOP_NS = (5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)


def load_common():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    spec = importlib.util.spec_from_file_location(
        "fig5b_v56_plot_common", CODE / "plot_tcga_v56_original80_fold20_common.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_data():
    path = DATA_FILE if DATA_FILE.is_file() else Path(str(DATA_FILE) + ".gz")
    if not path.is_file():
        raise FileNotFoundError("Prepare Fig5b data first")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        data = pd.read_csv(handle, low_memory=False)
    columns = ["Method", "TopN", "MedianCIndex", "Q25CIndex", "Q75CIndex"]
    if list(data.columns) != columns or data.empty:
        raise ValueError("Fig5b data has an unexpected format")
    data["Method"] = data["Method"].astype(str)
    for column in columns[1:]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    if data[columns[1:]].isna().any().any():
        raise ValueError("Fig5b data has missing numeric values")
    data["TopN"] = data["TopN"].astype(int)
    expected = {(method, n) for method in METHODS for n in TOP_NS}
    actual = set(zip(data["Method"], data["TopN"]))
    if actual != expected or data.duplicated(["Method", "TopN"]).any():
        raise ValueError("Fig5b data has an incomplete method/TopN grid")
    if not data[["MedianCIndex", "Q25CIndex", "Q75CIndex"]].apply(
            lambda column: column.between(0, 1).all()).all():
        raise ValueError("Fig5b C-index values must lie in [0, 1]")
    if (data["Q25CIndex"].gt(data["MedianCIndex"])
            | data["MedianCIndex"].gt(data["Q75CIndex"])).any():
        raise ValueError("Fig5b quartiles are not ordered around the median")
    return data


def run():
    common = load_common()
    style = common.load_style("plot_tcga_fig5b_pooled_nested_original_style.py")
    common.configure_style(style)
    style.configure_style()

    data = read_data()
    figure = style.plt.figure(
        figsize=(72 / 25.4, 48 / 25.4), facecolor="white",
        constrained_layout=False,
    )
    ax = figure.add_axes([.105, .155, .880, .755])
    for method in METHODS:
        part = data.loc[data["Method"].eq(method)].sort_values("TopN")
        x = part["TopN"].to_numpy(float)
        median = part["MedianCIndex"].to_numpy(float)
        lower = part["Q25CIndex"].to_numpy(float)
        upper = part["Q75CIndex"].to_numpy(float)
        color = style.MODEL_COLORS[method]
        ax.fill_between(x, lower, upper, color=color, alpha=.12,
                        linewidth=0, zorder=1)
        ax.plot(x, median, marker=style.MODEL_MARKERS[method],
                linestyle=style.MODEL_LINESTYLES[method], color=color,
                label=style.MODEL_LABELS[method],
                lw=style.MODEL_LINEWIDTHS[method],
                alpha=style.MODEL_ALPHAS[method],
                markersize=style.MODEL_MARKERSIZES[method],
                zorder=style.MODEL_ZORDERS[method], markeredgewidth=0)

    ax.axhline(.5, color="gray", linestyle="--", linewidth=style.LW_RANDOM,
               label="Random, 0.5", zorder=0, alpha=.75)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor(style.COLOR_AXIS)
        spine.set_linewidth(style.LW_AXIS)
    ax.tick_params(axis="both", which="major", direction="out", length=2,
                   width=style.LW_AXIS, colors=style.COLOR_AXIS,
                   labelsize=style.FS_TICK, pad=1.2)
    ax.set_xlim(-2, 105)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Number of features (top N)", fontsize=style.FS_AXIS,
                  fontweight="normal", labelpad=.8, color=style.COLOR_AXIS)
    q_values = data[["Q25CIndex", "Q75CIndex"]].to_numpy(float)
    y_min, y_max = common.y_limits(q_values, "5b")
    ax.set_ylim(y_min, y_max)
    ax.set_yticks([tick / 10 for tick in range(11)
                   if y_min <= tick / 10 <= y_max])
    ax.set_ylabel("Median C-index", fontsize=style.FS_AXIS,
                  fontweight="normal", labelpad=.8, color=style.COLOR_AXIS)
    figure.text(.55, .975, "Prognostic performance", ha="center", va="top",
                fontsize=style.FS_TITLE, fontweight="normal",
                color=style.COLOR_AXIS)
    ax.legend(loc="upper left", bbox_to_anchor=(.010, .985), frameon=False,
              fontsize=style.FS_LEGEND, handlelength=1.05,
              handletextpad=.30, labelspacing=.18, borderaxespad=0)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", facecolor="white",
                   bbox_inches=None, pad_inches=0)
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None and hasattr(get_ipython(), "kernel"):
            display(figure)
    except ImportError:
        pass
    style.plt.close(figure)


if __name__ == "__main__":
    run()
