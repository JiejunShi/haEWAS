"""Draw FigS19 from the prepared ROC coordinates."""

import gzip
import importlib
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA_FILE = ROOT / "data" / "FigS19.csv"
OUTPUT = ROOT / "figures" / "FigS19.pdf"
METHODS = ("Baseline", "EWAS", "CHALM", "CAMDA", "haEWAS", "Combined")
COLORS = {
    "Baseline": "#777777",
    "EWAS": "#0F4D92",
    "CHALM": "#E69F00",
    "CAMDA": "#009E73",
    "haEWAS": "#B64342",
    "Combined": "#6A3D9A",
}


def read_data():
    path = DATA_FILE if DATA_FILE.is_file() else Path(str(DATA_FILE) + ".gz")
    if not path.is_file():
        raise FileNotFoundError("Prepare FigS19 data first")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        data = pd.read_csv(handle, low_memory=False)
    expected = ["Method", "FPR", "TPR", "AUC"]
    if list(data.columns) != expected or data.empty:
        raise ValueError("FigS19 data has an unexpected format")
    if not set(METHODS).issubset(set(data["Method"].astype(str))):
        raise ValueError("FigS19 is missing a plotted method")
    return data


def run():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    importlib.invalidate_caches()
    from figure_workflow_common import style

    style()
    mpl.rcParams.update({
        "font.size": 6,
        "axes.labelsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    data = read_data()
    figure = plt.figure(figsize=(55 / 25.4, 55 / 25.4))
    ax = figure.add_axes([12 / 55, 24 / 55, 27 / 55, 27 / 55])
    legend_ax = figure.add_axes([1 / 55, 1 / 55, 53 / 55, 15 / 55])
    legend_ax.set_axis_off()
    ax.plot([0, 1], [0, 1], color="#B8B8B8", linestyle="--", linewidth=.45)

    handles = []
    for method in METHODS:
        part = data.loc[data["Method"].eq(method)].sort_values(
            ["FPR", "TPR"], kind="mergesort"
        )
        aucs = part["AUC"].dropna().unique()
        if len(aucs) != 1:
            raise ValueError("FigS19 has inconsistent AUC values for {}".format(method))
        ax.plot(part["FPR"].to_numpy(float), part["TPR"].to_numpy(float),
                color=COLORS[method], linewidth=1)
        handles.append(Line2D(
            [0], [0], color=COLORS[method], linewidth=1,
            label="{} (AUC={:.3f})".format(method, float(aucs[0])),
        ))

    ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[0, .25, .5, .75, 1],
           yticks=[0, .25, .5, .75, 1], xlabel="1 - specificity",
           ylabel="Sensitivity")
    ax.set_aspect("equal")
    ax.set_title("London blood", fontsize=6, pad=2.5)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.35)
    ax.tick_params(length=1.3, width=.35, pad=1.5)
    legend_ax.legend(handles=handles, loc="center", ncol=2, frameon=False,
                     fontsize=6, handlelength=.8, handletextpad=.25,
                     columnspacing=.6, labelspacing=.4)

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
