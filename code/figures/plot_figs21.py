"""Draw the per-cancer FigS21 Top-N C-index curves."""

import gzip
import importlib.util
import os
from pathlib import Path
import sys

import pandas as pd


ROOT = Path(os.environ.get(
    "REPRO_ROOT", str(Path(__file__).resolve().parents[2])
))
REVISION_ROOT = Path(os.environ.get("REVISION_ROOT", str(ROOT.parent)))
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA_FILE = ROOT / "data" / "FigS21.csv"
OUTPUT = ROOT / "figures" / "FigS21.pdf"
METHODS = ("m", "hm", "m+hm")
TOP_NS = (5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)


def load_common():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    spec = importlib.util.spec_from_file_location(
        "figs21_v56_plot_common", CODE / "plot_tcga_v56_original80_fold20_common.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_data():
    path = DATA_FILE if DATA_FILE.is_file() else Path(str(DATA_FILE) + ".gz")
    if not path.is_file():
        raise FileNotFoundError("Prepare FigS21 data first")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        data = pd.read_csv(handle, low_memory=False)
    columns = ["Cancer", "Method", "TopN", "CIndex"]
    if list(data.columns) != columns or data.empty:
        raise ValueError("FigS21 data has an unexpected format")
    data["Cancer"] = data["Cancer"].astype(str)
    data["Method"] = data["Method"].astype(str)
    data["TopN"] = pd.to_numeric(data["TopN"], errors="coerce")
    data["CIndex"] = pd.to_numeric(data["CIndex"], errors="coerce")
    if data[["TopN", "CIndex"]].isna().any().any():
        raise ValueError("FigS21 data has missing numeric values")
    data["TopN"] = data["TopN"].astype(int)
    if (not data["Method"].isin(METHODS).all()
            or not data["CIndex"].between(0, 1).all()
            or data.duplicated(["Cancer", "Method", "TopN"]).any()):
        raise ValueError("FigS21 data contains invalid or duplicate points")
    expected = {(method, top_n) for method in METHODS for top_n in TOP_NS}
    for cancer, group in data.groupby("Cancer", sort=False):
        actual = set(zip(group["Method"], group["TopN"]))
        if actual != expected:
            raise ValueError("Incomplete method/TopN grid for {}".format(cancer))
    return data


def run():
    common = load_common()
    data = read_data()
    curves = data.rename(columns={
        "Cancer": "project",
        "Method": "model",
        "TopN": "selected_top_n",
        "CIndex": "harrell_c",
    })
    curves["outer_fold"] = 0

    style = common.load_style("plot_tcga_figs7d_pooled_nested_original_style.py")
    common.configure_style(style, compact=True)
    figure, projects = common.draw_s7d_compat(style, curves)
    for label in list(figure.texts):
        if label.get_text() == "d":
            label.remove()
    figure.set_size_inches(168 / 25.4, 111 / 25.4, forward=True)
    limits = common.y_limits(curves["harrell_c"].to_numpy(float), "s7d")
    for ax in figure.axes[:len(projects)]:
        ax.set_ylim(*limits)
        ticks = sorted(set([limits[0], .5, .7, .9, 1.0]))
        ax.set_yticks([tick for tick in ticks if limits[0] <= tick <= limits[1]])
        ax.set_yticklabels(["{:.1f}".format(tick) for tick in ax.get_yticks()])

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
