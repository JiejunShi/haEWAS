"""Plot the four matched-group panels for FigS2."""

import contextlib
import importlib.util
import io
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = ROOT.resolve().parent
os.environ["REVISION_ROOT"] = str(REVISION_ROOT)
SOURCE = REVISION_ROOT / "A2H" / "code" / "plot_response_1_5_original_figs1_raw2_min2_mainonly_v1.py"
PANELS = (
    ("CAMDA", "FigS2a.csv", "FigS2c.csv"),
    ("CHALM", "FigS2b.csv", "FigS2d.csv"),
)
METHODS = (("A2H", "A2H"), ("m2", "M2"), ("m3", "M3"))


def read_panel(filename):
    data = pd.read_csv(ROOT / "data" / filename)
    if list(data.columns) != ["Group", "A2H", "M2", "M3"] or data.empty:
        raise ValueError("Unexpected data in " + filename)
    if data["Group"].duplicated().any():
        raise ValueError("Duplicate groups in " + filename)
    return data


def load_renderer():
    spec = importlib.util.spec_from_file_location("figs2_source", SOURCE)
    source = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        spec.loader.exec_module(source)
    return source


def matched_groups():
    rows = []
    group_order = None
    for target, pearson_file, rmse_file in PANELS:
        pearson = read_panel(pearson_file)
        rmse = read_panel(rmse_file)
        current = pearson["Group"].astype(str).tolist()
        if current != rmse["Group"].astype(str).tolist():
            raise ValueError("Group order differs between " + pearson_file + " and " + rmse_file)
        if group_order is not None and current != group_order:
            raise ValueError("The four panels do not use the same matched groups")
        group_order = current
        for method, column in METHODS:
            rows.append(pd.DataFrame({
                "group_id": pearson["Group"],
                "target": target,
                "method": method,
                "pearson_r": pearson[column],
                "rmse": rmse[column],
            }))
    return pd.concat(rows, ignore_index=True)


def run():
    source = load_renderer()
    engine = source.engine
    frame = matched_groups()
    group_count = frame["group_id"].nunique()
    output = ROOT / "figures" / "FigS2.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    old_targets, old_show = engine.TARGETS, engine.show_pdf

    def save(fig, path):
        source.renderer.fix_axes(fig)
        if fig._suptitle is not None:
            fig._suptitle.remove()
        fig.axes[0].set_title("hm(CAMDA) (N={:,} groups)".format(group_count), pad=5)
        fig.axes[1].set_title("hm(CHALM) (N={:,} groups)".format(group_count), pad=5)
        for ax in fig.axes[:2]:
            ax.set_ylim(-.5, 1)
            ax.set_yticks([-.5, 0, .5, 1])
            ax.spines["left"].set_bounds(-.5, 1)
        for ax in fig.axes[2:]:
            ax.set_ylim(0, 1)
            ax.set_yticks([0, .25, .5, .75, 1])
            ax.spines["left"].set_bounds(0, 1)
        fig.subplots_adjust(top=.92)
        fig.savefig(path, format="pdf", dpi=300, facecolor="white")
        try:
            from IPython import get_ipython
            from IPython.display import display
            if get_ipython() is not None:
                display(fig)
            else:
                plt.show()
        except ImportError:
            plt.show()
        plt.close(fig)

    engine.TARGETS = ("CAMDA", "CHALM")
    engine.show_pdf = save
    try:
        with mpl.rc_context({
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.labelsize": 6.5,
            "axes.titlesize": 7,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 6,
            "axes.linewidth": .55,
            "pdf.fonttype": 42,
        }):
            engine.distribution_figure(
                frame, output,
                "Matched 10-CpG groups; at least 2 valid participants/CpG",
                unit="groups", centered=False,
            )
    finally:
        engine.TARGETS = old_targets
        engine.show_pdf = old_show


if __name__ == "__main__":
    run()
