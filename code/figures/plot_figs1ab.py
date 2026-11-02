"""Plot FigS1a/b from the prepared original CpG-sample pairs."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize


ROOT = Path(os.environ["REPRO_ROOT"])
GRID_SIZE, SIGMA = 50, 1.0
TARGETS = (
    ("CAMDA", "a", "WGBS_CAMDA", "Predicted_CAMDA",
     "Raw CAMDA from WGBS", "A2H inferred CAMDA from array"),
    ("CHALM", "b", "WGBS_CHALM", "Predicted_CHALM",
     "Raw CHALM from WGBS", "A2H inferred CHALM from array"),
)
PALETTE = ("#F4F7F9", "#B9CDD9", "#6C8EAD", "#244960")


def style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 6,
        "axes.titlesize": 6.5,
        "axes.labelsize": 6,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 5.5,
        "axes.linewidth": 0.5,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def smooth_histogram(grid):
    radius = max(1, int(np.ceil(4 * SIGMA)))
    offsets = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (offsets / SIGMA) ** 2)
    kernel /= kernel.sum()
    result = np.asarray(grid, dtype=float).copy()
    for axis in (0, 1):
        pad = [(0, 0), (0, 0)]
        pad[axis] = (radius, radius)
        result = np.apply_along_axis(
            lambda values: np.convolve(values, kernel, mode="valid"),
            axis,
            np.pad(result, pad, mode="symmetric"),
        )
    return result


def count_to_color(count):
    count = np.clip(np.asarray(count, dtype=float), 0, 10000)
    result = np.empty_like(count)
    low = count <= 10
    tens = (count > 10) & (count <= 100)
    hundreds = (count > 100) & (count <= 1000)
    thousands = count > 1000
    result[low] = np.sqrt(count[low] / 10.0) / 4.0
    result[tens] = 0.25 + (count[tens] - 10) / 90 / 4
    result[hundreds] = 0.50 + (count[hundreds] - 100) / 900 / 4
    result[thousands] = 0.75 + (count[thousands] - 1000) / 9000 / 4
    return result


def _plot_target(target, panel, value_columns, x_label, y_label, metrics):
    data = pd.read_csv(
        ROOT / "data" / ("FigS1" + panel + ".csv"),
        dtype={"CpG": str, "Sample": str},
        float_precision="round_trip",
    )
    expected_columns = ["CpG", "Sample", value_columns[0], value_columns[1]]
    if list(data.columns) != expected_columns:
        raise ValueError("FigS1 pair table has an unexpected layout")
    if data.duplicated(["CpG", "Sample"]).any():
        raise ValueError("Duplicate CpG-sample pairs in FigS1 data")
    if not np.isfinite(data[list(value_columns)].to_numpy(float)).all():
        raise ValueError("Invalid original FigS1 pair values")

    target_metrics = metrics.loc[metrics["Target"].eq(target)]
    samples = target_metrics["Sample"].tolist()
    if len(samples) != 10 or target_metrics["Sample"].duplicated().any():
        raise ValueError("Expected ten participant metrics for FigS1")
    if set(data["Sample"]) != set(samples):
        raise ValueError("FigS1 pair and metric participants do not match")
    if not np.isfinite(target_metrics[["R", "RMSE", "NRMSE"]].to_numpy(float)).all():
        raise ValueError("Invalid FigS1 participant metrics")

    style()
    cmap = LinearSegmentedColormap.from_list("figs1_blue", PALETTE)
    fig, axes = plt.subplots(
        2, 5, figsize=(129 / 25.4, 88 / 25.4), sharex=True, sharey=True
    )
    fig.subplots_adjust(
        left=0.085, right=0.985, bottom=0.25, top=0.88,
        wspace=0.18, hspace=0.55,
    )
    image = None
    axis_pad = 0.02
    edges = np.linspace(0, 1, GRID_SIZE + 1)

    for ax, sample in zip(axes.flat, samples):
        panel_data = data.loc[data["Sample"].eq(sample)]
        x = panel_data[value_columns[0]].to_numpy(dtype=float)
        y = panel_data[value_columns[1]].to_numpy(dtype=float)
        histogram = np.histogram2d(
            np.clip(x, 0, 1), np.clip(y, 0, 1), bins=(edges, edges)
        )[0]
        grid = smooth_histogram(histogram)
        color_values = np.ma.array(count_to_color(grid.T), mask=grid.T <= 0)
        image = ax.imshow(
            color_values, origin="lower", extent=(0, 1, 0, 1),
            interpolation="bilinear", cmap=cmap, norm=Normalize(0, 1),
            aspect="auto", rasterized=True,
        )
        ax.plot(
            [-axis_pad, 1 + axis_pad], [-axis_pad, 1 + axis_pad],
            "--", color="#666666", linewidth=0.55,
        )
        ax.set(
            xlim=(-axis_pad, 1 + axis_pad),
            ylim=(-axis_pad, 1 + axis_pad),
            title=sample.replace("PGPUK_", ""),
        )
        ax.set_aspect("equal")
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda value, _: format(value, "g")))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda value, _: format(value, "g")))
        ax.tick_params(length=2, width=0.5, pad=1)
        for spine in ax.spines.values():
            spine.set_visible(True)

        row = target_metrics.loc[target_metrics["Sample"].eq(sample)].iloc[0]
        label = "r = {:.3f}\nRMSE = {:.3f}\nNRMSE = {:.3f}".format(
            row["R"], row["RMSE"], row["NRMSE"]
        )
        ax.text(
            0.035, 0.965, label, transform=ax.transAxes,
            va="top", fontsize=5.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1},
        )

    fig.suptitle("hm({})".format(target), fontsize=7, y=0.98)
    fig.text(0.535, 0.17, x_label, ha="center", va="center", fontsize=6)
    fig.text(0.018, 0.56, y_label, ha="center", va="center", rotation=90, fontsize=6)
    color_axis = fig.add_axes([0.40, 0.095, 0.27, 0.018])
    colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal")
    tick_counts = np.array([0, 10, 100, 1000, 10000], dtype=float)
    colorbar.set_ticks(count_to_color(tick_counts))
    colorbar.set_ticklabels(["0", "10", "100", "1,000", "10,000"])
    colorbar.set_label("CpGs", fontsize=6)
    colorbar.ax.tick_params(labelsize=5.5, length=2)
    colorbar.ax.xaxis.set_minor_locator(ticker.NullLocator())

    output = ROOT / "figures" / ("FigS1" + panel + ".pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="pdf", facecolor="white", dpi=600)
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


def run():
    metrics = pd.read_csv(
        ROOT / "data" / "FigS1_metrics.csv",
        dtype={"Target": str, "Sample": str},
        float_precision="round_trip",
    )
    if list(metrics.columns) != ["Target", "Sample", "R", "RMSE", "NRMSE"]:
        raise ValueError("FigS1 metric table has an unexpected layout")
    for target, panel, wgbs, predicted, x_label, y_label in TARGETS:
        _plot_target(target, panel, (wgbs, predicted), x_label, y_label, metrics)


if __name__ == "__main__":
    run()