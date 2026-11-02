"""Plot FigS6 from the two compact CpG-level files."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
PANELS = (
    ("a", "FigS6a.csv.gz", "hm(CAMDA)"),
    ("b", "FigS6b.csv.gz", "hm(CHALM)"),
)
MODELS = ("WholeAtlasCV", "MinusBloodImmune", "MinusEpithelial")
LABELS = ("Whole-atlas CV", "\u2212Blood/immune", "\u2212Epithelial")
COLORS = ("#335C75", "#8BA8BA", "#C4D5DF")
YMAX = 0.4


def density_curve(values, y):
    bandwidth = float(np.std(values, ddof=1) * len(values) ** (-0.2))
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError("Invalid violin bandwidth")
    low = min(float(values.min()), 0.0) - 4 * bandwidth
    high = max(float(values.max()), float(y[-1])) + 4 * bandwidth
    bins = min(100000, max(200, int(np.ceil((high - low) / (bandwidth / 4)))))
    edges = np.linspace(low, high, bins + 1)
    counts, _ = np.histogram(values, bins=edges)
    step = edges[1] - edges[0]
    radius = max(1, int(np.ceil(4 * bandwidth / step)))
    kernel = np.exp(
        -0.5 * (np.arange(-radius, radius + 1) * step / bandwidth) ** 2
    )
    kernel /= kernel.sum()
    estimate = np.interp(
        y, (edges[:-1] + edges[1:]) / 2,
        np.convolve(counts, kernel, mode="same"),
    )
    return estimate / estimate.max()


def run():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 6,
        "axes.labelsize": 6,
        "axes.titlesize": 6.5,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 5.5,
        "axes.linewidth": 0.55,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })

    fig, axes = plt.subplots(1, 2, figsize=(129 / 25.4, 85 / 25.4))
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.35, top=0.86, wspace=0.35)
    y = np.linspace(0, YMAX, 500)

    for ax, (panel, filename, title) in zip(axes, PANELS):
        data = pd.read_csv(ROOT / "data" / filename)
        if list(data.columns) != ["CpG"] + list(MODELS):
            raise ValueError("{} has unexpected columns".format(filename))
        if data.empty or data["CpG"].isna().any() or data["CpG"].duplicated().any():
            raise ValueError("{} has invalid CpG rows".format(filename))

        for x, (model, color) in enumerate(zip(MODELS, COLORS)):
            values = pd.to_numeric(data[model], errors="coerce").to_numpy(float)
            if not len(values) or not np.isfinite(values).all():
                raise ValueError("Invalid {} values".format(model))
            density = density_curve(values, y)
            ax.fill_betweenx(y, x - 0.32 * density, x + 0.32 * density,
                             color=color, alpha=0.9, linewidth=0)
            q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
            ax.plot([x, x], [q1, q3], color="#244960", linewidth=1.2, zorder=3)
            ax.scatter(x, median, facecolor="white", edgecolor="#244960",
                       s=12, linewidth=0.5, zorder=4)

        ax.set_xticks(range(3))
        ax.set_xticklabels(LABELS, rotation=30, ha="right", va="top", fontsize=5.5)
        ax.set(xlim=(-0.6, 2.6), ylim=(0, YMAX),
               ylabel="Per-CpG RMSE", title=title)
        ax.set_yticks(np.linspace(0, YMAX, 5))
        ax.tick_params(length=2, width=0.5, pad=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.17, 1.07, panel, transform=ax.transAxes,
                fontsize=7, fontweight="bold")

    output = ROOT / "figures" / "FigS6.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="pdf", facecolor="white", dpi=300)
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


if __name__ == "__main__":
    run()
