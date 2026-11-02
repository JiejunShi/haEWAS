"""Plot FigS4a-b from compact CpG-level tables."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator


ROOT = Path(os.environ["REPRO_ROOT"])
TARGETS = ("CAMDA", "CHALM")
PANELS = (
    ("a", "FigS4a.csv.gz", "Depth", (20, 10, 5), "Reference depth (\u00d7)"),
    ("b", "FigS4b.csv.gz", "ReadLength", (125, 100, 75), "Reference read length (nt)"),
)
BLUE = "#DCE6EC"
EDGE = "#335C75"


def box_stats(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        raise ValueError("No finite FigS4 values")
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    iqr = q3 - q1
    inside = values[(values >= q1 - 1.5 * iqr) & (values <= q3 + 1.5 * iqr)]
    return {
        "q1": q1,
        "med": median,
        "q3": q3,
        "whislo": float(inside.min()),
        "whishi": float(inside.max()),
    }


def format_delta(value, _):
    return "{:g}".format(value)


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

    fig = plt.figure(figsize=(129 / 25.4, 50 / 25.4))
    grid = fig.add_gridspec(
        1, 2, left=0.09, right=0.99, bottom=0.30, top=0.80,
        wspace=0.75,
    )

    for index, (panel, filename, name, levels, xlabel) in enumerate(PANELS):
        data = pd.read_csv(ROOT / "data" / filename)
        value_columns = ["{}{}".format(name, level) for level in levels]
        if list(data.columns) != ["Target", "CpG", "Locus"] + value_columns:
            raise ValueError("{} has unexpected columns".format(filename))
        if (data.empty or not data["Target"].isin(TARGETS).all()
                or data[["Target", "CpG", "Locus"]].isna().any().any()
                or data.duplicated(["Target", "CpG"]).any()):
            raise ValueError("{} has invalid CpG rows".format(filename))

        boxes_by_target = {}
        for target in TARGETS:
            target_data = data.loc[data["Target"].eq(target)]
            boxes_by_target[target] = []
            for column in value_columns:
                values = pd.to_numeric(target_data[column], errors="coerce").dropna().to_numpy(float)
                if not len(values) or not np.isfinite(values).all():
                    raise ValueError("Invalid {} data for {}".format(column, target))
                boxes_by_target[target].append(box_stats(values))

        all_boxes = [box for target in TARGETS for box in boxes_by_target[target]]
        low = min([0.0] + [box["whislo"] for box in all_boxes])
        high = max([0.0] + [box["whishi"] for box in all_boxes])
        pad = max((high - low) * 0.12, 0.00005)

        inner = grid[0, index].subgridspec(1, 2, wspace=0.30)
        axes = [fig.add_subplot(inner[0, j]) for j in range(2)]
        for j, (ax, target) in enumerate(zip(axes, TARGETS)):
            ax.bxp(
                boxes_by_target[target], positions=range(len(levels)),
                widths=0.5, showfliers=False, patch_artist=True,
                boxprops={"facecolor": BLUE, "edgecolor": EDGE, "linewidth": 0.6},
                medianprops={"color": "#244960", "linewidth": 0.8},
                whiskerprops={"color": EDGE, "linewidth": 0.6},
                capprops={"color": EDGE, "linewidth": 0.6},
            )
            ax.axhline(0, color="#777777", linewidth=0.6, linestyle="--")
            ax.set(
                xlim=(-0.5, len(levels) - 0.5),
                ylim=(low - pad, high + pad),
                title="hm({})".format(target),
            )
            ax.set_xticks(range(len(levels)))
            ax.set_xticklabels([str(value) for value in levels])
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_formatter(FuncFormatter(format_delta))
            if j == 0:
                ax.set_ylabel("Per-CpG \u0394RMSE")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(length=2, width=0.5, pad=2)

        left = axes[0].get_position()
        right = axes[1].get_position()
        fig.text((left.x0 + right.x1) / 2, 0.15, xlabel,
                 ha="center", va="top", fontsize=6)
        axes[0].text(-0.22, 1.12, panel, transform=axes[0].transAxes,
                     fontsize=7, fontweight="bold")

    output = ROOT / "figures" / "FigS4.pdf"
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
