"""Plot FigS3a-d from four panel-specific CSV files."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter


ROOT = Path(os.environ["REPRO_ROOT"])
TARGETS = ("CAMDA", "CHALM")
SUPPORT = ((50, 75), (75, 100), (100, 150), (150, 200), (200, 206))
BLUE = "#6C8EAD"
FILES = {
    "CAMDA": ("FigS3a.csv.gz", "FigS3b.csv.gz"),
    "CHALM": ("FigS3c.csv.gz", "FigS3d.csv.gz"),
}


def box_stats(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        raise ValueError("No finite values for the FigS3 box summary")
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    iqr = q3 - q1
    within = values[(values >= q1 - 1.5 * iqr) & (values <= q3 + 1.5 * iqr)]
    return {
        "q1": q1,
        "med": median,
        "q3": q3,
        "whislo": float(within.min()),
        "whishi": float(within.max()),
    }


def _draw(data):
    fig, axes = plt.subplots(
        2, 2, figsize=(129 / 25.4, 128 / 25.4),
        gridspec_kw={"width_ratios": [1, 2.4]},
    )
    fig.subplots_adjust(
        left=0.13, right=0.98, bottom=0.13, top=0.88,
        hspace=0.68, wspace=0.40,
    )

    for row, target in enumerate(TARGETS):
        hist_data, support_data = data[target]
        rmse = hist_data["RMSE"].to_numpy(dtype=float)
        edge_max = max(0.3, float(np.ceil(rmse.max() / 0.0025) * 0.0025))
        edges = np.linspace(0, edge_max, int(round(edge_max / 0.0025)) + 1)
        counts = np.histogram(rmse, bins=edges)[0]

        ax = axes[row, 0]
        ax.bar(edges[:-1], counts / 1000.0, width=np.diff(edges), align="edge",
               color=BLUE, edgecolor="none")
        ax.set(
            xlim=(0, 0.3), ylim=(0, 32),
            xlabel="Per-CpG RMSE", ylabel="CpGs (\u00d710\u00b3)",
            title="hm({})".format(target),
        )
        ax.set_xticks([0, 0.1, 0.2, 0.3])
        ax.set_yticks([0, 8, 16, 24, 32])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: "{:g}".format(value)))
        stats = box_stats(rmse)
        ax.axvline(stats["med"], color="#244960", linewidth=0.8,
                   linestyle="--", zorder=3)
        ax.text(
            0.97, 0.97,
            "N={:,.0f}\nMedian={:.3f}\nIQR: {:.3f}\u2013{:.3f}".format(
                len(rmse), stats["med"], stats["q1"], stats["q3"]),
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5,
        )

        ax = axes[row, 1]
        records, labels = [], []
        for lower, upper in SUPPORT:
            selected = support_data.loc[
                support_data["ReferenceN"].ge(lower)
                & support_data["ReferenceN"].lt(upper), "RMSE"
            ].to_numpy(dtype=float)
            if not len(selected):
                continue
            records.append(box_stats(selected))
            labels.append("{}\u2013{}\n({:,})".format(lower, upper - 1, len(selected)))
        if not records:
            raise ValueError("No reference-support groups for {}".format(target))
        ax.bxp(
            records, positions=range(len(records)), showfliers=False,
            patch_artist=True, widths=0.55,
            boxprops={"facecolor": "#B9CDD9", "edgecolor": "#335C75", "linewidth": 0.6},
            medianprops={"color": "#244960", "linewidth": 0.9},
            whiskerprops={"color": "#335C75", "linewidth": 0.6},
            capprops={"color": "#335C75", "linewidth": 0.6},
        )
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=5.5, rotation=0, ha="center", va="top")
        ax.set_xlim(-0.5, len(labels) - 0.5)
        ax.set_ylim(-0.018, 0.3)
        ax.set_yticks([0, 0.1, 0.2, 0.3])
        ax.set_xlabel("Reference samples")
        ax.set_ylabel("Per-CpG RMSE")
        ax.set_title("hm({})".format(target))
        ax.xaxis.labelpad = 5

        for col in range(2):
            axis = axes[row, col]
            axis.text(-0.23, 1.08, "abcd"[row * 2 + col],
                      transform=axis.transAxes, fontsize=7, fontweight="bold")
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.tick_params(length=2, width=0.5, pad=2)

    return fig


def run():
    data = {}
    for target, filenames in FILES.items():
        hist_data = pd.read_csv(ROOT / "data" / filenames[0])
        support_data = pd.read_csv(ROOT / "data" / filenames[1])
        if list(hist_data.columns) != ["CpG", "RMSE"]:
            raise ValueError("{} has unexpected columns".format(filenames[0]))
        if list(support_data.columns) != ["CpG", "RMSE", "ReferenceN"]:
            raise ValueError("{} has unexpected columns".format(filenames[1]))
        if not hist_data["CpG"].equals(support_data["CpG"]):
            raise ValueError("{} and {} contain different CpGs".format(*filenames))
        if hist_data.empty or support_data.empty:
            raise ValueError("FigS3 panel data is empty")
        if (not np.isfinite(hist_data["RMSE"].to_numpy(dtype=float)).all()
                or hist_data["RMSE"].lt(0).any()
                or not np.isfinite(support_data[["RMSE", "ReferenceN"]].to_numpy(dtype=float)).all()
                or support_data["RMSE"].lt(0).any()):
            raise ValueError("Invalid FigS3 panel data")
        data[target] = (hist_data, support_data)

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
    fig = _draw(data)
    output = ROOT / "figures" / "FigS3.pdf"
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
