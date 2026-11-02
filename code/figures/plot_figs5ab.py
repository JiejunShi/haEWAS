"""Plot FigS5a-b from the two prepared CpG-level CSV files."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
PANELS = (
    ("a", "FigS5a.csv", "hm(CAMDA)"),
    ("b", "FigS5b.csv", "hm(CHALM)"),
)
MODELS = ("Full205", "MinusBloodImmune161", "MinusEpithelial122")
COLORS = ("#335C75", "#8BA8BA", "#C4D5DF")
LABELS = ("Full\n205", "\u2212Blood/immune\n161", "\u2212Epithelial\n122")


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

    for ax, (panel, filename, title) in zip(axes, PANELS):
        data = pd.read_csv(ROOT / "data" / filename)
        if list(data.columns) != ["CpG", "Locus"] + list(MODELS):
            raise ValueError("{} has unexpected columns".format(filename))
        if data.empty or data[["CpG", "Locus"]].isna().any().any() or data["CpG"].duplicated().any():
            raise ValueError("{} has invalid CpG rows".format(filename))

        for index, (model, color) in enumerate(zip(MODELS, COLORS)):
            values = pd.to_numeric(data[model], errors="coerce").dropna().to_numpy(float)
            if not len(values) or not np.isfinite(values).all():
                raise ValueError("No finite values in {} for {}".format(filename, model))
            if len(values) > 1 and np.ptp(values) > 0:
                violin = ax.violinplot(
                    [values], positions=[index], widths=0.72, points=400,
                    bw_method="scott", showmeans=False, showmedians=False,
                    showextrema=False,
                )
                for body in violin["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor("none")
                    body.set_alpha(0.9)
            else:
                ax.plot([index - 0.25, index + 0.25], [values[0], values[0]],
                        color=color, linewidth=2)

            q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
            ax.plot([index, index], [q1, q3], color="#244960", linewidth=1.2, zorder=3)
            ax.scatter(index, median, facecolor="white", edgecolor="#244960",
                       s=12, linewidth=0.5, zorder=4)

        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(LABELS, rotation=30, ha="right", va="top",
                           rotation_mode="default", fontsize=5.5)
        ax.set(xlim=(-0.6, 2.6), ylim=(0, 0.3),
               ylabel="Per-CpG RMSE", title=title)
        ax.set_yticks(np.arange(0, 0.301, 0.1))
        ax.tick_params(length=2, width=0.5, pad=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.17, 1.07, panel, transform=ax.transAxes,
                fontsize=7, fontweight="bold")

    output = ROOT / "figures" / "FigS5.pdf"
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
