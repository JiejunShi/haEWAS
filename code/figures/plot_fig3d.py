"""Plot Fig3d using the prepared paired-bootstrap effect differences."""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
EXTERNAL = ("GSE110607", "GSE118144", "GSE96879")
WIDTH, HEIGHT = 46.0, 66.0


def draw_ci(ax, y, low, high, color, width=.62):
    ax.plot([low, high], [y, y], color=color, linewidth=width)
    ax.plot([low, low], [y - .04, y + .04], color=color, linewidth=.42)
    ax.plot([high, high], [y - .04, y + .04], color=color, linewidth=.42)


def run():
    data = pd.read_csv(ROOT / "data" / "Fig3d.csv", float_precision="round_trip")
    labels = [*EXTERNAL, "Mean across cohorts"]
    if (list(data.columns) != ["Cohort", "Delta", "LowerCI", "UpperCI"]
            or len(data) != 4 or set(data["Cohort"]) != set(labels)):
        raise ValueError("Fig3d.csv must contain three external cohorts and their equal-weight mean")
    data = data.set_index("Cohort").loc[labels].reset_index()
    values = data[["Delta", "LowerCI", "UpperCI"]].to_numpy(float)
    if not np.isfinite(values).all() or (data["LowerCI"] > data["UpperCI"]).any():
        raise ValueError("Fig3d.csv contains invalid effect differences or intervals")
    font_file = ROOT.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
    font_family = ["Arial", "DejaVu Sans"]
    if font_file.is_file():
        fm.fontManager.addfont(str(font_file))
        font_family = fm.FontProperties(fname=str(font_file)).get_name()
    plt.rcParams.update({
        "font.family": font_family, "font.size": 5.5,
        "axes.labelsize": 6.0, "axes.titlesize": 6.0, "xtick.labelsize": 5.0,
        "ytick.labelsize": 5.0, "legend.fontsize": 5.0, "axes.linewidth": .55,
        "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 300,
        "savefig.facecolor": "white", "axes.unicode_minus": False,
    })
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
    ax = fig.add_axes([15 / WIDTH, 10.5 / HEIGHT, 27 / WIDTH, 46 / HEIGHT])
    for y, row in zip((3, 2, 1), data.iloc[:3].itertuples(index=False)):
        draw_ci(ax, y, row.LowerCI, row.UpperCI, "#6A6A6A")
        ax.scatter(row.Delta, y, s=11, marker="s", facecolor="#B64342",
                   edgecolor="#272727", linewidth=.25, zorder=3)
    pooled = data.iloc[3]
    draw_ci(ax, 0, pooled["LowerCI"], pooled["UpperCI"], "#272727", width=.72)
    ax.scatter(pooled["Delta"], 0, s=10, marker="D", facecolor="#272727",
               edgecolor="#272727", linewidth=.25, zorder=4)
    ax.axhline(.5, color="#B8B8B8", linewidth=.35)
    ax.axvline(0, color="#858585", linewidth=.48, linestyle=(0, (2.1, 1.8)))
    ax.set_yticks([3, 2, 1, 0])
    ax.set_yticklabels(labels)
    limits = np.r_[data["LowerCI"].to_numpy(float), data["UpperCI"].to_numpy(float), 0]
    pad = max(float(np.ptp(limits)) * .13, .08)
    ax.set_xlim(float(np.min(limits) - pad), float(np.max(limits) + pad))
    ax.set_ylim(-.28, len(EXTERNAL) + .72)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
    ax.set_xlabel("Difference in Cohen's $d$\n(haEWAS $-$ EWAS)", labelpad=2)
    ax.tick_params(axis="both", length=1.8, width=.45, pad=1.4, labelsize=5.0, colors="#272727")
    for spine in ax.spines.values():
        spine.set_linewidth(.5)
        spine.set_color("#272727")
    ax.grid(False)
    output = ROOT / "figures" / "Fig3d.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="pdf", facecolor="white")
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    run()
