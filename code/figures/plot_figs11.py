"""Plot FigS11 from the compact cohort lambda table."""

import importlib.util
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA = ROOT / "data" / "FigS11.csv"
OUTPUT = ROOT / "figures" / "FigS11.pdf"
COLUMNS = ["Cohort", "Method", "PC", "Lambda"]


def load_source():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    path = CODE / "plot_R1_7_SLE_PC_lambda_by_cohort.py"
    spec = importlib.util.spec_from_file_location("figs11_plot_source", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run():
    if not DATA.is_file():
        raise FileNotFoundError("FigS11 data is missing; prepare it first")

    data = pd.read_csv(DATA, dtype={"Cohort": str, "Method": str})
    if list(data.columns) != COLUMNS:
        raise ValueError("FigS11.csv has unexpected columns")
    data = data.rename(columns={
        "Cohort": "cohort",
        "Method": "method",
        "PC": "pc",
        "Lambda": "lambda_gc",
    })

    source = load_source()
    base = source._lambda
    base._configure_style()
    cohorts = list(dict.fromkeys(data["cohort"].astype(str)))
    if not cohorts:
        raise ValueError("FigS11 has no cohort data")

    ncols = 3
    nrows = (len(cohorts) + ncols - 1) // ncols
    figure = plt.figure(
        figsize=(126 / 25.4, (nrows * 42 + 15) / 25.4),
        facecolor="white",
    )
    grid = figure.add_gridspec(
        nrows=nrows, ncols=ncols, hspace=0.58, wspace=0.22,
        left=0.10, right=0.99, top=0.98, bottom=0.17,
    )
    axes = []
    for index, cohort in enumerate(cohorts):
        row, col = divmod(index, ncols)
        ax = figure.add_subplot(grid[row, col])
        axes.append(ax)
        subset = data.loc[data["cohort"].astype(str).eq(cohort)]
        ticks = source._panel_y_ticks(float(subset["lambda_gc"].max()))
        minimum = float(subset["lambda_gc"].min())
        ymin = min(1.0, max(0.0, minimum - max(.05, .05 * max(1.0, abs(minimum)))))
        if set(subset["method"].astype(str)) != set(base.METHODS):
            raise ValueError("Incomplete method grid for {}".format(cohort))

        for method in base.METHODS:
            curve = subset.loc[subset["method"].eq(method)].sort_values("pc")
            if tuple(curve["pc"].astype(int)) != base.PC_LEVELS:
                raise ValueError("Incomplete PC curve for {}/{}".format(cohort, method))
            ax.plot(
                curve["pc"].to_numpy(float), curve["lambda_gc"].to_numpy(float),
                color=base.COLORS[method], marker="o", markersize=2.0,
                linewidth=0.85, label=base.LABELS[method],
            )
            below = curve.loc[curve["lambda_gc"].lt(1.3)]
            choice = (
                below.iloc[0] if not below.empty
                else curve.sort_values(["lambda_gc", "pc"], kind="stable").iloc[0]
            )
            ax.scatter(
                [choice["pc"]], [choice["lambda_gc"]], marker="*", s=30,
                facecolor=base.COLORS[method], edgecolor="black",
                linewidth=0.35, zorder=8, clip_on=True,
            )

        ax.set_title(cohort, fontsize=6.0, pad=1.8)
        ax.set_xlim(-1.5, 41.5)
        ax.set_ylim(ymin, float(ticks[-1]))
        ax.set_xticks(base.PC_LEVELS)
        ax.set_xticklabels([str(pc) for pc in base.PC_LEVELS])
        ax.set_yticks(ticks)
        ax.set_yticklabels(["{:.1f}".format(value) for value in ticks])
        ax.tick_params(direction="out", length=1.5, width=0.45, pad=0.8, labelsize=5.0)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.55)
            spine.set_color("black")
        if col == 0:
            ax.set_ylabel("Genomic inflation factor (λ)", fontsize=5.8, labelpad=1.0)
        if row == nrows - 1:
            ax.set_xlabel("Number of PCs", fontsize=5.7, labelpad=1.0)

    for index in range(len(cohorts), nrows * ncols):
        figure.add_subplot(grid[index // ncols, index % ncols]).axis("off")

    handles, labels = axes[0].get_legend_handles_labels()
    handles.append(Line2D(
        [0], [0], marker="*", linestyle="none", markerfacecolor="white",
        markeredgecolor="black", markeredgewidth=0.4, markersize=6,
    ))
    labels.append("selected PC")
    figure.legend(
        handles, labels, frameon=False, ncol=5, loc="lower center",
        bbox_to_anchor=(0.55, 0.018), fontsize=5.4, handlelength=0.95,
        columnspacing=0.7, handletextpad=0.25, borderaxespad=0.0,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", facecolor="white")
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None:
            display(figure)
    except ImportError:
        pass
    plt.close(figure)


if __name__ == "__main__":
    run()