"""Draw FigS16a and b-d together from compact input tables."""

import contextlib
import importlib
import importlib.util
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code" / "sle_mrs_four_methods"
DATA = ROOT / "data"
OUTPUT = ROOT / "figures" / "FigS16.pdf"
A_WIDTH_MM, A_HEIGHT_MM = 120.0, 117.6
ROC_HEIGHT_MM = 54.0
ROW_GAP_MM = 8.0
FIG_WIDTH_MM = 120.0
FIG_HEIGHT_MM = A_HEIGHT_MM + ROW_GAP_MM + ROC_HEIGHT_MM
A_BOTTOM_MM = ROC_HEIGHT_MM + ROW_GAP_MM
ROC_PANEL_GAP_MM = (FIG_WIDTH_MM - 3 * 39.0) / 2
ROC_X = (0.0, 39.0 + ROC_PANEL_GAP_MM, 2 * (39.0 + ROC_PANEL_GAP_MM))
A_COLUMNS = [
    "Cohort", "Group", "EWAS_MRS_Z", "CHALM_MRS_Z",
    "CAMDA_MRS_Z", "haEWAS_MRS_Z",
]


def reload_density_data():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    importlib.invalidate_caches()
    module = sys.modules.get("four_method_density_data")
    if module is None:
        return importlib.import_module("four_method_density_data")
    return importlib.reload(module)

def source_module(filename, name):
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    spec = importlib.util.spec_from_file_location(name, CODE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def data_file(stem):
    path = DATA / (stem + ".csv")
    if path.is_file():
        return path
    compressed = Path(str(path) + ".gz")
    if compressed.is_file():
        return compressed
    raise FileNotFoundError("{} data is missing; prepare it first".format(stem))


def add_rect_a(figure, source, x, y, width, height):
    return figure.add_axes([
        x / A_WIDTH_MM, (A_BOTTOM_MM + y) / FIG_HEIGHT_MM,
        width / A_WIDTH_MM, height / FIG_HEIGHT_MM,
    ])


def draw_density(figure, data, source):
    source.style()
    left, bottom, row_h, col_gap, row_gap = 10.5, 9.2, 22.0, 1.5, 1.6
    strip_w, strip_gap, top_h = 7.2, 0.7, 4.3
    grid_w = A_WIDTH_MM - left - 1.5 - strip_w - strip_gap
    col_w = (grid_w - 2 * col_gap) / 3
    top = bottom + 4 * row_h + 3 * row_gap + 0.7

    def add_strip(x, y, w, h, color, text, rotation=0):
        ax = add_rect_a(figure, source, x, y, w, h)
        ax.set_facecolor(color)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.35)
        ax.text(
            0.5, 0.5, text, ha="center", va="center", rotation=rotation,
            fontsize=6.0, color="white" if rotation else "#272727",
            transform=ax.transAxes,
        )

    for j, cohort in enumerate(source.EXTERNAL):
        add_strip(
            left + j * (col_w + col_gap), top, col_w, top_h,
            "#DCDCDC", cohort,
        )

    display_names = {
        "EWAS": "EWAS", "CHALM": "CHALM",
        "CAMDA": "CAMDA", "haEWAS": "haEWAS",
    }
    for i, method in enumerate(source.METHODS):
        y = bottom + (3 - i) * (row_h + row_gap)
        add_strip(
            left + grid_w + strip_gap, y, strip_w, row_h,
            source.METHOD_COLORS[method],
            display_names[method],
            -90,
        )
        for j, cohort in enumerate(source.EXTERNAL):
            subset = data.loc[data["Cohort"].eq(cohort)]
            control = subset.loc[
                subset["Group"].eq("Normal"), method + "_MRS_Z"
            ].to_numpy(float)
            case = subset.loc[
                subset["Group"].eq("SLE"), method + "_MRS_Z"
            ].to_numpy(float)
            grid, curves, overlap, p_value = source.density_pair(control, case)
            ax = add_rect_a(
                figure, source, left + j * (col_w + col_gap), y, col_w, row_h
            )
            for group, values in (("Normal", control), ("SLE", case)):
                color = source.DISEASE_COLORS[group]
                ax.fill_between(
                    grid, curves[group], color=color, alpha=0.52, linewidth=0
                )
                ax.plot(grid, curves[group], color=color, lw=0.75)
                ax.axvline(
                    np.median(values), color=color, ls=(0, (3, 2)), lw=0.65
                )
            p_text = (
                "P = {:.1e}".format(p_value) if p_value < 0.001
                else "P = {:.3f}".format(p_value)
            )
            ax.text(
                0.97, 0.96, "Overlap = {:.1f}%\n{}".format(
                    overlap * 100, p_text
                ),
                ha="right", va="top", transform=ax.transAxes,
                fontsize=6.0, linespacing=1.0,
            )
            ax.set_xlim(grid[0], grid[-1])
            ax.set_ylim(0, 1.04)
            ax.set_yticks([0, 0.5, 1])
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(
                labelsize=6.0, length=1.5, width=0.4, pad=1.0,
                labelbottom=(i == 3), labelleft=(j == 0),
            )
            for spine in ax.spines.values():
                spine.set_visible(True)

    figure.text(
        0.5, (A_BOTTOM_MM + 1.7) / FIG_HEIGHT_MM,
        "MRS (Z-normalized)", ha="center", va="bottom", fontsize=6.5,
    )
    figure.text(
        1.6 / FIG_WIDTH_MM,
        (A_BOTTOM_MM + A_HEIGHT_MM / 2) / FIG_HEIGHT_MM,
        "Relative density", rotation=90, ha="center", va="center",
        fontsize=6.5,
    )
    handles = [
        Line2D([0], [0], color=color, lw=1.7, label=label)
        for label, color in source.DISEASE_COLORS.items()
    ]
    figure.legend(
        handles=handles, loc="upper center",
        bbox_to_anchor=(0.5, (A_BOTTOM_MM + 0.992 * A_HEIGHT_MM)
                        / FIG_HEIGHT_MM),
        ncol=2, frameon=False, fontsize=6.5, handlelength=1.2,
        handletextpad=0.35, columnspacing=1.0, borderaxespad=0,
    )
    figure.text(
        1.6 / FIG_WIDTH_MM,
        (A_BOTTOM_MM + A_HEIGHT_MM - 1.5) / FIG_HEIGHT_MM,
        "a", fontsize=7.0, fontweight="bold", va="top",
    )


def draw_roc(figure, data, source):
    source.style()
    for i, cohort in enumerate(source.EXTERNAL):
        panel_x = ROC_X[i]
        ax = figure.add_axes([
            (panel_x + 7.0) / FIG_WIDTH_MM,
            8.1 / FIG_HEIGHT_MM,
            30.8 / FIG_WIDTH_MM,
            30.8 / FIG_HEIGHT_MM,
        ])
        ax.plot(
            [0, 1], [0, 1], linestyle="--", color="#9D9D9D",
            linewidth=0.4, zorder=1,
        )
        handles = []
        for method in source.PLOT_METHODS:
            curve = data.loc[
                data["Cohort"].eq(cohort) & data["Method"].eq(method)
            ].sort_values(["FPR", "TPR"])
            if curve.empty:
                raise ValueError("Missing ROC data for {}/{}".format(cohort, method))
            auc = float(curve["AUC"].iloc[0])
            ax.plot(
                curve["FPR"].to_numpy(float),
                curve["TPR"].to_numpy(float),
                color=source.COLORS[method], linewidth=1.0,
                solid_capstyle="round", zorder=2,
            )
            handles.append(Line2D(
                [0], [0], color=source.COLORS[method], linewidth=1.0,
                label="{} (AUC={:.3f})".format(source.LABELS[method], auc),
            ))
        ax.set_xlim(-0.025, 1.025)
        ax.set_ylim(-0.025, 1.025)
        ax.set_aspect("equal", adjustable="box")
        ticks = [0, 0.25, 0.50, 0.75, 1.0]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1"])
        ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1"])
        ax.set_xlabel("1 - specificity", labelpad=1.5)
        ax.set_ylabel("Sensitivity", labelpad=1.5)
        ax.legend(
            handles=handles, loc="lower right", bbox_to_anchor=(0.98, 0.04),
            frameon=False, fontsize=4.8, handlelength=1.0,
            handletextpad=0.3, borderaxespad=0.0, labelspacing=0.3,
        )
        ax.tick_params(
            axis="both", which="major", length=1.4, width=0.35,
            pad=1.0, labelsize=5.0, colors="#272727",
        )
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.4)
            spine.set_color("#272727")

        strip = figure.add_axes([
            (panel_x + 7.0) / FIG_WIDTH_MM,
            (8.1 + 30.8 + 0.22) / FIG_HEIGHT_MM,
            30.8 / FIG_WIDTH_MM,
            3.0 / FIG_HEIGHT_MM,
        ])
        strip.set_facecolor("#DCDCDC")
        strip.set_xticks([])
        strip.set_yticks([])
        for spine in strip.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.35)
            spine.set_color("#272727")
        strip.text(
            0.5, 0.5, cohort, ha="center", va="center", fontsize=6.0,
            fontweight="bold", color="#272727", transform=strip.transAxes,
        )
        figure.text(
            (panel_x + 1.5) / FIG_WIDTH_MM,
            (8.1 + 30.8 + 0.22 + 3.0 + 4.0) / FIG_HEIGHT_MM,
            "bcd"[i], fontsize=7.0, fontweight="bold", va="bottom",
        )


def run():
    reload_density_data()
    density = source_module(
        "plot_sle_four_method_mrs_density.py", "figs16_density_plot"
    )
    auc = source_module("plot_sle_four_method_mrs_auc.py", "figs16_auc_plot")
    scores = pd.read_csv(data_file("FigS16a"))
    roc = pd.read_csv(
        data_file("FigS16b-d"),
        dtype={"Cohort": str, "Method": str},
    )
    if list(scores.columns) != A_COLUMNS:
        raise ValueError("FigS16a data has unexpected columns")
    if list(roc.columns) != ["Cohort", "Method", "FPR", "TPR", "AUC"]:
        raise ValueError("FigS16b-d data has unexpected columns")

    density.style()
    auc.style()
    figure = plt.figure(
        figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4),
        facecolor="white",
    )
    draw_density(figure, scores, density)
    draw_roc(figure, roc, auc)

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