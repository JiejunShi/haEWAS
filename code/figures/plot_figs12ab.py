"""Draw FigS12a-b together from their compact data tables."""

import copy
import importlib.util
import itertools
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA = ROOT / "data"
OUTPUT = ROOT / "figures" / "FigS12ab.pdf"

B_WIDTH_MM = 129.0
B_HEIGHT_MM = 132.0
GAP_MM = 6.0
A_SECTION_MM = 43.0
FIG_HEIGHT_MM = B_HEIGHT_MM + GAP_MM + A_SECTION_MM
GRID_LEFT, GRID_RIGHT, GRID_BOTTOM, GRID_TOP = 0.12, 0.83, 0.25, 0.91
GRID_WSPACE, GRID_HSPACE = 0.28, 0.37
AX_WIDTH = (GRID_RIGHT - GRID_LEFT) / (2 + GRID_WSPACE)
AX_HEIGHT = (GRID_TOP - GRID_BOTTOM) / (2 + GRID_HSPACE)
A_AX_HEIGHT_MM = 31.0
A_AX_BOTTOM_MM = B_HEIGHT_MM + GAP_MM + 6.0
A_METHODS = ("m", "m2", "m3", "chalm_camda")
FAMILIES = ("m", "m2", "m3", "haEWAS")


def source_module(filename, name):
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    spec = importlib.util.spec_from_file_location(name, CODE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def draw_a(figure, data, source):
    source._configure_style()
    source.METHOD_ORDER = A_METHODS
    source.METHOD_LABELS["chalm_camda"] = "hm"

    ax = figure.add_axes([
        GRID_LEFT,
        A_AX_BOTTOM_MM / FIG_HEIGHT_MM,
        AX_WIDTH,
        A_AX_HEIGHT_MM / FIG_HEIGHT_MM,
    ])
    for method in A_METHODS:
        curve = data.loc[data["Method"].eq(method)].sort_values("PC")
        ax.plot(
            curve["PC"].to_numpy(float),
            curve["Significant CpGs"].to_numpy(float),
            color=source.COLORS[method], marker="o", markersize=2.1,
            linewidth=0.8, label=source.METHOD_LABELS[method],
        )
    ax.set_xticks(source.PC_LEVELS)
    ax.set_xlim(-1, 41)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4, integer=True))
    ax.set_xlabel("Number of PCs", labelpad=1.5)
    ax.set_ylabel("Significant CpGs", labelpad=1.7)
    ax.set_title("SLE whole blood", pad=2.5)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.55)
        spine.set_color("black")
    ax.tick_params(direction="out", length=1.7, width=0.45, pad=1.3)
    ax.legend(
        frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.30),
        handlelength=1.0, columnspacing=0.8, handletextpad=0.3,
        labelspacing=0.3, borderaxespad=0.0,
    )
    figure.text(
        0.012, (FIG_HEIGHT_MM - 1.5) / FIG_HEIGHT_MM, "a",
        ha="left", va="top", fontsize=8, fontweight="bold",
    )


def draw_b(figure, data, source):
    colors = {
        "m": "#0072B2", "m2": "#D55E00", "m3": "#CC79A7",
        "haEWAS": "#8E6C9E",
    }
    method_labels = {"m": "m", "m2": "m²", "m3": "m³", "haEWAS": "hm"}
    codes = {}
    for cohort, code in data[["Discovery", "DiscoveryCode"]].drop_duplicates().itertuples(index=False):
        codes[str(cohort)] = int(code)
    for cohort, code in data[["Validation", "ValidationCode"]].drop_duplicates().itertuples(index=False):
        cohort, code = str(cohort), int(code)
        if cohort in codes and codes[cohort] != code:
            raise ValueError("Inconsistent FigS12b cohort code")
        codes[cohort] = code
    cohorts = [name for name, _ in sorted(codes.items(), key=lambda item: item[1])]
    pairs = list(itertools.permutations(cohorts, 2))
    labels = ["{}→{}".format(codes[a], codes[b]) for a, b in pairs]
    bad_color = "#E5E7EB"
    cmap = copy.copy(plt.get_cmap("RdBu_r"))
    cmap.set_bad(bad_color)
    norm = Normalize(vmin=-1, vmax=1)

    style = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.2, "axes.labelsize": 6.2, "axes.titlesize": 7.2,
        "xtick.labelsize": 5.5, "ytick.labelsize": 4.7,
        "axes.linewidth": 0.5, "pdf.fonttype": 42, "svg.fonttype": "none",
    }
    ax_width = AX_WIDTH
    ax_height = AX_HEIGHT
    col_gap = GRID_WSPACE * ax_width
    row_gap = GRID_HSPACE * ax_height
    b_offset = 0.0

    with mpl.rc_context(style):
        for index, family in enumerate(FAMILIES):
            row, col = divmod(index, 2)
            x = GRID_LEFT + col * (ax_width + col_gap)
            y = GRID_BOTTOM + (1 - row) * (ax_height + row_gap)
            ax = figure.add_axes([
                x,
                (b_offset + y * B_HEIGHT_MM) / FIG_HEIGHT_MM,
                ax_width,
                ax_height * B_HEIGHT_MM / FIG_HEIGHT_MM,
            ])
            block = data.loc[data["Family"].eq(family)]
            matrix = np.full((len(pairs), len(source.PC_LEVELS)), np.nan)
            for i, pc in enumerate(source.PC_LEVELS):
                for j, (discovery, validation) in enumerate(pairs):
                    selected = block.loc[
                        block["PC"].eq(pc)
                        & block["Discovery"].eq(discovery)
                        & block["Validation"].eq(validation)
                    ]
                    if len(selected) != 1:
                        raise ValueError("Missing or duplicate FigS12b heatmap cell")
                    value = selected["PearsonR"].iloc[0]
                    if pd.notna(value):
                        matrix[j, i] = float(value)

            image = ax.imshow(
                np.ma.masked_invalid(matrix), cmap=cmap, norm=norm,
                interpolation="none", aspect="auto",
            )
            for pair_index in range(len(pairs)):
                for pc_index in range(len(source.PC_LEVELS)):
                    if not np.isfinite(matrix[pair_index, pc_index]):
                        ax.text(
                            pc_index, pair_index, "×", ha="center", va="center",
                            color="#687078", fontsize=5.2, zorder=3,
                        )
            ax.set_title(method_labels[family], pad=3)
            ax.set_xticks(np.arange(len(source.PC_LEVELS)))
            ax.set_xticklabels([str(pc) for pc in source.PC_LEVELS])
            ax.set_yticks(np.arange(len(pairs)))
            ax.set_yticklabels(labels)
            if col == 0:
                ax.set_ylabel("Discovery → validation cohort")
            else:
                ax.tick_params(axis="y", labelleft=False)
            if row == 1:
                ax.set_xlabel("Number of PCs")
            else:
                ax.tick_params(axis="x", labelbottom=False)
            ax.set_xticks(np.arange(-0.5, len(source.PC_LEVELS), 1), minor=True)
            ax.set_yticks(np.arange(-0.5, len(pairs), 1), minor=True)
            ax.grid(which="minor", color="white", linewidth=0.6)
            ax.tick_params(which="major", length=0, pad=1.5)
            ax.tick_params(which="minor", bottom=False, left=False)
            ax.set_xlim(-0.5, len(source.PC_LEVELS) - 0.5)
            ax.set_ylim(len(pairs) - 0.5, -0.5)
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.45)
                spine.set_color("#C7CCD1")

        scalar = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        scalar.set_array([])
        colorbar_height = 0.38 / 3
        cax = figure.add_axes([
            0.865,
            (b_offset + (0.57 - colorbar_height / 2) * B_HEIGHT_MM) / FIG_HEIGHT_MM,
            0.018,
            colorbar_height * B_HEIGHT_MM / FIG_HEIGHT_MM,
        ])
        colorbar = figure.colorbar(
            scalar, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1]
        )
        colorbar.set_label("Pearson's r", labelpad=3, fontsize=5.8)
        colorbar.ax.tick_params(labelsize=5, length=1.5, width=0.45, pad=1)
        colorbar.outline.set_linewidth(0.45)
        figure.legend(
            handles=[Patch(
                facecolor=bad_color, edgecolor="#D0D4D8", label="×  Undefined r"
            )],
            loc="center left",
            bbox_to_anchor=(0.84, (b_offset + 0.44 * B_HEIGHT_MM) / FIG_HEIGHT_MM),
            frameon=False, fontsize=5, handlelength=0.9, handleheight=0.7,
            handletextpad=0.35,
        )
        key = "; ".join("{} = {}".format(codes[c], c) for c in cohorts)
        figure.text(
            0.5, (b_offset + 0.075 * B_HEIGHT_MM) / FIG_HEIGHT_MM,
            key, ha="center", va="bottom", fontsize=5.4,
        )
        figure.text(
            0.012, (b_offset + 0.985 * B_HEIGHT_MM) / FIG_HEIGHT_MM, "b",
            ha="left", va="top", fontsize=8, fontweight="bold",
        )


def run():
    data_a = pd.read_csv(DATA / "FigS12a.csv", dtype={"Method": str})
    data_b = pd.read_csv(
        DATA / "FigS12b.csv",
        dtype={"Discovery": str, "Validation": str, "Family": str},
    )
    expected_a = ["PC", "Method", "Significant CpGs"]
    expected_b = [
        "PC", "Discovery", "Validation", "Family", "PearsonR",
        "DiscoveryCode", "ValidationCode",
    ]
    if list(data_a.columns) != expected_a or list(data_b.columns) != expected_b:
        raise ValueError("FigS12 input data has unexpected columns")

    source_a = source_module(
        "plot_R1_7_SLE_PC_meta_counts.py", "figs12a_plot_source"
    )
    source_b = source_module(
        "plot_R1_7_SLE_PC_pairwise_effect_correlation.py", "figs12b_plot_source"
    )
    figure = plt.figure(
        figsize=(B_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4), facecolor="white"
    )
    draw_a(figure, data_a, source_a)
    draw_b(figure, data_b, source_b)

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
