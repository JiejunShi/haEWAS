"""Plot the three cohort groups and connected GO panel for FigS15."""

import importlib.util
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter


REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
ROOT = Path(os.environ["REPRO_ROOT"])
MRS_CODE = REVISION_ROOT / "20_R2_7_SLE_EWAS_prediction_audit" / "code"
GO_CODE = REVISION_ROOT / "20_R2_7_SLE_EWAS_prediction_audit" / "code"
DATA = ROOT / "data"
OUTPUT = ROOT / "figures" / "FigS15.pdf"
FIG_W_MM, FIG_H_MM = 128.0, 112.0
GO_X, GO_Y, GO_W, GO_H = 71.0, 3.0, 49.0, 55.0
HELDOUTS = ("GSE110607", "GSE118144", "GSE96879")
POINTS_COLUMNS = [
    "CpG", "Cohort", "Method", "DiscoveryEffect", "ValidationEffect"
]
GO_COLUMNS = ["Term", "EWAS -log10(FDR)", "haEWAS -log10(FDR)"]


def source_module(code, filename, name):
    if str(code) not in sys.path:
        sys.path.insert(0, str(code))
    spec = importlib.util.spec_from_file_location(name, code / filename)
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


def quadrant_stats(frame):
    x = frame["DiscoveryEffect"].to_numpy(float)
    y = frame["ValidationEffect"].to_numpy(float)
    masks = {
        "I": (x > 0) & (y > 0),
        "II": (x < 0) & (y > 0),
        "III": (x < 0) & (y < 0),
        "IV": (x > 0) & (y < 0),
    }
    counts = {q: int(mask.sum()) for q, mask in masks.items()}
    total = sum(counts.values())
    stats = {}
    for q, n in counts.items():
        stats["quadrant_{}_n".format(q)] = n
        stats["quadrant_{}_pct".format(q)] = (
            100.0 * n / total if total else np.nan
        )
    return stats, total


def draw_mrs(figure, points, source):
    values = points[["DiscoveryEffect", "ValidationEffect"]].to_numpy(float).ravel()
    limit = source.nice_limit(values)
    width, _, panels = source.layout()
    letters = dict(zip(HELDOUTS, "abc"))
    group_methods = {
        "EWAS": ("beta",),
        "haEWAS": ("chalm", "camda"),
    }

    for cohort, group, x, y, size in panels:
        y = 75.0 if cohort in HELDOUTS[:2] else 32.0
        frame = points.loc[
            points["Cohort"].eq(cohort)
            & points["Method"].isin(group_methods[group])
        ]
        stats, total = quadrant_stats(frame)
        ax = figure.add_axes([
            x / width, y / FIG_H_MM, size / width, size / FIG_H_MM,
        ])
        if len(frame):
            ax.scatter(
                frame["DiscoveryEffect"], frame["ValidationEffect"],
                s=7, c=source.COLORS[group], marker="o",
                edgecolors="white", linewidths=0.15, alpha=0.85, zorder=2,
            )
        ax.axhline(0, color="#A5A5A5", linewidth=0.4, zorder=0)
        ax.axvline(0, color="#A5A5A5", linewidth=0.4, zorder=0)
        ax.plot(
            [-limit, limit], [-limit, limit],
            color="#666666", linestyle="--", linewidth=0.4, zorder=0,
        )
        source.annotate_quadrants(ax, stats)
        ax.set(
            xlim=(-limit, limit), ylim=(-limit, limit),
            xlabel="Discovery effect", ylabel="Validation effect",
        )
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([-limit, 0, limit])
        ax.set_yticks([-limit, 0, limit])
        for axis in (ax.xaxis, ax.yaxis):
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-2, 2))
            axis.set_major_formatter(formatter)
            axis.get_offset_text().set_fontsize(6)
        ax.tick_params(length=1.5, width=0.4, pad=1)
        ax.xaxis.labelpad = 1.5
        ax.yaxis.labelpad = 1.5
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.45)
        ax.set_title(
            "{}\n{} (N={})".format(cohort, group, total),
            fontsize=6.5, pad=3, linespacing=1.1,
        )
        if group == "EWAS":
            letter_y = y + 33.5 if cohort in HELDOUTS[:2] else y + size + 6.5
            figure.text(
                x / width, letter_y / FIG_H_MM,
                letters[cohort], fontsize=7, fontweight="bold", va="bottom",
            )


def draw_go(figure, data, source):
    source._style()
    max_value = max(
        float(data["EWAS -log10(FDR)"].max()),
        float(data["haEWAS -log10(FDR)"].max()),
    )
    x_max = max(5.0, float(np.ceil(max_value / 5.0) * 5.0))
    tick_step = max(5.0, float(np.ceil(x_max / 4.0 / 5.0) * 5.0))
    y = np.arange(len(data), dtype=float)

    left, right, bottom, top = 0.59, 0.985, 0.14, 0.90
    ax = figure.add_axes([
        (GO_X + left * GO_W) / FIG_W_MM,
        (GO_Y + bottom * GO_H) / FIG_H_MM,
        (right - left) * GO_W / FIG_W_MM,
        (top - bottom) * GO_H / FIG_H_MM,
    ])
    ewas = data["EWAS -log10(FDR)"].to_numpy(float)
    haewas = data["haEWAS -log10(FDR)"].to_numpy(float)
    for i in range(len(data)):
        ax.annotate(
            "", xy=(haewas[i], y[i]), xytext=(ewas[i], y[i]),
            arrowprops={
                "arrowstyle": "-|>", "color": source.ARROW_COLOR,
                "linewidth": 0.65, "mutation_scale": 5.5,
                "shrinkA": 4.5, "shrinkB": 4.5,
            },
            zorder=1,
        )
    ax.scatter(
        ewas, y, s=source.POINT_AREA, color=source.METHOD_COLORS["EWAS"],
        edgecolors="white", linewidths=0.55, alpha=0.95, zorder=3,
    )
    ax.scatter(
        haewas, y, s=source.POINT_AREA, color=source.METHOD_COLORS["haEWAS"],
        edgecolors="white", linewidths=0.55, alpha=0.95, zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(
        [source._wrap_term(term) for term in data["Term"]], fontsize=6.0
    )
    ax.set_ylim(len(data) - 0.48, -0.52)
    ax.set_xlabel(r"$-\log_{10}$(FDR)", labelpad=3)
    ax.axvline(
        -np.log10(source.FDR_THRESHOLD), color="#A6A6A6",
        linestyle="--", linewidth=0.55, zorder=0,
    )
    ax.set_xlim(0, x_max)
    ticks = np.unique(np.r_[np.arange(0.0, x_max, tick_step), x_max])
    ax.set_xticks(ticks)
    ax.set_xticklabels([
        str(int(value)) if float(value).is_integer() else str(float(value))
        for value in ticks
    ])
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.45, linestyle=":")
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", length=1.4, width=0.45, pad=1.2, labelsize=6.0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.55)

    handles = [
        Line2D([], [], linestyle="none", marker="o",
               markerfacecolor=source.METHOD_COLORS[method],
               markeredgecolor="white", markeredgewidth=0.45,
               markersize=6.0, label=method)
        for method in ("EWAS", "haEWAS")
    ]
    figure.legend(
        handles=handles, loc="upper center",
        bbox_to_anchor=((GO_X + GO_W * 0.5) / FIG_W_MM,
                        (GO_Y + GO_H * 0.99) / FIG_H_MM),
        ncol=2, frameon=False, fontsize=6.0,
        handletextpad=0.35, columnspacing=1.0, handlelength=0.9,
    )
    figure.text(
        (GO_X + 0.5) / FIG_W_MM, (GO_Y + GO_H - 1.5) / FIG_H_MM,
        "d", fontsize=7, fontweight="bold", va="top",
    )


def run():
    points = pd.read_csv(
        data_file("FigS15abc"),
        dtype={"CpG": str, "Cohort": str, "Method": str},
    )
    go_data = pd.read_csv(data_file("FigS15d"), dtype={"Term": str})
    if list(points.columns) != POINTS_COLUMNS or list(go_data.columns) != GO_COLUMNS:
        raise ValueError("FigS15 data has unexpected columns")
    if points.duplicated(["Cohort", "Method", "CpG"]).any():
        raise ValueError("FigS15 has duplicate cohort/method/CpG rows")
    if points.empty or go_data.empty:
        raise ValueError("FigS15 data is empty")

    mrs_source = source_module(
        MRS_CODE, "plot_sle_frozen_mrs_effect_direction.py",
        "figs15_mrs_plot_source",
    )
    go_source = source_module(
        GO_CODE, "plot_sle_pooled_model_go_terms_connected.py",
        "figs15_go_plot_source",
    )
    go_source._style()
    mpl.rcParams.update({
        "font.size": 6.5, "axes.labelsize": 6.5, "axes.titlesize": 6.5,
        "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
        "pdf.fonttype": 42, "axes.unicode_minus": False,
    })
    figure = plt.figure(
        figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4), facecolor="white"
    )
    draw_mrs(figure, points, mrs_source)
    draw_go(figure, go_data, go_source)

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