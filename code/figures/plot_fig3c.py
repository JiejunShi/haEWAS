"""Plot Fig3c using only the prepared KDE coordinates and annotations."""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
EXTERNAL = ("GSE110607", "GSE118144", "GSE96879")
METHODS = ("EWAS", "haEWAS")
METHOD_COLORS = {"EWAS": "#0F4D92", "haEWAS": "#B64342"}
COLORS = {"Normal": "#7EAFC5", "SLE": "#CB5F48"}
WIDTH, HEIGHT = 120.0, 70.4


def rect(left, bottom, width, height):
    return [left / WIDTH, bottom / HEIGHT, width / WIDTH, height / HEIGHT]


def run():
    data = pd.read_csv(ROOT / "data" / "Fig3c.csv", float_precision="round_trip")
    columns = ["Cohort", "Model", "X", "NormalDensity", "SLEDensity", "NormalMedian",
               "SLEMedian", "Overlap", "PValue"]
    if list(data.columns) != columns or set(zip(data["Cohort"], data["Model"])) != {
        (cohort, method) for cohort in EXTERNAL for method in METHODS
    }:
        raise ValueError("Fig3c.csv has an invalid schema or panel roster")
    font_file = ROOT.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
    font_family = ["Arial", "DejaVu Sans"]
    if font_file.is_file():
        fm.fontManager.addfont(str(font_file))
        font_family = fm.FontProperties(fname=str(font_file)).get_name()
    plt.rcParams.update({
        "font.family": font_family, "font.size": 5.5,
        "axes.labelsize": 6.0, "axes.titlesize": 6.0,
        "xtick.labelsize": 5.0, "ytick.labelsize": 5.0, "legend.fontsize": 5.0,
        "axes.linewidth": .55, "legend.frameon": False, "pdf.fonttype": 42,
        "ps.fonttype": 42, "figure.dpi": 300, "savefig.facecolor": "white",
        "axes.unicode_minus": False,
    })
    left, bottom, row_gap, col_gap, cell_h = 10.5, 9.2, 1.6, 1.5, 22.0
    plot_w = WIDTH - left - 1.5 - .7 - 7.2
    cell_w = (plot_w - 2 * col_gap) / len(EXTERNAL)
    rows = {method: bottom + (len(METHODS) - 1 - i) * (cell_h + row_gap)
            for i, method in enumerate(METHODS)}
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
    for col, cohort in enumerate(EXTERNAL):
        strip = fig.add_axes(rect(left + col * (cell_w + col_gap),
                                  bottom + 2 * cell_h + row_gap + .7, cell_w, 4.3))
        strip.set_facecolor("#DCDCDC")
        strip.set_xticks([])
        strip.set_yticks([])
        for spine in strip.spines.values():
            spine.set_linewidth(.35)
            spine.set_color("#272727")
        strip.text(.5, .5, cohort, ha="center", va="center", fontsize=5.4,
                   fontweight="bold", color="#272727", transform=strip.transAxes)
    for row, method in enumerate(METHODS):
        method_data = data.loc[data["Model"].eq(method)]
        strip = fig.add_axes(rect(left + plot_w + .7, rows[method], 7.2, cell_h))
        strip.set_facecolor(METHOD_COLORS[method])
        strip.set_xticks([])
        strip.set_yticks([])
        for spine in strip.spines.values():
            spine.set_linewidth(.35)
            spine.set_color("#272727")
        strip.text(.5, .5, method, ha="center", va="center",
                   rotation=-90, fontsize=5.1, fontweight="bold", color="white",
                   transform=strip.transAxes)
        for col, cohort in enumerate(EXTERNAL):
            panel = method_data.loc[method_data["Cohort"].eq(cohort)].sort_values("X")
            numeric = panel[["X", "NormalDensity", "SLEDensity", "NormalMedian", "SLEMedian"]].to_numpy(float)
            if len(panel) != 512 or not np.isfinite(numeric).all():
                raise ValueError("Fig3c.csv contains incomplete or nonfinite KDE coordinates")
            metadata = panel.iloc[0]
            ax = fig.add_axes(rect(left + col * (cell_w + col_gap), rows[method], cell_w, cell_h))
            x = panel["X"].to_numpy(dtype=float)
            for group in ("Normal", "SLE"):
                density = panel[group + "Density"].to_numpy(dtype=float)
                ax.fill_between(x, density, color=COLORS[group], alpha=.52, linewidth=0)
                ax.plot(x, density, color=COLORS[group], linewidth=.75)
                ax.axvline(float(metadata[group + "Median"]), color=COLORS[group],
                           linestyle=(0, (3.0, 2.0)), linewidth=.65, alpha=.9, zorder=4)
            overlap = metadata["Overlap"]
            p_value = metadata["PValue"]
            overlap_text = "NA" if not np.isfinite(overlap) else f"{100 * overlap:.1f}%"
            p_text = "P = NA" if not np.isfinite(p_value) else (
                f"P = {p_value:.1e}" if p_value < 1e-3 else f"P = {p_value:.3f}")
            ax.text(.96, .95, f"Overlap = {overlap_text}\n{p_text}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=5.0, color="#272727", linespacing=1.0)
            ax.set_xlim(float(panel["X"].iloc[0]), float(panel["X"].iloc[-1]))
            ax.set_ylim(0, 1.04)
            ax.set_yticks([0, .5, 1.0])
            if row != len(METHODS) - 1:
                ax.set_xticklabels([])
            if col != 0:
                ax.set_yticklabels([])
            for spine in ax.spines.values():
                spine.set_linewidth(.45)
                spine.set_color("#272727")
            ax.tick_params(axis="both", labelsize=5.0, length=1.5, width=.4, pad=1.0)
    fig.text(.50, 1.7 / HEIGHT, "MRS (Z-normalized)", ha="center", va="bottom", fontsize=6.0, color="#272727")
    fig.text(1.6 / WIDTH, .50, "Density", ha="center", va="center", rotation=90, fontsize=6.0, color="#272727")
    handles = [Line2D([0], [0], color=COLORS[group], linewidth=1.7, label=group)
               for group in ("Normal", "SLE")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .99), ncol=2,
               frameon=False, fontsize=5.1, handlelength=1.2, handletextpad=.35,
               columnspacing=1.0, borderaxespad=0)
    output = ROOT / "figures" / "Fig3c.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="pdf", facecolor="white")
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    run()

