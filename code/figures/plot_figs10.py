"""Draw the five FigS10 panels on one three-row page."""

import importlib.util
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.lines import Line2D


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA = ROOT / "data"
FIG_W_MM = 183.0


def panel_file(name):
    path = DATA / (name + ".csv")
    return path if path.is_file() else Path(str(path) + ".gz")


def source_module(filename, name):
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    spec = importlib.util.spec_from_file_location(name, CODE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def render(fig):
    fig.set_dpi(600)
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    image = np.asarray(canvas.buffer_rgba()).copy()
    plt.close(fig)
    return image


def upset_panel():
    source = source_module("plot_r3_7_sle_ad_upset.py", "figs10_upset_plot")
    counts = pd.read_csv(panel_file("FigS10a"))
    method_columns = ["EWAS", "haEWAS (CAMDA)", "haEWAS (CHALM)"]
    pattern_keys = {
        (1, 0, 0): "EWAS_only",
        (0, 1, 0): "CAMDA_only",
        (0, 0, 1): "CHALM_only",
        (1, 1, 0): "EWAS_CAMDA_only",
        (1, 0, 1): "EWAS_CHALM_only",
        (0, 1, 1): "CHALM_CAMDA_only",
        (1, 1, 1): "All_three",
    }
    set_sizes = {
        "EWAS set size": "n_significant_EWAS",
        "haEWAS (CAMDA) set size": "n_significant_CAMDA",
        "haEWAS (CHALM) set size": "n_significant_CHALM",
    }
    row = {"phenotype": "Systemic Lupus Erythematosus"}
    for _, item in counts.iterrows():
        measure = item["Measure"]
        if measure in set_sizes:
            column = set_sizes[measure]
        else:
            pattern = tuple(int(item[name]) for name in method_columns)
            if pattern not in pattern_keys:
                raise ValueError("FigS10a has an unknown method combination")
            column = pattern_keys[pattern]
        if column in row:
            raise ValueError("FigS10a has a repeated measure")
        row[column] = int(item["Count"])
    expected = set(pattern_keys.values()) | set(set_sizes.values())
    if expected.difference(row):
        raise ValueError("FigS10a is missing a measure")
    source.apply_style()
    fig = plt.figure(
        figsize=(82 / 25.4, 82 * 45 / 59 / 25.4),
        facecolor="white",
    )
    grid = fig.add_gridspec(1, 1, left=0.06, right=0.985, top=0.90, bottom=0.07)
    source._draw_panel(
        fig, grid[0, 0], row,
        source.INTERSECTION_Y_MAX["Systemic Lupus Erythematosus"],
    )
    fig.text(0.012, 0.985, "a", ha="left", va="top", fontsize=8,
             fontweight="bold", color=source.AXIS_COLOR)
    return render(fig)


def lambda_panel():
    source = source_module("plot_sle_figs4a_lambda.py", "figs10_lambda_plot")
    data = pd.read_csv(panel_file("FigS10b"))
    source.style()
    mpl.rcParams.update({
        "font.size": 5.5, "axes.labelsize": source.FS_AXIS,
        "axes.titlesize": source.FS_AXIS, "xtick.labelsize": source.FS_TICK,
        "ytick.labelsize": source.FS_TICK, "legend.fontsize": source.FS_LEGEND,
        "axes.linewidth": source.LW_AXIS, "figure.dpi": 300,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    beta = (data.loc[data["method"].eq("beta")]
            .sort_values(["lambda_gc", "cohort"], kind="mergesort")["cohort"]
            .drop_duplicates().tolist())
    cohorts = beta + sorted(set(data["cohort"]) - set(beta))
    y_pos = {cohort: i for i, cohort in enumerate(cohorts)}
    xmax = max(2.0, float(np.ceil((data["lambda_gc"].max() + 0.08) / 0.25) * 0.25))
    ticks = source.MaxNLocator(
        nbins=5, steps=[1, 2, 2.5, 5, 10], min_n_ticks=4
    ).tick_values(0, xmax)
    ticks = np.unique(np.append(ticks[(ticks >= 0) & (ticks <= xmax)], 1.0))

    fig = plt.figure(figsize=(source.PANEL_W_MM / 25.4,
                              source.PANEL_H_MM / 25.4), facecolor="white")
    ax = fig.add_axes([0.285, 0.220, 0.55317, 0.715])
    fig.text(0.012, 0.985, "b", ha="left", va="top", fontsize=source.FS_PANEL,
             fontweight="bold", color="#222222")
    ax.axvline(1.0, color="#6E6E6E", linewidth=source.LW_REF,
               linestyle=(0, (3, 2)), zorder=1)
    for method in source.METHOD_ORDER:
        part = data.loc[data["method"].eq(method)]
        y = np.asarray([y_pos[c] + source.OFFSETS[method] for c in part["cohort"]])
        ax.scatter(part["lambda_gc"].to_numpy(float), y, s=16, marker="o",
                   facecolor=source.COLORS[method], edgecolor="white",
                   linewidth=0.22, alpha=0.94, zorder=3)
    ax.set_yticks(np.arange(len(cohorts)))
    ax.set_yticklabels(cohorts)
    ax.set_ylim(-0.70, len(cohorts) - 0.30)
    ax.set_xlim(0, xmax)
    ax.xaxis.set_major_locator(mticker.FixedLocator(np.sort(ticks)))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda value, _: "{:.1f}".format(value)))
    ax.set_xlabel(r"Inflation factor ($\lambda$)", labelpad=2.0)
    ax.set_ylabel("")
    ax.set_axisbelow(True)
    ax.grid(axis="x", which="major", color="#E8E8E8", linewidth=source.LW_GRID)
    ax.grid(axis="y", which="major", color="#DDDDDD", linewidth=source.LW_GRID)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set(color="#222222", linewidth=source.LW_AXIS)
    ax.tick_params(axis="both", which="major", direction="out", width=source.LW_AXIS,
                   length=1.3, pad=1.5, colors="#222222", labelsize=source.FS_TICK)
    handles = [Line2D([0], [0], marker="o", linestyle="none",
                      markerfacecolor=source.COLORS[m], markeredgecolor="white",
                      markeredgewidth=0.22, markersize=4.2,
                      label=source.METHOD_LABELS[m]) for m in source.METHOD_ORDER]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.135),
              ncol=3, frameon=False, fontsize=source.FS_LEGEND,
              handlelength=0.8, handletextpad=0.35, columnspacing=0.95,
              borderaxespad=0.0, labelspacing=0.2)
    return render(fig)


def qq_panel():
    source = source_module("plot_sle_figs4b_meta_qq.py", "figs10_qq_plot")
    wide = pd.read_csv(panel_file("FigS10c"))
    method_names = {
        "beta": "EWAS",
        "camda": "haEWAS-CAMDA",
        "chalm": "haEWAS-CHALM",
    }
    parts = []
    for method, prefix in method_names.items():
        expected_column = prefix + "_Expected"
        observed_column = prefix + "_Observed"
        points = wide[["CpG", expected_column, observed_column]].dropna(
            subset=[expected_column, observed_column]
        ).sort_values(expected_column, kind="mergesort")
        expected = points[expected_column].to_numpy(float)
        observed = points[observed_column].to_numpy(float)
        lambda_values = wide[prefix + "_Lambda"].dropna()
        if len(expected) != len(observed) or len(lambda_values) != 1:
            raise ValueError("FigS10c has an incomplete method column")
        part = source.thin_qq_data(expected, observed, method)
        part["Lambda"] = float(lambda_values.iloc[0])
        parts.append(part)
    data = pd.concat(parts, ignore_index=True)
    source.style()
    source.mpl.rcParams.update({
        "font.size": 5.5, "axes.labelsize": source.FS_AXIS,
        "xtick.labelsize": source.FS_TICK, "ytick.labelsize": source.FS_TICK,
        "legend.fontsize": source.FS_LEGEND, "axes.linewidth": source.LW_AXIS,
        "figure.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    ymax, yticks = source.determine_y_axis(data)
    fig = plt.figure(figsize=(source.PANEL_W_MM / 25.4,
                              source.PANEL_H_MM / 25.4), facecolor="white")
    ax = fig.add_axes([0.165, 0.220, 0.810, 0.715])
    fig.text(0.012, 0.985, "c", ha="left", va="top", fontsize=source.FS_PANEL,
             fontweight="bold", color="#222222")
    ax.plot([0, source.X_AXIS_MAX], [0, source.X_AXIS_MAX], color="#666666",
            linewidth=source.LW_REF, linestyle=(0, (3, 2)), zorder=1)
    for method in source.METHOD_ORDER:
        part = data.loc[data["Method"].eq(method)]
        ax.scatter(part["Expected"].to_numpy(float), part["Observed"].to_numpy(float),
                   s=2.2, marker="o", facecolor=source.COLORS[method],
                   edgecolor="none", linewidth=0, alpha=0.52,
                   rasterized=True, zorder=2)
    ax.set_xlim(0, source.X_AXIS_MAX)
    ax.set_ylim(0, ymax)
    ax.set_xticks([0, 2, 4, 6])
    ax.set_yticks(yticks)
    ax.set_xlabel(r"Expected $-\log_{10}(P)$", labelpad=2.0)
    ax.set_ylabel(r"Observed $-\log_{10}(P)$", labelpad=2.0)
    ax.set_axisbelow(True)
    ax.grid(axis="both", which="major", color="#EEEEEE", linewidth=source.LW_GRID)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set(color="#222222", linewidth=source.LW_AXIS)
    ax.tick_params(axis="both", which="major", direction="out", width=source.LW_AXIS,
                   length=1.3, pad=1.5, colors="#222222", labelsize=source.FS_TICK)
    handles, labels = [], []
    for method in source.METHOD_ORDER:
        handles.append(Line2D([0], [0], marker="o", linestyle="none",
                              markerfacecolor=source.COLORS[method], markeredgecolor="none",
                              markersize=3.0, alpha=1.0))
        value = float(data.loc[data["Method"].eq(method), "Lambda"].iloc[0])
        labels.append("{}: λ={:.2f}".format(source.METHOD_LABELS[method], value))
    ax.legend(handles=handles, labels=labels, loc="upper left",
              bbox_to_anchor=(0.020, 0.980), frameon=False,
              fontsize=source.FS_LEGEND, handlelength=0.8,
              handletextpad=0.35, borderaxespad=0.0, labelspacing=0.25)
    return render(fig)


def go_panel():
    source = source_module("plot_sle_figs4c_gene_lists_pooled_v2.py", "figs10_go_plot")
    data = pd.read_csv(panel_file("FigS10d"))
    fig = source._draw(data)
    for label in fig.texts:
        if label.get_text() == "c":
            label.set_text("d")
    return render(fig)


def location_panel():
    source = source_module("plot_sle_figs4d_pooled_genomic_locations_v2.py", "figs10_location_plot")
    data = pd.read_csv(panel_file("FigS10e"))
    captured = []
    original_figure = plt.figure
    original_finish = source.finish

    def capture_figure(*args, **kwargs):
        fig = original_figure(*args, **kwargs)
        captured.append(fig)
        return fig

    plt.figure = capture_figure
    source.finish = lambda fig, path: None
    try:
        source._draw_figure(data)
    finally:
        plt.figure = original_figure
        source.finish = original_finish
    if not captured:
        raise RuntimeError("Could not capture the FigS10e source plot")
    fig = captured[-1]
    for label in fig.texts:
        if label.get_text() == "d":
            label.set_text("e")
    return render(fig)


def run():
    panels = {
        "a": upset_panel(),
        "b": lambda_panel(),
        "c": qq_panel(),
        "d": go_panel(),
        "e": location_panel(),
    }
    a_width = 82.0
    a_height = a_width * 45.0 / 59.0
    b_y = 4.0 + a_height + 8.0
    d_y = b_y + 58.0 + 8.0
    fig_height = d_y + 84.0 + 4.0
    fig = plt.figure(figsize=(FIG_W_MM / 25.4, fig_height / 25.4), facecolor="white")
    layout = {
        "a": (18.5, 4.0, a_width, a_height),
        "b": (18.5, b_y, 82.0, 58.0),
        "c": (108.5, b_y, 56.0, 58.0),
        "d": (5.5, d_y, 88.0, 84.0),
        "e": (101.5, d_y, 76.0, 84.0),
    }
    for key, image in panels.items():
        x, y, width, height = layout[key]
        ax = fig.add_axes([
            x / FIG_W_MM, 1 - (y + height) / fig_height,
            width / FIG_W_MM, height / fig_height,
        ])
        ax.imshow(image, interpolation="none", aspect="auto")
        ax.set_axis_off()

    out = ROOT / "figures" / "FigS10.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="pdf", facecolor="white", dpi=600)
    from IPython.display import display
    display(fig)
    plt.close(fig)


if __name__ == "__main__":
    run()
