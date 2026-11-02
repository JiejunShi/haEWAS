"""Plot FigS7a-d from the prepared source tables."""

import os
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu


ROOT = Path(os.environ["REPRO_ROOT"])
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
METHODS = ("Beta", "CAMDA", "CHALM")
METHOD_COLORS = {"Beta": "#4C78A8", "CAMDA": "#72A555", "CHALM": "#B85C5C"}
METHOD_TITLES = {"Beta": r"$\bar{m}$", "CAMDA": r"$h_m$(CAMDA)", "CHALM": r"$h_m$(CHALM)"}
PANEL_TITLES = {
    "a": "Medulloblastoma / Brain",
    "b": "Schizophrenia / Blood",
    "c": "Liver Cancer / Liver",
}
DISEASES = ("Medulloblastoma", "Schizophrenia", "Liver Cancer")
CATEGORIES = [
    "Beta (Raw)", "CAMDA (Raw)", "CHALM (Raw)",
    "Beta (M-val)", "CAMDA (M-val)", "CHALM (M-val)",
]
AXIS = "#272727"
M_VALUE = "#1F1F1F"
PANEL_W = 171


def to_m_value(values):
    values = np.asarray(values, dtype=float)
    values = np.clip(values, 1e-6, 1 - 1e-6)
    return np.log2(values / (1 - values))


def panel_file(name):
    path = DATA / (name + ".csv")
    return path if path.is_file() else Path(str(path) + ".gz")


def clean_values(values):
    values = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return values[np.isfinite(values)]


def draw_kde(ax, values, color, **kwargs):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3 or np.unique(values).size < 2:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                ha="center", va="center", fontsize=5, color=AXIS)
        return

    original_plot = ax.plot

    def plot_numpy(*args, **plot_kwargs):
        args = tuple(
            np.asarray(value) if isinstance(value, (pd.Index, pd.Series)) else value
            for value in args
        )
        return original_plot(*args, **plot_kwargs)

    ax.plot = plot_numpy
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sns.kdeplot(values, color=color, ax=ax, **kwargs)
    finally:
        ax.plot = original_plot


def style_raw_axis(ax, color, show_ylabel):
    ax.set(xlim=(0, 1), xticks=[0, 0.25, 0.5, 0.75, 1])
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda value, _: "0" if abs(value) < 1e-9 else
        ("1.00" if abs(value - 1) < 1e-9 else "{:.2f}".format(value))
    ))
    ax.set_xlabel("Raw value (0-1)", labelpad=1, color=color)
    ax.set_ylabel("Density (raw)" if show_ylabel else "", labelpad=1.4, color=color)
    ax.tick_params(axis="both", labelsize=5, width=0.45, length=1.8,
                   pad=1, colors=color)
    for side in ("left", "bottom"):
        ax.spines[side].set(color=color, linewidth=0.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    top = max(2, int(np.ceil(ax.get_ylim()[1] / 2) * 2))
    ticks = np.arange(0, top + 0.001, 2)
    ax.set(yticks=ticks, ylim=(0, ticks[-1]))
    ax.grid(False)


def style_m_axes(ax_right, ax_top, show_ylabel):
    ax_top.set(xlim=(-10, 10), xticks=[-10, -5, 0, 5, 10])
    ax_top.set_xlabel("M-value", labelpad=1, color=M_VALUE)
    ax_top.tick_params(axis="x", labelsize=5, width=0.45, length=1.8,
                       pad=0.8, colors=M_VALUE)
    ax_top.tick_params(axis="y", left=False, right=False,
                       labelleft=False, labelright=False)
    ax_top.spines["top"].set(visible=True, color=M_VALUE, linewidth=0.5)
    for side in ("bottom", "left", "right"):
        ax_top.spines[side].set_visible(False)
    ax_right.set_ylabel("Density (M-value)" if show_ylabel else "",
                        fontsize=5.4, rotation=270, labelpad=6, color=M_VALUE)
    ax_right.tick_params(axis="y", labelsize=5, width=0.45, length=1.8,
                         pad=1, colors=M_VALUE, right=True, labelright=True)
    top = ax_right.get_ylim()[1]
    step = 0.2 if top > 0.2 else 0.1
    top = np.ceil(top / step) * step
    ticks = np.arange(0, top + step * 0.5, step)
    ax_right.set(yticks=ticks, ylim=(0, top))
    ax_right.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax_right.spines["right"].set(visible=True, color=M_VALUE, linewidth=0.5)
    for side in ("left", "top", "bottom"):
        ax_right.spines[side].set_visible(False)
    ax_right.grid(False)
    ax_top.grid(False)


def plot_kde_panel(panel):
    data = pd.read_csv(panel_file("FigS7" + panel))
    if list(data.columns) != ["Method", "RawValue", "MValue"]:
        raise ValueError("Unexpected FigS7{} data columns".format(panel))
    fig = plt.figure(figsize=(171 / 25.4, 50 / 25.4), facecolor="white")
    fig.text(0.5, 0.982, PANEL_TITLES[panel], ha="center", va="top", fontsize=7)
    fig.text(0.004, 0.982, panel, ha="left", va="top", fontsize=8, fontweight="bold")
    for index, method in enumerate(METHODS):
        part = data.loc[data["Method"].eq(method)]
        raw, m_value = clean_values(part["RawValue"]), clean_values(part["MValue"])
        if not raw.size or not m_value.size:
            if not raw.size:
                raise ValueError("FigS7{} is missing {}".format(panel, method))
            m_value = to_m_value(raw)
        start = index * 60
        ax = fig.add_axes([(start + 9) / 171, 8.3 / 50, 36.5 / 171, 22 / 50])
        fig.text((start + 25.5) / 171, 41.8 / 50, METHOD_TITLES[method],
                 ha="center", va="top", fontsize=6)
        draw_kde(ax, raw, METHOD_COLORS[method], fill=True, alpha=0.24,
                 linewidth=0.72, cut=0, clip=(0, 1))
        style_raw_axis(ax, METHOD_COLORS[method], index == 0)
        right = ax.twinx()
        top = right.twiny()
        right.patch.set_alpha(0)
        top.patch.set_alpha(0)
        draw_kde(top, m_value, M_VALUE, linestyle=(0, (2, 2)),
                 linewidth=0.7, cut=0, clip=(-10, 10))
        style_m_axes(right, top, index == 2)
    fig.savefig(FIGURES / ("FigS7" + panel + ".pdf"), facecolor="white")
    from IPython.display import display
    display(fig)
    plt.close(fig)


def significance(p):
    return "ns" if p > 0.05 else ("*" if p > 0.01 else ("**" if p > 0.001 else "***"))


def y_axis(values, data):
    if {"Y_max", "Y_step"}.issubset(data.columns):
        ymax = pd.to_numeric(data["Y_max"], errors="coerce").dropna()
        step = pd.to_numeric(data["Y_step"], errors="coerce").dropna()
        if not ymax.empty and not step.empty and ymax.iloc[0] > 0 and step.iloc[0] > 0:
            return float(ymax.iloc[0]), float(step.iloc[0])
    maximum = max(float(value.max()) for value in values if len(value))
    if maximum <= 0.3:
        return 0.3, 0.1
    ymax = float(np.ceil(maximum * 1.2 * 10) / 10)
    return ymax, ymax / 3


def make_d_item(data):
    data["Category"] = data["Method"].astype(str) + " (" + data["Scale"].astype(str) + ")"
    pivot = data.pivot_table(index="Comparison", columns="Category",
                             values="KSStatistic", aggfunc="mean").dropna()
    if any(name not in pivot for name in CATEGORIES):
        raise ValueError("FigS7d is missing a method/scale category")
    values = [pivot[name].to_numpy(float) for name in CATEGORIES]
    pairs = (("Beta", "CAMDA"), ("CAMDA", "CHALM"), ("Beta", "CHALM"))
    def compare(scale):
        return [float(mannwhitneyu(pivot["{} ({})".format(left, scale)],
                                   pivot["{} ({})".format(right, scale)],
                                   alternative="two-sided")[1])
                for left, right in pairs]
    return {
        "disease": str(data["Disease"].iloc[0]),
        "cohort": str(data["Cohort"].iloc[0]),
        "values": values,
        "p_raw": compare("Raw"),
        "p_m": compare("M-val"),
        "axis": y_axis(values, data),
    }


def draw_d_item(ax, item, show_ylabel):
    raw = np.array([0.0, 0.7, 1.4])
    mval = raw + 2.55
    positions = np.r_[raw, mval]
    colors = ["#4C78A8", "#72A555", "#B85C5C", "#2E5F8A", "#4F7F3A", "#B85C5C"]
    boxes = ax.boxplot(item["values"], positions=positions, widths=0.28,
                       patch_artist=True, showfliers=False,
                       medianprops={"color": "#333333", "linewidth": 0.65},
                       whiskerprops={"color": "#333333", "linewidth": 0.45},
                       capprops={"color": "#333333", "linewidth": 0.45})
    for index, box in enumerate(boxes["boxes"]):
        box.set(facecolor=colors[index], alpha=0.55 if index < 3 else 0.82,
                edgecolor="#333333", linewidth=0.45,
                hatch="" if index < 3 else "///")
    ax.axvline((raw[-1] + mval[0]) / 2, color="#A8A8A8", linewidth=0.4,
               linestyle=(0, (2, 2)))
    ymax, step = item["axis"]
    ax.set(xlim=(positions[0] - 0.35, positions[-1] + 0.35), ylim=(0, ymax),
           xticks=positions, yticks=np.arange(0, ymax + 1e-9, step))
    ax.set_xticklabels(["Beta", "hm(CAMDA)", "hm(CHALM)"] * 2,
                       fontsize=5, rotation=28, ha="right")
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax.tick_params(width=0.5, length=1.7, pad=1, labelsize=5, colors=AXIS)
    ax.set_ylabel("KS statistic" if show_ylabel else "", labelpad=1.2)
    ax.text(raw.mean(), -0.45, "Raw", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=5, clip_on=False)
    ax.text(mval.mean(), -0.45, "M-value", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=5, clip_on=False)
    ax.set_title("{}\n{}".format(item["disease"], item["cohort"]), fontsize=6, pad=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_bounds(positions[0], positions[-1])
    ax.spines["left"].set_bounds(0, ymax)

    y, h, label_pad, offset = ymax * 0.70, ymax * 0.075, ymax * 0.018, 0.05
    for group, pvalues in ((raw, item["p_raw"]), (mval, item["p_m"])):
        for left, right, level, p in (
            (group[0], group[1], y + h, pvalues[0]),
            (group[1], group[2], y + h, pvalues[1]),
            (group[0], group[2], y + h * 2.6, pvalues[2]),
        ):
            ax.plot([left + offset, right - offset], [level, level],
                    linewidth=0.45, color="#333333")
            ax.text((left + right) / 2, level + label_pad, significance(p),
                    ha="center", va="bottom", fontsize=5, clip_on=True)


def plot_d():
    data = pd.read_csv(panel_file("FigS7d"))
    items = []
    for disease in DISEASES:
        part = data.loc[data["Disease"].eq(disease)].copy()
        if part.empty:
            raise ValueError("FigS7d is missing {}".format(disease))
        items.append(make_d_item(part))
    fig = plt.figure(figsize=(171 / 25.4, 46 / 25.4), facecolor="white")
    group = 171 / len(items)
    for index, item in enumerate(items):
        ax = fig.add_axes([(index * group + 4.6) / 171, 12.2 / 46,
                           (group - 8.2) / 171, 22 / 46])
        draw_d_item(ax, item, index == 0)
    fig.text(0.004, 0.985, "d", ha="left", va="top", fontsize=8, fontweight="bold")
    fig.savefig(FIGURES / "FigS7d.pdf", facecolor="white")
    from IPython.display import display
    display(fig)
    plt.close(fig)


def run():
    FIGURES.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
        "font.size": 5.5, "axes.labelsize": 5.4, "axes.titlesize": 6,
        "xtick.labelsize": 5, "ytick.labelsize": 5, "axes.linewidth": 0.5,
        "figure.dpi": 300, "savefig.dpi": 600, "savefig.facecolor": "white",
        "savefig.edgecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.unicode_minus": False, "xtick.direction": "out", "ytick.direction": "out",
        "hatch.linewidth": 0.35,
    })
    for panel in ("a", "b", "c"):
        plot_kde_panel(panel)
    plot_d()


if __name__ == "__main__":
    run()
