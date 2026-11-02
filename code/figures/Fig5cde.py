"""Plot Fig5c-e and calculate curves and annotations from prepared data."""

import os
from pathlib import Path

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

ROOT = Path(os.environ["REPRO_ROOT"])

PANEL_CONFIGS = (
    {"panel": "c", "project": "TCGA-BLCA", "label": "BLCA",
     "cpg": "cg25900902", "gene": "ACAP1"},
    {"panel": "d", "project": "TCGA-BRCA", "label": "BRCA",
     "cpg": "cg24375409", "gene": "EPHA10"},
    {"panel": "e", "project": "TCGA-LIHC", "label": "LIHC",
     "cpg": "cg00067414", "gene": "MTHFD1L"},
)

COLOR_AXIS = "#272727"
COLOR_M_HIGH, COLOR_M_LOW = "#0F4D92", "#8FAFD3"
COLOR_HM_HIGH, COLOR_HM_LOW = "#8E5EA7", "#CDB8DB"
COLOR_INT_HIGH, COLOR_INT_LOW = "#B64342", "#E5B0AA"
LW_CURVE, LW_AXIS, LW_TICK = 1.00, 0.55, 0.50
FS_PANEL, FS_TITLE = 7.0, 6.0
FS_AXIS, FS_TICK, FS_LEGEND, FS_ANNOT = 6.0, 5.1, 4.9, 4.8
MM = 1 / 25.4
PANEL_W_MM, PANEL_H_MM = 82, 68
MIN_GROUP_SIZE = 6
FIGSIZE = (PANEL_W_MM * MM, PANEL_H_MM * MM)


def _format_p(value):
    if pd.isna(value):
        return "NA"
    if value < 1e-300:
        return "<1e-300"
    if value < 0.001:
        return "{:.1e}".format(value)
    return "{:.3f}".format(value)


def _style_axis(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(LW_AXIS)
        spine.set_color(COLOR_AXIS)
    ax.grid(False)
    ax.set_facecolor("white")
    ax.tick_params(axis="both", which="major", labelsize=FS_TICK,
        width=LW_TICK, length=2.0, direction="out", pad=0.8, colors=COLOR_AXIS)
    ax.set_ylim(0.0, 1.03)
    ax.yaxis.set_major_locator(ticker.FixedLocator([0.0, 0.25, 0.50, 0.75, 1.0]))
    ax.yaxis.set_major_formatter(ticker.FixedFormatter(["0", "0.25", "0.5", "0.75", "1"]))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4, prune=None))
    ax.set_xlabel("")
    ax.set_ylabel("")


def _group_statistics(subset):
    high = subset.loc[subset["Group"].eq("High"), ["Time", "Status"]].copy()
    low = subset.loc[subset["Group"].eq("Low"), ["Time", "Status"]].copy()
    n_high, n_low = len(high), len(low)
    if min(n_high, n_low) < MIN_GROUP_SIZE:
        raise ValueError("A Fig5 median group is too small")
    p_value = float(logrank_test(
        high["Time"], low["Time"],
        event_observed_A=high["Status"], event_observed_B=low["Status"],
    ).p_value)
    groups = pd.concat([high.assign(Group=1), low.assign(Group=0)], ignore_index=True)
    try:
        cph = CoxPHFitter()
        cph.fit(groups, duration_col="Time", event_col="Status")
        hr = float(cph.hazard_ratios_["Group"])
        ci = cph.confidence_intervals_.loc["Group"]
        hr_text = "HR = {:.2f} ({:.2f}-{:.2f})".format(
            hr, np.exp(ci.iloc[0]), np.exp(ci.iloc[1]))
    except Exception:
        hr_text = "HR = N/A"
    return p_value, hr_text


def plot_panel(table, cfg):
    required = {"Model", "Group", "Time", "Status"}
    if not required.issubset(table.columns):
        raise ValueError("Prepared survival table lacks required columns")
    table = table.copy()
    table["Time"] = pd.to_numeric(table["Time"], errors="coerce")
    table["Status"] = pd.to_numeric(table["Status"], errors="coerce")
    if (table["Time"].isna().any() or not table["Status"].isin((0, 1)).all()
            or not table["Model"].isin(("m", "hm")).all()
            or not table["Group"].isin(("High", "Low")).all()):
        raise ValueError("Prepared survival data contain invalid values")
    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    axes = [fig.add_axes([0.095, 0.255, 0.355, 0.430]),
            fig.add_axes([0.550, 0.255, 0.355, 0.430])]
    x_max = table["Time"].max()
    for ax, model, high_color, low_color in zip(
            axes, ("m", "hm"), (COLOR_M_HIGH, COLOR_INT_HIGH),
            (COLOR_M_LOW, COLOR_INT_LOW)):
        ax.set_box_aspect(1)
        subset = table.loc[table["Model"] == model]
        if subset.empty:
            raise ValueError("Missing {} model data".format(model))
        p_value, hr_text = _group_statistics(subset)
        for group, color in (("High", high_color), ("Low", low_color)):
            patients = subset.loc[subset["Group"] == group]
            if patients.empty:
                raise ValueError("Missing {} {} group".format(model, group))
            kmf = KaplanMeierFitter()
            kmf.fit(patients["Time"], event_observed=patients["Status"])
            curve = kmf.survival_function_.iloc[:, 0]
            ax.step(curve.index.to_numpy(dtype=float), curve.to_numpy(dtype=float),
                    where="post", color=color, linewidth=LW_CURVE,
                    label="{} {}".format(group, model))
        _style_axis(ax)
        ax.set_xlim(0, x_max * 1.02 if x_max > 0 else 1)
        ax.legend(frameon=False, loc="upper right", fontsize=FS_LEGEND,
                  handlelength=0.78, handletextpad=0.26, labelspacing=0.10,
                  borderaxespad=0.06)
        annotation = "P = {}\n{}".format(_format_p(p_value), hr_text)
        ax.text(0.05, 0.06, annotation, transform=ax.transAxes, fontsize=FS_ANNOT,
                color=COLOR_AXIS, ha="left", va="bottom", linespacing=1.10, zorder=20)
    title = "{}: {} ({})".format(cfg["label"], cfg["cpg"], cfg["gene"])
    fig.text(0.50, 0.905, title, ha="center", va="center", fontsize=FS_TITLE,
             fontweight="normal", color=COLOR_AXIS)
    fig.text(0.014, 0.985, cfg["panel"], ha="left", va="top", fontsize=FS_PANEL,
             fontweight="bold", color=COLOR_AXIS)
    fig.text(0.50, 0.085, "Overall survival (years)", ha="center", va="center",
             fontsize=FS_AXIS, color=COLOR_AXIS)
    fig.text(0.025, 0.470, "Survival probability", ha="center", va="center",
             rotation=90, fontsize=FS_AXIS, color=COLOR_AXIS)
    return fig


def run():
    plt.rcParams.update({"font.family": "Arial", "font.size": 6,
                         "axes.labelsize": 6, "axes.titlesize": 7,
                         "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
                         "legend.fontsize": 5.5, "axes.linewidth": 0.5,
                         "pdf.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "white",
                         "savefig.facecolor": "white"})
    (ROOT / "figures").mkdir(parents=True, exist_ok=True)
    for cfg in PANEL_CONFIGS:
        stem = f"Fig5{cfg['panel']}"
        table = pd.read_csv(ROOT / "data" / f"{stem}.csv")
        fig = plot_panel(table, cfg)
        fig.savefig(ROOT / "figures" / f"{stem}.pdf", format="pdf",
                    bbox_inches="tight", pad_inches=0.1)
        try:
            from IPython import get_ipython
            from IPython.display import display
            if get_ipython() is not None:
                display(fig)
        except ImportError:
            pass
        plt.close(fig)


if __name__ == "__main__":
    run()


