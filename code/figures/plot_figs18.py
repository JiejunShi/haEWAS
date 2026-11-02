"""Draw the haEWAS-only FigS18 panels from the prepared data."""

import gzip
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
DATA_FILE = ROOT / "data" / "FigS18.csv"
OUTPUT = ROOT / "figures" / "FigS18.pdf"
HELDOUTS = ("ROSMAP", "MtSinai", "LONDON", "GASPARONI")
COLOR = "#C44E52"
AXIS = "#272727"
ZERO = "#A5A5A5"
DIAGONAL = "#666666"
EFFECT_MIN, EFFECT_MAX = -1.0, 1.0
WIDTH_MM, HEIGHT_MM = 125.0, 90.0


def read_data():
    path = DATA_FILE if DATA_FILE.is_file() else Path(str(DATA_FILE) + ".gz")
    if not path.is_file():
        raise FileNotFoundError("Prepare FigS18 data first")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        data = pd.read_csv(handle, low_memory=False)
    columns = ["CpG", "Cohort", "DiscoveryEffect", "HeldoutEffect"]
    if list(data.columns) != columns or data.empty:
        raise ValueError("FigS18 data has an unexpected format")
    return data


def draw_quadrant_label(ax, counts, total, quadrant, upper, right):
    n = counts[quadrant]
    pct = 100 * n / total if total else np.nan
    pct_text = "{:.1f}%".format(pct) if np.isfinite(pct) else "NA"
    x_center, x_left = (0.93, 0.82) if right else (0.07, 0.03)
    box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
           "boxstyle": "round,pad=0.12"}
    if upper:
        label_y, detail_y, va = 0.98, 0.93, "top"
    else:
        label_y, detail_y, va = 0.15, 0.03, "bottom"
    ax.text(x_center, label_y, quadrant, transform=ax.transAxes,
            ha="center", va=va, fontsize=6.0, color=AXIS, bbox=box, zorder=3)
    ax.text(x_left, detail_y, "n={}\n{}".format(n, pct_text),
            transform=ax.transAxes, ha="left", va=va, fontsize=6.0,
            linespacing=0.9, color=AXIS, bbox=box, zorder=3)


def quadrant_counts(x, y):
    masks = {
        "I": (x > 0) & (y > 0),
        "II": (x < 0) & (y > 0),
        "III": (x < 0) & (y < 0),
        "IV": (x > 0) & (y < 0),
    }
    counts = {name: int(mask.sum()) for name, mask in masks.items()}
    return counts, sum(counts.values())


def run():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.2,
        "axes.labelsize": 6.2,
        "axes.titlesize": 6.2,
        "xtick.labelsize": 6.0,
        "ytick.labelsize": 6.0,
        "axes.linewidth": 0.45,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    data = read_data()
    figure, axes = plt.subplots(
        2, 3, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
        squeeze=False,
    )
    letters = "abcd"
    for i, cohort in enumerate(HELDOUTS):
        row, col = divmod(i, 3)
        ax = axes[row, col]
        part = data.loc[data["Cohort"].astype(str).eq(cohort)]
        x = part["DiscoveryEffect"].to_numpy(float)
        y = part["HeldoutEffect"].to_numpy(float)
        counts, total = quadrant_counts(x, y)

        ax.axhline(0, color=ZERO, linewidth=0.45, zorder=0)
        ax.axvline(0, color=ZERO, linewidth=0.45, zorder=0)
        ax.plot([EFFECT_MIN, EFFECT_MAX], [EFFECT_MIN, EFFECT_MAX],
                color=DIAGONAL, linestyle=(0, (2.0, 2.0)),
                linewidth=0.55, zorder=0)
        if len(part):
            ax.scatter(x, y, s=8.5, marker="o", c=COLOR,
                       edgecolors="white", linewidths=0.22,
                       alpha=0.82, zorder=2)

        draw_quadrant_label(ax, counts, total, "I", True, True)
        draw_quadrant_label(ax, counts, total, "II", True, False)
        draw_quadrant_label(ax, counts, total, "III", False, False)
        draw_quadrant_label(ax, counts, total, "IV", False, True)
        ax.set_title("{} (N={})".format(cohort, total), loc="center", pad=2.4)
        ax.text(-0.12, 1.02, letters[i], transform=ax.transAxes,
                ha="left", va="bottom", fontsize=7.0,
                fontweight="bold", color=AXIS, clip_on=False)
        ax.set_xlim(EFFECT_MIN, EFFECT_MAX)
        ax.set_ylim(EFFECT_MIN, EFFECT_MAX)
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(length=2.0, width=0.4, pad=1.3)
        ax.set_xlabel("Discovery effect", labelpad=1.4)
        ax.set_ylabel("Held-out effect", labelpad=1.4)

    for ax in axes.ravel()[len(HELDOUTS):]:
        ax.set_axis_off()
    figure.subplots_adjust(left=0.09, right=0.995, bottom=0.04, top=0.94,
                           wspace=0.31, hspace=0.48)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", facecolor="white",
                   bbox_inches=None, pad_inches=0)
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
