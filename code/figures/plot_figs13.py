"""Plot the haEWAS-only FigS13 panels."""

import importlib.util
import os
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA = ROOT / "data"
OUTPUT = ROOT / "figures" / "FigS13.pdf"
HELDOUTS = ("GSE110607", "GSE118144", "GSE59250", "GSE82218", "GSE96879")
COLUMNS = ["CpG", "Heldout", "DiscoveryEffect", "HeldoutEffect"]


def load_helper():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    path = CODE / "r3_3_loco_direction_common.py"
    spec = importlib.util.spec_from_file_location("figs13_plot_helper", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def data_file():
    path = DATA / "FigS13.csv"
    if path.is_file():
        return path
    compressed = Path(str(path) + ".gz")
    if compressed.is_file():
        return compressed
    raise FileNotFoundError("FigS13 data is missing; prepare it first")


def panel_stats(frame):
    x = frame["DiscoveryEffect"].to_numpy(float)
    y = frame["HeldoutEffect"].to_numpy(float)
    masks = {
        "I": (x > 0) & (y > 0),
        "II": (x < 0) & (y > 0),
        "III": (x < 0) & (y < 0),
        "IV": (x > 0) & (y < 0),
    }
    counts = {name: int(mask.sum()) for name, mask in masks.items()}
    total = sum(counts.values())
    stats = {"panel_N": total}
    for name, count in counts.items():
        stats["quadrant_{}_n".format(name)] = count
        stats["quadrant_{}_pct".format(name)] = (
            100.0 * count / total if total else np.nan
        )
    return stats


def run():
    data = pd.read_csv(
        data_file(), dtype={"CpG": str, "Heldout": str}
    )
    if list(data.columns) != COLUMNS:
        raise ValueError("FigS13 data has unexpected columns")
    if data.duplicated(["Heldout", "CpG"]).any():
        raise ValueError("FigS13 has duplicate heldout/CpG rows")

    source = load_helper()
    source.style()
    mpl.rcParams.update({
        "font.size": 6.2, "axes.labelsize": 6.2, "axes.titlesize": 6.2,
        "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
        "axes.linewidth": 0.45,
    })

    ncols = 3
    nrows = int(np.ceil(len(HELDOUTS) / ncols))
    figure, axes = plt.subplots(
        nrows, ncols, figsize=(125.0 / 25.4, 86.0 / 25.4), squeeze=False
    )
    used = set()
    for index, heldout in enumerate(HELDOUTS):
        row, col = divmod(index, ncols)
        ax = axes[row, col]
        used.add(row * ncols + col)
        frame = data.loc[data["Heldout"].eq(heldout)]
        stats = panel_stats(frame)

        ax.set_box_aspect(1)
        ax.axhline(0, color=source.ZERO_COLOR, linewidth=0.45, zorder=0)
        ax.axvline(0, color=source.ZERO_COLOR, linewidth=0.45, zorder=0)
        ax.plot(
            [source.EFFECT_MIN, source.EFFECT_MAX],
            [source.EFFECT_MIN, source.EFFECT_MAX],
            color=source.DIAGONAL_COLOR, linestyle=(0, (2.0, 2.0)),
            linewidth=0.55, zorder=0,
        )
        if len(frame):
            ax.scatter(
                frame["DiscoveryEffect"], frame["HeldoutEffect"],
                s=8.5, marker="o", c=source.HAEWAS_COLOR,
                edgecolors="white", linewidths=0.22, alpha=0.82, zorder=2,
            )

        for quadrant in ("I", "II", "III", "IV"):
            source._draw_quadrant_label(
                ax, stats, quadrant,
                quadrant in ("I", "II"), quadrant in ("I", "IV"),
            )
        ax.set_title("{} (N={})".format(heldout, stats["panel_N"]),
                     loc="center", pad=2.4)
        ax.text(
            -0.12, 1.02, "abcde"[index], transform=ax.transAxes,
            ha="left", va="bottom", fontsize=7.0, fontweight="bold",
            color=source.AXIS_COLOR, clip_on=False,
        )
        ax.set_xlim(source.EFFECT_MIN, source.EFFECT_MAX)
        ax.set_ylim(source.EFFECT_MIN, source.EFFECT_MAX)
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(length=2.0, width=0.4, pad=1.3,
                       labelbottom=True, labelleft=True)
        ax.set_xlabel("Discovery effect", labelpad=1.4)
        ax.set_ylabel("Held-out effect", labelpad=1.4)

    for index, ax in enumerate(axes.ravel()):
        if index not in used:
            ax.set_axis_off()
    figure.subplots_adjust(
        left=0.09, right=0.995, bottom=0.05, top=0.98,
        wspace=0.31, hspace=0.48,
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", pad_inches=0)
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