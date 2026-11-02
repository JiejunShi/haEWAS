"""Plot Fig5a.pdf from the prepared Fig5a.csv only."""

import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ["REPRO_ROOT"])


def style():
    import matplotlib as mpl
    mpl.rcParams.update({"font.family": "Arial", "font.size": 6,
                         "pdf.fonttype": 42, "figure.facecolor": "white"})


SOURCE = ROOT / "data" / "Fig5a.csv"
PDF = ROOT / "figures" / "Fig5a.pdf"

COLUMNS = ("haEWAS-specific", "Common", "EWAS-specific")
LABELS = {
    "haEWAS-specific": "haEWAS-\nspecific",
    "Common": "Common",
    "EWAS-specific": "EWAS-\nspecific",
}
PLOT_COLUMNS = list(LABELS)

MM = 1 / 25.4
PANEL_W_MM = 54
PANEL_H_MM = 82.8
FS_TITLE = 7.0
FS_AXIS = 6.0
FS_TICK = 5.6
FS_VALUE = 5.2
FS_CBAR = 5.4
LW_TICK = 0.45
LW_CBAR = 0.40
COLOR_AXIS = "#272727"
COLOR_ZERO = "#EEEEEE"
CBAR_TITLE = "Number"

PALETTE = ("#EEF2F7", "#C8D6EA", "#8FAFD3", "#4F8FC0", "#0F5B8C")


def _display_and_save(fig):
    """Save the exact-size PDF and display the same figure in Jupyter."""
    import matplotlib.pyplot as plt

    PDF.parent.mkdir(parents=True, exist_ok=True)
    # Deliberately omit bbox_inches="tight" to retain the requested 54 x 82.8 mm.
    fig.savefig(PDF, format="pdf", facecolor="white")
    try:
        from IPython import get_ipython
        from IPython.display import display

        if get_ipython() is not None:
            display(fig)
    except ImportError:
        pass
    plt.close(fig)


def run():
    table = pd.read_csv(SOURCE)
    table = table.sort_values("Cancer").copy()
    if len(table) != 15 or table["Cancer"].duplicated().any():
        raise ValueError("Fig. 5a requires one row for each of 15 TCGA cancers")
    for name in COLUMNS:
        table[name] = np.asarray(table[name], dtype=float)
    values = table[PLOT_COLUMNS].to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Locus counts must be finite, non-negative values")

    # Match the original Fig. 5a: cancers with no significant loci are omitted.
    keep = values.sum(axis=1) > 0
    table = table.loc[keep].reset_index(drop=True)
    values = values[keep]
    if table.empty:
        raise ValueError("No cancer types with survival-associated loci were found")

    row_sums = values.sum(axis=1, keepdims=True)
    proportions = np.divide(
        values, row_sums, out=np.zeros_like(values), where=row_sums != 0
    )

    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.colors import LinearSegmentedColormap, PowerNorm

    style()
    plt.rcParams.update({
        "font.size": 5.5,
        "axes.labelsize": FS_AXIS,
        "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "axes.linewidth": 0.55,
        "patch.linewidth": 0.45,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    })
    cmap = LinearSegmentedColormap.from_list("nature_blue_scale", PALETTE, N=256)
    norm = PowerNorm(gamma=0.55, vmin=0, vmax=1.0)

    fig = plt.figure(
        figsize=(PANEL_W_MM * MM, PANEL_H_MM * MM),
        facecolor="white",
        constrained_layout=False,
    )
    ax = fig.add_axes([0.105, 0.090, 0.675, 0.820])

    for row in range(len(table)):
        for col in range(len(PLOT_COLUMNS)):
            raw_value = values[row, col]
            share = proportions[row, col]
            face_color = cmap(norm(share)) if raw_value > 0 else COLOR_ZERO
            rect = patches.FancyBboxPatch(
                (col - 0.46, row - 0.42), 0.92, 0.84,
                boxstyle="round,pad=0.025",
                facecolor=face_color,
                edgecolor="white",
                linewidth=0.20,
                zorder=3,
            )
            ax.add_patch(rect)
            ax.text(
                col, row, "{:,}".format(int(raw_value)),
                ha="center", va="center", fontsize=FS_VALUE,
                fontweight="normal",
                color="white" if share > 0.45 else COLOR_AXIS,
                zorder=4,
            )

    ax.set_xlim(-0.50, len(PLOT_COLUMNS) - 0.50)
    ax.set_ylim(len(table) - 0.50, -0.50)
    ax.set_yticks(np.arange(len(table)))
    ax.set_yticklabels(
        table["Cancer"].astype(str).str.replace("TCGA-", "", regex=False),
        fontsize=FS_TICK, fontweight="normal", color=COLOR_AXIS,
    )
    ax.set_xticks(np.arange(len(PLOT_COLUMNS)))
    ax.set_xticklabels(
        [LABELS[name] for name in PLOT_COLUMNS],
        rotation=45, ha="right", rotation_mode="anchor",
        fontsize=FS_TICK, fontweight="normal", color=COLOR_AXIS,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(
        axis="both", which="major", bottom=True, left=True,
        direction="out", length=2.0, width=LW_TICK,
        color=COLOR_AXIS, pad=1.5,
    )

    fig.text(
        0.455, 0.982, "Survival-associated loci",
        ha="center", va="top", fontsize=FS_TITLE,
        fontweight="normal", color=COLOR_AXIS,
    )

    cax = fig.add_axes([0.835, 0.610, 0.030, 0.205])
    scalar = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(scalar, cax=cax, orientation="vertical")
    cbar.set_ticks([])
    cbar.outline.set_linewidth(LW_CBAR)
    cbar.outline.set_edgecolor(COLOR_AXIS)
    cax.text(
        0.50, 1.18, CBAR_TITLE, va="bottom", ha="center",
        transform=cax.transAxes, fontsize=FS_CBAR,
        fontweight="normal", color=COLOR_AXIS,
    )
    cax.text(
        1.50, 1.00, "High", va="center", ha="left",
        transform=cax.transAxes, fontsize=FS_CBAR,
        fontweight="normal", color=COLOR_AXIS,
    )
    cax.text(
        1.50, 0.00, "Low", va="center", ha="left",
        transform=cax.transAxes, fontsize=FS_CBAR,
        fontweight="normal", color=COLOR_AXIS,
    )
    _display_and_save(fig)


if __name__ == "__main__":
    run()
