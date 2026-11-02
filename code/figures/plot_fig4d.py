"""Draw Fig4d using only its prepared CSV; export PDF and display inline."""
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(os.environ["REPRO_ROOT"])
PDF = ROOT / "figures" / "Fig4d.pdf"

def style():
    import matplotlib as mpl
    import matplotlib.font_manager as fm
    font_file = ROOT.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
    if font_file.is_file():
        fm.fontManager.addfont(str(font_file))
        family = fm.FontProperties(fname=str(font_file)).get_name()
    else:
        family = "Arial"
    mpl.rcParams.update({
        "font.family": family, "font.size": 6, "axes.labelsize": 6,
        "axes.titlesize": 7, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
        "legend.fontsize": 5.5, "axes.linewidth": 0.5, "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })

def finish(fig, bbox_inches=None):
    import matplotlib.pyplot as plt
    PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PDF, format="pdf", facecolor="white", bbox_inches=bbox_inches)
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None:
            display(fig)
    except ImportError:
        pass
    plt.close(fig)

from collections import OrderedDict
FDR_THRESHOLD = 0.05

MM_TO_INCH = 1.0 / 25.4

PANEL_WIDTH_MM = 101.0

PANEL_HEIGHT_MM = 90.0

COLOR_EWAS = "#0F4D92"

COLOR_HAEWAS = "#B64342"

COLOR_AXIS = "#272727"

COLOR_CONNECTOR = "#7A7A7A"

COLOR_THRESHOLD = "#6F6F6F"

PATHWAY_DATA = OrderedDict([
    ("GO: Postsynapse organization", {
        "color": "#D7E7D9",
        "genes": ["INSR", "TREM2", "SIGMAR1", "UBE3B", "PICK1", "NLGN2",
                  "ITSN1", "KIF1A", "WASF2", "FGFR1", "DBN1", "SLC12A5", "CNTNAP1"],
    }),
    ("GO: Vesicle-mediated transport", {
        "color": "#DCE7F2",
        "genes": ["CEBPE", "INSR", "TREM2", "IL10RA", "CHMP4A", "PIK3CB", "EXOC7",
                  "VAMP3", "UNC13D", "CLNK", "SSC5D", "PRKAR1B", "SLC2A4", "DGKQ",
                  "TRARG1", "INSIG1", "GNAI2", "WDR81", "PICK1", "NLGN2", "ITSN1",
                  "KIF1A", "WASF2", "LAT2", "PLEK", "S100A13", "KIF5A", "RAB5B",
                  "VPS51", "AP1G2", "RABIF", "DNAJC5", "CBARP", "ATP6V0A4", "CNN2"],
    }),
    ("GO: Regulation of immune system process", {
        "color": "#F2E1C5",
        "genes": ["IKBKB", "BRD4", "IRF4", "HOXA9", "ZNF580", "HOXA7", "SKAP1",
                  "SMARCD2", "ICAM1", "TREM2", "PLA2G2A", "ISG15", "TRIM11", "DTX4",
                  "PIK3CB", "EIF4E2", "VAMP3", "UNC13D", "CLNK", "CD2", "CFD",
                  "TREML4", "IFNL2", "UBR2", "PRKAR1B", "LAT2", "PRELID1", "PARP3",
                  "GPNMB", "CD27", "ZBTB46", "CD70", "SPPL2B", "SECTM1"],
    }),
    ("GO: Cell activation", {
        "color": "#E4DDEE",
        "genes": ["HDAC5", "SBNO2", "BRD4", "IRF4", "HOXA9", "SMARCD2", "ICAM1",
                  "TREM2", "PLA2G2A", "AICDA", "PIK3CB", "VAMP3", "UNC13D", "CLNK",
                  "CD2", "DGKQ", "ZP1", "CYGB", "LAT2", "PLEK", "S100A13",
                  "PRELID1", "PARP3", "GPNMB", "CD27", "ZBTB46", "CD70", "SKAP2"],
    }),
    ("KEGG: Insulin signaling pathway", {
        "color": "#E6DDD6",
        "genes": ["IKBKB", "INSR", "INPP5K", "PIK3CB", "EXOC7", "EIF4E2", "FASN",
                  "PRKAR1B", "SLC2A4", "MKNK2"],
    }),
])

FS_PANEL = 7.0

FS_AXIS = 6.0

FS_TICK = 5.2

FS_INLINE = 5.3

FS_PATHWAY = 4.8

FS_THRESHOLD = 5.0

LW_AXIS = 0.55

LW_TICK = 0.50

LW_CONNECTOR = 0.45

LW_THRESHOLD = 0.50

LW_SEPARATOR = 0.55

LW_POINT_EDGE = 0.25

POINT_SIZE = 18

BG_ALPHA = 0.68

def _escape_mathtext(value):
    return (str(value).replace("\\", r"\\").replace("_", r"\_")
            .replace("{", r"\{").replace("}", r"\}").replace("$", r"\$")
            .replace("#", r"\#").replace("%", r"\%").replace("&", r"\&"))

def _draw(selected):
    import matplotlib.pyplot as plt


    style()
    beta_x = selected["Log_beta"].to_numpy(dtype=float)
    ha_x = selected["Log_haewas"].to_numpy(dtype=float)
    y_pos = np.arange(len(selected))
    threshold = -np.log10(FDR_THRESHOLD)
    x_max_data = float(np.nanmax([np.nanmax(beta_x), np.nanmax(ha_x), threshold]))
    x_right = np.ceil((x_max_data + max(0.9, 0.16 * x_max_data)) * 2.0) / 2.0
    x_left = -0.15
    x_span = x_right - x_left

    fig = plt.figure(
        figsize=(PANEL_WIDTH_MM * MM_TO_INCH, PANEL_HEIGHT_MM * MM_TO_INCH),
        facecolor="white",
        constrained_layout=False,
    )
    ax = fig.add_axes([0.30, 0.10, 0.66, 0.82])

    for pathway, info in PATHWAY_DATA.items():
        idx = selected.index[selected["Pathway"].eq(pathway)].to_numpy()
        if not len(idx):
            continue
        y_min, y_max = idx.min() - 0.5, idx.max() + 0.5
        ax.axhspan(y_min, y_max, facecolor=info["color"], alpha=BG_ALPHA,
                   edgecolor="none", zorder=0)
        ax.text(x_right - 0.025 * x_span, (y_min + y_max) / 2.0, pathway,
                ha="right", va="center", fontsize=FS_PATHWAY, color=COLOR_AXIS,
                zorder=6, bbox={"facecolor": "white", "edgecolor": "none",
                                "alpha": 0.55, "pad": 0.6})

    for pathway in selected["Pathway"].drop_duplicates():
        idx = selected.index[selected["Pathway"].eq(pathway)].to_numpy()
        if len(idx):
            ax.axhline(idx.min() - 0.5, color="white", linewidth=LW_SEPARATOR, zorder=1)
            ax.axhline(idx.max() + 0.5, color="white", linewidth=LW_SEPARATOR, zorder=1)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(LW_AXIS)
        spine.set_color(COLOR_AXIS)
    ax.grid(False)
    ax.set_facecolor("white")
    ax.tick_params(axis="both", labelsize=FS_TICK, width=LW_TICK, length=2.0,
                   direction="out", pad=1.2, colors=COLOR_AXIS)
    ax.axvline(threshold, color=COLOR_THRESHOLD, linestyle=(0, (2.2, 2.2)),
               linewidth=LW_THRESHOLD, alpha=0.85, zorder=2)
    ax.text(threshold + 0.06, 1.01, "FDR = 0.05", transform=ax.get_xaxis_transform(),
            color=COLOR_THRESHOLD, fontsize=FS_THRESHOLD, va="bottom", ha="left",
            clip_on=False)

    for x, color, label in ((0.70, COLOR_EWAS, "EWAS"), (0.84, COLOR_HAEWAS, "haEWAS")):
        ax.scatter([x], [1.035], transform=ax.transAxes, s=POINT_SIZE, color=color,
                   edgecolors=COLOR_AXIS, linewidth=LW_POINT_EDGE, clip_on=False, zorder=8)
        ax.text(x + 0.025, 1.035, label, transform=ax.transAxes, ha="left", va="center",
                fontsize=FS_INLINE, color=color, clip_on=False)

    for i, (beta_value, ha_value) in enumerate(zip(beta_x, ha_x)):
        ax.annotate("", xy=(ha_value, i), xytext=(beta_value, i),
                    arrowprops={"arrowstyle": "-|>", "color": COLOR_CONNECTOR,
                                "lw": LW_CONNECTOR, "mutation_scale": 4.8,
                                "shrinkA": 1.2, "shrinkB": 1.2}, zorder=3)
        ax.scatter(beta_value, i, color=COLOR_EWAS, s=POINT_SIZE,
                   edgecolors=COLOR_AXIS, linewidth=LW_POINT_EDGE, zorder=4)
        ax.scatter(ha_value, i, color=COLOR_HAEWAS, s=POINT_SIZE,
                   edgecolors=COLOR_AXIS, linewidth=LW_POINT_EDGE, zorder=5)

    labels = [
        r"$\mathit{{{}}}$ ({})".format(_escape_mathtext(row.Gene),
                                        _escape_mathtext(row.CpG))
        for row in selected.itertuples(index=False)
    ]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=FS_TICK, color=COLOR_AXIS)
    ax.set_xlabel(r"$-\log_{10}(\mathrm{FDR})$", fontsize=FS_AXIS,
                  fontweight="normal", labelpad=1.5, color=COLOR_AXIS)
    ax.set_xlim(left=x_left, right=x_right)
    ax.set_ylim(len(selected) - 0.5, -0.5)
    fig.text(0.018, 0.972, "d", ha="left", va="top", fontsize=FS_PANEL,
             fontweight="bold", color=COLOR_AXIS)
    return fig
data = pd.read_csv(ROOT / "data" / "Fig4d.csv")
finish(_draw(data))
