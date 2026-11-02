"""Draw Fig4c using only its prepared CSV; export PDF and display inline."""
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(os.environ["REPRO_ROOT"])
PDF = ROOT / "figures" / "Fig4c.pdf"

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

FDR_THRESHOLD = 0.05

FIGURE_WIDTH_MM = 62.0

FIGURE_HEIGHT_MM = 90.0

MM_TO_INCH = 1.0 / 25.4

COLOR_EWAS_ONLY = "#0F4D92"

COLOR_COMMON = "#8F8F8F"

COLOR_HAEWAS_ONLY = "#B64342"

COLOR_AXIS = "#272727"

COLOR_CONNECTOR = "#6F6F6F"

FS_TITLE = 6.0

FS_VENN_SET = 6.0

FS_VENN_NUM = 6.0

FS_AXIS = 6.0

FS_TICK = 5.2

FS_SIG = 6.0

FS_PANEL = 7.0

LW_AXIS = 0.55

LW_BOX = 0.50

LW_MEDIAN = 0.60

LW_ARROW = 0.48

LW_VENN_EDGE = 0.70

def _format_p(p_value):
    if pd.isna(p_value):
        return "NA"
    if p_value < 2.2e-16 or p_value == 0.0:
        return "< 2.2e-16"
    return "{:.2e}".format(p_value)

def _get_star(p_value):
    if pd.isna(p_value):
        return ""
    if p_value < 1e-3:
        return "***"
    if p_value < 1e-2:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""

def _get_whisker_max(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0
    q1, q3 = np.percentile(values, [25, 75])
    upper = q3 + 1.5 * (q3 - q1)
    inside = values[values <= upper]
    return float(np.max(inside)) if inside.size else float(np.max(values))

def _draw_figure(table):
    from matplotlib.patches import FancyArrowPatch
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    ewas_only = table["Pair_class"].eq("EWAS-specific")
    common = table["Pair_class"].eq("Common")
    haewas_only = table["Pair_class"].eq("haEWAS-specific")

    n_ewas_only = int(ewas_only.sum())
    n_common = int(common.sum())
    n_haewas_only = int(haewas_only.sum())
    n_union = n_ewas_only + n_common + n_haewas_only
    if n_union == 0:
        raise ValueError("No significant eQTM pairs in the complete-case table")

    shared = table.loc[common].copy()
    if shared.empty:
        raise ValueError("No Common significant eQTM pairs for the effect-size boxplot")
    beta_abs = shared["Beta_abs"].to_numpy(dtype=float)
    haewas_abs = shared["Haewas_abs"].to_numpy(dtype=float)
    p_value = float(table["Wilcoxon_p"].iloc[0])

    style()
    fig = plt.figure(
        figsize=(FIGURE_WIDTH_MM * MM_TO_INCH, FIGURE_HEIGHT_MM * MM_TO_INCH),
        facecolor="white",
        constrained_layout=False,
    )
    ax_venn = fig.add_axes([0.08, 0.50, 0.84, 0.39])
    ax_box = fig.add_axes([0.25, 0.10, 0.50, 0.35])

    venn = venn2(
        subsets=(n_ewas_only, n_haewas_only, n_common),
        set_labels=("", ""),
        alpha=1.0,
        normalize_to=1.0,
        ax=ax_venn,
    )
    patch_config = {
        "10": (COLOR_EWAS_ONLY, n_ewas_only),
        "01": (COLOR_HAEWAS_ONLY, n_haewas_only),
        "11": (COLOR_COMMON, n_common),
    }
    for region, (color, count) in patch_config.items():
        patch = venn.get_patch_by_id(region)
        if patch is not None:
            patch.set_facecolor(color)
            patch.set_edgecolor("white")
            patch.set_linewidth(LW_VENN_EDGE)
            patch.set_alpha(0.88 if region != "11" else 0.78)
        label = venn.get_label_by_id(region)
        if label is not None:
            percentage = 100.0 * count / n_union
            label.set_text("{:,}\n({:.1f}%)".format(count, percentage))
            label.set_fontsize(FS_VENN_NUM)
            label.set_fontweight("bold")
            label.set_color("black")

    ax_venn.text(-0.42, -0.73, "EWAS", ha="center", va="center",
                 fontsize=FS_VENN_SET, fontweight="bold", color=COLOR_EWAS_ONLY,
                 clip_on=False)
    ax_venn.text(0.42, -0.73, "haEWAS", ha="center", va="center",
                 fontsize=FS_VENN_SET, fontweight="bold", color=COLOR_HAEWAS_ONLY,
                 clip_on=False)
    ax_venn.set_xlim(-1.04, 1.04)
    ax_venn.set_ylim(-0.86, 0.84)
    ax_venn.set_aspect("equal")
    ax_venn.set_axis_off()
    ax_venn.set_title("Significant eQTM pairs", fontsize=FS_TITLE,
                      fontweight="normal", color=COLOR_AXIS, pad=1.0)

    bp = ax_box.boxplot(
        [beta_abs, haewas_abs],
        positions=[1, 2],
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        zorder=3,
        medianprops={"color": COLOR_AXIS, "linewidth": LW_MEDIAN},
        boxprops={"linewidth": LW_BOX, "color": COLOR_AXIS},
        whiskerprops={"linewidth": LW_BOX, "color": COLOR_AXIS},
        capprops={"linewidth": LW_BOX, "color": COLOR_AXIS},
    )
    for patch, color in zip((bp["boxes"][0], bp["boxes"][1]),
                            (COLOR_EWAS_ONLY, COLOR_HAEWAS_ONLY)):
        patch.set_facecolor(color)
        patch.set_alpha(0.82)
        patch.set_edgecolor(COLOR_AXIS)
        patch.set_linewidth(LW_BOX)

    for spine_name, spine in ax_box.spines.items():
        spine.set_visible(spine_name in ("left", "bottom"))
        spine.set_linewidth(LW_AXIS)
        spine.set_color(COLOR_AXIS)
    ax_box.grid(False)
    ax_box.set_facecolor("white")
    ax_box.tick_params(axis="both", labelsize=FS_TICK, width=LW_AXIS,
                       length=2.0, pad=1.0, colors=COLOR_AXIS,
                       top=False, right=False)
    ax_box.set_xticks([1, 2])
    ax_box.set_xticklabels([r"$\bar{m}$", r"$h_{m}$"], fontsize=FS_TICK,
                           fontweight="normal", rotation=0, ha="center")
    ax_box.set_ylabel("|Effect size of eQTM|", fontsize=FS_AXIS,
                      fontweight="normal", labelpad=1.5, color=COLOR_AXIS)

    y_max = max(_get_whisker_max(beta_abs), _get_whisker_max(haewas_abs))
    if y_max <= 0:
        y_max = max(float(np.nanmax(beta_abs)), float(np.nanmax(haewas_abs)))
    y_sig = y_max * 1.08
    y_bar_top = y_sig + y_max * 0.035
    ax_box.plot([1, 1, 2, 2], [y_sig, y_bar_top, y_bar_top, y_sig],
                color=COLOR_AXIS, linewidth=LW_BOX)
    ax_box.text(1.5, y_bar_top + y_max * 0.025, _get_star(p_value),
                ha="center", va="bottom", fontsize=FS_SIG,
                fontweight="bold", color=COLOR_AXIS)
    y_upper = max(1.0, y_bar_top + y_max * 0.12)
    ax_box.set_ylim(0, y_upper)
    ax_box.set_yticks([0, 0.5, 1.0])
    ax_box.set_yticklabels(["0", "0.5", "1"])

    fig.canvas.draw()
    common_label = venn.get_label_by_id("11")
    if common_label is not None:
        common_x, common_y = common_label.get_position()
        connector_start = fig.transFigure.inverted().transform(
            ax_venn.transData.transform((common_x, common_y - 0.20))
        )
    else:
        connector_start = fig.transFigure.inverted().transform(
            ax_venn.transAxes.transform((0.50, 0.25))
        )
    connector_end = fig.transFigure.inverted().transform(
        ax_box.transAxes.transform((0.50, 1.015))
    )
    fig.add_artist(FancyArrowPatch(
        posA=connector_start,
        posB=connector_end,
        transform=fig.transFigure,
        arrowstyle="-|>",
        connectionstyle="arc3,rad=0",
        mutation_scale=5.4,
        linewidth=LW_ARROW,
        color=COLOR_CONNECTOR,
        shrinkA=1.0,
        shrinkB=1.0,
        zorder=10,
    ))

    return fig

data = pd.read_csv(ROOT / "data" / "Fig4c.csv")
finish(_draw_figure(data))
