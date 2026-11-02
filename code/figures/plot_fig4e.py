"""Draw Fig4e using only its prepared CSV; export PDF and display inline."""
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(os.environ["REPRO_ROOT"])
PDF = ROOT / "figures" / "Fig4e.pdf"

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

import matplotlib.pyplot as plt
TARGETS = [
    {"cpg": "cg10091752", "display_gene": "SLC22A23", "method": "camda", "order": 0},
    {"cpg": "cg03605404", "display_gene": "CPNE1", "method": "chalm", "order": 1},
]

STAGE_LABELS = {0: "0", 1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI"}

RED = "#B64342"

BLUE = "#0F4D92"

GRAY = "#666666"

Y_LIMITS = {
    ("cg10091752", "hm"): (0.32, 0.352),
    ("cg10091752", "beta"): (0.30, 0.60),
    ("cg03605404", "hm"): (0.065, 0.081),
    ("cg03605404", "beta"): (0.06, 0.10),
}

Y_TICKS = {
    ("cg10091752", "hm"): np.linspace(0.32, 0.35, 4),
    ("cg10091752", "beta"): np.arange(0.30, 0.6001, 0.10),
    ("cg03605404", "hm"): np.array([0.065, 0.070, 0.075, 0.080]),
    ("cg03605404", "beta"): np.arange(0.06, 0.1001, 0.02),
}

def _format_p(p_value):
    if not np.isfinite(p_value):
        return "NA"
    if p_value < 1e-4:
        return "{:.1e}".format(p_value)
    return "{:.3f}".format(p_value)

def _add_box_points(ax, table, cpg, method):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    from scipy import stats

    # CHALM and CAMDA are separate data sources, but both use the h_m axis
    # scale. Keep their method labels for audit/statistics and normalize only
    # the key used to look up the figure's y-axis limits and ticks.
    axis_method = "hm" if method in {"chalm", "camda"} else method
    axis_key = (cpg, axis_method)
    if axis_key not in Y_LIMITS or axis_key not in Y_TICKS:
        raise ValueError("No y-axis configuration for CpG/method {}".format(axis_key))

    labels = [STAGE_LABELS[i] for i in range(7)]
    observed_values = pd.to_numeric(table["Observed_value"], errors="coerce")
    value_median = observed_values.median()
    if not np.isfinite(value_median):
        raise ValueError("No observed blood values for {} / {}".format(cpg, method))
    plot_values = observed_values.fillna(value_median)
    groups = [
        plot_values.loc[table["Ordered_stage"].eq(i)].dropna().to_numpy(dtype=float)
        for i in range(7)
    ]
    boxes = ax.boxplot(
        groups, positions=np.arange(7, dtype=float), widths=0.56,
        patch_artist=True, showfliers=False,
        medianprops={"color": GRAY, "linewidth": 0.55},
        boxprops={"facecolor": "#CFCFCF", "edgecolor": "none", "linewidth": 0},
        whiskerprops={"color": "#CFCFCF", "linewidth": 0.45},
        capprops={"color": "#CFCFCF", "linewidth": 0.45},
    )
    for patch in boxes["boxes"]:
        patch.set_alpha(0.62)

    rng = np.random.default_rng(20260723)
    stage_medians = []
    for i, values in enumerate(groups):
        if not len(values):
            stage_medians.append(np.nan)
            continue
        jitter = rng.normal(0.0, 0.05, size=len(values))
        ax.scatter(np.full(len(values), i) + jitter, values, s=4.2,
                   color=GRAY, alpha=0.62, linewidths=0, zorder=3)
        stage_medians.append(float(np.nanmedian(values)))
        ax.text(i, -0.18, "({})".format(len(values)), transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=4.4, clip_on=False)

    stage_medians = np.asarray(stage_medians, dtype=float)
    positions = np.arange(7, dtype=float)
    valid = np.isfinite(stage_medians)
    if valid.sum() >= 2:
        coef = np.polyfit(positions[valid], stage_medians[valid], 1)
        x_smooth = np.linspace(positions[valid].min(), positions[valid].max(), 200)
        line_color = RED if method != "beta" else BLUE
        ax.plot(x_smooth, np.polyval(coef, x_smooth), color=line_color, linewidth=0.95, zorder=5)

    ax.set_xlim(-0.55, 6.45)
    ax.set_ylim(*Y_LIMITS[axis_key])
    ax.margins(y=0)
    ax.set_yticks(Y_TICKS[axis_key])
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=0, ha="center", va="top", fontsize=5.2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.35)
        ax.spines[side].set_color("#333333")
    ax.tick_params(width=0.35, length=1.4, pad=1.0, labelsize=5.2)
    tick_step = float(abs(Y_TICKS[axis_key][1] - Y_TICKS[axis_key][0]))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f" if tick_step < 0.01 else "%.2f"))
    ax.set_box_aspect(1)
    # The original figure annotates Spearman statistics computed from
    # non-imputed observations, while boxplots/points/fit lines use the
    # median-imputed plotting table.
    # Match the original statistic using non-imputed observations.
    observed = table.loc[table["Observed_value"].notna()]
    medians = observed.groupby("Ordered_stage")["Observed_value"].median().reindex(range(7)).dropna()
    if len(medians) < 3 or medians.nunique() < 2:
        return np.nan, np.nan
    result = stats.spearmanr(medians.index.to_numpy(dtype=float), medians.to_numpy(dtype=float))
    return float(result[0]), float(result[1])

def _italic_gene_title(gene, cpg):
    safe_gene = str(gene).replace("_", r"\_")
    return r"$\it{{{}}}$ ({})".format(safe_gene, cpg)
data = pd.read_csv(ROOT / "data" / "Fig4e.csv")
panel_values = []
for target in TARGETS:
    part = data.loc[data["CpG"].eq(target["cpg"])]
    panel_values.append((target, part.loc[part["Plot_method"].eq(target["method"])], part.loc[part["Plot_method"].eq("beta")]))
style()
fig, axes = plt.subplots(2, 2, figsize=(68 / 25.4, 60 / 25.4), squeeze=False)
for col, (target, hm_values, beta_values) in enumerate(panel_values):
    axes[0, col].set_title(_italic_gene_title(target["display_gene"], target["cpg"]),
                           fontsize=6.2, pad=2.0)
    for row, (values, method) in enumerate(((hm_values, target["method"]), (beta_values, "beta"))):
        rho, p_value = _add_box_points(axes[row, col], values, target["cpg"], method)
        axes[row, col].text(
            0.04, 0.05,
            r"$\rho$ = {:.2f}".format(rho) + "\nP = " + _format_p(p_value),
            transform=axes[row, col].transAxes,
            ha="left", va="bottom", fontsize=4.8, linespacing=1.35, zorder=10,
        )

for ax in axes.flat:
    ax.set_xlabel("")
    ax.set_ylabel("")
axes[0, 0].set_ylabel(r"$h_m$", fontsize=6.2)
axes[1, 0].set_ylabel(r"$\bar{m}$", fontsize=6.2)
fig.text(0.51, 0.018, "Braak stage", ha="center", va="bottom", fontsize=6.2)
fig.subplots_adjust(left=0.125, right=0.990, top=0.925, bottom=0.195,
                    wspace=0.36, hspace=0.50)

finish(fig)
