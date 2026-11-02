import os
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
import seaborn as sns


root = Path(os.environ["REPRO_ROOT"])
input_file = root / "data" / "Fig2a.csv"
output_file = root / "figures" / "Fig2a.pdf"

panel_w_mm = 55.7
panel_h_mm = 66.0
mm_to_inch = 1 / 25.4
max_logp = 30
groups = ["haEWAS-specific", "EWAS-specific", "Common"]
colors = {
    "haEWAS-specific": "#B64342",
    "EWAS-specific": "#0F4D92",
    "Common": "#8F8F8F",
}
alphas = {"haEWAS-specific": 0.72, "EWAS-specific": 0.68, "Common": 0.38}
zorder = {"Common": 1, "EWAS-specific": 2, "haEWAS-specific": 3}

data = pd.read_csv(input_file)
required = {"P_Beta", "P_haEWAS", "Group"}
if not required.issubset(data.columns):
    raise ValueError("The Fig2a input table is missing required columns")

for column in ("P_Beta", "P_haEWAS"):
    data[column] = pd.to_numeric(data[column], errors="coerce")
    values = data[column].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values <= 0).any() or (values > 1).any():
        raise ValueError("The Fig2a input table has invalid P values")
if not data["Group"].isin(groups).all():
    raise ValueError("The Fig2a input table has an unknown group label")

font_file = root.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
if font_file.is_file():
    fm.fontManager.addfont(str(font_file))
    font_family = fm.FontProperties(fname=str(font_file)).get_name()
else:
    font_family = "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_family,
    "font.size": 5.6,
    "axes.labelsize": 6.0,
    "axes.titlesize": 6.0,
    "xtick.labelsize": 5.2,
    "ytick.labelsize": 5.2,
    "legend.fontsize": 5.1,
    "axes.linewidth": 0.55,
    "lines.linewidth": 0.7,
    "patch.linewidth": 0.5,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

data["logP_Beta"] = -np.log10(data["P_Beta"])
data["logP_haEWAS"] = -np.log10(data["P_haEWAS"])
data["Group"] = pd.Categorical(data["Group"], categories=groups, ordered=True)
data = data.sort_values("Group")
counts = data["Group"].value_counts().reindex(groups).fillna(0).astype(int)


def rect_mm(left, bottom, width, height):
    return [
        left / panel_w_mm,
        bottom / panel_h_mm,
        width / panel_w_mm,
        height / panel_h_mm,
    ]


fig = plt.figure(
    figsize=(panel_w_mm * mm_to_inch, panel_h_mm * mm_to_inch),
    facecolor="white",
)
main_left, main_bottom, main_size = 9.9, 11.0, 38.0
gap, top_height, right_width = 1.2, 6.5, 5.8
ax_main = fig.add_axes(rect_mm(main_left, main_bottom, main_size, main_size))
ax_top = fig.add_axes(
    rect_mm(main_left, main_bottom + main_size + gap, main_size, top_height),
    sharex=ax_main,
)
ax_right = fig.add_axes(
    rect_mm(main_left + main_size + gap, main_bottom, right_width, main_size),
    sharey=ax_main,
)

for group in ("Common", "EWAS-specific", "haEWAS-specific"):
    subset = data[data["Group"] == group]
    if not subset.empty:
        ax_main.scatter(
            subset["logP_Beta"],
            subset["logP_haEWAS"],
            c=colors[group],
            s=1.0,
            alpha=alphas[group],
            edgecolors="none",
            linewidths=0,
            rasterized=True,
            zorder=zorder[group],
        )

ax_main.plot(
    [0, max_logp],
    [0, max_logp],
    color="#606060",
    linestyle=(0, (2.2, 2.2)),
    linewidth=0.5,
    zorder=10,
)
ax_main.set_xlim(0, max_logp)
ax_main.set_ylim(0, max_logp)
ax_main.set_aspect("equal", adjustable="box")
ax_main.xaxis.set_major_locator(MultipleLocator(10))
ax_main.yaxis.set_major_locator(MultipleLocator(10))
ax_main.set_xlabel(r"EWAS $-\log_{10}(P)$", labelpad=1.8)
ax_main.set_ylabel(r"haEWAS $-\log_{10}(P)$", labelpad=1.8)
ax_main.tick_params(
    axis="both", which="major", length=2.0, width=0.5, pad=1.4, colors="#272727"
)
for spine in ax_main.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.55)
    spine.set_color("#272727")

for group in groups:
    subset = data[data["Group"] == group]
    alpha = 0.18 if group == "Common" else 0.22
    if len(subset) >= 3 and subset["logP_Beta"].nunique() > 1:
        sns.kdeplot(
            data=subset,
            x="logP_Beta",
            ax=ax_top,
            color=colors[group],
            fill=True,
            alpha=alpha,
            linewidth=0.55,
            bw_adjust=0.85,
            cut=0,
            clip=(0, max_logp),
            warn_singular=False,
            zorder=zorder[group],
        )
    if len(subset) >= 3 and subset["logP_haEWAS"].nunique() > 1:
        sns.kdeplot(
            data=subset,
            y="logP_haEWAS",
            ax=ax_right,
            color=colors[group],
            fill=True,
            alpha=alpha,
            linewidth=0.55,
            bw_adjust=0.85,
            cut=0,
            clip=(0, max_logp),
            warn_singular=False,
            zorder=zorder[group],
        )

for axis in (ax_top, ax_right):
    axis.set_facecolor("white")
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for spine in axis.spines.values():
        spine.set_visible(False)
ax_top.set_xlim(0, max_logp)
ax_right.set_ylim(0, max_logp)

handles = []
for group, label in (
    ("haEWAS-specific", "haEWAS-spec."),
    ("EWAS-specific", "EWAS-spec."),
    ("Common", "Common"),
):
    handles.append(Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor=colors[group],
        markeredgecolor="none",
        markersize=3.0,
        alpha=0.72 if group != "Common" else 0.65,
        label="{} (n={:,})".format(label, counts[group]),
    ))
fig.legend(
    handles=handles,
    loc="upper center",
    bbox_to_anchor=(0.57, 0.956),
    ncol=1,
    frameon=False,
    handlelength=0.8,
    handletextpad=0.22,
    borderaxespad=0.0,
    columnspacing=0.62,
    fontsize=5.0,
)

output_file.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output_file, format="pdf", facecolor="white")


try:
    from IPython import get_ipython
    from IPython.display import display

    shell = get_ipython()
    if shell is not None and hasattr(shell, "kernel"):
        display(fig)
except ImportError:
    pass

plt.close(fig)



