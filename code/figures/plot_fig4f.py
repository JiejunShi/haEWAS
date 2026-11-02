"""Draw Fig4f using only its prepared CSV; export PDF and display inline."""
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(os.environ["REPRO_ROOT"])
PDF = ROOT / "figures" / "Fig4f.pdf"

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
from matplotlib.lines import Line2D
MODELS = ("Baseline", "EWAS", "haEWAS", "Combined_component_joint")

MODEL_LABELS = {
    "Baseline": "Baseline",
    "EWAS": "EWAS",
    "haEWAS": "haEWAS",
    "Combined_component_joint": "Combined",
}

COLORS = {
    "Baseline": "#777777",
    "EWAS": "#0F4D92",
    "haEWAS": "#B64342",
    "Combined_component_joint": "#6A3D9A",
}
frame = pd.read_csv(ROOT / "data" / "Fig4f.csv")
style()
# Restore the original Fig. 4f typography and axes settings. Keep the
# pooled-BH prediction input; only the visual formatting follows the paper.
plt.rcParams.update({
    "font.size": 5.5,
    "axes.labelsize": 5.8,
    "xtick.labelsize": 5.2,
    "ytick.labelsize": 5.2,
    "legend.fontsize": 5.2,
    "axes.linewidth": 0.35,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.unicode_minus": False,
})
fig, ax = plt.subplots(figsize=(45 / 25.4, 42 / 25.4))
fig.subplots_adjust(left=.18, right=.98, top=.90, bottom=.18)
ax.plot(
    [0, 1], [0, 1], color="#B8B8B8", linewidth=.45,
    linestyle=(0, (3, 2)), zorder=1,
)
handles = []
for model in MODELS:
    part = frame.loc[frame["Model"].eq(model)]
    fpr, tpr = part["Fpr"].to_numpy(float), part["Tpr"].to_numpy(float)
    auc = float(part["Auc"].iloc[0])
    ax.plot(fpr, tpr, color=COLORS[model], linewidth=1.0, zorder=3)
    handles.append(Line2D(
        [0], [0], color=COLORS[model], linewidth=1.0,
        label="{} (AUC = {:.3f})".format(MODEL_LABELS[model], auc),
    ))

ax.set_xticks([0, 0.25, 0.50, 0.75, 1])
ax.set_yticks([0, 0.25, 0.50, 0.75, 1])
ax.set_xlabel("1 - specificity")
ax.set_ylabel("Sensitivity")
ax.set_title("London blood", fontsize=5.8, pad=2.5)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.35)
    spine.set_color("#222222")
ax.tick_params(
    axis="both", width=0.35, length=1.3, pad=1.5, colors="#222222",
)
ax.legend(
    handles=handles, loc="lower right", frameon=False,
    handlelength=1.5, handletextpad=0.45, borderpad=0.1,
)

finish(fig)
