"""Plot Fig3e, Fig3f and Fig3g from their shared prepared ROC table."""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
PANELS = {"GSE110607": "Fig3e", "GSE118144": "Fig3f", "GSE96879": "Fig3g"}
MODELS = ("EWAS", "haEWAS", "Joint_MRS_Logistic")
LABELS = {"EWAS": "EWAS", "haEWAS": "haEWAS", "Joint_MRS_Logistic": "Combined"}
COLORS = {"EWAS": "#0F4D92", "haEWAS": "#B64342", "Joint_MRS_Logistic": "#5A3E85"}
WIDTH, HEIGHT = 39.0, 54.0


def plot_panel(cohort, data):
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
    left, bottom, size = 7.0, 8.1, 30.8
    ax = fig.add_axes([left / WIDTH, bottom / HEIGHT, size / WIDTH, size / HEIGHT])
    ax.plot([0, 1], [0, 1], linestyle="--", color="#9D9D9D", linewidth=.4, zorder=1)
    handles = []
    for model in MODELS:
        curve = data.loc[data["Cohort"].eq(cohort) & data["Model"].eq(model)].sort_values(["FPR", "TPR"])
        auc_values = curve["AUC"].unique()
        if len(curve) < 2 or len(auc_values) != 1:
            raise ValueError(f"Incomplete ROC curve or inconsistent AUC for {cohort}/{model}")
        auc = float(auc_values[0])
        fpr = curve["FPR"].to_numpy(dtype=float)
        tpr = curve["TPR"].to_numpy(dtype=float)
        ax.plot(fpr, tpr, color=COLORS[model], linewidth=1.0,
                solid_capstyle="round", zorder=2)
        handles.append(Line2D([0], [0], color=COLORS[model], linewidth=1.0,
                              label=f"{LABELS[model]} (AUC={auc:.3f})"))
    ax.set_xlim(-.025, 1.025)
    ax.set_ylim(-.025, 1.025)
    ax.set_aspect("equal", adjustable="box")
    ticks = [0, .25, .50, .75, 1.0]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1"])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1"])
    ax.set_xlabel("1 - specificity", labelpad=1.5)
    ax.set_ylabel("Sensitivity", labelpad=1.5)
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(.98, .04),
              frameon=False, fontsize=5.0, handlelength=1.0,
              handletextpad=.35, borderaxespad=0.0, labelspacing=.35)
    ax.tick_params(axis="both", which="major", length=1.4, width=.35,
                   pad=1.0, labelsize=5.0, colors="#272727")
    for spine in ax.spines.values():
        spine.set_linewidth(.4)
        spine.set_color("#272727")
    strip = fig.add_axes([left / WIDTH, (bottom + size + .22) / HEIGHT,
                          size / WIDTH, 3.0 / HEIGHT])
    strip.set_facecolor("#DCDCDC")
    strip.set_xticks([])
    strip.set_yticks([])
    for spine in strip.spines.values():
        spine.set_linewidth(.35)
        spine.set_color("#272727")
    strip.text(.5, .5, cohort, ha="center", va="center", fontsize=6.0,
               fontweight="bold", color="#272727", transform=strip.transAxes)
    output = ROOT / "figures" / (PANELS[cohort] + ".pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="pdf", facecolor="white")
    plt.show()
    plt.close(fig)


def run():
    data = pd.read_csv(ROOT / "data" / "Fig3e_g.csv", float_precision="round_trip")
    if (list(data.columns) != ["Cohort", "Model", "FPR", "TPR", "AUC"]
            or set(zip(data["Cohort"], data["Model"])) != {
                (cohort, model) for cohort in PANELS for model in MODELS}):
        raise ValueError("Fig3e_g.csv has an invalid schema or cohort/model roster")
    if (not np.isfinite(data[["FPR", "TPR", "AUC"]].to_numpy(float)).all()
            or not data[["FPR", "TPR", "AUC"]].apply(lambda column: column.between(0, 1).all()).all()):
        raise ValueError("ROC coordinates and AUCs must be finite and between zero and one")
    font_file = ROOT.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
    font_family = ["Arial", "DejaVu Sans"]
    if font_file.is_file():
        fm.fontManager.addfont(str(font_file))
        font_family = fm.FontProperties(fname=str(font_file)).get_name()
    plt.rcParams.update({
        "font.family": font_family, "font.size": 5.5,
        "axes.labelsize": 5.8, "axes.titlesize": 6.0, "xtick.labelsize": 5.0,
        "ytick.labelsize": 5.0, "legend.fontsize": 5.0, "axes.linewidth": .5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": False,
        "figure.facecolor": "white", "savefig.facecolor": "white", "legend.frameon": False,
    })
    for cohort in PANELS:
        plot_panel(cohort, data)


if __name__ == "__main__":
    run()

