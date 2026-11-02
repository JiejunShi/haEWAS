import os
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


root = Path(os.environ["REPRO_ROOT"])
input_file = root / "data" / "Fig2c.csv"
output_file = root / "figures" / "Fig2c.pdf"

panel_w_mm = 55.7
panel_h_mm = 66.0
mm_to_inch = 1 / 25.4
category_order = [
    "Neurological", "Psychiatric", "Immune system disease", "Infection",
    "Respiratory", "Cancer", "Metabolic disease", "Cardiovascular disease",
    "Body weights and measures", "Reproduction", "Nutritional", "Environment",
    "Other trait",
]
category_colors = {
    "Neurological": "#4C78A8", "Psychiatric": "#72A7C7",
    "Immune system disease": "#B64342", "Infection": "#E6A0A2",
    "Respiratory": "#D6A23A", "Cancer": "#8E6A9E",
    "Metabolic disease": "#D9822B", "Cardiovascular disease": "#9C5A36",
    "Body weights and measures": "#5E9F5A", "Reproduction": "#4BA3C7",
    "Nutritional": "#8C6D62", "Environment": "#B8A15A",
    "Other trait": "#9D9D9D",
}
category_map = {
    "Huntington's Disease": "Neurological",
    "Parkinson's Disease": "Neurological",
    "Sporadic Creutzfeldt-Jakob Disease": "Neurological",
    "Braak Stage": "Neurological",
    "Schizophrenia": "Psychiatric",
    "First Episode Psychosis": "Psychiatric",
    "Intellectual Disability": "Psychiatric",
    "Behavior Score Index": "Psychiatric",
    "Ulcerative Colitis": "Immune system disease",
    "Inflammatory Bowel Disease": "Immune system disease",
    "Arthritis": "Immune system disease",
    "Psoriasis": "Immune system disease",
    "Allergy": "Immune system disease",
    "Systemic Lupus Erythematosus": "Immune system disease",
    "SLEDAI Score": "Immune system disease",
    "COVID19": "Infection",
    "HCV Status": "Infection",
    "HIV": "Infection",
    "Viral Load(log10)": "Infection",
    "Asthma": "Respiratory",
    "Idiopathic Pulmonary Fibrosis": "Respiratory",
    "Breast Cancer": "Cancer",
    "Colorectal Cancer": "Cancer",
    "Liver Cancer": "Cancer",
    "PSA Amount": "Cancer",
    "HOMA-IR": "Metabolic disease",
    "Diabetes": "Metabolic disease",
    "Gestational Diabetes Mellitus": "Metabolic disease",
    "Edematous Severe Acute Malnutrition": "Nutritional",
    "Atherosclerosis": "Cardiovascular disease",
    "Preeclampsia": "Cardiovascular disease",
    "BMI": "Body weights and measures",
    "Waist Circumference(CM)": "Body weights and measures",
    "Birth Weight(Kg)": "Body weights and measures",
    "Gestational Age(Weeks)": "Reproduction",
    "Postmenstrual Age(Weeks)": "Reproduction",
    "Assisted Reproductive Technologies(ART)": "Reproduction",
    "Birth Mode": "Reproduction",
    "Parity": "Reproduction",
    "Alcohol Use in Pregnancy": "Environment",
    "Smoking Status": "Environment",
    "Total PBB Exposure": "Environment",
    "Chronic Stress": "Environment",
    "War Trauma": "Environment",
    "Adolescent Victimization": "Environment",
    "Zygosity": "Other trait",
}
type_map = {
    "Adolescent Victimization": "Binary",
    "Alcohol Use in Pregnancy": "Binary",
    "Allergy": "Binary",
    "Arthritis": "Binary",
    "Assisted Reproductive Technologies(ART)": "Binary",
    "Asthma": "Binary",
    "Atherosclerosis": "Binary",
    "Behavior Score Index": "Continuous",
    "Birth Mode": "Binary",
    "Birth Weight(Kg)": "Continuous",
    "BMI": "Continuous",
    "Braak Stage": "Continuous",
    "Breast Cancer": "Binary",
    "Chronic Stress": "Continuous",
    "Colorectal Cancer": "Binary",
    "COVID19": "Binary",
    "Diabetes": "Binary",
    "Edematous Severe Acute Malnutrition": "Binary",
    "First Episode Psychosis": "Binary",
    "Gestational Age(Weeks)": "Continuous",
    "Gestational Diabetes Mellitus": "Binary",
    "HCV Status": "Binary",
    "HIV": "Binary",
    "HOMA-IR": "Continuous",
    "Huntington's Disease": "Binary",
    "Idiopathic Pulmonary Fibrosis": "Binary",
    "Inflammatory Bowel Disease": "Binary",
    "Intellectual Disability": "Binary",
    "Liver Cancer": "Binary",
    "Parity": "Binary",
    "Parkinson's Disease": "Binary",
    "Postmenstrual Age(Weeks)": "Continuous",
    "Preeclampsia": "Binary",
    "PSA Amount": "Continuous",
    "Psoriasis": "Binary",
    "Schizophrenia": "Binary",
    "SLEDAI Score": "Continuous",
    "Smoking Status": "Binary",
    "Sporadic Creutzfeldt-Jakob Disease": "Binary",
    "Systemic Lupus Erythematosus": "Binary",
    "Total PBB Exposure": "Continuous",
    "Ulcerative Colitis": "Binary",
    "Viral Load(log10)": "Continuous",
    "Waist Circumference(CM)": "Continuous",
    "War Trauma": "Continuous",
    "Zygosity": "Binary",
}
markers = {"Binary": "o", "Continuous": "^"}
side_colors = {"EWAS-favored": "#0F4D92", "haEWAS-favored": "#B64342"}
background_colors = {"EWAS-favored": "#D8E6F5", "haEWAS-favored": "#F4D7D7"}

data = pd.read_csv(input_file)
columns = ["Phenotype", "EWAS-favored", "HaEWAS-favored"]
if list(data.columns) != columns or data.empty:
    raise ValueError("Fig2c.csv must contain 46 phenotypes and the two favored columns")
if data["Phenotype"].duplicated().any():
    raise ValueError("Fig2c.csv contains duplicate phenotypes")

for column in ("EWAS-favored", "HaEWAS-favored"):
    data[column] = pd.to_numeric(data[column], errors="coerce")
    values = data[column].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Fig2c.csv contains an invalid favored score")
if ((data["EWAS-favored"] > 0) & (data["HaEWAS-favored"] > 0)).any():
    raise ValueError("A phenotype cannot favor both methods at once")

data["haEWAS_gain_score"] = data["HaEWAS-favored"] - data["EWAS-favored"]
data["Category"] = data["Phenotype"].map(category_map)
data["Type"] = data["Phenotype"].map(type_map)
if data["Category"].isna().any() or data["Type"].isna().any():
    raise ValueError("A phenotype is missing its plot classification")

traits = data.copy()
if traits.empty:
    raise ValueError("No nonzero gain scores were found")
traits["Category"] = pd.Categorical(
    traits["Category"], categories=category_order, ordered=True
)
order = (
    traits.groupby("Category", observed=True)["haEWAS_gain_score"]
    .median()
    .sort_values(ascending=False)
    .index.astype(str)
    .tolist()
)

font_file = root.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
if font_file.is_file():
    fm.fontManager.addfont(str(font_file))
    font_family = fm.FontProperties(fname=str(font_file)).get_name()
else:
    font_family = "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_family,
    "font.size": 5.5,
    "axes.labelsize": 5.8,
    "axes.titlesize": 6.0,
    "xtick.labelsize": 5.0,
    "ytick.labelsize": 5.0,
    "legend.fontsize": 5.0,
    "axes.linewidth": 0.55,
    "lines.linewidth": 0.7,
    "patch.linewidth": 0.45,
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
ax = fig.add_axes(rect_mm(3.4, 12.4, 27.0, 46.0))
rng = np.random.default_rng(20250304)
x_min = np.floor(traits["haEWAS_gain_score"].min() - 0.5)
x_max = np.ceil(traits["haEWAS_gain_score"].max() + 0.5)
limit = max(abs(x_min), abs(x_max), 2)
x_min, x_max = -limit, limit
ax.axvspan(x_min, 0, color=background_colors["EWAS-favored"], zorder=0)
ax.axvspan(0, x_max, color=background_colors["haEWAS-favored"], zorder=0)
ax.axvline(
    0, color="#272727", linestyle=(0, (2.2, 2.2)), linewidth=0.55, zorder=1
)

for row_index, category in enumerate(order):
    subset = traits[traits["Category"].astype(str).eq(category)]
    values = subset["haEWAS_gain_score"].dropna().to_numpy(dtype=float)
    if not len(values):
        continue
    color = category_colors[category]
    box = ax.boxplot(
        values,
        positions=[row_index],
        vert=False,
        widths=0.46,
        patch_artist=True,
        showfliers=False,
        zorder=2,
        medianprops={"color": "#222222", "linewidth": 0.65},
        whiskerprops={"color": "#444444", "linewidth": 0.50},
        capprops={"color": "#444444", "linewidth": 0.50},
    )
    for patch in box["boxes"]:
        patch.set_facecolor(mcolors.to_rgba(color, alpha=0.20))
        patch.set_edgecolor("#444444")
        patch.set_linewidth(0.50)
    for trait_type in ("Binary", "Continuous"):
        part = subset[subset["Type"].eq(trait_type)]
        if part.empty:
            continue
        jitter = rng.uniform(-0.13, 0.13, size=len(part))
        ax.scatter(
            part["haEWAS_gain_score"],
            np.full(len(part), row_index) + jitter,
            s=9.5,
            marker=markers[trait_type],
            color=color,
            edgecolor="#272727",
            linewidth=0.22,
            alpha=0.88,
            zorder=3,
        )

ax.set_yticks(np.arange(len(order)))
ax.set_yticklabels(order, fontsize=5.0)
ax.invert_yaxis()
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.tick_params(
    axis="y", labelleft=False, labelright=True, left=False,
    right=True, length=2.0, width=0.5, pad=1.0,
)
for tick, category in zip(ax.get_yticklabels(), order):
    tick.set_clip_on(False)
    tick.set_horizontalalignment("left")
    tick.set_color(category_colors[category])
ax.set_xlim(x_min, x_max)
ax.xaxis.set_major_locator(MultipleLocator(2 if limit <= 4 else 3))
ax.tick_params(
    axis="x", which="major", length=2.0, width=0.5,
    pad=1.4, labelsize=5.0, colors="#272727",
)
ax.set_xlabel(
    "Gain score",
    fontsize=5.7,
    labelpad=1.8,
)
ax.text(
    x_min, 1.02, "EWAS-favored", transform=ax.get_xaxis_transform(),
    ha="left", va="bottom", fontsize=5.0,
    color=side_colors["EWAS-favored"], fontweight="bold",
)
ax.text(
    x_max, 1.02, "haEWAS-favored", transform=ax.get_xaxis_transform(),
    ha="right", va="bottom", fontsize=5.0,
    color=side_colors["haEWAS-favored"], fontweight="bold",
)
fig.text(
    0.018, 0.972, "c", ha="left", va="top",
    fontsize=7.0, fontweight="bold", color="#272727",
)

handles = [
    Line2D(
        [0], [0], marker=markers[name], linestyle="none",
        markerfacecolor="white", markeredgecolor="#272727",
        markeredgewidth=0.55, markersize=3.2 if name == "Binary" else 3.6,
        label=name,
    )
    for name in ("Binary", "Continuous")
]
fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.28, 0.026),
    ncol=2,
    frameon=False,
    fontsize=5.0,
    handlelength=0.8,
    handletextpad=0.25,
    columnspacing=0.8,
    borderaxespad=0.0,
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.55)
    spine.set_color("#272727")
ax.grid(axis="x", linestyle=(0, (2.0, 2.0)), linewidth=0.32, alpha=0.22)
ax.set_axisbelow(True)

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





