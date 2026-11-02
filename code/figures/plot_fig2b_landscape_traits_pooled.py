import logging
import os
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None


class _HideAdjustTextVersionLog(logging.Filter):
    def filter(self, record):
        return "Looks like you are using an old matplotlib version." not in record.getMessage()


logging.getLogger("adjustText").addFilter(_HideAdjustTextVersionLog())


root = Path(os.environ["REPRO_ROOT"])
input_file = root / "data" / "Fig2b.csv"
output_file = root / "figures" / "Fig2b.pdf"

panel_w_mm = 55.7
panel_h_mm = 66.0
mm_to_inch = 1 / 25.4
max_log_count = 5
category_order = [
    "Neurological", "Psychiatric", "Immune system disease", "Infection", "Respiratory",
    "Cancer", "Metabolic disease", "Cardiovascular disease", "Body weights and measures",
    "Reproduction", "Nutritional", "Environment", "Other trait",
]
category_colors = {
    "Neurological": "#4C78A8", "Psychiatric": "#72A7C7",
    "Immune system disease": "#B64342", "Infection": "#E6A0A2",
    "Respiratory": "#D6A23A", "Cancer": "#8E6A9E",
    "Metabolic disease": "#D9822B", "Cardiovascular disease": "#9C5A36",
    "Body weights and measures": "#5E9F5A", "Reproduction": "#4BA3C7",
    "Nutritional": "#8C6D62", "Environment": "#B8A15A", "Other trait": "#9D9D9D",
}
markers = {"Binary": "o", "Continuous": "^"}
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
category_map = {
    "Adolescent Victimization": "Environment",
    "Alcohol Use in Pregnancy": "Environment",
    "Allergy": "Immune system disease",
    "Arthritis": "Immune system disease",
    "Assisted Reproductive Technologies(ART)": "Reproduction",
    "Asthma": "Respiratory",
    "Atherosclerosis": "Cardiovascular disease",
    "Behavior Score Index": "Psychiatric",
    "Birth Mode": "Reproduction",
    "Birth Weight(Kg)": "Body weights and measures",
    "BMI": "Body weights and measures",
    "Braak Stage": "Neurological",
    "Breast Cancer": "Cancer",
    "Chronic Stress": "Environment",
    "Colorectal Cancer": "Cancer",
    "COVID19": "Infection",
    "Diabetes": "Metabolic disease",
    "Edematous Severe Acute Malnutrition": "Nutritional",
    "First Episode Psychosis": "Psychiatric",
    "Gestational Age(Weeks)": "Reproduction",
    "Gestational Diabetes Mellitus": "Metabolic disease",
    "HCV Status": "Infection",
    "HIV": "Infection",
    "HOMA-IR": "Metabolic disease",
    "Huntington's Disease": "Neurological",
    "Idiopathic Pulmonary Fibrosis": "Respiratory",
    "Inflammatory Bowel Disease": "Immune system disease",
    "Intellectual Disability": "Psychiatric",
    "Liver Cancer": "Cancer",
    "Parity": "Reproduction",
    "Parkinson's Disease": "Neurological",
    "Postmenstrual Age(Weeks)": "Reproduction",
    "Preeclampsia": "Cardiovascular disease",
    "PSA Amount": "Cancer",
    "Psoriasis": "Immune system disease",
    "Schizophrenia": "Psychiatric",
    "SLEDAI Score": "Immune system disease",
    "Smoking Status": "Environment",
    "Sporadic Creutzfeldt-Jakob Disease": "Neurological",
    "Systemic Lupus Erythematosus": "Immune system disease",
    "Total PBB Exposure": "Environment",
    "Ulcerative Colitis": "Immune system disease",
    "Viral Load(log10)": "Infection",
    "Waist Circumference(CM)": "Body weights and measures",
    "War Trauma": "Environment",
    "Zygosity": "Other trait",
}

data = pd.read_csv(input_file)
columns = ["Phenotype", "EWAS-specific loci", "HaEWAS-specific loci"]
if list(data.columns) != columns or len(data) != 46:
    raise ValueError("Fig2b.csv must contain the 46 plotted traits and three count columns")
if data["Phenotype"].duplicated().any():
    raise ValueError("Fig2b.csv contains duplicate phenotypes")

data = data.rename(columns={
    "Phenotype": "phenotype",
    "EWAS-specific loci": "EWAS_count",
    "HaEWAS-specific loci": "haEWAS_count",
})
for column in ("EWAS_count", "haEWAS_count"):
    data[column] = pd.to_numeric(data[column], errors="coerce")
    values = data[column].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Fig2b.csv contains an invalid locus count")
    data[column] = data[column].astype(int)

data["Type"] = data["phenotype"].map(type_map)
data["Category"] = data["phenotype"].map(category_map)
if data["Type"].isna().any() or data["Category"].isna().any():
    raise ValueError("A phenotype is missing its plot category")
data["Category"] = pd.Categorical(
    data["Category"], categories=category_order, ordered=True
)

data["x_log_raw"] = np.log10(data["EWAS_count"] + 1)
data["y_log_raw"] = np.log10(data["haEWAS_count"] + 1)
data["x_plot"] = data["x_log_raw"].astype(float)
data["y_plot"] = data["y_log_raw"].astype(float)

for _, indices in data.groupby(
    ["x_log_raw", "y_log_raw"], sort=False, dropna=False
).groups.items():
    if len(indices) <= 1:
        continue
    ordered = data.loc[list(indices)].sort_values(
        ["Category", "Type", "phenotype"], kind="mergesort"
    ).index.tolist()
    radius = min(0.075, 0.030 + 0.008 * (len(ordered) - 1))
    angles = np.linspace(0, 2 * np.pi, len(ordered), endpoint=False) + np.pi / 4
    for position, index in enumerate(ordered):
        data.at[index, "x_plot"] += radius * np.cos(angles[position])
        data.at[index, "y_plot"] += radius * np.sin(angles[position])

total = data["EWAS_count"] + data["haEWAS_count"]
share = data["haEWAS_count"] / total.replace(0, np.nan)
max_total = max(float(total.max()), 1.0)
data["label_score"] = (
    0.72 * share.fillna(0)
    + 0.28 * np.log10(total + 1) / np.log10(max_total + 1)
)
data["label_rank"] = 0
label_indices = (
    data[total > 0]
    .assign(_score=data.loc[total > 0, "label_score"])
    .sort_values("_score", ascending=False)
    .head(10)
    .index
)
data.loc[label_indices, "label_rank"] = range(1, len(label_indices) + 1)
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
    "xtick.labelsize": 5.1,
    "ytick.labelsize": 5.1,
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


def add_annotation(ax, row, xy, xytext, ha, va, coords):
    annotation = ax.annotate(
        row["phenotype"],
        xy=xy,
        xytext=xytext,
        textcoords=coords,
        fontsize=5.0,
        ha=ha,
        va=va,
        color=category_colors[str(row["Category"])],
        zorder=6,
        clip_on=False,
        arrowprops=dict(
            arrowstyle="-",
            color="#8A8A8A",
            lw=0.24,
            shrinkA=0,
            shrinkB=0,
            connectionstyle="arc3,rad=0",
        ),
    )
    annotation.set_path_effects([
        pe.Stroke(linewidth=0.85, foreground="white"),
        pe.Normal(),
    ])


def annotate(ax, labels):
    if adjust_text is None:
        for _, row in labels.iterrows():
            above = row["y_log_raw"] >= row["x_log_raw"]
            add_annotation(
                ax,
                row,
                (row["x_plot"], row["y_plot"]),
                (6, 6 if above else -6),
                "left",
                "bottom" if above else "top",
                "offset points",
            )
        return

    upper_offsets = [
        (0.18, 0.26), (-0.18, 0.26), (0.28, 0.16), (-0.28, 0.16),
        (0.06, 0.34), (-0.06, 0.34), (0.34, 0.08), (-0.34, 0.08),
    ]
    lower_offsets = [
        (0.18, -0.26), (-0.18, -0.26), (0.28, -0.16), (-0.28, -0.16),
        (0.06, -0.34), (-0.06, -0.34), (0.34, -0.08), (-0.34, -0.08),
    ]
    temporary = []
    points = []
    upper_i = lower_i = 0

    for _, row in labels.iterrows():
        x, y = row["x_plot"], row["y_plot"]
        points.append((x, y))
        if row["y_log_raw"] >= row["x_log_raw"]:
            dx, dy = upper_offsets[upper_i % len(upper_offsets)]
            upper_i += 1
            va = "bottom"
        else:
            dx, dy = lower_offsets[lower_i % len(lower_offsets)]
            lower_i += 1
            va = "top"
        temporary.append(ax.text(
            x + dx,
            y + dy,
            row["phenotype"],
            fontsize=5.0,
            ha="left" if dx >= 0 else "right",
            va=va,
            color=category_colors[str(row["Category"])],
            zorder=6,
            clip_on=False,
            alpha=0,
        ))

    adjust_text(
        temporary,
        x=[point[0] for point in points],
        y=[point[1] for point in points],
        ax=ax,
        expand=(1.12, 1.22),
        force_text=(0.55, 0.70),
        force_static=(0.38, 0.50),
        max_move=(28, 28),
        only_move={"text": "xy", "static": "xy", "explode": "xy", "pull": "xy"},
        ensure_inside_axes=True,
        expand_axes=False,
        iter_lim=800,
    )
    locations = [
        (text.get_position(), text.get_ha(), text.get_va()) for text in temporary
    ]
    for text in temporary:
        text.remove()
    for (_, row), point, location in zip(labels.iterrows(), points, locations):
        (tx, ty), ha, va = location
        add_annotation(ax, row, point, (tx, ty), ha, va, "data")


fig = plt.figure(
    figsize=(panel_w_mm * mm_to_inch, panel_h_mm * mm_to_inch),
    facecolor="white",
)
ax = fig.add_axes(rect_mm(8.8, 28.5, 28.5, 28.5))

for _, row in data.iterrows():
    ax.scatter(
        row["x_plot"],
        row["y_plot"],
        s=10.5,
        marker=markers[row["Type"]],
        color=category_colors[str(row["Category"])],
        edgecolor="#272727",
        linewidth=0.25,
        alpha=0.88,
        zorder=3,
    )

ax.plot(
    [0, max_log_count],
    [0, max_log_count],
    linestyle=(0, (2.2, 2.2)),
    linewidth=0.5,
    color="#707070",
    zorder=1,
)
ax.set_xlim(-0.20, max_log_count + 0.45)
ax.set_ylim(-0.20, max_log_count + 0.45)
ax.set_aspect("equal", adjustable="box")
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_xticks(range(6))
ax.set_yticks(range(6))
ax.set_xlabel("EWAS-specific loci (log scale)", labelpad=1.4)
ax.set_ylabel("haEWAS-specific (log scale)", labelpad=1.4)
ax.tick_params(
    axis="both", which="major", length=2.0, width=0.5, pad=1.2, colors="#272727"
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.55)
    spine.set_color("#272727")

labels = data[data["label_rank"] > 0].sort_values("label_rank")
annotate(ax, labels)
fig.text(
    0.018, 0.972, "b", ha="left", va="top",
    fontsize=7.0, fontweight="bold", color="#272727",
)

type_handles = [
    Line2D(
        [0], [0], marker="o", linestyle="none", markerfacecolor="white",
        markeredgecolor="#272727", markeredgewidth=0.55, markersize=3.1,
        label="Binary",
    ),
    Line2D(
        [0], [0], marker="^", linestyle="none", markerfacecolor="white",
        markeredgecolor="#272727", markeredgewidth=0.55, markersize=3.5,
        label="Continuous",
    ),
]
fig.legend(
    handles=type_handles, loc="lower center", bbox_to_anchor=(0.50, 0.005),
    ncol=2, frameon=False, fontsize=5.0, handlelength=0.8,
    handletextpad=0.25, columnspacing=0.9, borderaxespad=0.0,
)
category_handles = [
    Line2D(
        [0], [0], marker="o", linestyle="none",
        markerfacecolor=category_colors[category], markeredgecolor="#272727",
        markeredgewidth=0.22, markersize=2.7, label=category,
    )
    for category in category_order
]
fig.legend(
    handles=category_handles, loc="lower center", bbox_to_anchor=(0.50, 0.092),
    ncol=2, frameon=False, fontsize=5.0, handlelength=0.7,
    handletextpad=0.20, columnspacing=0.55, labelspacing=0.13,
    borderaxespad=0.0,
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







