import logging
import os
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke
from matplotlib.ticker import FixedFormatter, FixedLocator, MaxNLocator
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
input_file = root / "data" / "Fig3b.csv"
output_file = root / "figures" / "Fig3b.pdf"

panel_w_mm = 116.0
panel_h_mm = 56.0
mm_to_inch = 1 / 25.4
fs_panel = 6.5
fs_axis = 6.0
fs_tick = 6.0
fs_legend = 6.0
fs_gene = 6.2
y_headroom_factor = 1.30
label_category_order = ["Beta_Shared", "Hetero_Unique", "Beta_Only"]
label_side_order = [1, -1]
legend_map = {
    "Hetero_Unique": "haEWAS-specific",
    "Beta_Only": "EWAS-specific",
    "Beta_Shared": "Common",
}
plot_colors = {
    "Beta_Only": "#0F4D92",
    "Beta_Shared": "#8F8F8F",
    "Hetero_Unique": "#B64342",
    "Background_Odd": "#B8B8B8",
    "Background_Even": "#D1D1D1",
    "Bg_Strip_Odd": "#F6F7F8",
    "Bg_Strip_Even": "#FFFFFF",
    "Axis": "#272727",
}

data = pd.read_csv(
    input_file,
    dtype={"CpG": str, "Gene": str, "Direction": str, "Group": str},
    low_memory=False,
)
columns = ["Chr", "CpG", "Gene", "P", "Direction", "Group", "Pos"]
if list(data.columns) != columns or data.empty:
    raise ValueError("Fig3b.csv must contain genomic points, labels, P values and groups")

for column in ("Chr", "P", "Pos"):
    data[column] = pd.to_numeric(data[column], errors="coerce")
values = data[["Chr", "P", "Pos"]].to_numpy(dtype=float)
if not np.isfinite(values).all():
    raise ValueError("Fig3b.csv contains invalid numeric values")
if (
    (data["Chr"] < 1).any()
    or (data["Chr"] > 22).any()
    or (data["Chr"] != np.floor(data["Chr"])).any()
    or (data["Pos"] < 0).any()
    or (data["P"] <= 0).any()
    or (data["P"] > 1).any()
):
    raise ValueError("Fig3b.csv contains invalid chromosome, position or P values")

group_map = {
    "haEWAS-specific": "Hetero_Unique",
    "EWAS-specific": "Beta_Only",
    "Common": "Beta_Shared",
    "Non_Sig": "Non_Sig",
}
if not data["Group"].isin(group_map).all():
    raise ValueError("Fig3b.csv contains an unknown group")
if not {"haEWAS-specific", "EWAS-specific", "Common"}.issubset(set(data["Group"])):
    raise ValueError("Fig3b.csv is missing a significant group")

if not data["Direction"].isin(["+", "-"]).all():
    raise ValueError("Fig3b.csv contains an invalid direction")

data["Category"] = data["Group"].map(group_map)
data["Score"] = -np.log10(data["P"]) * data["Direction"].map({"+": 1, "-": -1})
data["Label"] = data["Gene"].fillna("").astype(str)

chromosome_min = data.groupby("Chr")["Pos"].transform("min")
chromosome_max = data.groupby("Chr")["Pos"].transform("max")
data["Position"] = (
    data["Chr"] - 1 + (data["Pos"] - chromosome_min) / (chromosome_max - chromosome_min + 1)
)
significant_indices = data.index[data["Category"].ne("Non_Sig")]
rng = np.random.default_rng(42)
data.loc[significant_indices, "Position"] += rng.uniform(
    -0.03, 0.03, size=len(significant_indices)
)

def italic_gene(label):
    escaped = str(label)
    for char, replacement in (
        ("\\", r"\backslash "),
        ("_", r"\_"),
        ("%", r"\%"),
        ("&", r"\&"),
        ("#", r"\#"),
        ("{", r"\{"),
        ("}", r"\}"),
    ):
        escaped = escaped.replace(char, replacement)
    return r"$\mathit{" + escaped + r"}$"


def draw_labels(ax, significant):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    span = ymax - ymin
    margin_y = 0.035 * span
    offsets_x = {1: 0, 2: -0.10, 3: 0.10, 4: -0.20, 5: 0.20}
    grouped = {1: [], -1: []}

    for side in label_side_order:
        for category in label_category_order:
            selected = significant[
                significant["Category"].eq(category)
                & (np.sign(significant["Score"]) == side)
                & significant["Label"].ne("")
            ].copy()
            selected["Abs_Score"] = selected["Score"].abs()
            selected = selected.sort_values(
                ["Abs_Score", "Label"], ascending=[False, True], kind="mergesort"
            )
            for rank, row in enumerate(selected.itertuples(index=False), start=1):
                rank5 = ((rank - 1) % 5) + 1
                x = np.clip(row.Position + offsets_x[rank5], xmin + 0.16, xmax - 0.16)
                y = np.clip(
                    row.Score + side * (0.050 * span + (rank5 - 1) * 0.010 * span),
                    ymin + margin_y,
                    ymax - margin_y,
                )
                text = ax.text(
                    x,
                    y,
                    italic_gene(row.Label),
                    fontsize=fs_gene,
                    ha="center",
                    va="center",
                    color=plot_colors[category],
                    zorder=10,
                    clip_on=True,
                )
                text.set_path_effects([withStroke(linewidth=1.05, foreground="white")])
                grouped[side].append((row.Position, row.Score, text))

    for side, items in grouped.items():
        if not items:
            continue
        texts = [text for _, _, text in items]
        if adjust_text is not None:
            options = dict(
                ax=ax,
                expand=(1.04, 1.10),
                force_text=(0.14, 0.22),
                force_static=(0.06, 0.10),
                max_move=(7, 7),
                iter_lim=200,
            )
            try:
                adjust_text(
                    texts,
                    target_x=[x for x, _, _ in items],
                    target_y=[y for _, y, _ in items],
                    ensure_inside_axes=True,
                    expand_axes=False,
                    prevent_crossings=True,
                    only_move={
                        "text": "xy",
                        "static": "xy",
                        "explode": "xy",
                        "pull": "xy",
                    },
                    **options,
                )
            except TypeError:
                adjust_text(texts, only_move={"text": "xy"}, **options)
        for anchor_x, anchor_y, text in items:
            x, y = text.get_position()
            x = float(np.clip(x, xmin + 0.16, xmax - 0.16))
            if side > 0:
                y = float(np.clip(y, max(0.30, anchor_y + 0.35), ymax - margin_y))
            else:
                y = float(np.clip(y, ymin + margin_y, min(-0.30, anchor_y - 0.35)))
            text.set_position((x, y))
            ax.annotate(
                "",
                xy=(anchor_x, anchor_y),
                xytext=(x, y),
                annotation_clip=True,
                zorder=6,
                arrowprops=dict(
                    arrowstyle="-",
                    color=text.get_color(),
                    lw=0.32,
                    alpha=0.90,
                    shrinkA=2.0,
                    shrinkB=1.5,
                ),
            )


font_file = root.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
if font_file.is_file():
    fm.fontManager.addfont(str(font_file))
    font_family = fm.FontProperties(fname=str(font_file)).get_name()
else:
    font_family = "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_family,
    "font.size": 6.0,
    "axes.labelsize": fs_axis,
    "axes.titlesize": fs_axis,
    "xtick.labelsize": fs_tick,
    "ytick.labelsize": fs_tick,
    "legend.fontsize": fs_legend,
    "axes.linewidth": 0.55,
    "lines.linewidth": 0.7,
    "patch.linewidth": 0.5,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
})

data["Display_Score"] = data["Score"]
y_limit = max(6.0, data["Display_Score"].abs().max() * y_headroom_factor)
fig = plt.figure(
    figsize=(panel_w_mm * mm_to_inch, panel_h_mm * mm_to_inch),
    facecolor="white",
)
ax = fig.add_axes([8.5 / panel_w_mm, 8.0 / panel_h_mm, 105.0 / panel_w_mm, 39.5 / panel_h_mm])
ax.set_xlim(0, 22)
ax.set_ylim(-y_limit, y_limit)

for chromosome in range(1, 23):
    ax.axvspan(
        chromosome - 1,
        chromosome,
        facecolor=plot_colors["Bg_Strip_Even" if chromosome % 2 == 0 else "Bg_Strip_Odd"],
        edgecolor="none",
        zorder=0,
    )

background = data[data["Category"].eq("Non_Sig")].copy()
background["Chromosome"] = np.floor(background["Position"]).astype(int) + 1
for chromosome in range(1, 23):
    part = background[background["Chromosome"].eq(chromosome)]
    if not part.empty:
        ax.scatter(
            part["Position"],
            part["Display_Score"],
            marker=".",
            s=0.72,
            color=plot_colors["Background_Even" if chromosome % 2 == 0 else "Background_Odd"],
            alpha=0.40,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

significant = data[data["Category"].ne("Non_Sig")].copy()
if not significant.empty:
    points = significant.sample(frac=1, random_state=42)
    ax.scatter(
        points["Position"],
        points["Display_Score"],
        s=2.8,
        color=points["Category"].map(plot_colors),
        edgecolors="white",
        linewidths=0.16,
        alpha=0.98,
        zorder=3,
    )
    draw_labels(ax, significant)

ax.axhline(0, color=plot_colors["Axis"], linestyle=(0, (2.2, 2.2)), lw=0.5, zorder=2)
ax.yaxis.set_major_locator(MaxNLocator(nbins=5, symmetric=True))
ax.xaxis.set_major_locator(FixedLocator([chromosome - 0.5 for chromosome in range(1, 23)]))
ax.xaxis.set_major_formatter(FixedFormatter([str(chromosome) for chromosome in range(1, 23)]))
ax.tick_params(
    which="major",
    length=2.0,
    width=0.5,
    labelsize=fs_tick,
    pad=1.5,
    colors=plot_colors["Axis"],
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.55)
    spine.set_color(plot_colors["Axis"])
ax.set_ylabel(r"Signed $-\log_{10}(P)$", fontsize=fs_axis, labelpad=1.8)
ax.set_xlabel("Chromosome", fontsize=fs_axis, labelpad=1.8)

legend_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        label=legend_map[category],
        markerfacecolor=plot_colors[category],
        markeredgecolor="white",
        markeredgewidth=0.25,
        markersize=3.0,
    )
    for category in ("Hetero_Unique", "Beta_Only", "Beta_Shared")
]
fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.56, 0.975),
    ncol=3,
    fontsize=fs_legend,
    handlelength=0.8,
    handletextpad=0.25,
    columnspacing=0.75,
    borderaxespad=0,
)
fig.text(
    0.010, 0.972, "b", ha="left", va="top",
    fontsize=fs_panel, fontweight="bold", color=plot_colors["Axis"],
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








