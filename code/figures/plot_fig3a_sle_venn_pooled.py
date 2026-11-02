import os
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
import numpy as np
import pandas as pd
from matplotlib_venn import venn2


root = Path(os.environ["REPRO_ROOT"])
input_file = root / "data" / "Fig3a.csv"
output_file = root / "figures" / "Fig3a.pdf"

panel_w_mm = 50.0
panel_h_mm = 56.0
mm_to_inch = 1 / 25.4
font_size_number = 6.0
font_size_set_label = 6.5
color_ewas_only = "#0F4D92"
color_shared = "#8F8F8F"
color_haewas_only = "#B64342"
color_text = "#272727"
color_shared_line = "#666666"

data = pd.read_csv(input_file)
segments = ["EWAS-only", "Shared", "haEWAS-only"]
if list(data.columns) != ["Segment", "Count"] or set(data["Segment"]) != set(segments):
    raise ValueError("Fig3a.csv must contain the three Venn segments and their counts")

data["Count"] = pd.to_numeric(data["Count"], errors="coerce")
values = data.set_index("Segment").loc[segments, "Count"].to_numpy(dtype=float)
if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
    raise ValueError("Fig3a.csv contains an invalid count")

n_ewas_only, n_shared, n_haewas_only = values.astype(int)
n_total = n_ewas_only + n_shared + n_haewas_only
if n_total == 0:
    raise ValueError("SLE has no reported significant loci")

font_file = root.parents[2] / "0.substitution_supplementary_data_code" / "arial.ttf"
if font_file.is_file():
    fm.fontManager.addfont(str(font_file))
    font_family = fm.FontProperties(fname=str(font_file)).get_name()
else:
    font_family = "DejaVu Sans"

plt.rcParams.update({
    "font.family": font_family,
    "font.size": 5.5,
    "axes.labelsize": 6.0,
    "axes.titlesize": 6.0,
    "xtick.labelsize": 5.0,
    "ytick.labelsize": 5.0,
    "legend.fontsize": 5.0,
    "axes.linewidth": 0.55,
    "lines.linewidth": 0.7,
    "patch.linewidth": 0.55,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

fig = plt.figure(
    figsize=(panel_w_mm * mm_to_inch, panel_h_mm * mm_to_inch),
    facecolor="white",
)
ax = fig.add_axes([
    2.5 / panel_w_mm,
    4.0 / panel_h_mm,
    45.0 / panel_w_mm,
    48.5 / panel_h_mm,
])
venn = venn2(
    subsets=(n_ewas_only, n_haewas_only, n_shared),
    set_labels=("", ""),
    alpha=1.0,
    ax=ax,
)

for region, color, alpha in (
    ("10", color_ewas_only, 0.88),
    ("11", color_shared, 0.78),
    ("01", color_haewas_only, 0.88),
):
    patch = venn.get_patch_by_id(region)
    if patch is not None:
        patch.set_facecolor(color)
        patch.set_edgecolor("white")
        patch.set_linewidth(0.70)
        patch.set_alpha(alpha)

anchors = {}
for region in ("10", "11", "01"):
    label = venn.get_label_by_id(region)
    if label is not None:
        anchors[region] = tuple(float(value) for value in label.get_position())
        label.set_visible(False)

centers = [(float(center.x), float(center.y)) for center in venn.centers]
radii = [float(radius) for radius in venn.radii]
circle_left = min(x - radius for (x, _), radius in zip(centers, radii))
circle_right = max(x + radius for (x, _), radius in zip(centers, radii))
circle_bottom = min(y - radius for (_, y), radius in zip(centers, radii))
circle_top = max(y + radius for (_, y), radius in zip(centers, radii))

if n_ewas_only > 0 and "10" in anchors:
    x, y = anchors["10"]
    ax.annotate(
        "{:,}\n({:.1f}%)".format(n_ewas_only, 100 * n_ewas_only / n_total),
        xy=(x, y),
        xytext=(x - 0.20, y + 0.08),
        ha="right",
        va="center",
        fontsize=font_size_number,
        fontweight="bold",
        color=color_ewas_only,
        linespacing=1.10,
        annotation_clip=False,
        bbox=dict(boxstyle="square,pad=0.04", facecolor="white", edgecolor="none", alpha=0.98),
        arrowprops=dict(
            arrowstyle="-",
            color=color_ewas_only,
            linewidth=0.55,
            shrinkA=0.8,
            shrinkB=0.4,
            connectionstyle="arc3,rad=0.0",
        ),
    )

if n_shared > 0 and "11" in anchors:
    x, y = anchors["11"]
    ax.annotate(
        "{:,}\n({:.1f}%)".format(n_shared, 100 * n_shared / n_total),
        xy=(x, y),
        xytext=(x - 0.01, y + 0.48),
        ha="center",
        va="bottom",
        fontsize=font_size_number,
        fontweight="bold",
        color=color_text,
        linespacing=1.10,
        annotation_clip=False,
        bbox=dict(boxstyle="square,pad=0.04", facecolor="white", edgecolor="none", alpha=0.98),
        arrowprops=dict(
            arrowstyle="-",
            color=color_shared_line,
            linewidth=0.55,
            shrinkA=0.8,
            shrinkB=0.4,
            connectionstyle="arc3,rad=0.0",
        ),
    )

if n_haewas_only > 0 and "01" in anchors:
    x, y = anchors["01"]
    ax.text(
        x,
        y,
        "{:,}\n({:.1f}%)".format(n_haewas_only, 100 * n_haewas_only / n_total),
        ha="center",
        va="center",
        fontsize=font_size_number,
        fontweight="bold",
        color="white",
        linespacing=1.10,
        path_effects=[withStroke(linewidth=0.75, foreground="#6F2929")],
        clip_on=False,
    )

label_y = circle_bottom - 0.12
ax.text(
    centers[0][0], label_y, "EWAS",
    ha="center", va="center", fontsize=font_size_set_label,
    fontweight="bold", color=color_ewas_only, clip_on=False,
)
ax.text(
    centers[1][0], label_y, "haEWAS",
    ha="center", va="center", fontsize=font_size_set_label,
    fontweight="bold", color=color_haewas_only, clip_on=False,
)
for label in venn.set_labels:
    if label is not None:
        label.set_visible(False)

ax.set_xlim(circle_left - 0.22, circle_right + 0.12)
ax.set_ylim(circle_bottom - 0.22, circle_top + 0.22)
ax.set_aspect("equal")
ax.axis("off")
fig.text(
    0.018, 0.972, "a", ha="left", va="top",
    fontsize=7.0, fontweight="bold", color=color_text,
)

output_file.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output_file, format="pdf", facecolor="white", edgecolor="white")

try:
    from IPython import get_ipython
    from IPython.display import display

    shell = get_ipython()
    if shell is not None and hasattr(shell, "kernel"):
        display(fig)
except ImportError:
    pass

plt.close(fig)
