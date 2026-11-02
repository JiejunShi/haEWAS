"""Draw FigS17 from its prepared data tables."""

import contextlib
import gzip
import io
import os
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke
from matplotlib.ticker import FixedFormatter, FixedLocator, MaxNLocator
from matplotlib_venn import venn2


ROOT = Path(os.environ["REPRO_ROOT"])
DATA = ROOT / "data"
OUTPUT = ROOT / "figures" / "FigS17.pdf"
WIDTH_MM = 135.0
HEIGHT_MM = 189.0
MM = 1 / 25.4
ROW2_Y = 53.0
ROW3_Y = 111.0
A_SIZE = (68.0, 42.0)
B_SIZE = (59.0, 45.0)
CD_SIZE = (56.0, 50.0)
A_X, B_X = 0.0, 76.0
CD_X = 7.5
C_X, D_X = CD_X, CD_X + CD_SIZE[0] + 8.0
E_X = C_X + .28 * CD_SIZE[0] - 8.5
LEFT_LABEL_X = C_X + .28 * CD_SIZE[0] - 1.8
RIGHT_LABEL_X = D_X + .28 * CD_SIZE[0] - 1.8

METHODS = ("CAMDA", "CHALM", "EWAS")
METHOD_LABELS = {
    "CAMDA": "haEWAS (CAMDA)",
    "CHALM": "haEWAS (CHALM)",
    "EWAS": "EWAS",
}
METHOD_COLORS = {
    "EWAS": "#5F7FA1",
    "CHALM": "#D9893D",
    "CAMDA": "#3D9AA8",
}
COLORS = {
    "EWAS-specific": "#5F7FA1",
    "Common": "#C9B48A",
    "haEWAS-specific": "#CA6A6D",
    "Background_Odd": "#B9B9B9",
    "Background_Even": "#D0D0D0",
    "Bg_Strip_Odd": "#F6F6F6",
    "Bg_Strip_Even": "#FFFFFF",
    "Axis": "#222222",
}
LAMBDA_METHODS = ("EWAS", "haEWAS(CAMDA)", "haEWAS(CHALM)")
QQ_METHODS = ("EWAS", "haEWAS-CAMDA", "haEWAS-CHALM")
QQ_LABELS = {
    "EWAS": "EWAS",
    "haEWAS-CAMDA": "haEWAS (CAMDA)",
    "haEWAS-CHALM": "haEWAS (CHALM)",
}
QQ_COLORS = {
    "EWAS": "#0072B2",
    "haEWAS-CAMDA": "#009E73",
    "haEWAS-CHALM": "#E69F00",
}
E_COLUMNS = ["Chr", "CpG", "Gene", "P", "Direction", "Group", "Pos"]
LAMBDA_COLORS = {
    "EWAS": QQ_COLORS["EWAS"],
    "haEWAS(CAMDA)": QQ_COLORS["haEWAS-CAMDA"],
    "haEWAS(CHALM)": QQ_COLORS["haEWAS-CHALM"],
}


def read_data(name):
    path = DATA / (name + ".csv")
    if not path.is_file():
        path = Path(str(path) + ".gz")
    if not path.is_file():
        raise FileNotFoundError("Prepare {} data first".format(name))
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return pd.read_csv(handle, low_memory=False)


def add_axes_mm(fig, x, y_top, width, height):
    return fig.add_axes([
        x / WIDTH_MM, 1 - (y_top + height) / HEIGHT_MM,
        width / WIDTH_MM, height / HEIGHT_MM,
    ])


def add_panel_label(fig, text, x, y_top, ha="left", va="bottom"):
    fig.text(x / WIDTH_MM, 1 - y_top / HEIGHT_MM, text, ha=ha, va=va,
             fontsize=6.5, fontweight="bold", color=COLORS["Axis"])


def panel_a(ax, data):
    counts = data.set_index("Category")["Count"].astype(int)
    ewas_only = int(counts["EWAS-specific"])
    common = int(counts["Common"])
    haewas_only = int(counts["haEWAS-specific"])
    venn = venn2(
        subsets=(ewas_only, haewas_only, common),
        set_labels=("EWAS", "pooled haEWAS"), ax=ax,
    )
    for region, color in (("10", COLORS["EWAS-specific"]),
                          ("01", COLORS["haEWAS-specific"]),
                          ("11", COLORS["Common"])):
        patch = venn.get_patch_by_id(region)
        if patch is not None:
            patch.set_facecolor(color)
            patch.set_edgecolor("white")
            patch.set_linewidth(.55)
            patch.set_alpha(.96)
        text = venn.get_label_by_id(region)
        if text is not None:
            value = {"10": ewas_only, "01": haewas_only, "11": common}[region]
            total = ewas_only + haewas_only + common
            text.set_text("{:,}\n({:.1f}%)".format(value, 100 * value / total))
            text.set_fontsize(5.4)
            text.set_color("#111111")
            text.set_path_effects([withStroke(linewidth=.9, foreground="white")])
    ax.set_axis_off()


def panel_b(fig, x, y_top, data):
    left = (x + .06 * B_SIZE[0]) / WIDTH_MM
    right = (x + .985 * B_SIZE[0]) / WIDTH_MM
    top = 1 - (y_top + .10 * B_SIZE[1]) / HEIGHT_MM
    bottom = 1 - (y_top + .93 * B_SIZE[1]) / HEIGHT_MM
    layout = fig.add_gridspec(
        2, 3, left=left, right=right, top=top, bottom=bottom,
        width_ratios=(3.5, 4.4, 2.1), height_ratios=(3.0, .90),
        wspace=.04, hspace=.025,
    )
    intersection = data.loc[data["Type"].eq("Intersection")].copy()
    sizes = data.loc[data["Type"].eq("SetSize")].set_index("Combination")
    counts = intersection["Count"].to_numpy(int)
    if len(counts) and counts.max() > 1500:
        raise ValueError("FigS17b intersection exceeds the source panel scale")
    x = np.arange(len(intersection))
    y = np.arange(3)

    blank_left = fig.add_subplot(layout[0, 0])
    blank_right = fig.add_subplot(layout[0, 2])
    blank_left.set_axis_off()
    blank_right.set_axis_off()
    ax_count = fig.add_subplot(layout[0, 1])
    ax_method = fig.add_subplot(layout[1, 0])
    ax_matrix = fig.add_subplot(layout[1, 1], sharex=ax_count)
    ax_size = fig.add_subplot(layout[1, 2], sharey=ax_matrix)

    ax_count.bar(x, counts, width=.46, color="#666666", edgecolor="none", zorder=2)
    ax_count.set_xlim(-.5, len(intersection) - .5)
    ax_count.set_ylim(0, 1500)
    ax_count.set_yticks(np.linspace(0, 1500, 6))
    ax_count.set_ylabel("Intersection size", labelpad=2, fontsize=6)
    ax_count.set_xticks([])
    ax_count.tick_params(axis="x", bottom=False, top=False, labelbottom=False)
    ax_count.tick_params(axis="y", length=1.7, width=.4, pad=1.2,
                         labelsize=6, colors=COLORS["Axis"])
    ax_count.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax_count.spines["left"].set_color(COLORS["Axis"])
    ax_count.spines["left"].set_linewidth(.45)
    ax_count.spines["left"].set_bounds(0, 1500)
    for side in ("top", "right", "bottom"):
        ax_count.spines[side].set_visible(False)
    for xi, count in zip(x, counts):
        ax_count.annotate("{:,}".format(count), xy=(xi, count), xytext=(.5, 3),
                          textcoords="offset points", ha="center", va="bottom",
                          rotation=90, fontsize=5.5, color=COLORS["Axis"],
                          annotation_clip=False)

    ax_method.set_xlim(0, 1)
    ax_method.set_ylim(2.5, -.5)
    ax_method.set_axis_off()
    for yi, method in zip(y, METHODS):
        ax_method.text(.98, yi, METHOD_LABELS[method], ha="right", va="center",
                       fontsize=5.8, color=METHOD_COLORS[method])

    ax_matrix.set_xlim(-.5, len(intersection) - .5)
    ax_matrix.set_ylim(2.5, -.5)
    ax_matrix.set_xticks([])
    ax_matrix.set_yticks([])
    ax_matrix.set_axis_off()
    for xi, row in enumerate(intersection.itertuples(index=False)):
        active = [i for i, method in enumerate(METHODS)
                  if bool(getattr(row, method))]
        if len(active) > 1:
            ax_matrix.plot([xi, xi], [min(active), max(active)], color="#858585",
                           linewidth=.55, zorder=1)
        ax_matrix.scatter([xi] * 3, y, s=13, color="#E2E2E2", zorder=2)
        if active:
            ax_matrix.scatter(
                [xi] * len(active), active, s=23,
                color=[METHOD_COLORS[METHODS[i]] for i in active], zorder=3,
            )

    set_sizes = np.asarray([int(sizes.loc[method, "Count"]) for method in METHODS])
    max_size = int(set_sizes.max())
    ax_size.barh(y, set_sizes, height=.43,
                 color=[METHOD_COLORS[method] for method in METHODS], edgecolor="none")
    ax_size.set_xlim(0, max_size * 2 if max_size else 1)
    ax_size.set_ylim(2.5, -.5)
    ax_size.set_axis_off()
    for yi, value in zip(y, set_sizes):
        ax_size.text(value + max(max_size * .035, .25), yi, "{:,}".format(value),
                     ha="left", va="center", fontsize=5.8, color=COLORS["Axis"])


def panel_c(ax, data):
    offset = {"EWAS": -.20, "haEWAS(CAMDA)": 0, "haEWAS(CHALM)": .20}
    cohorts = (data.loc[data["Method"].eq("EWAS")]
               .sort_values("Lambda")["Cohort"].astype(str).tolist())
    cohorts += [x for x in sorted(data["Cohort"].astype(str).unique()) if x not in cohorts]
    y_pos = {cohort: i for i, cohort in enumerate(cohorts)}
    upper = max(2.0, float(np.ceil((data["Lambda"].max() + .08) / .25) * .25))
    ax.axvline(1, color="#666666", linewidth=.35, linestyle=(0, (3, 2)), zorder=1)
    for method in LAMBDA_METHODS:
        subset = data.loc[data["Method"].eq(method)]
        y = [y_pos[str(cohort)] + offset[method] for cohort in subset["Cohort"]]
        ax.scatter(subset["Lambda"], y, s=12, facecolor=LAMBDA_COLORS[method],
                   edgecolor="white", linewidth=.2, alpha=.94, zorder=3)
    ax.set_yticks(range(len(cohorts)))
    ax.set_yticklabels(cohorts)
    ax.set_ylim(-.6, len(cohorts) - .4)
    ax.set_xlim(0, upper)
    ticks = sorted(set(x for x in (0, .5, 1, 1.5, 2, upper) if x <= upper))
    ax.set_xticks(ticks)
    ax.set_xlabel(r"Inflation factor ($\lambda$)")
    ax.set_ylabel("")
    ax.grid(axis="x", color="#E9E9E9", linewidth=.2)
    ax.grid(axis="y", color="#DDDDDD", linewidth=.2)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(.35)
        ax.spines[side].set_color(COLORS["Axis"])
    ax.tick_params(axis="both", width=.35, length=1.2, pad=1.2,
                   colors=COLORS["Axis"], labelsize=5)
    handles = [Line2D([0], [0], marker="o", color="none",
                      markerfacecolor=LAMBDA_COLORS[method], markeredgecolor="white",
                      markeredgewidth=.2, markersize=3.4, label=method)
               for method in LAMBDA_METHODS]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.16),
              ncol=2, frameon=False, fontsize=4.8, handletextpad=.35,
              columnspacing=.65)


def panel_d(ax, data):
    ax.plot([0, 6], [0, 6], color="#666666", linewidth=.35,
            linestyle=(0, (3, 2)), zorder=1)
    lambdas = {}
    for method in QQ_METHODS:
        prefix = method
        x_col = prefix + "_Expected"
        y_col = prefix + "_Observed"
        lambda_col = prefix + "_Lambda"
        points = data[["CpG", x_col, y_col]].dropna(subset=[x_col, y_col])
        if points.empty:
            raise ValueError("Missing QQ data for {}".format(method))
        expected = points[x_col].to_numpy(float)
        observed = points[y_col].to_numpy(float)
        tail = observed > 4.0
        body = pd.DataFrame({"Expected": expected[~tail],
                             "Observed": observed[~tail]})
        body["Expected_bin"] = np.round(body["Expected"] * 22) / 22
        body["Observed_bin"] = np.round(body["Observed"] * 22) / 22
        body = body.drop_duplicates(["Expected_bin", "Observed_bin"])
        plot_points = pd.concat([
            pd.DataFrame({"Expected": expected[tail], "Observed": observed[tail]}),
            body[["Expected", "Observed"]],
        ], ignore_index=True)
        ax.scatter(plot_points["Expected"], plot_points["Observed"], s=2,
                   color=QQ_COLORS[method], alpha=.52, linewidths=0,
                   rasterized=True, zorder=2)
        values = data[lambda_col].dropna()
        if len(values) != 1:
            raise ValueError("Missing lambda value for {}".format(method))
        lambdas[method] = float(values.iloc[0])
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 25)
    ax.set_xticks([0, 2, 4, 6])
    ax.set_yticks([0, 5, 10, 15, 20, 25])
    ax.set_xlabel(r"Expected -log$_{10}$($P$)")
    ax.set_ylabel(r"Observed -log$_{10}$($P$)")
    ax.grid(True, color="#EEEEEE", linewidth=.2)
    ax.set_axisbelow(True)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(.35)
        ax.spines[side].set_color(COLORS["Axis"])
    ax.tick_params(axis="both", width=.35, length=1.2, pad=1.2,
                   colors=COLORS["Axis"], labelsize=5)
    handles = [Line2D([0], [0], marker="o", linestyle="", color=QQ_COLORS[m],
                      markersize=2.8, label="{}: $\\lambda$={:.2f}".format(
                          QQ_LABELS[m], lambdas[m])) for m in QQ_METHODS]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(.03, .98),
              frameon=False, fontsize=4.6, handletextpad=.3,
              borderaxespad=0, labelspacing=.2)


def gene_label(value):
    if pd.isna(value):
        return None
    for gene in str(value).split(";"):
        gene = gene.strip()
        if gene and gene.upper() not in {"NA", "NAN", "NONE"}:
            return gene
    return None


def italic_gene(value):
    text = str(value)
    for char, replacement in (("\\", r"\backslash "), ("_", r"\_"),
                              ("%", r"\%"), ("&", r"\&"), ("#", r"\#"),
                              ("{", r"\{"), ("}", r"\}")):
        text = text.replace(char, replacement)
    return r"$\mathit{" + text + r"}$"


def panel_e(fig, data):
    ax = add_axes_mm(fig, E_X + 8.5, ROW3_Y + 15.5, 105.0, 54.5)
    frame = data.copy()
    frame["Chr"] = pd.to_numeric(frame["Chr"], errors="coerce")
    frame["Pos"] = pd.to_numeric(frame["Pos"], errors="coerce")
    frame["P"] = pd.to_numeric(frame["P"], errors="coerce")
    frame = frame.dropna(subset=["Chr", "Pos", "P"])
    frame["Chr"] = frame["Chr"].astype(int)
    frame["Group"] = frame["Group"].astype(str)
    groups = frame["Group"]
    if (not frame["Direction"].isin(("+", "-")).all()
            or not groups.isin(("EWAS-specific", "Common",
                                "haEWAS-specific", "Non_Sig")).all()):
        raise ValueError("FigS17e data contains an unknown group")
    significant_groups = groups.loc[groups.ne("Non_Sig")]
    if not {"EWAS-specific", "Common", "haEWAS-specific"}.issubset(
            set(significant_groups)):
        raise ValueError("FigS17e data is missing a significant group")
    frame = frame.loc[frame["Chr"].between(1, 22) & frame["Pos"].gt(0)
                      & frame["P"].gt(0) & frame["P"].le(1)].copy()
    frame["SignedScore"] = -np.log10(frame["P"]) * frame["Direction"].map(
        {"+": 1, "-": -1}
    )
    frame["Category"] = frame["Group"]
    bounds = frame.groupby("Chr")["Pos"].agg(["min", "max"])
    minimum = frame["Chr"].map(bounds["min"])
    span = frame["Chr"].map(bounds["max"] - bounds["min"] + 1)
    frame["X"] = frame["Chr"] - 1 + (frame["Pos"] - minimum) / span
    frame["PlotScore"] = frame["SignedScore"].clip(-30, 30)
    y_limit = max(5.0, float(frame["PlotScore"].abs().max()) * 1.12)
    ax.set_xlim(0, 22)
    ax.set_ylim(-y_limit, y_limit)

    for chrom in range(1, 23):
        ax.axvspan(chrom - 1, chrom,
                   facecolor=COLORS["Bg_Strip_Even"] if chrom % 2 == 0
                   else COLORS["Bg_Strip_Odd"], edgecolor="none", zorder=0)
    background = frame.loc[frame["Group"].eq("Non_Sig")]
    if len(background) > 80000:
        background = background.sample(80000, random_state=42)
    for chrom in range(1, 23):
        subset = background.loc[background["Chr"].eq(chrom)]
        if not subset.empty:
            color = COLORS["Background_Even"] if chrom % 2 == 0 else COLORS["Background_Odd"]
            ax.scatter(subset["X"], subset["PlotScore"], marker=".", s=.72,
                       color=color, alpha=.40, linewidths=0, rasterized=True, zorder=1)

    significant = frame.loc[frame["Group"].ne("Non_Sig")].copy()
    if not significant.empty:
        rng = np.random.default_rng(42)
        significant["Xj"] = significant["X"] + rng.uniform(-.03, .03, len(significant))
        plotted = significant.sample(frac=1, random_state=42).reset_index(drop=True)
        ax.scatter(plotted["Xj"], plotted["PlotScore"], s=2.8,
                   color=plotted["Category"].map(COLORS), edgecolors="white",
                   linewidths=.16, alpha=.98, zorder=3)
        clipped = plotted.loc[plotted["SignedScore"].abs().gt(30)]
        if not clipped.empty:
            ax.scatter(clipped["Xj"], clipped["PlotScore"], s=8, marker="^",
                       color=clipped["Category"].map(COLORS), edgecolors="white",
                       linewidths=.2, zorder=4)

        labels, anchor_x, anchor_y, label_colors = [], [], [], []
        assigned = set()
        step = 2 * y_limit * .095
        for category in ("Common", "haEWAS-specific", "EWAS-specific"):
            for sign in (1, -1):
                candidates = significant.loc[
                    significant["Category"].eq(category)
                    & np.sign(significant["SignedScore"]).eq(sign)
                ].copy()
                candidates = candidates.reindex(
                    candidates["SignedScore"].abs().sort_values(ascending=False).index
                )
                shown = 0
                taken = 0
                for _, row in candidates.iterrows():
                    gene = gene_label(row["Gene"])
                    if not gene or gene in assigned:
                        continue
                    assigned.add(gene)
                    taken += 1
                    if shown < 5:
                        color = COLORS[category]
                        y_text = np.clip(row["PlotScore"] + sign * step,
                                         -y_limit * .94, y_limit * .94)
                        text = ax.text(row["Xj"], y_text, italic_gene(gene),
                                       fontsize=5.8, ha="center", va="center",
                                       color=color, clip_on=False, zorder=10)
                        text.set_path_effects([withStroke(linewidth=1.05,
                                                          foreground="white")])
                        labels.append(text)
                        anchor_x.append(float(row["Xj"]))
                        anchor_y.append(float(row["PlotScore"]))
                        label_colors.append(color)
                        shown += 1
                    if taken >= 8:
                        break

        adjust_text = None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                try:
                    from adjustText import adjust_text as _adjust_text
                    adjust_text = _adjust_text
                except Exception:
                    pass
        if labels and adjust_text is not None:
            adjust_text(labels, x=anchor_x, y=anchor_y, ax=ax,
                        only_move={"text": "xy"}, expand=(1.08, 1.20),
                        force_text=(.45, .65), force_static=(.30, .42),
                        max_move=(28, 28), iter_lim=300)
        for text, x0, y0, color in zip(labels, anchor_x, anchor_y, label_colors):
            x1, y1 = text.get_position()
            ax.plot([x0, x1], [y0, y1], color=color, linewidth=.32,
                    alpha=.95, solid_capstyle="round", clip_on=False, zorder=6)

    ax.axhline(0, color=COLORS["Axis"], linestyle=(0, (2.2, 2.2)), lw=.5, zorder=2)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, symmetric=True))
    ax.xaxis.set_major_locator(FixedLocator([chrom - .5 for chrom in range(1, 23)]))
    ax.xaxis.set_major_formatter(FixedFormatter([str(chrom) for chrom in range(1, 23)]))
    ax.tick_params(which="major", direction="out", length=2, width=.5,
                   bottom=True, left=True, labelsize=6, pad=1.5,
                   colors=COLORS["Axis"])
    ax.tick_params(which="minor", direction="out", length=1.4, width=.35,
                   bottom=True, left=True, colors=COLORS["Axis"])
    for spine in ax.spines.values():
        spine.set_linewidth(.55)
        spine.set_color(COLORS["Axis"])
        spine.set_visible(True)
    ax.set_ylabel(r"Signed $-\log_{10}(P)$", fontsize=6, labelpad=1.8)
    ax.set_xlabel("Chromosome", fontsize=6, labelpad=1.8)
    handles = [Line2D([0], [0], marker="o", linestyle="none", label=category,
                      markerfacecolor=COLORS[category], markeredgecolor="white",
                      markeredgewidth=.25, markersize=3)
               for category in ("haEWAS-specific", "Common", "EWAS-specific")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.56, 1.12),
              ncol=3, frameon=False, fontsize=6, handlelength=.8,
              handletextpad=.25, columnspacing=.75, borderaxespad=0)


def run():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 5.5,
        "axes.linewidth": .4,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    data = {
        "a": read_data("FigS17a"),
        "b": read_data("FigS17b"),
        "c": read_data("FigS17c"),
        "d": read_data("FigS17d"),
        "e": read_data("FigS17e"),
    }
    if list(data["e"].columns) != E_COLUMNS or data["e"].empty:
        raise ValueError("FigS17e data must use the Fig3b column format")
    figure = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), facecolor="white")
    ax_a = add_axes_mm(
        figure, A_X + .08 * A_SIZE[0], .06 * A_SIZE[1],
        .88 * A_SIZE[0], .89 * A_SIZE[1],
    )
    panel_a(ax_a, data["a"])
    add_panel_label(figure, "a", LEFT_LABEL_X, 2.0)
    panel_b(figure, B_X, 0.0, data["b"])
    add_panel_label(figure, "b", RIGHT_LABEL_X, 2.0)

    ax_c = add_axes_mm(
        figure, C_X + .28 * CD_SIZE[0], ROW2_Y + .06 * CD_SIZE[1],
        .70 * CD_SIZE[0], .69 * CD_SIZE[1],
    )
    ax_d = add_axes_mm(
        figure, D_X + .28 * CD_SIZE[0], ROW2_Y + .06 * CD_SIZE[1],
        .70 * CD_SIZE[0], .69 * CD_SIZE[1],
    )
    panel_c(ax_c, data["c"])
    panel_d(ax_d, data["d"])
    axis_y = ROW2_Y + .06 * CD_SIZE[1]
    add_panel_label(figure, "c", LEFT_LABEL_X, axis_y - .8)
    add_panel_label(figure, "d", RIGHT_LABEL_X, axis_y - .8)

    panel_e(figure, data["e"])
    add_panel_label(
        figure, "e", LEFT_LABEL_X, ROW3_Y + 15.5 - .8
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", facecolor="white", transparent=False,
                   bbox_inches=None, pad_inches=0)
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None and hasattr(get_ipython(), "kernel"):
            display(figure)
    except ImportError:
        pass
    plt.close(figure)


if __name__ == "__main__":
    run()
