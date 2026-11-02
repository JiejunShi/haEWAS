"""Render prepared GO heatmap CSVs only."""
from _fig6_plot_common import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
FDR_CUTOFF = 0.05
FDR_CAP = 10.0
PANEL_W_MM, PANEL_H_MM = 82.0, 66.0
FS_TITLE, FS_ROW, FS_COL, FS_LEGEND = 7.0, 5.0, 5.0, 5.0
CELL_CONFIGS = {
    "Hep": {
        "panel": "e", "title": "Hep", "cell_size": 4.6,
        "highlight_first_column_rows": 4,
        "terms": (
            "GOBP_LIPID_METABOLIC_PROCESS",
            "GOBP_SMALL_MOLECULE_METABOLIC_PROCESS",
            "GOBP_MEMBRANE_ORGANIZATION",
            "GOBP_SECRETION",
            "GOBP_IMPORT_INTO_CELL",
            "GOBP_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION",
            "GOBP_ORGANIC_ACID_METABOLIC_PROCESS",
            "GOBP_MITOTIC_NUCLEAR_DIVISION",
            "GOBP_CHROMOSOME_SEGREGATION",
            "GOBP_SPINDLE_ORGANIZATION",
        ),
    },
    "Chol": {
        "panel": "f", "title": "Chol", "cell_size": 4.35,
        "highlight_first_column_rows": 4,
        "terms": (
            "GOBP_CELL_JUNCTION_ORGANIZATION",
            "GOBP_REGULATION_OF_CELL_CELL_ADHESION",
            "GOBP_TUBE_DEVELOPMENT",
            "GOBP_CELL_MOTILITY",
            "GOBP_CELL_SURFACE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY",
            "GOBP_IMPORT_INTO_CELL",
            "GOBP_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION",
            "GOBP_MITOTIC_CELL_CYCLE",
            "GOBP_CHROMOSOME_SEGREGATION",
            "GOBP_SPINDLE_ORGANIZATION",
        ),
    },
}
PLOT_COLUMNS = ("haEWAS-specific", "Common", "EWAS-specific")
FILE_BY_CATEGORY = {
    "haEWAS-specific": "haewas_specific.tsv",
    "Common": "common.tsv",
    "EWAS-specific": "ewas_specific.tsv",
}
COLOR_MAP = LinearSegmentedColormap.from_list(
    "fig6_go_fdr",
    ("#FFFFFF", "#F5E4E1", "#E3A19A", "#BD5B55", "#7F1D1D"),
)




def _pretty_term(term_id):
    text = term_id[len("GOBP_"):] if term_id.startswith("GOBP_") else term_id
    text = text.replace("_", " ").lower()
    return text[:1].upper() + text[1:]

def _draw_cell(ax, x_mm, y_bottom_mm, cell_size, value):
    fill = "#FFFFFF" if not np.isfinite(value) else COLOR_MAP(
        min(FDR_CAP, max(0.0, float(value))) / FDR_CAP)
    ax.add_patch(Rectangle(
        (x_mm, y_bottom_mm), cell_size, cell_size,
        facecolor=fill, edgecolor="#D7D7D7", linewidth=0.25,
    ))

def _draw_panel(cell, matrix, font_family):
    cfg = CELL_CONFIGS[cell]
    cell_size = cfg["cell_size"]
    n_rows, n_cols = matrix.shape
    hm_x, hm_y = 6.0, 13.5
    hm_w, hm_h = cell_size * n_cols, cell_size * n_rows
    row_label_x = hm_x + hm_w + 2.2
    fig = plt.figure(figsize=(PANEL_W_MM / 25.4, PANEL_H_MM / 25.4), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, PANEL_W_MM)
    ax.set_ylim(0, PANEL_H_MM)
    ax.axis("off")

    ax.text(PANEL_W_MM / 2, 62.5, cfg["title"], ha="center", va="center",
            fontsize=FS_TITLE, color="#111111", fontfamily=font_family)
    ax.text(8.0, 62.5, cfg["panel"], ha="left", va="center",
            fontsize=8.0, fontweight="bold", color="#111111", fontfamily=font_family)

    for i, term_id in enumerate(cfg["terms"]):
        y_bottom = hm_y + hm_h - (i + 1) * cell_size
        for j in range(n_cols):
            _draw_cell(ax, hm_x + j * cell_size, y_bottom, cell_size, matrix[i, j])
        y_mid = y_bottom + cell_size / 2
        ax.text(row_label_x, y_mid, _pretty_term(term_id), ha="left", va="center",
                fontsize=FS_ROW, color="#111111", fontfamily=font_family)
        ax.plot([hm_x + hm_w, hm_x + hm_w + 1.1], [y_mid, y_mid],
                color="#6F6F6F", linewidth=0.25)

    ax.add_patch(Rectangle((hm_x, hm_y), hm_w, hm_h, facecolor="none",
                           edgecolor="#6F6F6F", linewidth=0.35))
    highlight_h = cfg["highlight_first_column_rows"] * cell_size
    ax.add_patch(Rectangle(
        (hm_x, hm_y + hm_h - highlight_h), cell_size, highlight_h,
        facecolor="none", edgecolor="#111111", linewidth=0.8,
        linestyle=(0, (2.2, 1.6)), zorder=6,
    ))

    x_labels = ("haEWAS-\nspecific", "Common", "EWAS-\nspecific")
    for j, label in enumerate(x_labels):
        x_mid = hm_x + (j + 0.5) * cell_size
        ax.plot([x_mid, x_mid], [hm_y, hm_y - 1.1], color="#6F6F6F", linewidth=0.25)
        ax.text(x_mid, hm_y - 2.6, label, ha="right", va="center", rotation=45,
                rotation_mode="anchor", fontsize=FS_COL, linespacing=0.86,
                color="#111111", fontfamily=font_family)

    # Original compact red gradient legend with -log10(FDR) ticks at 0, 5, 10.
    legend_x, legend_y, legend_w, legend_h = 74.5, 40.0, 3.0, 18.0
    n_steps = 100
    for k in range(n_steps):
        y = legend_y + legend_h * k / n_steps
        ax.add_patch(Rectangle(
            (legend_x, y), legend_w, legend_h / n_steps + 0.01,
            facecolor=COLOR_MAP((k + 0.5) / n_steps), edgecolor="none"))
    ax.add_patch(Rectangle((legend_x, legend_y), legend_w, legend_h,
                           facecolor="none", edgecolor="#6F6F6F", linewidth=0.3))
    ax.text(legend_x + legend_w / 2, legend_y + legend_h + 5.0,
            r"$-\log_{10}$ FDR", ha="center", va="center",
            fontsize=FS_LEGEND, color="#111111", fontfamily=font_family)
    for tick in (0, 5, 10):
        y = legend_y + legend_h * tick / FDR_CAP
        ax.plot([legend_x + legend_w, legend_x + legend_w + 1.2], [y, y],
                color="#6F6F6F", linewidth=0.25)
        ax.text(legend_x + legend_w + 2.0, y, str(tick), ha="left", va="center",
                fontsize=FS_LEGEND, color="#111111", fontfamily=font_family)
    return fig

def run(cell, panel):
    frame = pd.read_csv(ROOT / 'data' / ('Fig6' + panel + '.csv'))
    expected = {(term, category) for term in CELL_CONFIGS[cell]['terms'] for category in PLOT_COLUMNS}
    if (set(frame[['Term', 'Category']].itertuples(index=False, name=None)) != expected
            or frame.duplicated(['Term', 'Category']).any()):
        raise ValueError('Incomplete or duplicate GO heatmap grid')
    values = frame['MinusLog10FDR'].dropna()
    if not np.isfinite(values).all() or not values.between(0, FDR_CAP).all():
        raise ValueError('Invalid prepared GO heatmap values')
    matrix = frame.pivot(index='Term', columns='Category', values='MinusLog10FDR')
    matrix = matrix.reindex(index=CELL_CONFIGS[cell]['terms'], columns=PLOT_COLUMNS).to_numpy(float)
    style()
    family = plt.rcParams['font.family'][0]
    finish(_draw_panel(cell, matrix, family), panel)
