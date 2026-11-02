"""Plot FigS9 from its compact reproduction table."""

import contextlib
import importlib.util
import io
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
CODE = REVISION_ROOT / "Fig_new" / "code"
DATA = ROOT / "data" / "FigS9.csv"
OUTPUT = ROOT / "figures" / "FigS9.pdf"
COLUMNS = [
    "Trait", "Category", "CAMDA", "CHALM", "Shared",
    "haEWAS-specific", "haEWAS (CAMDA)", "haEWAS (CHALM)", "EWAS",
]


def load_source():
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    path = CODE / "plot_r3_7_figs1_contribution.py"
    spec = importlib.util.spec_from_file_location("figs9_plot_source", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run():
    if not DATA.is_file():
        raise FileNotFoundError("FigS9 data is missing; prepare it first")

    data = pd.read_csv(DATA, dtype={"Trait": str, "Category": str})
    if list(data.columns) != COLUMNS:
        raise ValueError("FigS9.csv has unexpected columns")

    data = data.rename(columns={
        "Trait": "phenotype",
        "CAMDA": "CAMDA_only",
        "CHALM": "CHALM_only",
        "Shared": "CHALM_CAMDA_only",
        "haEWAS-specific": "haEWAS_specific_n",
        "haEWAS (CAMDA)": "camda_n",
        "haEWAS (CHALM)": "chalm_n",
        "EWAS": "beta_n",
    })

    source = load_source()
    source.apply_style()
    label_height = source._label_panel_height_mm(data["phenotype"].astype(str))
    figure = plt.figure(
        figsize=(source.FIG_W_MM / source.MM_PER_INCH,
                 source.FIG_H_MM / source.MM_PER_INCH),
        facecolor="white",
    )
    ax = figure.add_axes(source._rect_mm(
        source.PLOT_LEFT_MM, source.BAR_TOP_MM,
        source.PLOT_WIDTH_MM, source.BAR_HEIGHT_MM,
    ))
    source._draw_stacked_bar(ax, data)
    source._draw_method_heatmap(figure, data)
    source._add_x_strip(figure, data, label_height)

    handles = [
        source.Patch(facecolor=color, edgecolor="none", label=label)
        for _, label, color in source.CONTRIBUTION_COLUMNS
    ]
    figure.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.50, 0.985),
        ncol=3, frameon=False, fontsize=6.0, handlelength=0.9,
        handletextpad=0.3, columnspacing=1.0, borderaxespad=0,
    )
    with contextlib.redirect_stdout(io.StringIO()):
        source.finish_fixed(figure, OUTPUT)


if __name__ == "__main__":
    run()