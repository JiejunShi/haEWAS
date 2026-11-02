"""Plot FigS8 from its compact combined data file."""

import contextlib
import importlib.util
import io
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure


ROOT = Path(os.environ["REPRO_ROOT"])
REVISION_ROOT = Path(os.environ["REVISION_ROOT"])
DATA = ROOT / "data" / "FigS8.csv"
SOURCE_CODE = REVISION_ROOT / "Fig_new" / "code"
SOURCE_SCRIPT = SOURCE_CODE / "plot_figs3_pooled_bh_method_counts.py"
STAR_LABEL = "haEWAS-specific > EWAS-specific"


def run():
    if not SOURCE_SCRIPT.is_file():
        raise FileNotFoundError("Original FigS8 plotting script is missing")
    if not DATA.is_file():
        raise FileNotFoundError("FigS8 data is missing; prepare it first")

    data = pd.read_csv(DATA, dtype={"Trait": str, "Trait_category": str})
    expected = {
        "Trait", "Trait_category", "CaseN", "ControlN", "SampleN",
        "EWAS-specific", "Common", "haEWAS-specific",
        "Total_Discovery",
    }
    if not expected.issubset(data.columns):
        raise ValueError("FigS8 data has unexpected columns")
    data = data.rename(columns={
        "EWAS-specific": "Narrative_Beta_Unique",
        "Common": "Narrative_Common",
        "haEWAS-specific": "Narrative_Hetero_Unique",
        "Trait": "Y_label",
    })

    sys.path.insert(0, str(SOURCE_CODE))
    try:
        spec = importlib.util.spec_from_file_location("figs8_original_plot", SOURCE_SCRIPT)
        source_plot = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(source_plot)
    finally:
        sys.path.remove(str(SOURCE_CODE))

    source_plot.BINARY_SUMMARY = DATA
    source_plot.CONTINUOUS_SUMMARY = DATA
    source_plot.PDF = ROOT / "figures" / "FigS8.pdf"
    source_plot._load_pooled_counts = lambda: data.copy()

    original_legend = Figure.legend

    def without_star_legend(fig, *args, **kwargs):
        handles = kwargs.get("handles", ())
        if any(getattr(handle, "get_label", lambda: "")() == STAR_LABEL
               for handle in handles):
            return None
        return original_legend(fig, *args, **kwargs)

    Figure.legend = without_star_legend
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            source_plot.run()
    finally:
        Figure.legend = original_legend


if __name__ == "__main__":
    run()
