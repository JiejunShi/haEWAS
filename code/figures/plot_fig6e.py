"""Plot prepared Fig6e.csv."""

import importlib
import _fig6_go_heatmap

if __name__ == "__main__":
    importlib.reload(_fig6_go_heatmap).run("Hep", "e")