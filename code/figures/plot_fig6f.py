"""Plot prepared Fig6f.csv."""

import importlib
import _fig6_go_heatmap

if __name__ == "__main__":
    importlib.reload(_fig6_go_heatmap).run("Chol", "f")