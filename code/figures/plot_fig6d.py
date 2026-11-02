"""Plot prepared Fig6d.csv."""

import importlib
import _fig6_landscape

if __name__ == "__main__":
    importlib.reload(_fig6_landscape).run("Chol", "d")