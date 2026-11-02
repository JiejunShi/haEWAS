import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT = Path(os.environ['REPRO_ROOT'])
ARIAL_TTF = ROOT.parents[2] / '0.substitution_supplementary_data_code/arial.ttf'

def style():
    import matplotlib.font_manager as fm
    family = 'DejaVu Sans'
    if ARIAL_TTF.is_file():
        fm.fontManager.addfont(str(ARIAL_TTF))
        family = fm.FontProperties(fname=str(ARIAL_TTF)).get_name()
    plt.rcParams.update({'font.family': family, 'font.size': 6,
        'axes.labelsize': 6, 'axes.titlesize': 7, 'xtick.labelsize': 5.5,
        'ytick.labelsize': 5.5, 'legend.fontsize': 5.5, 'axes.linewidth': .5,
        'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white'})

def finish(fig, panel, tight=False):
    output = ROOT / 'figures' / ('Fig6' + panel + '.pdf')
    output.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {'bbox_inches': 'tight', 'pad_inches': .1} if tight else {}
    fig.savefig(output, format='pdf', facecolor='white', **kwargs)
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None:
            display(fig)
    except ImportError:
        pass
    plt.close(fig)
