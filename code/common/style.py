"""Publication plotting style (matplotlib/seaborn). Import and call apply().

Reflects the original user settings: all text black, pdf.fonttype=42 (Illustrator-editable).
Size scheme (A4 basis): base 6 / tick 5 / panel letter bold 12.
The panel letter size is exposed as PANEL_LETTER.
"""
import matplotlib as mpl
import seaborn as sns

BASE = 6            # base text (titles, axis labels, legends, annotations)
TICK = 5            # axis tick labels
PANEL_LETTER = 12   # panel letter titles (bold)


def apply():
    mpl.rcParams["pdf.fonttype"] = 42       # editable text in Illustrator
    mpl.rcParams["ps.fonttype"] = 42
    # --- all text/lines black (original user settings) ---
    for k in ["text.color", "axes.edgecolor", "axes.labelcolor",
              "axes.titlecolor", "xtick.color", "ytick.color", "lines.color"]:
        mpl.rcParams[k] = "black"
    mpl.rcParams["patch.edgecolor"] = "black"
    # --- uniform font sizes ---
    mpl.rcParams["font.size"] = BASE
    mpl.rcParams["axes.titlesize"] = BASE
    mpl.rcParams["axes.labelsize"] = BASE
    mpl.rcParams["legend.fontsize"] = BASE
    mpl.rcParams["figure.titlesize"] = BASE
    mpl.rcParams["xtick.labelsize"] = TICK
    mpl.rcParams["ytick.labelsize"] = TICK
    # --- white background ---
    mpl.rcParams["figure.facecolor"] = "white"
    mpl.rcParams["axes.facecolor"] = "white"
    mpl.rcParams["savefig.facecolor"] = "white"
    sns.set_style("white")
    sns.set_style("ticks")
