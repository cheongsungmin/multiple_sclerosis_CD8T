"""Supplementary Fig. 4 panel b (flow, age/sex robustness) — individual (3,4) panels to cut out.

Each panel is a standalone figsize=(3,4) figure using the user's exact barstrip style
(sns.barplot width=.6 SE capsize=.2 lightgrey; sns.stripplot size=8; Set2 green/orange;
R_Pre/N_Pre labels), saved one-per-page to a multi-page PDF so proportions match exactly.

Panels (all two-sided Student t), assembled into Supplementary Fig. 4 b/c/d:
  b  KIR+/CD8+, CD4 CTL/CD4+, GZMK+/CD8+, CD8+/live  R vs NR  -> all n.s. (signal Temra-specific)
  c  Temra by sex (all)           CD29+KIR- by sex (all)      -> no sex difference
  d  Temra age-matched R vs NR    CD29+KIR- age-matched R vs NR -> difference holds, age-independent
Age-matched subset (10 R / 10 NR) reconstructed from the cowork Supplementary Table by
(Age, Sex, Responder); reproduces cowork numbers.

Source (READ-ONLY): <MS_SRC>/DATA/24NOV/FACS/250824_FACS_ODS.csv
Output: figures/supplementary/Supplementary_Fig4.pdf  (8 pages: b marker-specificity ×4, c sex ×2, d age-matched ×2)
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl, seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import ttest_ind

sns.set_style("white"); sns.set_style("ticks")
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({"text.color": "black", "axes.edgecolor": "black", "axes.labelcolor": "black",
                     "xtick.color": "black", "ytick.color": "black"})
mpl.rcParams.update({"lines.color": "black", "patch.edgecolor": "black"})

SRC = os.path.join(MS_SRC, "DATA/24NOV/FACS/250824_FACS_ODS.csv")
OUT = os.path.join(MS_ROOT, "figures/supplementary/Supplementary_Fig4.pdf")
GREEN, ORANGE = "#66C2A5", "#FC8D62"

ff = pd.read_csv(SRC).rename(columns={"Temra/CD8": "Temra", "ZNF683/CD8": "CD29posKIRneg",
                                      "Age(Sampling)": "Age", "KIR/CD8": "KIRpos",
                                      "CD4 CTL/CD4": "CD4CTL", "GZMK/CD8": "GZMKpos",
                                      "CD8 T cells/Live": "CD8live"})
for _c in ["Temra", "CD29posKIRneg", "KIRpos", "CD4CTL", "GZMKpos", "CD8live"]:
    ff[_c] = pd.to_numeric(ff[_c], errors="coerce")
cc = ff[ff.Responder_1yr.isin([0, 1])].copy()
cc["group"] = np.where(cc.Responder_1yr == 1, "Responder", "Non-responder")

# ---- age-matched subset (cowork Supplementary Table, matched by Age/Sex/Responder) ----
NR = [(24, 'F'), (25, 'F'), (25, 'F'), (26, 'F'), (29, 'F'), (30, 'F'), (32, 'M'), (34, 'M'), (36, 'M'), (41, 'M')]
R  = [(25, 'F'), (28, 'F'), (28, 'F'), (31, 'F'), (32, 'F'), (33, 'F'), (37, 'F'), (39, 'F'), (40, 'M'), (40, 'F')]


def match(df, tuples, resp):
    sub = df[df.Responder_1yr == resp]; used = set(); idx = []
    for age, sex in tuples:
        cand = sub[(sub.Age == age) & (sub.Sex == sex) & (~sub.index.isin(used))]
        if len(cand):
            used.add(cand.index[0]); idx.append(cand.index[0])
    return df.loc[idx]


am = pd.concat([match(cc, NR, 0), match(cc, R, 1)])
am["group"] = np.where(am.Responder_1yr == 1, "Responder", "Non-responder")


def barstrip(df, col, xcol, order, xlabels, pal, ylabel, title, figsize=(3, 4)):
    """Exact replica of the user's barstrip, generalized over df/x/order. Returns the figure."""
    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(data=df, x=xcol, y=col, order=order, ax=ax, width=.6,
                errorbar="se", capsize=.2, err_kws={"linewidth": .8},
                linewidth=.8, edgecolor="black", color="lightgrey", zorder=1)
    sns.stripplot(data=df, x=xcol, y=col, order=order, ax=ax, palette=pal,
                  size=8, linewidth=.5, edgecolor="black", zorder=2)
    a = df[df[xcol] == order[0]][col].dropna(); b = df[df[xcol] == order[1]][col].dropna()
    p = ttest_ind(a, b, equal_var=True).pvalue
    yt = df[col].max() * 1.12
    ax.plot([0, 0, 1, 1], [yt * 0.98, yt, yt, yt * 0.98], color="black", lw=0.7)
    ax.text(0.5, yt * 1.01, f"p = {p:.3f}", ha="center", va="bottom")
    ax.set_xticks([0, 1]); ax.set_xticklabels(xlabels)
    ax.set_ylabel(ylabel); ax.set_xlabel(""); ax.set_title(title)
    ax.set_ylim(0, yt * 1.2); sns.despine(ax=ax); plt.tight_layout()
    print(f"{title!r:32s} p={p:.4f}")
    return fig


GREY = {"F": "0.6", "M": "0.6"}
RNR = {"Responder": GREEN, "Non-responder": ORANGE}
SPEC = ["Responder", "Non-responder"]
panels = [
    # panel b — marker specificity (full R-vs-NR cohort; all non-significant -> signal
    #           specific to the naive-to-Temra shift). Same barstrip style as c/d.
    (cc, "KIRpos",  "group", SPEC, ["R_Pre", "N_Pre"], RNR, "fraction", "KIR⁺ / CD8⁺"),
    (cc, "CD4CTL",  "group", SPEC, ["R_Pre", "N_Pre"], RNR, "fraction", "CD4 CTL / CD4⁺"),
    (cc, "GZMKpos", "group", SPEC, ["R_Pre", "N_Pre"], RNR, "fraction", "GZMK⁺ / CD8⁺"),
    (cc, "CD8live", "group", SPEC, ["R_Pre", "N_Pre"], RNR, "fraction", "CD8⁺ / live"),
    # panel c — sex robustness
    (cc, "Temra",         "Sex",   ["F", "M"], ["Female", "Male"], GREY, "Temra / CD8⁺",      "Temra (by sex)"),
    (cc, "CD29posKIRneg", "Sex",   ["F", "M"], ["Female", "Male"], GREY, "CD29⁺KIR⁻ / CD8⁺",  "CD29⁺KIR⁻ (by sex)"),
    # panel d — age-matched R vs NR
    (am, "Temra",         "group", ["Responder", "Non-responder"], ["R_Pre", "N_Pre"], RNR, "Temra / CD8⁺",     "Temra (age-matched)"),
    (am, "CD29posKIRneg", "group", ["Responder", "Non-responder"], ["R_Pre", "N_Pre"], RNR, "CD29⁺KIR⁻ / CD8⁺", "CD29⁺KIR⁻ (age-matched)"),
]

with PdfPages(OUT) as pdf:
    for args in panels:
        fig = barstrip(*args)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

print(f"\nage-matched n={len(am)} (R {int((am.Responder_1yr==1).sum())}/NR {int((am.Responder_1yr==0).sum())})")
print("wrote (8 pages = 4 specificity[b] + 2 sex[c] + 2 age-matched[d], each figsize=(3,4)):", OUT)
