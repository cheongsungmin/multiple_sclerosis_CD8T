"""Supplementary Figure 6 (supports Main Figure 5) — one supple per main figure.

Main Fig. 5 (original) is already dense with panels a–j, and the monocyte UMAP and
markers overlap with Fig. 1h/5c and the dedicated marker figure (Supplementary Fig. 2).
This supplement therefore contains the single piece of evidence not shown in the main
figure: **quantification of the U-ISGF3 / P-ISGF3 module score across the 4 groups
(R/NR × Pre/Post)**. Main Fig. 5i shows pre-only scores and Fig. 5j is a gene dotplot,
so only here do we show, as a score, that "U-ISGF3 stays high in non-responders even
after treatment".

  a  U-ISGF3 monocyte module score, R/NR × Pre/Post (paired samples)
  b  P-ISGF3 monocyte module score, R/NR × Pre/Post (paired samples)

The module score and gene set use the same definition as the original
(03_Figures_NR/260201_2154_IFN.ipynb); recomputation is deterministic from the
self-contained data/objects.
Output: figures/supplementary/Supplementary_Fig6.pdf (single vector, dpi300, Illustrator-editable)
        data/results/supp_fig6_isgf3_scores.csv (per-sample module score, source data)
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import sys, warnings
sys.path.insert(0, os.path.join(MS_ROOT, "code"))
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from common import io, style, palette
style.apply()
BASE, TICK, PLET = style.BASE, style.TICK, style.PANEL_LETTER
OBJ = os.path.join(MS_ROOT, "data/objects")
OUT = os.path.join(MS_ROOT, "figures/supplementary")
RES = os.path.join(MS_ROOT, "data/results")
GCOL = palette.GROUP   # canonical Responder/Nonresponder × Before/After colors

# ------------------------------------------------------------------ gene sets
# U-ISGF3 / P-ISGF3 = original definition as-is (03_Figures_NR/260201_2154_IFN.ipynb, cells 19–20).
# only legacy symbols are remapped to current ones (matching panel j): G1P2->ISG15, BIRC4BP->XAF1.
_ALIAS = {"G1P2": "ISG15", "BIRC4BP": "XAF1"}
_mp = lambda g: _ALIAS.get(g, g)
PISGF3_RAW = [
    "ARHGAP27", "B3GALT3", "BLZF1", "C1orf80", "CUTA", "CXCL16", "CXorf38", "DLL1",
    "EHD4", "EIF2AK2", "FAM46A", "FLJ11286", "HES4", "HLA-F", "IFI16", "INSIG1", "LAMA1",
    "MASTL", "MGC19764", "MT1M", "MYD88", "NCOA7", "NFE2L3", "NUDCD1", "OGFR", "PANX1",
    "PARP10", "PCGF5", "PHACTR4", "PI4K2B", "PLEKHA4", "PNPT1", "PPM1K", "PRKD2", "PRP2",
    "RAB20", "RASGRP3", "RNF36", "RP3-473B4.1", "SAMD9", "SOX9", "STARD5", "TDRD7", "TGM2",
    "TNFAIP6", "TRIM5", "UBE1L", "UNC93B1", "ZFYVE26",
]
UISGF3_RAW = sorted(set([
    "BATF2", "BIRC4BP", "BST2", "DDX58", "DDX60", "EPSTI1", "G1P2", "HERC5", "HERC6",
    "IFI27", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT3", "IFITM1", "IRF7",
    "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "PLSCR1", "PARRES3", "RTP4", "STAT1",
    "TMEM140", "FLJ20035", "IPSTI1",
]))

# ------------------------------------------------------------------ load object
ad = sc.read_h5ad(f"{OBJ}/MS_scRNA_anno2_myeloid.h5ad")   # X = log-norm, metadata embedded
ad.obs["anno2"] = ad.obs["anno2"].astype(str)
UISGF3 = sorted({_mp(g) for g in UISGF3_RAW if _mp(g) in ad.var_names})
PISGF3 = sorted({_mp(g) for g in PISGF3_RAW if _mp(g) in ad.var_names})
print(f"U-ISGF3 genes scored: {len(UISGF3)} | P-ISGF3 genes scored: {len(PISGF3)}")

# ------------------------------------------------------------------ module scores (original method)
# original: scale all myeloid cells (max_value=10) then score_genes, then per-sample mean over M_Mono cells.
scaled = ad.copy()
sc.pp.scale(scaled, max_value=10)
sc.tl.score_genes(scaled, UISGF3, score_name="UISG")
sc.tl.score_genes(scaled, PISGF3, score_name="PISG")
mo = scaled.obs.copy()
mo["anno2"] = mo["anno2"].astype(str)
mono_obs = mo[mo["anno2"].str.startswith("M_Mono")].copy()


def per_sample(df):
    return df.groupby("SampleID", observed=True).agg(
        UISG=("UISG", "mean"), PISG=("PISG", "mean"),
        Group1=("Group1", "first"), NEDA3_1yr=("NEDA3_1yr", "first"),
        Timepoint=("Timepoint", "first"), Paired=("Paired", "first")).reset_index()


ps = per_sample(mono_obs[mono_obs["Paired"] == "Yes"])   # paired R/NR × Pre/Post
ps.to_csv(f"{RES}/supp_fig6_isgf3_scores.csv", index=False)

# ================================================================== figure
G4 = ["Responder_Before", "Nonresponder_Before", "Responder_After", "Nonresponder_After"]
G4LAB = ["R\nPre", "NR\nPre", "R\nPost", "NR\nPost"]
G4COL = [GCOL["Responder_Before"], GCOL["Nonresponder_Before"],
         GCOL["Responder_After"], GCOL["Nonresponder_After"]]

fig, axes = plt.subplots(1, 2, figsize=(4.6, 3.0))
plt.subplots_adjust(wspace=0.55, bottom=0.18, top=0.82, left=0.13, right=0.97)


def score_box(ax, score, ylab, letter):
    data = [ps.loc[ps.Group1 == g, score].values for g in G4]
    bp = ax.boxplot(data, positions=range(len(G4)), widths=0.6, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for patch, c in zip(bp["boxes"], G4COL):
        patch.set_facecolor(c); patch.set_alpha(0.35); patch.set_edgecolor("black")
    rng = np.random.default_rng(0)
    for i, d in enumerate(data):
        jit = (rng.random(len(d)) - 0.5) * 0.28
        ax.scatter(np.full(len(d), i) + jit, d, s=9, c=G4COL[i],
                   edgecolors="black", linewidths=0.4, zorder=3)
    ax.set_xticks(range(len(G4))); ax.set_xticklabels(G4LAB, fontsize=TICK)
    ax.set_ylabel(ylab, fontsize=BASE)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2)
    ax.set_title(ylab.replace(" score", ""), fontsize=BASE)
    ax.text(-0.28, 1.16, letter, transform=ax.transAxes, fontsize=PLET,
            fontweight="bold", va="top", ha="left")
    ppre = mannwhitneyu(data[0], data[1], alternative="two-sided")[1]
    ppost = mannwhitneyu(data[2], data[3], alternative="two-sided")[1]
    # pre R vs NR bracket
    y = ax.get_ylim()[1]
    ax.plot([0, 0, 1, 1], [y, y * 1.03 if y > 0 else y + 0.01,
            y * 1.03 if y > 0 else y + 0.01, y], lw=0.6, c="black")
    txt = f"P = {ppre:.3f}" if ppre >= 0.001 else "P < 0.001"
    ax.text(0.5, y * 1.05 if y > 0 else y + 0.015, txt, ha="center", va="bottom", fontsize=4.8)
    return ppre, ppost


pU_pre, pU_post = score_box(axes[0], "UISG", "U-ISGF3 score", "a")
pP_pre, pP_post = score_box(axes[1], "PISG", "P-ISGF3 score", "b")

fig.savefig(f"{OUT}/Supplementary_Fig6.pdf", dpi=300, bbox_inches="tight")
print("saved:", f"{OUT}/Supplementary_Fig6.pdf")
print(f"[stats] U-ISGF3 pre R vs NR P={pU_pre:.4f}, post P={pU_post:.4f}; "
      f"P-ISGF3 pre P={pP_pre:.4f}, post P={pP_post:.4f}")
print("[n per group]", {g: int((ps.Group1 == g).sum()) for g in G4})
