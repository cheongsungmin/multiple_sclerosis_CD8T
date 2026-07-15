"""Draft schematic for Main Figure 6 panel d — how the subset-proportion PCA
(Fig. 6e, f) is constructed. Text-rich process diagram, single standalone PDF.

RULES R0/R6: main figures are the user's originals; this is a DRAFT placement
asset for Illustrator, not an edit of the original figure PDF. Every number shown
is real (anno1 lineage composition of example pre-DMT samples for the stacked
bars; the PCA scatter and PC1/PC2 variance are the real embedding of all 23
pre-DMT samples). Deterministic (RULES R4); env = scanpy (RULES R2).

Pipeline: sample -> within-lineage subset proportions -> samples x subsets matrix
-> z-score per subset -> PCA -> 2D embedding (one point = one sample), coloured in
Fig. 6e, f.

Output: figures/panels/fig6d_pca_schematic.pdf   (single image, not the full figure)
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
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from common import io, style, naming, palette
style.apply()
BASE, TICK, PLET = style.BASE, style.TICK, style.PANEL_LETTER
OUTP = os.path.join(MS_ROOT, "figures/panels")

# ---- real data ---------------------------------------------------------------
meta = io.build_meta(); meta["Timepoint"] = meta["Timepoint"].astype(str)
pre_ids = meta.loc[meta["Timepoint"] == "Before", "SampleID"].astype(str).tolist()
prop = {}
for lin in ["T_CD8", "T_CD4", "B", "NK", "Myeloid"]:
    ad = io.load_lineage(lin)
    sid = ad.obs["SampleID"].astype(str); lab = ad.obs["anno2"].astype(str)
    tot = sid.value_counts()
    for sub in sorted(lab.unique()):
        prop[sub] = (sid[lab == sub].value_counts() / tot).reindex(tot.index).fillna(0.0)
SUB = pd.DataFrame(prop).reindex(pre_ids).dropna(how="all")
SUB = SUB.loc[:, SUB.std() > 0]
pre_ids = SUB.index.tolist()
n_samp, n_sub = SUB.shape
Zc = StandardScaler().fit_transform(SUB.values)
pca = PCA(n_components=2, random_state=0).fit(Zc)
scores = pca.transform(Zc)
evr = pca.explained_variance_ratio_

# per-sample responder status (for the colour note) + a subset for colouring
smap = meta.set_index("SampleID")
status = pd.Series(pre_ids, index=pre_ids).map(smap["NEDA3_1yr"])
colsub = "T_CD8_EM_ZNF683(hi)" if "T_CD8_EM_ZNF683(hi)" in SUB.columns else SUB.columns[0]

# anno1 composition for legible stacked bars of 3 example samples
a1 = io.load_anno1()
comp = pd.crosstab(a1.obs["SampleID"].astype(str), a1.obs["anno1"].astype(str),
                   normalize="index")
lineages = [c for c in naming.LEVEL1 if c in comp.columns]
comp = comp[lineages]
pal = palette.load_colors("anno1")
lin_cols = [pal.get(l, plt.cm.tab20(i / len(lineages))) for i, l in enumerate(lineages)]
examples = pre_ids[:3]

# =============================================================================
fig = plt.figure(figsize=(9.2, 5.0))
cv = fig.add_axes([0, 0, 1, 1]); cv.axis("off"); cv.set_xlim(0, 1); cv.set_ylim(0, 1)

# title + panel letter
cv.text(0.02, 0.965, "d", fontsize=PLET, fontweight="bold", va="top")
cv.text(0.065, 0.955, "Construction of the subset-proportion PCA (Fig. 6e, f)",
        fontsize=8, fontweight="bold", va="top")
cv.text(0.065, 0.915,
        "Each pre-treatment sample is reduced to its immune-composition vector, "
        "and samples are embedded by PCA so that one point represents one patient.",
        fontsize=5.6, va="top", style="italic", color="#333333")

GLY_B, GLY_T = 0.50, 0.86          # glyph band
ARR = dict(arrowstyle="-|>", color="#555555", lw=1.3, mutation_scale=11)


def hbox(x, w, title):
    cv.text(x + w / 2, GLY_T + 0.055, title, fontsize=6.0, fontweight="bold",
            ha="center", va="bottom")

# ---- step 1: samples ---------------------------------------------------------
axS = fig.add_axes([0.02, GLY_B, 0.135, GLY_T - GLY_B]); axS.axis("off")
axS.set_xlim(0, 1); axS.set_ylim(0, 1)
for i, y in enumerate([0.74, 0.45, 0.16]):
    axS.add_patch(Rectangle((0.10, y), 0.12, 0.16, facecolor="#e8e8e8",
                  edgecolor="black", lw=0.6))
    axS.add_patch(Rectangle((0.10, y), 0.12, 0.05, facecolor="#c0392b",
                  edgecolor="black", lw=0.6))
    axS.text(0.30, y + 0.08, f"Sample {i+1}", fontsize=5.2, va="center")
axS.text(0.5, -0.02, f"n = {n_samp} pre-DMT\nPBMC samples", fontsize=5.0,
         ha="center", va="top")
hbox(0.02, 0.135, "1  Patients")

# ---- step 2: subset proportions (stacked bars) -------------------------------
axB = fig.add_axes([0.205, GLY_B, 0.15, GLY_T - GLY_B]); axB.axis("off")
axB.set_xlim(0, 1); axB.set_ylim(0, 1)
for i, y in enumerate([0.74, 0.45, 0.16]):
    x0 = 0.0
    for l, c in zip(lineages, lin_cols):
        w = comp.loc[examples[i], l] * 1.0
        axB.add_patch(Rectangle((x0, y), w, 0.15, facecolor=c,
                      edgecolor="white", lw=0.15))
        x0 += w
axB.text(0.5, -0.02,
         r"$p_i=(p_{i,1},\,\ldots,\,p_{i,%d})$" % n_sub + "\nfraction within parent lineage",
         fontsize=5.0, ha="center", va="top")
hbox(0.205, 0.15, "2  Cell-type proportions")

# ---- step 3: samples x subsets matrix ----------------------------------------
axM = fig.add_axes([0.40, GLY_B + 0.02, 0.115, GLY_T - GLY_B - 0.02])
axM.imshow(Zc[:14, :14], cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
axM.set_xticks([]); axM.set_yticks([])
axM.set_xlabel("subsets", fontsize=5.0, labelpad=1)
axM.set_ylabel("samples", fontsize=5.0, labelpad=1)
cv.text(0.4575, GLY_B - 0.055,
        f"{n_samp} × {n_sub} matrix\nz-scored per subset (column)",
        fontsize=5.0, ha="center", va="top")
hbox(0.40, 0.115, "3  Compose & scale")

# ---- step 4: PCA operation ---------------------------------------------------
axP = fig.add_axes([0.55, GLY_B, 0.11, GLY_T - GLY_B]); axP.axis("off")
axP.set_xlim(0, 1); axP.set_ylim(0, 1)
axP.text(0.5, 0.74, "PCA", fontsize=9, fontweight="bold", ha="center")
axP.text(0.5, 0.50, "eigen-decomposition\nof the sample\ncovariance", fontsize=5.0,
         ha="center", va="center")
axP.text(0.5, 0.18, f"PC1  {evr[0]*100:.0f}%\nPC2  {evr[1]*100:.0f}%\nvariance",
         fontsize=5.2, ha="center", va="center", fontweight="bold")
hbox(0.55, 0.11, "4  Decompose")

# ---- step 5: embedding -------------------------------------------------------
axE = fig.add_axes([0.685, GLY_B + 0.01, 0.15, GLY_T - GLY_B - 0.01])
axE.scatter(scores[:, 0], scores[:, 1], s=16, c="#8a8a8a",
            edgecolor="black", linewidth=0.25, zorder=2)
for i in range(3):
    p = pre_ids.index(examples[i])
    axE.scatter(scores[p, 0], scores[p, 1], s=30,
                c=[plt.cm.Dark2(i / 6)], edgecolor="black", linewidth=0.5, zorder=3)
    axE.annotate(f"S{i+1}", (scores[p, 0], scores[p, 1]), fontsize=4.6,
                 xytext=(3, 2), textcoords="offset points")
axE.set_xlabel(f"PC1 ({evr[0]*100:.0f}%)", fontsize=5.4)
axE.set_ylabel(f"PC2 ({evr[1]*100:.0f}%)", fontsize=5.4)
axE.tick_params(labelsize=4.2)
for sp in ["top", "right"]:
    axE.spines[sp].set_visible(False)
cv.text(0.76, GLY_B - 0.045, "one point = one sample", fontsize=5.0,
        ha="center", va="top", style="italic")
hbox(0.685, 0.15, "5  PCA embedding")

# ---- flow arrows between glyph steps -----------------------------------------
for x0, x1 in [(0.165, 0.20), (0.36, 0.40), (0.525, 0.55), (0.665, 0.685)]:
    cv.add_patch(FancyArrowPatch((x0, 0.66), (x1, 0.66),
                 transform=cv.transData, **ARR))

# ---- colour note (-> Fig. 6e, f) --------------------------------------------
cv.add_patch(FancyArrowPatch((0.855, 0.66), (0.885, 0.66),
             transform=cv.transData, **ARR))
cv.add_patch(FancyBboxPatch((0.865, 0.50), 0.125, 0.30,
             boxstyle="round,pad=0.006", facecolor="#f4f4f4",
             edgecolor="#999999", lw=0.6))
cv.text(0.9275, 0.775, "then colour\neach point by", fontsize=5.2, ha="center",
        va="top", fontweight="bold")
cv.text(0.872, 0.665, "•  subset proportion", fontsize=4.9, va="center")
cv.text(0.986, 0.635, "(Fig. 6e)", fontsize=4.4, va="center", ha="right",
        style="italic", color="#666666")
cv.text(0.872, 0.575, "•  clinical variable /", fontsize=4.9, va="center")
cv.text(0.882, 0.535, "responder status", fontsize=4.9, va="center")
cv.text(0.986, 0.515, "(Fig. 6f)", fontsize=4.4, va="center", ha="right",
        style="italic", color="#666666")

# =============================================================================
# bottom: numbered step captions (text-rich)
# =============================================================================
cv.plot([0.02, 0.98], [0.40, 0.40], color="#cccccc", lw=0.5)
steps = [
    ("1.  Patients.",
     f"Each of the {n_samp} pre-treatment (baseline) PBMC samples enters the analysis once; "
     "post-treatment samples are excluded so that the embedding reflects the untreated immune state."),
    ("2.  Cell-type proportions.",
     "Within each parent lineage (CD8+ T, CD4+ T, B, NK, myeloid) every fine-grained subset is expressed "
     f"as its fraction of that lineage, giving each sample a {n_sub}-dimensional composition vector "
     "pᵢ. This is the same within-lineage proportion used for the correlations in Fig. 6a–c."),
    ("3.  Compose & scale.",
     f"The vectors are stacked into a {n_samp}×{n_sub} samples×subsets matrix, and each subset "
     "(column) is centred and scaled to unit variance (z-score) so that rare and abundant subsets "
     "contribute comparably rather than the analysis being dominated by the largest populations."),
    ("4.  Decompose.",
     "Principal-component analysis of this matrix finds orthogonal axes (principal components) that "
     f"capture the greatest sample-to-sample variation in composition; PC1 and PC2 explain "
     f"{evr[0]*100:.0f}% and {evr[1]*100:.0f}% of the variance and are dominated by the "
     "naive↔effector-memory CD8+ contrast (loadings in Supplementary Fig. 7b)."),
    ("5.  Embedding.",
     "Each sample is projected onto PC1 and PC2, so that one point represents one patient and nearby points "
     "have similar immune composition. The same embedding is then coloured by an individual subset "
     "proportion (Fig. 6e) or by a clinical variable and baseline responder status (Fig. 6f)."),
]
ytop = 0.365
for i, (h, b) in enumerate(steps):
    y = ytop - i * 0.072
    cv.text(0.03, y, h, fontsize=5.6, fontweight="bold", va="top")
    cv.text(0.145, y, b, fontsize=5.2, va="top", wrap=True)
    # manual wrap width via textwrap is unreliable in data coords; rely on短 lines
cv.text(0.03, ytop - 5 * 0.072 - 0.010,
        "Draft schematic for placement in Fig. 6d; stacked bars and PCA scatter use real data "
        "(code/supple/fig6d_pca_schematic.py).",
        fontsize=4.3, va="top", style="italic", color="#888888")

fig.savefig(f"{OUTP}/fig6d_pca_schematic.pdf", dpi=300, bbox_inches="tight")
print(f"[fig] wrote {OUTP}/fig6d_pca_schematic.pdf  (n={n_samp}, subsets={n_sub}, "
      f"PC1={evr[0]*100:.1f}%, PC2={evr[1]*100:.1f}%)")
