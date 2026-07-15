"""Supplementary Figure 1 (supports Main Figure 1) — one supple per main figure.

Combines into a single vector, Illustrator-editable multi-panel PDF:
  a-d  QC (genes/cell, UMIs/cell, %mito, per-sample cell yield)
  e-f  batch/integration evidence (UMAP by sample, per-lineage sample entropy)
  g    baseline R-vs-NR differences robust to confounders (adjusted forest)
패널 순서 = 본문 등장 순서(QC → batch → confounder).
Canonical lineage/subset markers live in the dedicated marker figure
(Supplementary Fig. 2); baseline 임상특성표는 Supplementary Data 2(원본).
출력: figures/supplementary/Supplementary_Fig1.pdf
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
import seaborn as sns
import statsmodels.formula.api as smf
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.sparse import issparse

from common import io, style, naming, palette
style.apply()
BASE, TICK, PLET = style.BASE, style.TICK, style.PANEL_LETTER   # 6 / 5 / 12(bold)
OBJ = os.path.join(MS_ROOT, "data/objects")
OUT = os.path.join(MS_ROOT, "figures/supplementary")

# ---- load anno1 (clean object: QC in obs, X log-norm, anno1 ordered + colors) ----
a1 = sc.read_h5ad(f"{OBJ}/MS_scRNA_anno1_total.h5ad")
io.join_meta(a1)
order = [c for c in naming.LEVEL1 if c in a1.obs["anno1"].astype(str).unique()]
pal = palette.load_colors("anno1")
df = a1.obs.copy(); df["anno1"] = df["anno1"].astype(str)
for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# ---- figure scaffold ----
fig = plt.figure(figsize=(12, 13.5))
gs = fig.add_gridspec(4, 4, height_ratios=[1, 1, 1.25, 1.15], hspace=0.55, wspace=0.5)


def vln(ax, col, ylab, logy=False):
    data = [df.loc[df.anno1 == c, col].values for c in order]
    parts = ax.violinplot(data, showmedians=True, widths=0.85)
    for pc, c in zip(parts["bodies"], order):
        pc.set_facecolor(pal.get(c, "#999")); pc.set_edgecolor("black"); pc.set_alpha(.9); pc.set_linewidth(.4)
    for k in ["cmedians", "cbars", "cmins", "cmaxes"]:
        if k in parts: parts[k].set_color("black"); parts[k].set_linewidth(.6)
    ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=90, fontsize=TICK)
    ax.set_ylabel(ylab, fontsize=BASE)
    if logy: ax.set_yscale("log")
    sns.despine(ax=ax)


# a-d QC
axa = fig.add_subplot(gs[0, 0:2]); vln(axa, "n_genes_by_counts", "Genes per cell")
axb = fig.add_subplot(gs[0, 2:4]); vln(axb, "total_counts", "UMIs per cell (log)", logy=True)
axc = fig.add_subplot(gs[1, 0:2]); vln(axc, "pct_counts_mt", "% mitochondrial")
axd = fig.add_subplot(gs[1, 2:4])
df["Group1"] = df["Group1"].astype(str); df["SampleID"] = df["SampleID"].astype(str)
sc_n = df.groupby(["SampleID", "Group1"], observed=True).size().rename("n").reset_index()
gcol = {"Responder_Before": "#66c2a5", "Nonresponder_Before": "#fc8d62",
        "Responder_After": "#b3e2cd", "Nonresponder_After": "#fdcdac"}
sc_n["c"] = sc_n["Group1"].map(gcol).fillna("#cccccc"); sc_n = sc_n.sort_values("n")
axd.bar(range(len(sc_n)), sc_n["n"], color=sc_n["c"], edgecolor="black", linewidth=.3)
axd.set_xlabel(f"samples (n={len(sc_n)})", fontsize=BASE); axd.set_ylabel("cells recovered", fontsize=BASE)
axd.set_title(f"Total {len(df):,} cells", fontsize=BASE)
axd.legend(handles=[Patch(facecolor=v, edgecolor="black", label=k.replace("_", " "))
                    for k, v in gcol.items()], frameon=False, fontsize=BASE, loc="upper left")
sns.despine(ax=axd)

# e-f batch (Harmony) integration evidence: lineages are not sample/batch artifacts
um = a1.obsm["X_umap"]; sids = a1.obs["SampleID"].astype(str)
uniq = sorted(sids.unique())
# 41개 구분되는 색(흰색 회피): tab20+tab20b+tab20c(=60) 앞부분 사용
disc = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
scol = {s: disc[i] for i, s in enumerate(uniq)}
axf = fig.add_subplot(gs[2, 0:2])
axf.scatter(um[:, 0], um[:, 1], s=1.2, c=[scol[s] for s in sids],
            linewidth=0, rasterized=True)
axf.set_title(f"UMAP colored by sample (n={len(uniq)}) — no per-sample islands", fontsize=BASE)
axf.axis("off")
sleg = [Line2D([0], [0], marker="o", linestyle="none", markersize=2.2,
               markerfacecolor=scol[s], markeredgecolor="none", label=s) for s in uniq]
axf.legend(handles=sleg, ncol=7, fontsize=4, frameon=False, loc="upper center",
           bbox_to_anchor=(0.5, -0.02), handletextpad=0.15, columnspacing=0.5,
           labelspacing=0.25, title="SampleID (n=41)", title_fontsize=TICK)
axg = fig.add_subplot(gs[2, 2:4])
ct = pd.crosstab(df["anno1"], a1.obs["SampleID"].astype(str))
p = ct.div(ct.sum(1), axis=0)
H = -(p.where(p > 0) * np.log(p.where(p > 0))).sum(1)
ent = (H / np.log(ct.shape[1])).reindex(order)
axg.barh(range(len(order)), ent.values[::-1],
         color=[pal.get(c, "#999") for c in order[::-1]], edgecolor="black", linewidth=.4)
axg.set_yticks(range(len(order))); axg.set_yticklabels(order[::-1], fontsize=TICK)
axg.axvline(1, color="grey", lw=.5, ls="--")
axg.set_xlim(0, 1.05); axg.set_xlabel("normalized sample entropy\n(1 = evenly mixed across samples)", fontsize=BASE)
axg.set_title("Every lineage draws from most samples", fontsize=BASE); sns.despine(ax=axg)

# h confounder-adjusted forest
cd8 = sc.read_h5ad(f"{OBJ}/MS_scRNA_anno2_T_CD8.h5ad", backed="r")
cd4 = sc.read_h5ad(f"{OBJ}/MS_scRNA_anno2_T_CD4.h5ad", backed="r")
lin_prop = pd.crosstab(a1.obs["SampleID"].astype(str), df["anno1"], normalize=0)
lin_n = pd.crosstab(a1.obs["SampleID"].astype(str), df["anno1"])
cd8p = pd.crosstab(cd8.obs["SampleID"].astype(str), cd8.obs["anno2"].astype(str), normalize=0)
ctl_n = cd8.obs[cd8.obs["anno2"].astype(str) == "T_CD4_CTL"].groupby(
    cd8.obs["SampleID"].astype(str), observed=True).size()
ctl_prop = (ctl_n / lin_n["T_CD4"]).reindex(lin_n.index).fillna(0.0)
FEAT = {"CD8 lineage": lin_prop["T_CD8"], "NK lineage": lin_prop["NK"],
        "CD8 EM ZNF683(hi)": cd8p["T_CD8_EM_ZNF683(hi)"], "CD8 EM GZMK(hi)": cd8p["T_CD8_EM_GZMK(hi)"],
        "CD8 EM KIR": cd8p["T_CD8_EM_KIR"], "CD8 Naive/CM": cd8p["T_CD8_Naive/CM"], "CD4 CTL": ctl_prop}
feat = pd.DataFrame(FEAT); feat.index.name = "SampleID"
m = io.build_meta(); m["dur"] = m["AgeSampling"] - m["AgeOnset"]
cov = m.set_index("SampleID")[["NEDA3_1yr", "Timepoint", "AgeSampling", "Sex", "DMTint", "dur"]]
fd = feat.join(cov)
fd = fd[(fd.Timepoint == "Before") & fd.NEDA3_1yr.isin(["Responder", "Nonresponder"])].copy()
fd["Resp"] = (fd.NEDA3_1yr == "Responder").astype(int)
rows = []
for name in FEAT:
    d = fd.rename(columns={name: "y"}).dropna(subset=["y", "AgeSampling", "Sex", "DMTint", "dur"])
    mm = smf.ols("y ~ Resp + AgeSampling + C(Sex) + C(DMTint) + dur", data=d).fit()
    ci = mm.conf_int().loc["Resp"]
    rows.append((name, mm.params["Resp"], mm.pvalues["Resp"], ci[0], ci[1]))
res = pd.DataFrame(rows, columns=["feature", "beta", "p", "lo", "hi"])
axh = fig.add_subplot(gs[3, 1:3])
yy = np.arange(len(res))[::-1]
axh.axvline(0, color="grey", lw=.6, ls="--")
for yi, (_, r) in zip(yy, res.iterrows()):
    c = "#2c7fb8" if r.beta > 0 else "#d95f0e"; sig = r.p < 0.05
    axh.plot([r.lo, r.hi], [yi, yi], color=c, lw=1.4, alpha=.9)
    axh.scatter(r.beta, yi, s=40, color=c, edgecolor="black", lw=.5, zorder=5, marker="o" if sig else "x")
    star = "**" if r.p < 0.01 else "*" if r.p < 0.05 else "ns"
    axh.text(r.hi, yi, f"  p={r.p:.3f} {star}", va="center", fontsize=BASE)
axh.set_yticks(yy); axh.set_yticklabels(res.feature, fontsize=TICK)
axh.set_xlabel("Adjusted Responder effect (β, 95% CI;\n+Age +Sex +DMT +disease duration)", fontsize=BASE)
axh.set_title("Baseline R-vs-NR robust to confounders", fontsize=BASE)
sns.despine(ax=axh)

# panel labels (볼드 12)
for ax, lab in [(axa, "a"), (axb, "b"), (axc, "c"), (axd, "d"),
                (axf, "e"), (axg, "f"), (axh, "g")]:
    ax.text(-0.06, 1.04, lab, transform=ax.transAxes, fontsize=PLET, fontweight="bold",
            color="black", va="bottom", ha="right")

plt.savefig(f"{OUT}/Supplementary_Fig1.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("saved: figures/supplementary/Supplementary_Fig1.pdf (panels a-g)")
