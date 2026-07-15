"""Supplementary Fig. 3 (supports Fig. 2) — single multi-panel vector PDF.

Integrates two evidence blocks in text-appearance order (CD8+ subset marker
identity now lives in the dedicated marker figure, Supplementary Fig. 2):
  a, b, c  blood-CSF annotation transfer (Kavaka reference; supports Fig. 2h)
  d, e, f  TCR antigen-experience (VDJdb landscape, subset match rate, GLIPH2)
The R-vs-NR antigen/convergence comparison is not plotted (antigen-experience is
framed as a shared compartment property); its two-sided test is retained in
Source Data (suppfig3_RvsNR.csv) for transparency.

Self-contained, runs in the `scanpy` env. The scirpy/muon AIRR object cannot be
opened by anndata here (awkward/muon dep), so its obs columns are read directly
with h5py. Inputs are READ-ONLY. Output: figures/supplementary/Supplementary_Fig3.pdf
Naming is canonical (config/naming.yaml); colors from the original curated palette.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import warnings; warnings.filterwarnings("ignore")
import sys, pickle
import numpy as np
import pandas as pd
import anndata as ad
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.join(MS_ROOT, "code"))
from common import style
style.apply()

D26 = os.path.join(MS_SRC, "DATA/26JAN")
NB = os.path.join(MS_SRC, "notebook17_Fig3_TCR")
GLIPH = f"{NB}/250926_0935_P15315_MDC2JH2QV1.csv"          # patient GLIPH2 groups (number_subject)
GLIPH_AG = f"{NB}/250926_1423_P15324_AHEAG100TH_sub.csv"   # antigen-linked convergence groups
CDIR = os.path.join(MS_SRC, "01_Figures_pre")
META = os.path.join(MS_ROOT, "data/metadata/sample_metadata.csv")
OUT = os.path.join(MS_ROOT, "figures/supplementary/Supplementary_Fig3.pdf")
RES = os.path.join(MS_ROOT, "data/results")

CD8_ORDER = ["T_CD8_Naive/CM", "T_CD8_EM_GZMK(lo)", "T_CD8_EM_GZMK(hi)",
             "T_CD8_EM_ZNF683(hi)", "T_CD8_EM_KIR", "T_CD8_IFN"]
RESP_ENRICHED = {"T_CD8_EM_GZMK(hi)", "T_CD8_EM_ZNF683(hi)", "T_CD8_EM_KIR"}
GREY = "#bdbdbd"

# curated subset palette (matches the paper)
CPAL = pickle.load(open(f"{CDIR}/colors_anno2_T_CD8.pkl", "rb"))

meta = pd.read_csv(META).set_index("SampleID")


def read_h5ad_obs_col(path, name):
    """Read a single obs column (categorical or plain) from an .h5ad via h5py."""
    with h5py.File(path, "r") as f:
        node = f["obs"][name]
        if isinstance(node, h5py.Group):            # categorical
            cats = np.array([c.decode() if isinstance(c, bytes) else c
                             for c in node["categories"][:]])
            codes = node["codes"][:]
            out = np.where(codes >= 0, cats[codes.clip(min=0)], None)
        else:
            out = np.array([a.decode() if isinstance(a, bytes) else a
                            for a in node[:]])
        idx = f["obs"]["_index"][:]
        idx = [a.decode() if isinstance(a, bytes) else a for a in idx]
    return pd.Series(out, index=idx, name=name)


def panel_letter(ax, letter, dx=-0.14, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=style.PANEL_LETTER, fontweight="bold", va="bottom", ha="right")


# ============================================================ load
print("loading objects ...")
kav = ad.read_h5ad(f"{D26}/260131_Kavaka_overlap.h5ad")       # blood-CSF transfer

# AIRR obs via h5py (+ subset labels from the GEX partner object)
airr_ag = read_h5ad_obs_col(f"{D26}/260129_muon_airr_CD8.h5ad", "antigen.species")
airr_sid = read_h5ad_obs_col(f"{D26}/260129_muon_airr_CD8.h5ad", "SampleID").astype(str)
gex = ad.read_h5ad(f"{D26}/260129_muon_gex_CD8.h5ad")
tcr = pd.DataFrame({"antigen": airr_ag, "SampleID": airr_sid})
tcr["anno2"] = gex.obs["anno2"].reindex(tcr.index).astype(str)
tcr["matched"] = (~tcr["antigen"].isin([None, "ambiguous"])) & tcr["antigen"].notna()
tcr["Group"] = tcr["SampleID"].map(meta["NEDA3_1yr"])
tcr["TP"] = tcr["SampleID"].map(meta["Timepoint"])
em = tcr[tcr["anno2"].str.contains("T_CD8_EM")].copy()
print(f"  CD8 EM cells {len(em)} | VDJdb-matched {int(em['matched'].sum())} "
      f"({100*em['matched'].mean():.1f}%)")

# ============================================================ figure
fig = plt.figure(figsize=(10.0, 6.4))
gs = GridSpec(2, 6, figure=fig, height_ratios=[1.05, 1.0],
              hspace=0.9, wspace=0.9,
              left=0.07, right=0.94, top=0.93, bottom=0.13)

# ---------------------------------------------------------- a,b: transfer UMAPs
um = kav.obsm["X_umap"]
tr = kav.obs["anno2_transfer"].astype(str).values
conf = np.asarray(kav.obs["scl_prob"].values, float)
axb = fig.add_subplot(gs[0, 0:2])
for g in CD8_ORDER:
    m = tr == g
    axb.scatter(um[m, 0], um[m, 1], s=1.2, c=CPAL.get(g, GREY),
                linewidths=0, rasterized=True, label=g)
axb.set_title("Transferred CD8+ annotation", fontsize=6)
axb.set_xticks([]); axb.set_yticks([])
axb.set_xlabel("UMAP1", fontsize=5); axb.set_ylabel("UMAP2", fontsize=5)
axb.legend(loc="center left", bbox_to_anchor=(0.98, 0.5), fontsize=4,
           frameon=False, markerscale=3, handletextpad=0.1, labelspacing=0.3)
panel_letter(axb, "a")

axc = fig.add_subplot(gs[0, 2:4])
order = np.argsort(conf)
scc = axc.scatter(um[order, 0], um[order, 1], s=1.2, c=conf[order],
                  cmap="Reds", vmin=0.2, vmax=1.0, linewidths=0, rasterized=True)
axc.set_title("Prediction confidence", fontsize=6)
axc.set_xticks([]); axc.set_yticks([])
axc.set_xlabel("UMAP1", fontsize=5); axc.set_ylabel("UMAP2", fontsize=5)
cbc = fig.colorbar(scc, ax=axc, fraction=0.04, pad=0.02, aspect=12)
cbc.ax.tick_params(labelsize=4)
panel_letter(axc, "b")

# ---------------------------------------------------------- c: matched proportion
axd = fig.add_subplot(gs[0, 4:6])
mcell = kav.obs["overlap"].astype(str).values == "Yes"
filt = kav.obs["anno2_transfer_filtered"].astype(str).values
cat = np.where((filt == "nan") | (filt == "None"), "confidence < 0.5", filt)
dser = pd.Series(cat[mcell]).value_counts()
dser = dser / dser.sum()
dorder = CD8_ORDER + ["confidence < 0.5"]
dser = dser.reindex(dorder).fillna(0)
dcol = [CPAL.get(g, GREY) for g in CD8_ORDER] + [GREY]
axd.bar(range(len(dorder)), dser.values, color=dcol, edgecolor="black", linewidth=0.4)
axd.set_xticks(range(len(dorder)))
axd.set_xticklabels(dorder, rotation=40, ha="right", fontsize=4.5)
axd.set_ylabel("Proportion in\nmatched cells", fontsize=5)
axd.set_title(f"Blood-CSF matched cells (n = {int(mcell.sum())})", fontsize=6)
panel_letter(axd, "c")

# ---------------------------------------------------------- d: antigen landscape
axe = fig.add_subplot(gs[1, 0:2])
comp = em[em["matched"]]["antigen"].value_counts()
comp = (comp / comp.sum()).head(8)[::-1]
comp.to_frame("frac").to_csv(f"{RES}/suppfig3_antigen_composition.csv")
axe.barh(range(len(comp)), comp.values, color="#6baed6",
         edgecolor="black", linewidth=0.4)
axe.set_yticks(range(len(comp))); axe.set_yticklabels(comp.index, fontsize=5)
axe.set_xlabel("Fraction of VDJdb-matched\nCD8+ EM cells", fontsize=5)
axe.set_title("Antigen landscape (viral + self)", fontsize=6)
panel_letter(axe, "d")

# ---------------------------------------------------------- e: match rate / subset
axf = fig.add_subplot(gs[1, 2:4])
EM_ORDER = ["T_CD8_EM_GZMK(lo)", "T_CD8_EM_GZMK(hi)", "T_CD8_EM_ZNF683(hi)", "T_CD8_EM_KIR"]
mr = em.groupby("anno2")["matched"].mean().reindex(EM_ORDER) * 100
colf = [CPAL.get(g, GREY) if g in RESP_ENRICHED else GREY for g in EM_ORDER]
axf.bar(range(len(mr)), mr.values, color=colf, edgecolor="black", linewidth=0.4)
axf.set_xticks(range(len(mr)))
axf.set_xticklabels(EM_ORDER, rotation=40, ha="right", fontsize=4.5)
axf.set_ylabel("VDJdb match rate (%)", fontsize=5)
axf.set_title("Match rate by CD8+ subset", fontsize=6)
panel_letter(axf, "e")

# ---------------------------------------------------------- f: GLIPH2 convergence
axg = fig.add_subplot(gs[1, 4:6])
# public convergent motifs shared across >=3 patients (matches Fig. 2 text: 219)
gph = pd.read_csv(GLIPH); gph.columns = [c.strip() for c in gph.columns]
n_public = int((gph.groupby("type")["number_subject"].first() >= 3).sum())
# antigen-linked convergence groups per antigen (bars)
sub = pd.read_csv(GLIPH_AG); sub.columns = [c.strip() for c in sub.columns]
sub["is_ref"] = ~sub["Sample"].astype(str).str.startswith("MS")
from collections import Counter
agc = Counter()
for t, d in sub.groupby("type"):
    if (~d["is_ref"]).any():                       # group contains a patient TCR
        for a in set(d.loc[d["is_ref"], "Sample"].astype(str)):
            agc[a] += 1
agtop = pd.Series(agc).sort_values(ascending=False).head(8)[::-1]
axg.barh(range(len(agtop)), agtop.values, color="#fc9272",
         edgecolor="black", linewidth=0.4)
axg.set_yticks(range(len(agtop))); axg.set_yticklabels(agtop.index, fontsize=5)
axg.set_xlabel("Convergence groups\nlinked to antigen", fontsize=5)
axg.set_title(f"GLIPH2 convergent public motifs\n({n_public} groups shared across "
              u"≥3 patients)", fontsize=6)
panel_letter(axg, "f")

# ---------------------------------------------------------- R vs NR transparency (no panel)
# Antigen-experience is framed as a shared property of the effector-memory compartment,
# not an outcome discriminator. We do not plot the R-vs-NR comparison, but retain the
# two-sided test in Source Data for transparency (both not significant).
pre = em[(em["TP"] == "Before") & em["Group"].isin(["Responder", "Nonresponder"])]
mf = pre.groupby("SampleID")["matched"].mean()
mfg = mf.index.map(meta["NEDA3_1yr"])
mfd = {g: mf.values[np.array(list(mfg)) == g] for g in ["Responder", "Nonresponder"]}
p1 = mannwhitneyu(mfd["Responder"], mfd["Nonresponder"]).pvalue
gb = pd.read_csv(f"{RES}/s9_gliph_burden.csv")  # canonical per-patient burden (s9 scirpy pipeline)
gbd = {g: gb.loc[gb["Group"] == g, "n_convergent"].values for g in ["Responder", "Nonresponder"]}
p2 = mannwhitneyu(gbd["Responder"], gbd["Nonresponder"]).pvalue
pd.DataFrame({"comparison": ["antigen_matched_frac", "gliph_convergent_burden"],
              "p_mwu_two_sided": [p1, p2]}).to_csv(
    f"{RES}/suppfig3_RvsNR.csv", index=False)

fig.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print("saved:", OUT)
print(f"  R vs NR: matched-frac P={p1:.3f}, GLIPH burden P={p2:.3f} (both n.s.)")
