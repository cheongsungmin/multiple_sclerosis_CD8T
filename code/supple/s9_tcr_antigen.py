"""Supple S9: the responder-enriched CD8+ effector-memory repertoire is antigen-experienced.

Purpose (narrative, completes Fig2 "antigen-driven"): CD8+ EM cells carry TCRs
that (i) match defined VDJdb antigens (chronic viruses + self) and (ii) form
GLIPH2 convergent public motifs shared across patients and linked to those
antigens. This is a property of the differentiated/expanded compartment; it is
NOT presented as a responder-vs-non-responder discriminator (that difference is
in abundance, shown in Fig 1-2, and antigen fractions here do not differ by
outcome — reported in data/results/s9_antigen_RvsNR.csv for transparency).
Inputs (READ-ONLY): 260129_muon_airr_CD8.h5ad (precomputed scirpy VDJdb identity
match) + 260129_muon_gex_CD8.h5ad (anno2) + GLIPH2 patient/antigen outputs.
Runs in scirpy env. Output: figures/supplementary/S5_tcr_antigen.pdf, data/results/s9_*.csv
(component of Supplementary Fig. 2 — to be integrated with the CD8 marker and CSF panels in the Fig2 session)
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.stats import mannwhitneyu

mpl.rcParams["pdf.fonttype"] = 42
for k in ["text.color", "axes.edgecolor", "axes.labelcolor", "xtick.color", "ytick.color"]:
    mpl.rcParams[k] = "black"
sns.set_style("white"); sns.set_style("ticks")

D26 = os.path.join(MS_SRC, "DATA/26JAN")
NB = os.path.join(MS_SRC, "notebook17_Fig3_TCR")
GLIPH = f"{NB}/250926_0935_P15315_MDC2JH2QV1.csv"       # patient (Before, SID:Group)
GLIPH_AG = f"{NB}/250926_1423_P15324_AHEAG100TH_sub.csv"  # antigen-linked groups
META = os.path.join(MS_ROOT, "data/metadata/sample_metadata.csv")
OUT = os.path.join(MS_ROOT, "figures/supplementary")   # Supplementary Fig. S5
RES = os.path.join(MS_ROOT, "data/results")

RESP_ENRICHED = {"T_CD8_EM_ZNF683(hi)", "T_CD8_EM_GZMK(hi)", "T_CD8_EM_KIR"}
CD8_ORDER = ["T_CD8_Naive/CM", "T_CD8_EM_GZMK(lo)", "T_CD8_EM_GZMK(hi)",
             "T_CD8_EM_ZNF683(hi)", "T_CD8_EM_KIR", "T_CD8_IFN"]

# ---- load + join ----
airr = ad.read_h5ad(f"{D26}/260129_muon_airr_CD8.h5ad")
gex = ad.read_h5ad(f"{D26}/260129_muon_gex_CD8.h5ad")
meta = pd.read_csv(META).set_index("SampleID")
o = airr.obs.copy()
o["anno2"] = gex.obs["anno2"].reindex(o.index).astype(str)
o["SampleID"] = o["SampleID"].astype(str)
o["matched"] = o["antigen.species"].notna() & (o["antigen.species"] != "ambiguous")
em = o[o["anno2"].str.contains("T_CD8_EM")].copy()
print("CD8 EM cells:", len(em), "| VDJdb-matched:", int(em['matched'].sum()),
      f"({100*em['matched'].mean():.1f}%)")

# ---- (transparency) R vs NR antigen fractions at baseline: report only ----
em["Group"] = em["SampleID"].map(meta["NEDA3_1yr"]); em["TP"] = em["SampleID"].map(meta["Timepoint"])
pre = em[(em.TP == "Before") & em.Group.isin(["Responder", "Nonresponder"])]
tot = pre.groupby("SampleID").size()
cc = pre[pre.matched].groupby(["SampleID", "antigen.species"]).size().rename("n").reset_index()
cc["prop"] = cc["n"] / cc["SampleID"].map(tot); cc["Group"] = cc.SampleID.map(meta["NEDA3_1yr"])
rr = []
for a in sorted(cc["antigen.species"].unique()):
    x = cc[(cc.Group == "Responder") & (cc["antigen.species"] == a)]["prop"]
    y = cc[(cc.Group == "Nonresponder") & (cc["antigen.species"] == a)]["prop"]
    if len(x) > 1 and len(y) > 1:
        rr.append((a, x.mean(), y.mean(), mannwhitneyu(x, y).pvalue))
pd.DataFrame(rr, columns=["antigen", "mean_R", "mean_NR", "p_mwu"]).to_csv(
    f"{RES}/s9_antigen_RvsNR.csv", index=False)

# ---- Panel A: antigen-species composition of matched CD8 EM cells ----
comp = (em[em.matched]["antigen.species"].value_counts()
        .rename_axis("antigen").reset_index(name="cells"))
comp["frac"] = comp["cells"] / comp["cells"].sum()
comp.to_csv(f"{RES}/s9_antigen_composition.csv", index=False)
topA = comp.head(8)[::-1]

# ---- Panel B: VDJdb match rate per CD8 subset ----
mr = em.groupby("anno2")["matched"].mean().reindex(CD8_ORDER).dropna()

# ---- Panel C: GLIPH2 convergence publicness + antigen links ----
g = pd.read_csv(GLIPH); g.columns = [c.strip() for c in g.columns]
nsub = g.groupby("type")["number_subject"].first()
sub = pd.read_csv(GLIPH_AG); sub.columns = [c.strip() for c in sub.columns]
sub["is_ref"] = ~sub["Sample"].astype(str).str.startswith("MS")
# antigen-linked convergence groups: contain both a patient TCR and an antigen ref
ag_link = (sub.groupby("type")
           .apply(lambda d: pd.Series({"has_pat": (~d.is_ref).any(),
                                       "ags": ",".join(sorted(set(d.loc[d.is_ref, "Sample"].astype(str))))}))
           .reset_index())
ag_link = ag_link[ag_link["has_pat"] & (ag_link["ags"] != "")]
# count groups per antigen (split multi)
from collections import Counter
agc = Counter()
for s in ag_link["ags"]:
    for a in set(s.split(",")):
        agc[a] += 1
agtop = pd.Series(agc).sort_values(ascending=False).head(8)[::-1]
print(f"GLIPH: {g['type'].nunique()} groups, {(nsub>=3).sum()} public(>=3 subj); "
      f"{len(ag_link)} antigen-linked groups")

# ---- figure ----
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))
# A
axes[0].barh(topA["antigen"], topA["frac"], color="#6baed6", edgecolor="black", linewidth=.5)
axes[0].set_xlabel("Fraction of VDJdb-matched\nCD8+ EM cells")
axes[0].set_title("Antigen landscape (viral + self)", fontsize=8); sns.despine(ax=axes[0])
# B
colB = ["#66c2a5" if s in RESP_ENRICHED else "#bdbdbd" for s in mr.index]
axes[1].bar(range(len(mr)), mr.values * 100, color=colB, edgecolor="black", linewidth=.5)
axes[1].set_xticks(range(len(mr))); axes[1].set_xticklabels(mr.index, rotation=40, ha="right", fontsize=6)
axes[1].set_ylabel("VDJdb match rate (%)")
axes[1].set_title("Match rate by CD8 subset\n(green = responder-enriched)", fontsize=8); sns.despine(ax=axes[1])
# C
axes[2].barh(agtop.index, agtop.values, color="#fc9272", edgecolor="black", linewidth=.5)
axes[2].set_xlabel("GLIPH2 convergence groups\nlinked to antigen")
axes[2].set_title(f"Convergent public motifs\n({(nsub>=3).sum()} groups ≥3 subjects)", fontsize=8)
sns.despine(ax=axes[2])
plt.tight_layout()
plt.savefig(f"{OUT}/S5_tcr_antigen.pdf", dpi=300, bbox_inches="tight")
plt.close()
print("saved: S5_tcr_antigen.pdf + s9_antigen_composition.csv, s9_antigen_RvsNR.csv")
