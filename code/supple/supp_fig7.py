"""Supplementary Figure 7 (supports Main Figure 6) — one supple per main figure.

Main Figure 6 organises baseline immune composition along two axes and relates it
to clinical heterogeneity (correlation heatmap, force-directed network, selected
pairwise correlations, PCA of subset proportions). This supplement provides the
quantitative underpinning and honesty checks for that analysis, as a single
Illustrator-editable vector PDF:

  a  Clinical/laboratory variable completeness across the 23 pre-DMT samples
  b  PCA of baseline subset proportions: scree + PC1/PC2 loadings (which subsets
     build the axes coloured in Fig. 6d,e)
  c  Bootstrap 95% CIs for the two highlighted pairwise correlations (Fig. 6c)
  d  Multiple-testing calibration of the full subset x clinical correlation matrix
     (Fig. 6a is exploratory / uncorrected)

Correlations reproduce Fig. 6a,c exactly under: within-parent-lineage subset
proportions + two-sided Spearman (M_Mono_C1QA - B_Activated r = -0.516, P = 0.012).
Deterministic (seed fixed, RULES R4). Env = scanpy (RULES R2).

Outputs:
  figures/supplementary/Supplementary_Fig7.pdf
  data/supplementary_tables/SupplementaryData5_Fig6_correlations.csv (+ .xlsx, 2 sheets)
  data/results/supp_fig7_pca_loadings.csv
  data/results/supp_fig7_bootstrap.csv
  data/results/supp_fig7_fdr_summary.csv
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import sys, warnings, glob
sys.path.insert(0, os.path.join(MS_ROOT, "code"))
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from common import io, style, palette
style.apply()
BASE, TICK, PLET = style.BASE, style.TICK, style.PANEL_LETTER   # 6 / 5 / 12(bold)
RNG = np.random.default_rng(0)                                  # deterministic

OUT   = os.path.join(MS_ROOT, "figures/supplementary")
RES   = os.path.join(MS_ROOT, "data/results")
STAB  = os.path.join(MS_ROOT, "data/supplementary_tables")

# --- clinical / laboratory variables shown in Fig. 6a (blue labels) ------------
CLIN = ["AgeSampling", "AgeOnset", "onset_to_PBMC", "attack_to_sampling",
        "Onset_to_OCB", "EDSS", "No.ofT2lesionsbrain", "No.ofGdlesionsbrain",
        "No.ofinfratentoriallesion", "No.ofSClesion", "No.ofT1blackholes",
        "WBC", "Prot", "IgG_index", "CSF_Alb", "CSF_IgG", "Serum_Alb",
        "Serum_IgG", "time_to_fail_NEDA3_or_last_fu", "WBCtest",
        "Segmentedneutrophil", "Lymphocyte", "lymph_count", "V66"]
# pretty CSF/serum protein label used in Fig. 6c
NICE = {"Prot": "CSF protein", "IgG_index": "IgG index",
        "No.ofT2lesionsbrain": "No. T2 lesions (brain)",
        "time_to_fail_NEDA3_or_last_fu": "Time to NEDA-3 failure"}


# =============================================================================
# 1. per-sample feature matrices
# =============================================================================
meta = io.build_meta()
meta["Timepoint"] = meta["Timepoint"].astype(str)
pre_ids = meta.loc[meta["Timepoint"] == "Before", "SampleID"].astype(str).tolist()

# within-parent-lineage subset proportions (validated against Fig. 6a,c)
LIN = {"T_CD8": "T_CD8", "T_CD4": "T_CD4", "B": "B", "NK": "NK", "Myeloid": "Myeloid"}
prop = {}
for lin in LIN:
    ad = io.load_lineage(lin)
    sid = ad.obs["SampleID"].astype(str)
    lab = ad.obs["anno2"].astype(str)
    tot = sid.value_counts()
    for sub in sorted(lab.unique()):
        cnt = sid[lab == sub].value_counts()
        prop[sub] = (cnt / tot).reindex(tot.index).fillna(0.0)
SUB = pd.DataFrame(prop)                          # samples x subsets
SUB = SUB.reindex(pre_ids).dropna(how="all")      # 23 pre-DMT samples
# drop all-zero / near-constant subsets (uninformative for PCA/correlation)
SUB = SUB.loc[:, SUB.std() > 0]
pre_ids = SUB.index.tolist()
print(f"[data] pre-DMT samples n={len(pre_ids)}, subsets n={SUB.shape[1]}")

# clinical matrix (numeric coercion), same sample index
d1 = pd.read_excel(glob.glob(
    os.path.join(MS_ROOT, "originals/supplementary_data/*Data1*"))[0])
d1["SampleID"] = d1["SampleID"].astype(str)
d1 = d1.set_index("SampleID")
CLINm = d1[CLIN].apply(pd.to_numeric, errors="coerce").reindex(pre_ids)


# =============================================================================
# 2. full subset x clinical Spearman correlation table (Fig. 6a support)
# =============================================================================
rows = []
for sub in SUB.columns:
    for cv in CLIN:
        x = SUB[sub]
        y = CLINm[cv]
        d = pd.concat([x, y], axis=1).dropna()
        if len(d) < 5:                            # too few pairs to correlate
            continue
        r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
        rows.append({"subset": sub, "clinical": cv, "n": len(d),
                     "spearman_r": r, "p_two_sided": p})
corr = pd.DataFrame(rows)
corr["p_adj_BH"] = multipletests(corr["p_two_sided"], method="fdr_bh")[1]
corr = corr.sort_values("p_two_sided").reset_index(drop=True)
corr.to_csv(f"{STAB}/SupplementaryData5_Fig6_correlations.csv", index=False)

# subset x subset Spearman table (the immune-immune part of the Fig. 6a matrix;
# includes the M_Mono_C1QA - B_Activated headline of Fig. 6c). Supplementary Data 5
# is delivered as a 2-sheet xlsx: subset_x_clinical (above) + subset_x_subset (here).
import itertools
ss_rows = []
for a, b in itertools.combinations(SUB.columns, 2):
    d = pd.concat([SUB[a], SUB[b]], axis=1).dropna()
    if len(d) < 5:
        continue
    r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    ss_rows.append({"subset_1": a, "subset_2": b, "n": len(d),
                    "spearman_r": r, "p_two_sided": p})
ss = pd.DataFrame(ss_rows)
ss["p_adj_BH"] = multipletests(ss["p_two_sided"], method="fdr_bh")[1]
ss = ss.sort_values("p_two_sided").reset_index(drop=True)
with pd.ExcelWriter(f"{STAB}/SupplementaryData5_Fig6_correlations.xlsx",
                    engine="openpyxl") as _xw:
    corr.to_excel(_xw, sheet_name="subset_x_clinical", index=False)
    ss.to_excel(_xw, sheet_name="subset_x_subset", index=False)
print(f"[corr] subset x subset pairs={len(ss)} "
      f"(Data 5 = 2-sheet xlsx: subset_x_clinical + subset_x_subset)")

n_raw = int((corr["p_two_sided"] < 0.05).sum())
n_q05 = int((corr["p_adj_BH"] < 0.05).sum())
n_q10 = int((corr["p_adj_BH"] < 0.10).sum())
print(f"[corr] pairs={len(corr)}  raw p<0.05={n_raw}  BH q<0.05={n_q05}  q<0.10={n_q10}")
pd.DataFrame([{"n_pairs": len(corr), "raw_p_lt_0.05": n_raw,
               "BH_q_lt_0.05": n_q05, "BH_q_lt_0.10": n_q10}]).to_csv(
    f"{RES}/supp_fig7_fdr_summary.csv", index=False)

# headline pairs (Fig. 6c) — one subset-subset, one subset-clinical
HEAD = [("M_Mono_C1QA", "B_Activated", -0.516, 0.012),
        ("T_CD4_CXCR5", "Prot", 0.633, 0.004)]


def series(name):
    """Resolve a Fig. 6c variable to its per-sample vector (subset or clinical)."""
    return SUB[name] if name in SUB.columns else CLINm[name]


# =============================================================================
# 3. PCA of baseline subset proportions (Fig. 6d,e support)
# =============================================================================
Z = StandardScaler().fit_transform(SUB.values)
pca = PCA(n_components=min(8, Z.shape[0] - 1, Z.shape[1]), random_state=0).fit(Z)
scores = pca.transform(Z)
evr = pca.explained_variance_ratio_
load = pd.DataFrame(pca.components_[:2].T, index=SUB.columns, columns=["PC1", "PC2"])
load.to_csv(f"{RES}/supp_fig7_pca_loadings.csv")
print(f"[pca] PC1={evr[0]*100:.1f}%  PC2={evr[1]*100:.1f}%")


# =============================================================================
# 4. bootstrap 95% CI for the two headline correlations (Fig. 6c robustness)
# =============================================================================
def boot_ci(sub, cv, nboot=5000):
    d = pd.concat([series(sub).rename("x"), series(cv).rename("y")],
                  axis=1).dropna()
    a, b = d.iloc[:, 0].values, d.iloc[:, 1].values
    obs = stats.spearmanr(a, b)[0]
    n = len(a)
    bs = np.empty(nboot)
    for i in range(nboot):
        idx = RNG.integers(0, n, n)
        bs[i] = stats.spearmanr(a[idx], b[idx])[0]
    bs = bs[np.isfinite(bs)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return obs, lo, hi, n, bs

boot = {}
brows = []
for sub, cv, _, _ in HEAD:
    obs, lo, hi, n, bs = boot_ci(sub, cv)
    boot[(sub, cv)] = bs
    brows.append({"pair": f"{sub} - {NICE.get(cv, cv)}", "subset": sub,
                  "clinical": cv, "n": n, "spearman_r": obs,
                  "ci95_low": lo, "ci95_high": hi})
    print(f"[boot] {sub}-{cv}: r={obs:.3f} 95%CI[{lo:.3f},{hi:.3f}] n={n}")
pd.DataFrame(brows).to_csv(f"{RES}/supp_fig7_bootstrap.csv", index=False)


# =============================================================================
# 5. figure — single multi-panel vector PDF (panels a-d)
#    a completeness · b FDR calibration · c bootstrap CI · d scree
#    Panel e (PCA embedding coloured by clinical variables) is added by the
#    authors in Illustrator; the embedding coordinates match main Fig. 6e, f.
# =============================================================================
fig = plt.figure(figsize=(7.2, 5.4))
gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1.35, 1.0],
                       hspace=0.55, wspace=0.55,
                       left=0.12, right=0.97, top=0.94, bottom=0.11)


def panel_letter(ax, s, dx=-0.09, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=PLET,
            fontweight="bold", va="top", ha="right")


# --- a  clinical variable completeness ---------------------------------------
ax_a = fig.add_subplot(gs[0, :])
avail = CLINm.notna().sum().reindex(CLIN)
labels_a = [NICE.get(c, c) for c in avail.index]
ypos = np.arange(len(avail))[::-1]
colors_a = ["#4d4d4d" if v == len(pre_ids) else "#bdbdbd" for v in avail.values]
ax_a.barh(ypos, avail.values, color=colors_a, edgecolor="black", linewidth=0.3)
ax_a.axvline(len(pre_ids), color="black", ls=":", lw=0.6)
ax_a.set_yticks(ypos); ax_a.set_yticklabels(labels_a, fontsize=4.3)
ax_a.set_xlim(0, len(pre_ids) + 1)
ax_a.set_xlabel(f"Pre-DMT samples with measurement (of {len(pre_ids)})", fontsize=BASE)
ax_a.set_title("Clinical / laboratory variable completeness", fontsize=BASE, pad=3)
ax_a.tick_params(axis="x", labelsize=TICK)
for sp in ["top", "right"]:
    ax_a.spines[sp].set_visible(False)
panel_letter(ax_a, "a", dx=-0.02)

# --- c  bootstrap CI of headline correlations --------------------------------
ax_c = fig.add_subplot(gs[1, 1])
parts = ax_c.violinplot([boot[(s, cv)] for s, cv, _, _ in HEAD],
                        positions=[0, 1], vert=False, showextrema=False,
                        widths=0.7)
for pc in parts["bodies"]:
    pc.set_facecolor("#cfe2f3"); pc.set_edgecolor("black")
    pc.set_linewidth(0.3); pc.set_alpha(0.9)
for i, (sub, cv, rfig, pfig) in enumerate(HEAD):
    d = pd.concat([series(sub), series(cv)], axis=1).dropna()
    obs = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0]
    lo, hi = np.percentile(boot[(sub, cv)], [2.5, 97.5])
    ax_c.plot([lo, hi], [i, i], color="black", lw=1.0, zorder=4)
    ax_c.plot(obs, i, "o", color="#d73027", ms=4, zorder=5)
ax_c.axvline(0, color="grey", ls="--", lw=0.5)
ax_c.set_yticks([0, 1])
ax_c.set_yticklabels([f"{s}\n- {NICE.get(cv, cv)}" for s, cv, _, _ in HEAD],
                     fontsize=4.3)
ax_c.set_xlim(-1, 1)
ax_c.set_xlabel("Spearman r (5000-sample bootstrap)", fontsize=BASE)
ax_c.set_title("Robustness of highlighted correlations", fontsize=BASE, pad=3)
ax_c.tick_params(axis="x", labelsize=TICK)
for sp in ["top", "right"]:
    ax_c.spines[sp].set_visible(False)
panel_letter(ax_c, "c")

# --- b  multiple-testing calibration -----------------------------------------
ax_d = fig.add_subplot(gs[1, 0])
ax_d.hist(corr["p_two_sided"], bins=20, range=(0, 1), color="#bdbdbd",
          edgecolor="black", linewidth=0.3)
exp = len(corr) / 20.0                              # uniform-null expectation
ax_d.axhline(exp, color="black", ls=":", lw=0.6)
ax_d.axvline(0.05, color="#d73027", ls="--", lw=0.7)
ax_d.set_xlabel("Two-sided Spearman P value", fontsize=BASE)
ax_d.set_ylabel("Subset x clinical pairs", fontsize=BASE)
ax_d.set_title("Multiple-testing calibration", fontsize=BASE, pad=3)
ax_d.tick_params(labelsize=TICK)
txt = (f"{len(corr)} pairs\nraw P<0.05: {n_raw}\n"
       f"BH q<0.05: {n_q05}\nBH q<0.10: {n_q10}")
ax_d.text(0.97, 0.95, txt, transform=ax_d.transAxes, fontsize=4.6,
          ha="right", va="top")
for sp in ["top", "right"]:
    ax_d.spines[sp].set_visible(False)
panel_letter(ax_d, "b")

# --- e  PCA scree ------------------------------------------------------------
ax_e = fig.add_subplot(gs[1, 2])
k = len(evr)
ax_e.bar(np.arange(1, k + 1), evr * 100, color="#4d4d4d",
         edgecolor="black", linewidth=0.3)
ax_e.plot(np.arange(1, k + 1), np.cumsum(evr) * 100, "-o", color="#d73027",
          ms=2.5, lw=0.8)
ax_e.set_xlabel("Principal component", fontsize=BASE)
ax_e.set_ylabel("Variance explained (%)", fontsize=BASE)
ax_e.set_title("PCA scree", fontsize=BASE, pad=3)
ax_e.set_xticks(np.arange(1, k + 1))
ax_e.tick_params(labelsize=TICK)
for sp in ["top", "right"]:
    ax_e.spines[sp].set_visible(False)
panel_letter(ax_e, "d")

fig.savefig(f"{OUT}/Supplementary_Fig7.pdf", dpi=300, bbox_inches="tight")
print(f"[fig] wrote {OUT}/Supplementary_Fig7.pdf")

# reproduction check (anchor correlation, subset-subset)
d = pd.concat([series("M_Mono_C1QA"), series("B_Activated")], axis=1).dropna()
ar, ap = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
print(f"[check] M_Mono_C1QA-B_Activated: r={ar:.3f} p={ap:.3f} (figure -0.516 / 0.012)")
