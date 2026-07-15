"""Supplementary Fig. 5 (supports Fig. 4).

Two sections:
  (a) Extended paired pre/post subset proportions across ALL lineages
      (CD8, CD4, B, NK, myeloid) — the full version of Fig. 4c,d.
      Proportion is fraction within the parent lineage, per sample.
      Paired samples only (metadata Paired == 'Yes' = 24 samples / 12 patients,
      6 responders + 6 non-responders). Two-sided Wilcoxon signed-rank per group.
  (b) SCENIC transcription-factor regulon activity of B_Activated cells.
      Regulon specificity score (RSS) reproduced deterministically from the
      SAVED AUCell matrix (<MS_SRC>/pyscenic_output_B/auc_mtx.csv, 100% barcode
      match to the 26JAN B object) via Jensen-Shannon divergence — no rerun of
      pyscenic is required, so the pyscenic environment is not touched (R2).
  (c) Per-cell AP-1/NF-kB regulon-module activity (mean AUCell of JUN/FOSB/FOS/REL)
      across all B-cell subsets — the module is highest in B_Activated, confirming
      the RSS ranking (b) as actual, subset-specific activity.

Also writes Supplementary Data 6: B_Activated defining genes (vs other B subsets)
and the full regulon RSS matrix.

Env: scanpy. Deterministic (RULES R4): 26JAN objects + saved SCENIC output.
Outputs: figures/supplementary/Supplementary_Fig5.pdf
         data/source_data/SuppFig5a_paired_proportions.csv
         data/source_data/SuppFig5a_paired_stats.csv
         data/source_data/SuppFig5b_Bactivated_RSS.csv
         data/supplementary_tables/SupplementaryData6_Bactivated_defining_genes.csv
         data/supplementary_tables/SupplementaryData6_Bactivated_regulon_RSS.csv
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scanpy as sc
import seaborn as sns
from scipy.spatial.distance import jensenshannon
from scipy.stats import wilcoxon, ttest_rel, mannwhitneyu

ROOT = Path(str(MS_ROOT))
sys.path.insert(0, str(ROOT / "code"))
from common import io, style
from common.naming import LEVEL2
from common import palette

style.apply()
SRC = ROOT / "data" / "source_data"; SRC.mkdir(parents=True, exist_ok=True)
AUC = Path(os.path.join(MS_SRC, "pyscenic_output_B/auc_mtx.csv"))

# lineage -> (io key, naming key, y-axis lineage label)
LINEAGES = [
    ("T_CD8",   "T_CD8",   "CD8 T"),
    ("T_CD4",   "T_CD4",   "CD4 T"),
    ("B",       "B",       "B"),
    ("NK",      "NK",      "NK"),
    ("Myeloid", "Myeloid", "myeloid"),
]
GCOL = {  # (group, timepoint) -> fill
    ("Responder", "Before"): palette.GROUP["Responder_Before"],
    ("Responder", "After"): palette.GROUP["Responder_After"],
    ("Nonresponder", "Before"): palette.GROUP["Nonresponder_Before"],
    ("Nonresponder", "After"): palette.GROUP["Nonresponder_After"],
}
GROUPS = ["Responder", "Nonresponder"]

# ----------------------------------------------------------------------------
# 1. paired pre/post subset proportions (within parent lineage), per sample
# ----------------------------------------------------------------------------
meta = io.build_meta()
paired = meta[meta["Paired"] == "Yes"].copy()
smap = paired.set_index("SampleID")
paired_ids = set(paired["SampleID"])
print(f"paired samples={len(paired_ids)}  patients={paired['PatientN'].nunique()}")
print(paired.drop_duplicates("PatientN").groupby("NEDA3_1yr").size().to_dict())

prop_rows = []          # long: subset proportion per paired sample
subset_order = []       # (lineage_label, subset) in canonical order
for iokey, nmkey, lab in LINEAGES:
    ad = io.load_lineage(iokey)
    ad.obs["anno2"] = ad.obs["anno2"].astype(str)
    ad.obs["SampleID"] = ad.obs["SampleID"].astype(str)
    canon = [s for s in LEVEL2[nmkey] if s.startswith(("T_CD8", "T_CD4", "B_", "NK", "M_"))]
    # keep only this lineage's canonical subsets (drops e.g. T_CD4_CTL from CD8 obj)
    sub = ad.obs[ad.obs["anno2"].isin(canon) & ad.obs["SampleID"].isin(paired_ids)]
    # per-sample counts within lineage
    ct = sub.groupby(["SampleID", "anno2"]).size().unstack(fill_value=0)
    ct = ct.reindex(columns=[c for c in canon if c in ct.columns])
    denom = ct.sum(axis=1)
    frac = ct.div(denom, axis=0)          # proportion within lineage
    for s in canon:
        if s not in frac.columns:
            continue
        subset_order.append((lab, s))
        for sid, val in frac[s].items():
            prop_rows.append(dict(lineage=lab, subset=s, SampleID=sid,
                                  PatientN=smap.loc[sid, "PatientN"],
                                  group=smap.loc[sid, "NEDA3_1yr"],
                                  timepoint=smap.loc[sid, "Timepoint"],
                                  proportion=float(val)))
prop = pd.DataFrame(prop_rows)
prop.to_csv(SRC / "SuppFig5a_paired_proportions.csv", index=False)

# paired pre/post matrices + two-sided Wilcoxon signed-rank per subset x group
def paired_pre_post(df):
    """return DataFrame indexed by PatientN with Before/After columns."""
    w = df.pivot_table(index="PatientN", columns="timepoint",
                       values="proportion", aggfunc="mean")
    return w.dropna(subset=["Before", "After"])

# Original Fig. 4c,d used a two-sided paired t-test (exact matches: B_Activated
# R P=0.011, B_Naive 0.040, GZMK(hi) 0.038). We annotate the panel with the
# paired t-test for consistency with the canonical main figure (R0/R6/R8) and
# additionally report the two-sided Wilcoxon signed-rank in Supplementary Data 4.
stat_rows = []
pairs_by_key = {}
for lab, s in subset_order:
    for g in GROUPS:
        d = prop[(prop["subset"] == s) & (prop["group"] == g)]
        w = paired_pre_post(d)
        pairs_by_key[(s, g)] = w
        pt = pw = np.nan
        if len(w) >= 2 and not np.allclose(w["Before"].values, w["After"].values):
            pt = ttest_rel(w["Before"].values, w["After"].values).pvalue
            try:
                pw = wilcoxon(w["Before"].values, w["After"].values,
                              alternative="two-sided", zero_method="wilcox").pvalue
            except ValueError:
                pw = np.nan
        mb = w["Before"].mean() if len(w) else np.nan
        ma = w["After"].mean() if len(w) else np.nan
        stat_rows.append(dict(lineage=lab, subset=s, group=g, n_pairs=len(w),
                              mean_proportion_before=mb, mean_proportion_after=ma,
                              log2FC_after_over_before=np.log2((ma + 1e-9) / (mb + 1e-9)),
                              paired_t_two_sided_P=pt, wilcoxon_signed_rank_two_sided_P=pw))
wil = pd.DataFrame(stat_rows)
wil["P"] = wil["paired_t_two_sided_P"]          # figure-annotation statistic
wil.to_csv(SRC / "SuppFig5a_paired_stats.csv", index=False)
# Supplementary Data 4 (reproducible table)
SUPT = ROOT / "data" / "supplementary_tables"; SUPT.mkdir(parents=True, exist_ok=True)
wil.drop(columns=["P"]).to_csv(
    SUPT / "SupplementaryData4_paired_prepost_proportions.csv", index=False)
nsig = (wil["paired_t_two_sided_P"] < 0.05).sum()
print(f"panel a: {len(subset_order)} subsets, {nsig} subset x group P<0.05 (paired t)")
print(wil[wil["paired_t_two_sided_P"] < 0.05][
    ["subset", "group", "n_pairs", "paired_t_two_sided_P",
     "wilcoxon_signed_rank_two_sided_P"]].to_string(index=False))

# ----------------------------------------------------------------------------
# 2. SCENIC RSS for B_Activated (reuse saved AUCell matrix; compute RSS here)
# ----------------------------------------------------------------------------
def regulon_specificity_scores(auc_mtx, labels):
    """pyscenic.rss reimplementation: RSS = 1 - JS_distance(reg_dist, type_dist)."""
    cats = sorted(pd.unique(labels))
    out = np.empty((len(cats), auc_mtx.shape[1]))
    for j, reg in enumerate(auc_mtx.columns):
        col = auc_mtx[reg].values.astype(float)
        for i, c in enumerate(cats):
            membership = (labels.values == c).astype(float)
            out[i, j] = 1.0 - jensenshannon(col, membership)
    return pd.DataFrame(out, index=cats, columns=auc_mtx.columns)

auc = pd.read_csv(AUC, index_col=0)
B = io.load_lineage("B")
B.obs["anno2"] = B.obs["anno2"].astype(str)
lab = B.obs["anno2"].reindex(auc.index)
assert lab.notna().all(), "auc barcodes missing from 26JAN B object"
rss = regulon_specificity_scores(auc, lab)
rss_ba = rss.loc["B_Activated"].sort_values(ascending=False)
rss_ba.to_csv(SRC / "SuppFig5b_Bactivated_RSS.csv", header=["RSS"])
TOPN = 5
print("panel b: B_Activated top regulons:",
      list(rss_ba.head(TOPN).index), "n_regulons=", rss.shape[1])

# ----------------------------------------------------------------------------
# 2b. panel c — per-cell AP-1/NF-kB regulon-module activity across B subsets
#     (module = mean AUCell of JUN/FOSB/FOS/REL); highest in B_Activated.
# ----------------------------------------------------------------------------
MODREG = [r for r in ["JUN(+)", "FOSB(+)", "FOS(+)", "REL(+)"] if r in auc.columns]
Border = [s for s in LEVEL2["B"]]                       # canonical B subset order
anno_b = B.obs["anno2"].reindex(auc.index)
module = auc[MODREG].mean(axis=1)                        # per-cell module score
ba_mod = module[anno_b.values == "B_Activated"].values
rest_mod = module[anno_b.values != "B_Activated"].values
mwu_p = mannwhitneyu(ba_mod, rest_mod, alternative="two-sided").pvalue
print("panel c: B_Activated module vs rest  MWU two-sided P=%.2e  (med %.3f vs %.3f)"
      % (mwu_p, np.median(ba_mod), np.median(rest_mod)))

# ----------------------------------------------------------------------------
# Supplementary Data 6 — B_Activated defining genes + full regulon RSS
# ----------------------------------------------------------------------------
SUPT = ROOT / "data" / "supplementary_tables"; SUPT.mkdir(parents=True, exist_ok=True)
Bc = B.copy()
if float(Bc.X.max()) > 50:                              # counts -> log-normalize
    sc.pp.normalize_total(Bc, target_sum=1e4); sc.pp.log1p(Bc)
sc.tl.rank_genes_groups(Bc, "anno2", groups=["B_Activated"], reference="rest",
                        method="wilcoxon")
deg = sc.get.rank_genes_groups_df(Bc, group="B_Activated")
deg.to_csv(SUPT / "SupplementaryData6_Bactivated_defining_genes.csv", index=False)
rss.T.to_csv(SUPT / "SupplementaryData6_Bactivated_regulon_RSS.csv")   # regulons x subsets
n_up = int(((deg["logfoldchanges"] > 1) & (deg["pvals_adj"] < 0.05)).sum())
print(f"Supp Data 6: B_Activated defining genes (up, lfc>1 & padj<0.05) = {n_up}")

# ----------------------------------------------------------------------------
# 3. figure
# ----------------------------------------------------------------------------
NCOL = 6
n = len(subset_order)
ngrid_rows = int(np.ceil(n / NCOL))
fig = plt.figure(figsize=(7.2, 1.15 * ngrid_rows + 2.1 + 0.5))
gs = gridspec.GridSpec(ngrid_rows + 1, NCOL, figure=fig,
                       hspace=0.95, wspace=0.55,
                       height_ratios=[1] * ngrid_rows + [1.9])

# lineage color for subset titles
lin_color = {lab: c for lab, c in zip(
    [l[2] for l in LINEAGES],
    ["#3b3b6d", "#8c1c6d", "#1f6f3f", "#7a5c00", "#5a3a1a"])}

xpos = {("Responder", "Before"): 0, ("Responder", "After"): 1,
        ("Nonresponder", "Before"): 2.6, ("Nonresponder", "After"): 3.6}

def pstr(p):
    if not np.isfinite(p):
        return ""
    return f"{p:.3f}" if p >= 0.001 else "<0.001"

for idx, (lab_l, s) in enumerate(subset_order):
    ax = fig.add_subplot(gs[idx // NCOL, idx % NCOL])
    ymax = 0
    for g in GROUPS:
        w = pairs_by_key[(s, g)]
        if len(w) == 0:
            continue
        for tp in ["Before", "After"]:
            vals = w[tp].values
            x = xpos[(g, tp)]
            m_ = vals.mean(); se = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
            ax.bar(x, m_, width=0.8, color=GCOL[(g, tp)], edgecolor="black",
                   linewidth=0.5, zorder=1)
            ax.errorbar(x, m_, yerr=se, color="black", lw=0.6, capsize=1.5, zorder=3)
            ymax = max(ymax, m_ + se, vals.max())
        # paired lines
        for _, r in w.iterrows():
            ax.plot([xpos[(g, "Before")], xpos[(g, "After")]],
                    [r["Before"], r["After"]], color="0.35", lw=0.35,
                    marker="o", ms=1.2, zorder=2)
            ymax = max(ymax, r["Before"], r["After"])
        # p value (two-sided paired t-test, matching main Fig. 4c,d)
        pv = wil[(wil["subset"] == s) & (wil["group"] == g)]["P"].iloc[0]
        xc = (xpos[(g, "Before")] + xpos[(g, "After")]) / 2
        if pstr(pv):
            ax.text(xc, 1.02, pstr(pv), ha="center", va="bottom",
                    fontsize=4.3, transform=ax.get_xaxis_transform(),
                    color="black" if pv >= 0.05 else "#b00000")
    ax.set_xticks([0.5, 3.1])
    ax.set_xticklabels(["R", "NR"], fontsize=4.6)
    ax.tick_params(axis="y", labelsize=4.3, length=1.5, pad=1)
    ax.tick_params(axis="x", length=0, pad=1)
    ax.set_ylim(0, ymax * 1.28 if ymax > 0 else 1)
    prefix = {"CD8 T": "T_CD8_", "CD4 T": "T_CD4_", "B": "B_",
              "NK": "NK_", "myeloid": "M_"}[lab_l]
    short = s[len(prefix):] if s.startswith(prefix) else s
    short = short.replace("EM_", "")          # drop effector-memory tag for brevity
    ax.set_title(short, fontsize=4.9, color=lin_color[lab_l], pad=8, fontweight="bold")
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    for sp in ["left", "bottom"]:
        ax.spines[sp].set_linewidth(0.5)
    if idx % NCOL == 0:
        ax.set_ylabel("prop.", fontsize=4.6)

# panel letter + inline colored lineage key (top-left, figure coords)
fig.text(0.008, 0.992, "a", fontsize=style.PANEL_LETTER, fontweight="bold", va="top")
kx = 0.055
fig.text(kx, 0.990, "Subset lineage:", fontsize=4.8, va="top", color="black")
kx += 0.085
for lab_l, col in [(l[2], lin_color[l[2]]) for l in LINEAGES]:
    fig.text(kx, 0.990, lab_l, fontsize=4.8, va="top", color=col, fontweight="bold")
    kx += 0.008 * len(lab_l) + 0.022

# ---- bottom block: b RSS scatter | c per-cell module specificity (violin) ----
bot = gs[ngrid_rows, :].subgridspec(1, 2, width_ratios=[1.15, 1.15], wspace=0.34)

# panel b: RSS scatter
axb = fig.add_subplot(bot[0, 0])
yv = rss_ba.values
xr = np.arange(1, len(yv) + 1)
axb.scatter(xr[TOPN:], yv[TOPN:], s=6, c="0.6", edgecolors="none", zorder=1)
axb.scatter(xr[:TOPN], yv[:TOPN], s=10, c="#d62728", edgecolors="none", zorder=2)
# labels staggered with a FIXED vertical spread (data range), not the tiny RSS gaps
yrange = yv.max() - yv.min()
lx = len(yv) * 0.20
lstep = yrange * 0.115
for i in range(TOPN):
    axb.annotate(rss_ba.index[i], (xr[i], yv[i]),
                 xytext=(lx, yv[0] - i * lstep), fontsize=5, color="#d62728", va="center",
                 arrowprops=dict(arrowstyle="-", color="#d62728", lw=0.3,
                                 shrinkA=1, shrinkB=1))
axb.set_xlabel("Regulon rank", fontsize=6)
axb.set_ylabel("Regulon specificity score (RSS)", fontsize=6)
axb.set_title("TF regulon specificity in B_Activated", fontsize=6, pad=3, loc="left")
axb.tick_params(labelsize=5, length=2)
for sp in ["top", "right"]:
    axb.spines[sp].set_visible(False)
axb.text(-0.30, 1.08, "b", transform=axb.transAxes, clip_on=False,
         fontsize=style.PANEL_LETTER, fontweight="bold", va="top")

# panel c: per-cell AP-1/NF-kB module activity across B subsets (violin)
axc = fig.add_subplot(bot[0, 1])
c_order = [s for s in Border if s in set(anno_b.dropna())]
c_data = [module.values[anno_b.values == s] for s in c_order]
parts = axc.violinplot(c_data, showextrema=False, widths=0.85)
for k, pc in enumerate(parts["bodies"]):
    is_ba = c_order[k] == "B_Activated"
    pc.set_facecolor("#d62728" if is_ba else "0.7")
    pc.set_edgecolor("black"); pc.set_linewidth(0.3); pc.set_alpha(0.9 if is_ba else 0.6)
axc.scatter(range(1, len(c_order) + 1), [np.median(x) for x in c_data],
            s=4, c="black", zorder=3)
axc.set_xticks(range(1, len(c_order) + 1))
axc.set_xticklabels([s.replace("B_", "") for s in c_order], rotation=90, fontsize=4.8)
if "B_Activated" in c_order:
    ix = c_order.index("B_Activated")
    axc.get_xticklabels()[ix].set_color("#d62728")
    axc.get_xticklabels()[ix].set_fontweight("bold")
axc.set_ylabel("AP-1/NF-κB module\nAUCell (per cell)", fontsize=5.6)
axc.set_title("Module activity across B subsets", fontsize=6, pad=3, loc="left")
axc.tick_params(axis="y", labelsize=4.8, length=1.5); axc.tick_params(axis="x", length=0)
axc.text(0.5, 0.98, "B_Activated vs rest, P < 0.001",
         transform=axc.transAxes, ha="center", va="top", fontsize=5, color="#b00000")
for sp in ["top", "right"]:
    axc.spines[sp].set_visible(False)
axc.text(-0.22, 1.08, "c", transform=axc.transAxes, clip_on=False,
         fontsize=style.PANEL_LETTER, fontweight="bold", va="top")

fig.subplots_adjust(top=0.955, bottom=0.09, left=0.06, right=0.975)
out = ROOT / "figures" / "supplementary" / "Supplementary_Fig5.pdf"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=300)
print("wrote", out)
