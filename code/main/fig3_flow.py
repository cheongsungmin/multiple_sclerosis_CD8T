"""Main Figure 3 (flow cytometry) — built in-house (user-directed; not the original PDF).

Temra-forward, honest gating: the headline is the CD45RA-defined Temra population
(CD8+ CCR7- CD45RA+), a marker phenotype (not a transcriptomic label). Gating strategy
and the finer CD45RA/CD29/KIR populations live in Supplementary Fig. 4.

Panels:
  a  study design (13 responders, 13 non-responders; flow cytometry)
  b  representative CCR7 vs CD45RA plots (one NR, one R) with the Temra gate  [PENDING FCS/FlowKit]
  c  CD8+ CCR7- CD45RA+ (Temra) proportion, responder vs non-responder (two-sided t-test)
  d  NEDA-3 failure Kaplan-Meier by Temra high/low (log-rank)

Source (READ-ONLY): <MS_SRC>/DATA/24NOV/FACS/250824_FACS_ODS.csv
Panel b needs event-level FCS (no reader in permitted envs -> placeholder until a flowkit env exists).
Output: figures/main/Figure3.pdf  (vector, pdf.fonttype=42, dpi=300)
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import sys, os, warnings
sys.path.insert(0, os.path.join(MS_ROOT, "code"))
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.stats import ttest_ind
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from common import style
style.apply()
BASE, TICK, PLET = style.BASE, style.TICK, style.PANEL_LETTER

SRC = os.path.join(MS_SRC, "DATA/24NOV/FACS/250824_FACS_ODS.csv")
OUT = os.path.join(MS_ROOT, "figures/main")
os.makedirs(OUT, exist_ok=True)
R_COL, NR_COL = "#3B6FB6", "#C0392B"

ff = pd.read_csv(SRC).rename(columns={"Age(Sampling)": "Age", "Temra/CD8": "Temra",
                                      "time_to_fail_NEDA3_or_last_fu_month": "time"})
ff["event"] = ff["NEDA3_fail_date"].notna().astype(int)
cc = ff[ff["Responder_1yr"].isin([0, 1])].copy()
cc["group"] = np.where(cc["Responder_1yr"] == 1, "Responder", "Non-responder")

fig = plt.figure(figsize=(13.5, 3.6))
gs = fig.add_gridspec(1, 4, width_ratios=[1.05, 1.15, 0.8, 1.25], wspace=0.5)


def plet(ax, letter, dx=-0.18, dy=1.12):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=PLET, fontweight="bold",
            va="bottom", ha="right")


# ---- a: study design -------------------------------------------------------------
axa = fig.add_subplot(gs[0]); axa.axis("off"); plet(axa, "a", dx=-0.02)
axa.set_xlim(0, 1); axa.set_ylim(0, 1)
for y, lab, col, n in [(0.72, "Responders", R_COL, 13), (0.42, "Non-responders", NR_COL, 13)]:
    box = FancyBboxPatch((0.02, y - 0.10), 0.42, 0.20, boxstyle="round,pad=0.01",
                         linewidth=0.8, edgecolor=col, facecolor=col, alpha=0.15)
    axa.add_patch(box)
    axa.text(0.23, y, f"{lab}\n(n = {n})", ha="center", va="center", fontsize=TICK, color=col)
axa.add_patch(FancyArrowPatch((0.46, 0.57), (0.60, 0.57), arrowstyle="-|>", mutation_scale=9,
                              lw=0.8, color="black"))
axa.text(0.80, 0.72, "Pre-DMT\nPBMCs", ha="center", va="center", fontsize=TICK)
axa.add_patch(FancyArrowPatch((0.80, 0.63), (0.80, 0.52), arrowstyle="-|>", mutation_scale=9,
                              lw=0.8, color="black"))
axa.text(0.80, 0.42, "Flow\ncytometry", ha="center", va="center", fontsize=TICK, fontweight="bold")
axa.text(0.5, 0.10, "13 vs 13; largely overlapping the single-cell cohort",
         ha="center", va="center", fontsize=TICK - 1, color="#555", transform=axa.transAxes)

# ---- b: representative plots (placeholder pending FCS/FlowKit) --------------------
axb = fig.add_subplot(gs[1]); axb.axis("off"); plet(axb, "b", dx=-0.05)
for i, (lab, col) in enumerate([("Non-responder", NR_COL), ("Responder", R_COL)]):
    x0 = 0.02 + i * 0.50
    box = FancyBboxPatch((x0, 0.18), 0.46, 0.64, boxstyle="round,pad=0.005",
                         linewidth=0.8, edgecolor="#999", facecolor="#f5f5f5", linestyle="--")
    axb.add_patch(box)
    axb.text(x0 + 0.23, 0.74, lab, ha="center", fontsize=TICK, color=col, fontweight="bold")
    axb.text(x0 + 0.23, 0.50, "CCR7 vs CD45RA\n(gated on CD8+)\nTemra gate", ha="center",
             va="center", fontsize=TICK - 1, color="#666")
axb.text(0.5, 0.06, "[representative FCS plots — pending FlowKit env]", ha="center",
         fontsize=TICK - 1, style="italic", color="#888")

# ---- c: Temra proportion R vs NR -------------------------------------------------
axc = fig.add_subplot(gs[2]); plet(axc, "c", dx=-0.42)
order = ["Non-responder", "Responder"]; pal = {"Non-responder": NR_COL, "Responder": R_COL}
sns.boxplot(data=cc, x="group", y="Temra", order=order, ax=axc, width=0.6, showfliers=False,
            palette=pal, boxprops=dict(alpha=.35), medianprops=dict(color="black", lw=1))
sns.stripplot(data=cc, x="group", y="Temra", order=order, ax=axc, size=3,
              palette=pal, edgecolor="black", linewidth=.3, jitter=.18)
p = ttest_ind(cc[cc.group == "Responder"].Temra, cc[cc.group == "Non-responder"].Temra, equal_var=True).pvalue
y = cc.Temra.max() * 1.08
axc.plot([0, 1], [y, y], color="black", lw=.6)
axc.text(0.5, y * 1.01, f"P = {p:.3f}", ha="center", va="bottom", fontsize=TICK)
axc.set_xlabel(""); axc.set_ylabel("CD8+ CCR7− CD45RA+ (Temra)\nfraction of CD8+", fontsize=TICK)
axc.set_title("Temra", fontsize=BASE, fontweight="bold")
axc.set_xticklabels(["NR", "R"], fontsize=TICK); axc.margins(y=.18); sns.despine(ax=axc)

# ---- d: NEDA-3 failure KM by Temra high/low --------------------------------------
axd = fig.add_subplot(gs[3]); plet(axd, "d", dx=-0.16)
s = ff.dropna(subset=["Temra_high_mean", "time"])
hi = s[s.Temra_high_mean == 1]; lo = s[s.Temra_high_mean == 0]
lr = logrank_test(hi.time, lo.time, hi.event, lo.event)
for d, col, lab in [(hi, R_COL, f"Temra high (n = {len(hi)})"), (lo, NR_COL, f"Temra low (n = {len(lo)})")]:
    km = KaplanMeierFitter().fit(d.time, d.event)
    ci = 1.0 - km.survival_function_.values[:, 0]
    t = km.survival_function_.index.values
    axd.step(t, ci, where="post", color=col, lw=1.3, label=lab)
    cens = d[d.event == 0].time.values
    if len(cens):
        cy = np.interp(cens, t, ci)
        axd.plot(cens, cy, "|", color=col, ms=5, mew=0.8)
axd.set_xlim(0, min(60, s.time.max()))
axd.set_ylim(-0.02, 1.02)
axd.set_xlabel("Time from DMT initiation (months)", fontsize=TICK)
axd.set_ylabel("Cumulative probability\nof NEDA-3 failure", fontsize=TICK)
axd.text(0.04, 0.94, f"log-rank P = {lr.p_value:.3f}", transform=axd.transAxes, fontsize=TICK, va="top")
axd.legend(loc="lower right", fontsize=TICK - 0.5, frameon=False)
sns.despine(ax=axd)

fig.savefig(f"{OUT}/Figure3.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Temra R vs NR t P = {p:.4f} | KM log-rank P = {lr.p_value:.4f} (high {len(hi)}/low {len(lo)})")
print("wrote:", f"{OUT}/Figure3.pdf")
