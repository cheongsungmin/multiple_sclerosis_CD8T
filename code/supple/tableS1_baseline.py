"""Reproduction source for Supplementary Data 2: baseline clinical/demographic characteristics (R vs NR).

The original Supplementary Data 2 (originals/supplementary_data/...Data2.xlsx) is canonical;
this script only builds the source CSV that reproduces/validates those values in code (not a
separate figure). At the pre-DMT sample level (one baseline per patient). Numeric: median (range)
+ MWU; categorical: n (%) + Fisher.
Output: data/supplementary_tables/TableS1_baseline_characteristics.csv
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import sys, os
sys.path.insert(0, os.path.join(MS_ROOT, "code"))
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import mannwhitneyu, fisher_exact

from common import io
mpl.rcParams["pdf.fonttype"] = 42
TAB = os.path.join(MS_ROOT, "data/supplementary_tables")
FIG = os.path.join(MS_ROOT, "figures/supplementary")
os.makedirs(TAB, exist_ok=True)

m = io.build_meta()
m["dur"] = m["AgeSampling"] - m["AgeOnset"]
base = m[(m.Timepoint == "Before") & (m.NEDA3_1yr.isin(["Responder", "Nonresponder"]))].copy()
R = base[base.NEDA3_1yr == "Responder"]; NR = base[base.NEDA3_1yr == "Nonresponder"]

def num(col):
    x = pd.to_numeric(R[col], errors="coerce").dropna()
    y = pd.to_numeric(NR[col], errors="coerce").dropna()
    p = mannwhitneyu(x, y).pvalue if len(x) > 1 and len(y) > 1 else np.nan
    fmt = lambda v: f"{v.median():.1f} ({v.min():.1f}–{v.max():.1f})"
    return fmt(x), fmt(y), f"{p:.2f}"

def cat_female():
    def femp(d):
        f = (d.Sex.astype(str).str.upper().str[0] == "F").sum()
        return f, len(d), f"{f} ({100*f/len(d):.0f}%)"
    fr, nr_, sr = femp(R); fn, nn_, sn = femp(NR)
    p = fisher_exact([[fr, nr_ - fr], [fn, nn_ - fn]])[1]
    return sr, sn, f"{p:.2f}"

def cat_ocb():
    def pos(d):
        s = d["OCBpos"].astype(str)
        p_ = s.isin(["1", "1.0", "(+)", "+", "Y", "Yes"]).sum()
        return f"{p_} ({100*p_/len(d):.0f}%)"
    return pos(R), pos(NR), "-"

rows = [("n (baseline samples/patients)", str(len(R)), str(len(NR)), "-"),
        ("Female, n (%)", *cat_female()),
        ("Age at sampling, median (range)", *num("AgeSampling")),
        ("Age at onset, median (range)", *num("AgeOnset")),
        ("Disease duration, yr, median (range)", *num("dur")),
        ("EDSS, median (range)", *num("EDSS")),
        ("Prior attacks, median (range)", *num("Cummulative_No_attack_before_sampling")),
        ("OCB positive, n (%)", *cat_ocb()),
        ("IgG index, median (range)", *num("IgG_index"))]
# DMT distribution rows
for dmt in ["DMF", "IFN", "TRF"]:
    rows.append((f"DMT: {dmt}, n", str((R.DMTint == dmt).sum()), str((NR.DMTint == dmt).sum()), "-"))

tab = pd.DataFrame(rows, columns=["Characteristic", "Responder", "Non-responder", "P"])
tab.to_csv(f"{TAB}/TableS1_baseline_characteristics.csv", index=False)
print(tab.to_string(index=False))

print(f"\nBaseline: {len(R)} responders vs {len(NR)} non-responders "
      "(28 patients, 41 samples; 23 pre / 18 post; 24 paired from 12).")
print("saved: TableS1_baseline_characteristics.csv  (= Supplementary Data 2 reproduction source)")
