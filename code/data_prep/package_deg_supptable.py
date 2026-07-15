"""Package baseline DEG lists as a Supplementary Data table.

Uses the ORIGINAL DEG tables (matches the published Fig 1e counts, e.g. CD8 = 979),
not the 26JAN re-computation. Significant DEG = |log2FoldChange| > 0.5 & pvalue < 0.05
(non-responder vs responder contrast as stored). One sheet per lineage + a summary.
출력: data/supplementary_tables/SupplementaryData3_DEG_baseline.xlsx (+ summary CSV)
번호: 원본 Supplementary Data 1(메타)·2(baseline) 다음이므로 DEG = Supplementary Data 3.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import os, glob
import pandas as pd

SRC = os.path.join(MS_SRC, "DATA/DEG/251110_0857_Total_DEG")      # READ-ONLY original
OUT = os.path.join(MS_ROOT, "data/supplementary_tables")
os.makedirs(OUT, exist_ok=True)
LFC, PCUT = 0.5, 0.05
ORDER = ["T_CD8", "B", "M_cDC", "M_Mono", "T_CD4", "NK", "T_CD4_Treg", "T_MAIT", "T_gd"]

xlsx = f"{OUT}/SupplementaryData3_DEG_baseline.xlsx"
summary = []
with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
    for ct in ORDER:
        f = f"{SRC}/Rpre_vs_NRpre_{ct}.csv"
        if not os.path.exists(f):
            print("missing:", ct); continue
        d = pd.read_csv(f, index_col=0)
        d.index.name = "gene"
        sig = d[(d["log2FoldChange"].abs() > LFC) & (d["pvalue"] < PCUT)].copy()
        sig = sig.sort_values("pvalue")
        # direction relative to responder: stored contrast is Rpre vs NRpre? label both
        up_R = (sig["log2FoldChange"] > 0).sum()
        up_NR = (sig["log2FoldChange"] < 0).sum()
        sig.to_excel(xw, sheet_name=ct[:31])
        summary.append((ct, len(d), len(sig), int(up_R), int(up_NR)))
sm = pd.DataFrame(summary, columns=["lineage", "n_tested", "n_sig", "up_in_R", "up_in_NR"])
sm.to_csv(f"{OUT}/SupplementaryData_DEG_summary.csv", index=False)
print(sm.to_string(index=False))
print("\nsaved:", xlsx)
print("(CD8 n_sig should match Fig 1e ~979)")
