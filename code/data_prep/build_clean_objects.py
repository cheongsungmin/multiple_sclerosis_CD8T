"""Build self-contained, analysis-ready RNA objects for the data-availability deposit.

For each 26JAN object (anno1 total + 5 anno2 lineages), produce a clean copy in
data/objects/ that needs NO manual relabeling/recoloring downstream:
  - X            : log-normalized expression (as-is, 26JAN fixed point)
  - layers['counts'] : raw UMI counts (injected from the count object, same 36,601 genes)
  - obs          : + clinical metadata (Group1, NEDA3_1yr, Timepoint, Paired, DMT...) via SampleID
  - anno* cats   : canonical order (from the naming module)
  - uns[*_colors]: curated palette (colors_*.pkl), 6-digit hex, aligned to category order
Originals in <MS_SRC> are READ-ONLY (R1); we only write into <MS_ROOT>/data (R5).
Output: data/objects/MS_scRNA_{anno1_total, anno2_<lineage>}.h5ad
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import sys, os, gc
sys.path.insert(0, os.path.join(MS_ROOT, "code"))
import warnings; warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc

from common import io, naming, palette

OUTDIR = os.path.join(MS_ROOT, "data/objects")
os.makedirs(OUTDIR, exist_ok=True)
COUNTOBJ = os.path.join(MS_SRC, "DATA/63_241025_MS_Control_count.h5ad")

# object spec: (short, source h5ad, anno cols, color pkl name for anno2)
D26 = io.D26
SPEC = [
    ("anno1_total",   D26 / "260122_anno1.h5ad",         ["anno1"],          "anno1"),
    ("anno2_T_CD8",   D26 / "260122_anno2_T_CD8.h5ad",   ["anno1", "anno2"], "anno2_T_CD8"),
    ("anno2_T_CD4",   D26 / "260122_anno2_T_CD4.h5ad",   ["anno1", "anno2"], "anno2_T_CD4"),
    ("anno2_B",       D26 / "260122_anno2_B.h5ad",       ["anno1", "anno2"], "anno2_B"),
    ("anno2_NK",      D26 / "260122_anno2_NK.h5ad",      ["anno1", "anno2"], "anno2_NK"),
    ("anno2_myeloid", D26 / "260122_anno2_myeloid.h5ad", ["anno1", "anno2"], "anno2_myeloid"),
]
META_ADD = ["PatientN", "Sex", "AgeSampling", "Timepoint", "NEDA3_1yr",
            "Paired", "Group1", "DMT1", "DMTint"]


def strip_alpha(hexc):
    """#rrggbbaa -> #rrggbb ; passthrough otherwise."""
    if isinstance(hexc, str) and hexc.startswith("#") and len(hexc) == 9:
        return hexc[:7]
    return hexc


def ordered_cats(present, key):
    """Canonical order for a lineage's subsets; append any extras (e.g. T_CD4_CTL in CD8)."""
    canon = naming.LEVEL1 if key == "anno1" else naming.LEVEL2.get(key, [])
    head = [c for c in canon if c in present]
    tail = [c for c in present if c not in head]
    return head + tail


# ---- load raw counts once (shared) ----
print("loading count object ...")
cnt = sc.read_h5ad(COUNTOBJ)
cnt_counts = cnt.layers["count"]
cnt_index = {b: i for i, b in enumerate(cnt.obs_names)}
cnt_var = list(cnt.var_names)
del cnt; gc.collect()
meta = io.build_meta()

for short, src, annocols, pkl2 in SPEC:
    print(f"\n=== {short} ===")
    ad = sc.read_h5ad(src)

    # 1) inject raw counts (align rows by barcode; identical 36,601-gene var order)
    assert list(ad.var_names) == cnt_var, f"{short}: var order mismatch"
    rows = np.array([cnt_index[b] for b in ad.obs_names])
    ad.layers["counts"] = cnt_counts[rows, :].copy()

    # 2) join clinical metadata via SampleID
    io.join_meta(ad, meta=meta, cols=[c for c in META_ADD if c in meta.columns])

    # 3) canonical category order + curated colors for each anno col
    lineage_key = short.replace("anno2_", "").replace("myeloid", "Myeloid")  # T_CD8/T_CD4/B/NK/Myeloid
    for col in annocols:
        s = ad.obs[col].astype("category")
        order_key = "anno1" if col == "anno1" else lineage_key
        cats = ordered_cats(list(s.cat.categories), order_key)
        ad.obs[col] = s.cat.reorder_categories(cats)
        cols = palette.load_colors("anno1" if col == "anno1" else pkl2)
        ad.uns[f"{col}_colors"] = [strip_alpha(cols.get(c, "#999999")) for c in cats]

    # 4) tidy: drop intermediate clustering artifacts (keep final leiden + anno + std scanpy)
    import re
    drop_obs = [c for c in ad.obs.columns if re.match(r"^leiden_R", c)]
    ad.obs.drop(columns=drop_obs, inplace=True)
    for k in list(ad.uns.keys()):
        if (re.match(r"^leiden_R", k) or k.startswith("cdm_")
                or re.match(r"^dendrogram_leiden_R", k)
                or k in {"pct_counts_mt_colors"}):
            del ad.uns[k]

    out = f"{OUTDIR}/MS_scRNA_{short}.h5ad"
    ad.write_h5ad(out)
    sz = os.path.getsize(out) / 1e9
    print(f"  n={ad.n_obs} vars={ad.n_vars} layers={list(ad.layers)} "
          f"metacols+={[c for c in META_ADD if c in ad.obs][:4]}... -> {out} ({sz:.2f}GB)")
    print(f"  {annocols[-1]} order: {list(ad.obs[annocols[-1]].cat.categories)}")
    del ad; gc.collect()

print("\nDONE. clean objects in", OUTDIR)
