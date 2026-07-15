"""IO + metadata for the MS reproduction. 26JAN annotated objects are the
fixed point (RULES R4). <MS_SRC> is READ-ONLY (R1)."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

from pathlib import Path
import pandas as pd
import scanpy as sc

ROOT = Path(str(MS_ROOT))
SRC = Path(str(MS_SRC))               # READ-ONLY source
D26 = SRC / "DATA" / "26JAN"
SUPP1 = ROOT / "originals" / "supplementary_data" / "260710_1232_MS_Supplementary_Data1.xlsx"

ANNO1 = D26 / "260122_anno1.h5ad"
LINEAGE = {
    "T_CD8":   D26 / "260122_anno2_T_CD8.h5ad",
    "T_CD4":   D26 / "260122_anno2_T_CD4.h5ad",
    "B":       D26 / "260122_anno2_B.h5ad",
    "NK":      D26 / "260122_anno2_NK.h5ad",
    "Myeloid": D26 / "260122_anno2_myeloid.h5ad",
}

META_COLS = ["PatientN", "Timepoint", "NEDA3_1yr", "Paired",
             "Group1", "DMT1", "DMTint"]


def load_anno1():
    return sc.read_h5ad(ANNO1)


def load_lineage(name):
    return sc.read_h5ad(LINEAGE[name])


def build_meta():
    """Canonical per-sample metadata from Supplementary Data 1.
    Group1 = NEDA3_1yr + '_' + Timepoint  (e.g. 'Responder_Before')."""
    m = pd.read_excel(SUPP1, "Sheet1")
    m["SampleID"] = m["SampleID"].astype(str)
    m["Group1"] = m["NEDA3_1yr"].astype(str) + "_" + m["Timepoint"].astype(str)
    return m


def join_meta(adata, meta=None, cols=META_COLS):
    """Attach per-sample clinical metadata onto adata.obs via SampleID."""
    if meta is None:
        meta = build_meta()
    smap = meta.set_index("SampleID")
    sid = adata.obs["SampleID"].astype(str)
    for c in cols:
        adata.obs[c] = sid.map(smap[c]).values
    return adata
