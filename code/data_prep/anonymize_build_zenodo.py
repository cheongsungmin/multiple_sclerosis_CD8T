#!/usr/bin/env python
"""Build <MS_ROOT>/zenodo/ : anonymized, upload-ready deposit.
PatientN (real reg numbers) -> P01..P28 in h5ad objects + sample_metadata.csv.
Private real<->anon map kept OUTSIDE zenodo/ (in ../private/).
Run in scanpy env.
"""
import os, sys, shutil, csv
import numpy as np, pandas as pd, h5py

ROOT = os.path.expanduser("<MS_ROOT>")
SRC_OBJ = os.path.join(ROOT, "data/objects")
SRC_META = os.path.join(ROOT, "data/metadata")
ZEN = os.path.join(ROOT, "zenodo")
PRIV = os.path.join(ROOT, "private")

OBJECTS = [
    "MS_scRNA_anno1_total.h5ad",
    "MS_scRNA_anno2_T_CD8.h5ad",
    "MS_scRNA_anno2_T_CD4.h5ad",
    "MS_scRNA_anno2_B.h5ad",
    "MS_scRNA_anno2_NK.h5ad",
    "MS_scRNA_anno2_myeloid.h5ad",
]

def log(*a): print(*a, flush=True)

# ---- 0. build deterministic map from sample_metadata reading order ----
meta = pd.read_csv(os.path.join(SRC_META, "sample_metadata.csv"), dtype=str)
seen, order = set(), []
for v in meta["PatientN"].astype(str):
    if v not in seen:
        seen.add(v); order.append(v)
assert len(order) == 28, f"expected 28 patients, got {len(order)}"
MAP = {reg: f"P{ i+1:02d}" for i, reg in enumerate(order)}   # real -> Pxx
CATS = [f"P{i+1:02d}" for i in range(28)]                     # canonical category order
log(f"[map] {len(MAP)} patients -> P01..P28")

# ---- write PRIVATE mapping (outside zenodo) ----
os.makedirs(PRIV, exist_ok=True)
with open(os.path.join(PRIV, "patient_id_map.PRIVATE.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["real_PatientN", "anon_id"])
    for reg, pid in MAP.items(): w.writerow([reg, pid])
with open(os.path.join(PRIV, "README.txt"), "w") as f:
    f.write("DO NOT UPLOAD. Real registration number <-> anonymized ID mapping (PHI).\n"
            "Keep private. Never place in zenodo/ or any public deposit.\n")
log(f"[private] wrote {PRIV}/patient_id_map.PRIVATE.csv")

# ---- 1. zenodo dirs ----
os.makedirs(os.path.join(ZEN, "objects"), exist_ok=True)
os.makedirs(os.path.join(ZEN, "metadata"), exist_ok=True)

# ---- 2. anonymize sample_metadata.csv ----
meta_anon = meta.copy()
meta_anon["PatientN"] = meta_anon["PatientN"].astype(str).map(MAP)
assert meta_anon["PatientN"].isna().sum() == 0
meta_anon.to_csv(os.path.join(ZEN, "metadata/sample_metadata.csv"), index=False)
log("[meta] wrote anonymized sample_metadata.csv")
# data_dictionary copied as-is (already says 'anonymized ID'); ensured below in README step.
shutil.copy2(os.path.join(SRC_META, "data_dictionary.md"),
             os.path.join(ZEN, "metadata/data_dictionary.md"))

# ---- 3. anonymize h5ad objects (cp then edit obs/PatientN only) ----
str_dt = h5py.special_dtype(vlen=str)
for name in OBJECTS:
    src = os.path.join(SRC_OBJ, name)
    dst = os.path.join(ZEN, "objects", name)
    log(f"[obj] copy {name} ({os.path.getsize(src)/1e9:.2f} GB) ...")
    shutil.copy2(src, dst)
    with h5py.File(dst, "r+") as h:
        node = h["obs/PatientN"]
        vals = node[...].astype(np.int64)
        keys = [str(int(v)) for v in vals]
        missing = set(keys) - set(MAP)
        assert not missing, f"{name}: reg numbers not in map: {missing}"
        codes = np.array([CATS.index(MAP[k]) for k in keys], dtype=np.int8)
        del h["obs/PatientN"]
        g = h.create_group("obs/PatientN")
        g.attrs["encoding-type"] = "categorical"
        g.attrs["encoding-version"] = "0.2.0"
        g.attrs["ordered"] = np.bool_(False)
        cds = g.create_dataset("categories", data=np.array(CATS, dtype=object), dtype=str_dt)
        cds.attrs["encoding-type"] = "string-array"; cds.attrs["encoding-version"] = "0.2.0"
        cod = g.create_dataset("codes", data=codes, dtype=np.int8)
        cod.attrs["encoding-type"] = "array"; cod.attrs["encoding-version"] = "0.2.0"
    # verify via anndata
    import anndata as ad
    a = ad.read_h5ad(dst, backed="r")
    uniq = set(a.obs["PatientN"].astype(str).unique())
    bad = uniq - set(CATS)
    assert not bad, f"{name}: unexpected PatientN values {bad}"
    assert not (uniq & set(MAP.keys())), f"{name}: real reg numbers still present!"
    log(f"[obj] {name}: OK n_obs={a.n_obs} patients={len(uniq)} -> {sorted(uniq)[:3]}...")

log("BUILD DONE")
