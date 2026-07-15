"""Canonical subset names and ordering for anno1 (12 lineages) and anno2
(per-lineage subsets). The deposited objects already store these categories in
canonical order; these lists reproduce that order for figures/tables."""

# anno1: 12 lineages, canonical order
LEVEL1 = [
    "T_CD4", "T_CD4_Treg", "T_CD8", "T_IFN", "T_MAIT", "T_Prolif",
    "T_gd", "NK", "B", "M_Mono", "M_cDC", "M_pDC",
]

# anno2: subsets per lineage, canonical order
LEVEL2 = {
    "T_CD8": [
        "T_CD8_Naive/CM", "T_CD8_EM_GZMK(lo)", "T_CD8_EM_GZMK(hi)",
        "T_CD8_EM_ZNF683(hi)", "T_CD8_EM_KIR", "T_CD8_IFN",
    ],
    "T_CD4": [
        "T_CD4_SOX4", "T_CD4_Naive/CM", "T_CD4_CXCR5", "T_CD4_Th1",
        "T_CD4_Th2", "T_CD4_Th17", "T_CD4_Treg", "T_CD4_CTL", "T_CD4_IFN",
    ],
    "B": [
        "B_Naive_SOX4", "B_Naive", "B_Activated", "B_Memory_Unswitched",
        "B_Memory_Switched", "B_Memory_Atypical", "B_Plasma",
        "B_Plasmablast", "B_IFN",
    ],
    "NK": ["NK_CD56(hi)", "NK_CD56(lo)", "NK_IFN", "NK_Prolif"],
    "Myeloid": [
        "M_Mono_CD14", "M_Mono_CD16", "M_Mono_C1QA", "M_Mono_IFN",
        "M_cDC1", "M_cDC2", "M_pDC",
    ],
}

# full subset order for total plots: B, Myeloid, CD8, CD4, NK
ORDER_TOTAL = (LEVEL2["B"] + LEVEL2["Myeloid"]
               + [c for c in LEVEL2["T_CD8"] if c.startswith("T_CD8")]
               + LEVEL2["T_CD4"] + LEVEL2["NK"])
