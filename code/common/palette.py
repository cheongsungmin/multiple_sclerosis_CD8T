"""Color palettes. Reuse original curated colors_*.pkl to match the paper.
Files live in <MS_SRC>/01_Figures_pre (READ-ONLY)."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

from pathlib import Path
import pickle

_CDIR = Path(os.path.join(MS_SRC, "01_Figures_pre"))

# fixed group colors used throughout (Set2-derived, from original notebooks)
GROUP = {
    "Responder_Before": "#66c2a5",
    "Nonresponder_Before": "#fc8d62",
    "Responder_After": "#b3e2cd",
    "Nonresponder_After": "#fdcdac",
    "Responder": "#66c2a5",
    "Nonresponder": "#fc8d62",
}


def load_colors(name):
    """name in {anno1, anno2_total, anno2_B, anno2_myeloid, anno2_NK,
    anno2_T_CD4, anno2_T_CD8}. Returns whatever was pickled (usually dict)."""
    p = _CDIR / f"colors_{name}.pkl"
    if not p.exists():
        return {}
    with open(p, "rb") as fh:
        return pickle.load(fh)


def palette_for(categories, name):
    """Ordered color list aligned to `categories`; None if incomplete."""
    cols = load_colors(name)
    if isinstance(cols, dict) and all(c in cols for c in categories):
        return [cols[c] for c in categories]
    return None
