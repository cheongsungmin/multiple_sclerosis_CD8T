"""Canonical subset naming — single source of truth = config/naming.yaml.
Do NOT hardcode subset names elsewhere; import from here (RULES R3)."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

from pathlib import Path
import yaml

_CFG = Path(os.path.join(MS_ROOT, "config/naming.yaml"))
_d = yaml.safe_load(open(_CFG))

LEVEL1 = list(_d["level1"])            # 12 lineages, canonical order
LEVEL2 = {k: list(v) for k, v in _d["level2"].items()}
ALIASES = dict(_d.get("aliases", {}))

# full subset order (for total dotplots/barplots): B, Myeloid, CD8, CD4, NK
ORDER_TOTAL = (LEVEL2["B"] + LEVEL2["Myeloid"]
               + [c for c in LEVEL2["T_CD8"] if c.startswith("T_CD8")]
               + LEVEL2["T_CD4"] + LEVEL2["NK"])


def canon(x):
    """Map a deviant/legacy label to its canonical form."""
    return ALIASES.get(x, x)


def canon_series(s):
    return s.map(lambda v: ALIASES.get(v, v))
