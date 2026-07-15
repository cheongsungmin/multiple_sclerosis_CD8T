# multiple_sclerosis_CD8T — reproduction code

Single-cell multi-omics of PBMCs in multiple sclerosis.

Analysis and figure-generation code for the manuscript
*[TITLE — to be filled]* ([AUTHORS — to be filled]).

The analysis-ready single-cell objects (`.h5ad`) and per-sample metadata are
deposited separately on **Zenodo: [DOI — to be filled]**. This repository holds
only the code that turns those objects into the manuscript's figures and tables.

## Layout

```
config/
└── naming.yaml            canonical subset names (loaded at runtime; do not hardcode)
code/
├── common/                shared utilities reused by every script
│   ├── io.py                 object loaders + path constants
│   ├── naming.py             naming.yaml loader (alias → canonical)
│   ├── style.py              global figure style
│   └── palette.py            canonical colors
├── data_prep/             build clean objects, export metadata, package deposits
├── main/                  main-figure scripts (fig3_flow.py)
└── supple/                supplementary figure / table scripts
```

See [code/README.md](code/README.md) for per-script details.

## Requirements

- Python ≥ 3.10 with `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`.
- Stochastic steps use a fixed seed; embeddings in the deposited objects are
  precomputed and must **not** be recomputed (results are deterministic).

## Data & paths

All paths are resolved from two environment variables (no hardcoded paths):

| Variable | Meaning | Default |
|---|---|---|
| `MS_ROOT` | project root — outputs (`figures/`, `data/`), `config/`, `originals/` | repository root |
| `MS_SRC`  | source objects downloaded from Zenodo + raw inputs | `$MS_ROOT/source_data` |

1. Download the analysis-ready objects from Zenodo ([DOI — to be filled]).
2. Point the variables at your local copies, then run any script:

```bash
export MS_SRC=/your/path/to/zenodo-objects
export MS_ROOT="$PWD"        # or wherever outputs should be written
python code/supple/supp_fig1.py
```

If unset, `MS_ROOT` defaults to the repository root and `MS_SRC` to
`$MS_ROOT/source_data`, so placing the downloaded objects under `source_data/`
works with no configuration.

## License

Code is released under the license in [LICENSE](LICENSE).
Deposited data are released separately under CC-BY-4.0 (see the Zenodo record).

## Citation

[Preferred citation — to be filled]
