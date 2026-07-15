"""Export a clean per-sample metadata table + data dictionary for the deposit.

Pulls the canonical clinical metadata (Supplementary Data 1) and writes a tidy,
publication-ready sample sheet plus a column dictionary.
Output: data/metadata/sample_metadata.csv, data/metadata/data_dictionary.md
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
MS_ROOT = os.environ.get("MS_ROOT", os.path.dirname(os.path.dirname(_HERE)))
MS_SRC = os.environ.get("MS_SRC", os.path.join(MS_ROOT, "source_data"))

import sys, os
sys.path.insert(0, os.path.join(MS_ROOT, "code"))
from common import io

OUT = os.path.join(MS_ROOT, "data/metadata")
os.makedirs(OUT, exist_ok=True)

m = io.build_meta()
KEEP = ["SampleID", "PatientN", "Sex", "AgeSampling", "AgeOnset",
        "Timepoint", "Paired", "NEDA3_1yr", "NEDA3_2yr", "DMT1", "DMTint",
        "EDSS", "Cummulative_No_attack_before_sampling", "OCBpos", "IgG_index",
        "Group1"]
keep = [c for c in KEEP if c in m.columns]
tab = m[keep].copy().sort_values(["Timepoint", "NEDA3_1yr", "SampleID"])
tab.to_csv(f"{OUT}/sample_metadata.csv", index=False)

DICT = {
    "SampleID": "Unique PBMC sample ID",
    "PatientN": "Anonymized patient number (same value for paired samples)",
    "Sex": "Sex",
    "AgeSampling": "Age at sampling (years)",
    "AgeOnset": "Age at onset (years)",
    "Timepoint": "Sampling timepoint (Before = pre-treatment / After = post-treatment)",
    "Paired": "Whether before+after are paired for the same patient (Yes)",
    "NEDA3_1yr": "NEDA-3 response at 1 year (Responder / Nonresponder / nan = undetermined)",
    "NEDA3_2yr": "NEDA-3 response at 2 years",
    "DMT1": "Initial prescribed DMT (IFN/DMF/TRF/RTX/GA)",
    "DMTint": "DMT group for analysis (intention-to-treat: DMF/IFN/TRF)",
    "EDSS": "EDSS score at sampling",
    "Cummulative_No_attack_before_sampling": "Cumulative number of relapses before sampling",
    "OCBpos": "Cerebrospinal fluid oligoclonal band positivity",
    "IgG_index": "CSF/serum IgG index",
    "Group1": "NEDA3_1yr + '_' + Timepoint (analysis group label)",
}
lines = ["# Sample metadata — data dictionary\n",
         f"`sample_metadata.csv`  ({len(tab)} samples · {tab['PatientN'].nunique()} patients)\n",
         "| Column | Description |", "|---|---|"]
for c in keep:
    lines.append(f"| `{c}` | {DICT.get(c, '')} |")
lines += ["", "## Cohort summary",
          f"- Samples {len(tab)} · Patients {tab['PatientN'].nunique()}",
          f"- Timepoint: " + ", ".join(f"{k} {v}" for k, v in tab['Timepoint'].value_counts().items()),
          f"- NEDA3_1yr: " + ", ".join(f"{k} {v}" for k, v in tab['NEDA3_1yr'].value_counts(dropna=False).items()),
          f"- DMTint: " + ", ".join(f"{k} {v}" for k, v in tab['DMTint'].value_counts().items())]
with open(f"{OUT}/data_dictionary.md", "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("saved:", f"{OUT}/sample_metadata.csv", "+ data_dictionary.md")
print(tab.head(6).to_string(index=False))
