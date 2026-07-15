"""Export a clean per-sample metadata table + data dictionary for the deposit.

Pulls the canonical clinical metadata (Supplementary Data 1) and writes a tidy,
publication-ready sample sheet plus a column dictionary.
출력: data/metadata/sample_metadata.csv, data/metadata/data_dictionary.md
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
    "SampleID": "PBMC 샘플 고유 ID",
    "PatientN": "환자 익명 번호 (paired 샘플은 동일 값)",
    "Sex": "성별",
    "AgeSampling": "채혈 시 나이(년)",
    "AgeOnset": "발병 나이(년)",
    "Timepoint": "채혈 시점 (Before=치료 전 / After=치료 후)",
    "Paired": "동일 환자 before+after 짝지음 여부 (Yes)",
    "NEDA3_1yr": "1년 시점 NEDA-3 반응 (Responder / Nonresponder / nan=미확정)",
    "NEDA3_2yr": "2년 시점 NEDA-3 반응",
    "DMT1": "초기 처방 DMT (IFN/DMF/TRF/RTX/GA)",
    "DMTint": "분석용 DMT 그룹 (intention-to-treat: DMF/IFN/TRF)",
    "EDSS": "채혈 시 EDSS 점수",
    "Cummulative_No_attack_before_sampling": "채혈 전 누적 재발 횟수",
    "OCBpos": "뇌척수액 oligoclonal band 양성 여부",
    "IgG_index": "CSF/serum IgG index",
    "Group1": "NEDA3_1yr + '_' + Timepoint (분석 그룹 라벨)",
}
lines = ["# Sample metadata — 데이터 사전\n",
         f"`sample_metadata.csv`  ({len(tab)} samples · {tab['PatientN'].nunique()} patients)\n",
         "| 컬럼 | 설명 |", "|---|---|"]
for c in keep:
    lines.append(f"| `{c}` | {DICT.get(c, '')} |")
lines += ["", "## 코호트 요약",
          f"- 샘플 {len(tab)} · 환자 {tab['PatientN'].nunique()}",
          f"- Timepoint: " + ", ".join(f"{k} {v}" for k, v in tab['Timepoint'].value_counts().items()),
          f"- NEDA3_1yr: " + ", ".join(f"{k} {v}" for k, v in tab['NEDA3_1yr'].value_counts(dropna=False).items()),
          f"- DMTint: " + ", ".join(f"{k} {v}" for k, v in tab['DMTint'].value_counts().items())]
with open(f"{OUT}/data_dictionary.md", "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("saved:", f"{OUT}/sample_metadata.csv", "+ data_dictionary.md")
print(tab.head(6).to_string(index=False))
