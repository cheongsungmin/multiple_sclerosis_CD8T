# code/ — 분석/그림 스크립트

## 규칙
- 실행: 기본 `scanpy` env의 절대경로 python. 단계별 env는 `../config/environments.md`.
- **환경 수정 금지** (RULES R2). naming은 `../config/naming.yaml`에서 로드 (하드코딩 금지, R3).
- 원본 `<MS_SRC>`·`../originals/`는 READ-ONLY 참조. 산출물은 `<MS_ROOT>` 안에만.
- stochastic 단계는 seed 고정 (결정론적, R4). 그림 스타일은 `common/style.py` 전역.

## 레이아웃 (2026-07-13 작전변경 후)
```
code/
├── common/        # 공용 유틸 (모든 그림이 재사용)
│   ├── io.py         # 객체 로더, 경로 상수(Supplementary Data는 originals/ 참조)
│   ├── naming.py     # naming.yaml 로더 + alias→canonical 매핑
│   ├── style.py      # 전역 스타일: 글씨 검정, 기본6/축5/패널문자 볼드12, pdf.fonttype=42
│   └── palette.py    # 색 정본(정비된 data/objects의 uns 색)
├── data_prep/     # 객체 정비·metadata·Supplementary Data 패키징
└── supple/        # supplementary figure 스크립트 (main 1개 ↔ supple 1개)
     supp_fig1.py         # Supplementary Fig. 1 (supports Fig. 1) — 확정
     s9_tcr_antigen.py    # Supp Fig. 2 구성요소(TCR 항원)
     tableS1_baseline.py  # Supplementary Data 2 재현 소스 CSV
```
- 메인 피겨(Fig 1–6)는 원본 사용(재현 X) → `../originals/published_figures/`.
- 옛 재현 파이프라인은 `../archive/fig1_overview/`로 이동(검증 기록).

## 산출 규약
- supplementary figure → `../figures/supplementary/<name>.pdf` (벡터, 일러 편집, dpi=300).
- 표/통계 → `../data/results/`, Supplementary Data → `../data/supplementary_tables/`.
- 패널 순서 = 본문 등장 순서. legend는 `../manuscript/legends.md`.
