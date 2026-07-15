"""Publication plotting style (matplotlib/seaborn). Import and call apply().

사용자 원본 세팅 반영: 모든 글씨 검은색, pdf.fonttype=42(일러 편집).
크기 체계(A4 기준): 기본 6 / 축(tick) 5 / 패널 알파벳 볼드 12.
패널 문자 크기는 PANEL_LETTER 로 노출.
"""
import matplotlib as mpl
import seaborn as sns

BASE = 6            # 기본 글씨(제목·축라벨·범례·주석)
TICK = 5            # 축에 있는 눈금 글씨
PANEL_LETTER = 12   # 패널 알파벳 제목 (볼드)


def apply():
    mpl.rcParams["pdf.fonttype"] = 42       # editable text in Illustrator
    mpl.rcParams["ps.fonttype"] = 42
    # --- 모든 글씨/선 검은색 (사용자 원본 세팅) ---
    for k in ["text.color", "axes.edgecolor", "axes.labelcolor",
              "axes.titlecolor", "xtick.color", "ytick.color", "lines.color"]:
        mpl.rcParams[k] = "black"
    mpl.rcParams["patch.edgecolor"] = "black"
    # --- 폰트 크기 통일 ---
    mpl.rcParams["font.size"] = BASE
    mpl.rcParams["axes.titlesize"] = BASE
    mpl.rcParams["axes.labelsize"] = BASE
    mpl.rcParams["legend.fontsize"] = BASE
    mpl.rcParams["figure.titlesize"] = BASE
    mpl.rcParams["xtick.labelsize"] = TICK
    mpl.rcParams["ytick.labelsize"] = TICK
    # --- 배경 흰색 ---
    mpl.rcParams["figure.facecolor"] = "white"
    mpl.rcParams["axes.facecolor"] = "white"
    mpl.rcParams["savefig.facecolor"] = "white"
    sns.set_style("white")
    sns.set_style("ticks")
