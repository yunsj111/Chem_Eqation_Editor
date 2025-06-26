import math
from typing import Tuple, List, Dict, Any
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io
import base64

def get_atom_color(element: str) -> str:
    """원소별 색상 반환"""
    colors = {
        'C': '#000000',  # 검정
        'N': '#0000FF',  # 파랑
        'O': '#FF0000',  # 빨강
        'S': '#FFFF00',  # 노랑
        'P': '#FFA500',  # 주황
        'F': '#00FF00',  # 초록
        'Cl': '#00FF00', # 초록
        'Br': '#8B4513', # 갈색
        'I': '#800080',  # 보라
        'H': '#FFFFFF',  # 흰색
        '*': '#808080'   # 회색 (wildcard)
    }
    return colors.get(element, '#808080')

def calculate_distance(x1: float, y1: float, x2: float, y2: float) -> float:
    """두 점 사이의 거리 계산"""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def point_on_line(px: float, py: float, x1: float, y1: float, 
                  x2: float, y2: float, threshold: float = 5) -> bool:
    """점이 선 위에 있는지 확인"""
    # 점과 선 사이의 최단 거리 계산
    A = px - x1
    B = py - y1
    C = x2 - x1
    D = y2 - y1

    dot = A * C + B * D
    len_sq = C * C + D * D

    if len_sq == 0:
        return False

    param = dot / len_sq

    if param < 0 or param > 1:
        return False

    xx = x1 + param * C
    yy = y1 + param * D

    dx = px - xx
    dy = py - yy

    return math.sqrt(dx * dx + dy * dy) <= threshold

def mol_to_image(mol, size=(300, 300)) -> Image.Image:
    """RDKit 분자를 PIL 이미지로 변환"""
    if mol is None:
        # 빈 이미지 생성
        img = Image.new('RGB', size, 'white')
        return img

    try:
        img = Draw.MolToImage(mol, size=size)
        return img
    except:
        # 오류 시 빈 이미지 반환
        img = Image.new('RGB', size, 'white')
        return img

def image_to_base64(img: Image.Image) -> str:
    """PIL 이미지를 base64 문자열로 변환"""
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode()
    return img_str

def validate_smarts(smarts: str) -> Tuple[bool, str]:
    """SMARTS 패턴 유효성 검사"""
    try:
        mol = Chem.MolFromSmarts(smarts)
        if mol is None:
            return False, "유효하지 않은 SMARTS 패턴입니다."
        return True, "유효한 SMARTS 패턴입니다."
    except Exception as e:
        return False, f"SMARTS 패턴 오류: {str(e)}"

def normalize_coordinates(atoms: List[Dict], canvas_width: int = 600, 
                         canvas_height: int = 400) -> List[Dict]:
    """좌표를 캔버스 크기에 맞게 정규화"""
    if not atoms:
        return atoms

    # 현재 좌표의 범위 계산
    min_x = min(atom['x'] for atom in atoms)
    max_x = max(atom['x'] for atom in atoms)
    min_y = min(atom['y'] for atom in atoms)
    max_y = max(atom['y'] for atom in atoms)

    # 범위가 0인 경우 처리
    range_x = max_x - min_x if max_x != min_x else 1
    range_y = max_y - min_y if max_y != min_y else 1

    # 여백을 고려한 스케일링
    margin = 50
    scale_x = (canvas_width - 2 * margin) / range_x
    scale_y = (canvas_height - 2 * margin) / range_y
    scale = min(scale_x, scale_y)

    # 정규화된 좌표로 변환
    normalized_atoms = []
    for atom in atoms:
        new_x = (atom['x'] - min_x) * scale + margin
        new_y = (atom['y'] - min_y) * scale + margin

        normalized_atom = atom.copy()
        normalized_atom['x'] = new_x
        normalized_atom['y'] = new_y
        normalized_atoms.append(normalized_atom)

    return normalized_atoms

@st.cache_data
def get_common_smarts_patterns() -> Dict[str, str]:
    """자주 사용되는 SMARTS 패턴들"""
    return {
        "방향족 고리": "c1ccccc1",
        "카르보닐": "C=O",
        "히드록실": "OH",
        "아미노": "N",
        "카르복실": "C(=O)O",
        "에스테르": "C(=O)O[!H]",
        "아미드": "C(=O)N",
        "니트로": "N(=O)=O",
        "술폰산": "S(=O)(=O)O",
        "인산": "P(=O)(O)(O)O",
        "알데히드": "C=O[H]",
        "케톤": "C(=O)[!H]",
        "에테르": "O[!H][!H]",
        "티올": "SH",
        "할로겐": "[F,Cl,Br,I]",
        "임의 원자": "*",
        "비수소": "[!H]",
        "방향족 탄소": "c",
        "지방족 탄소": "C",
        "양전하": "[+]",
        "음전하": "[-]"
    }

def get_element_info() -> Dict[str, Dict[str, Any]]:
    """원소 정보 반환"""
    return {
        'C': {'name': '탄소', 'color': '#000000', 'common': True},
        'N': {'name': '질소', 'color': '#0000FF', 'common': True},
        'O': {'name': '산소', 'color': '#FF0000', 'common': True},
        'S': {'name': '황', 'color': '#FFFF00', 'common': True},
        'P': {'name': '인', 'color': '#FFA500', 'common': True},
        'F': {'name': '플루오린', 'color': '#00FF00', 'common': True},
        'Cl': {'name': '염소', 'color': '#00FF00', 'common': True},
        'Br': {'name': '브롬', 'color': '#8B4513', 'common': True},
        'I': {'name': '요오드', 'color': '#800080', 'common': True},
        'H': {'name': '수소', 'color': '#FFFFFF', 'common': False},
        '*': {'name': '임의원자', 'color': '#808080', 'common': True}
    }