"""
Drawing 모듈 패키지 초기화
화학 구조 그리기 관련 기능들을 모아놓은 패키지
"""

from .molecule_drawer import MoleculeDrawer
from .js_molecule_canvas import create_interactive_canvas, get_canvas_info
from .utils import (
    get_atom_color, 
    calculate_distance, 
    point_on_line,
    get_common_smarts_patterns,
    get_element_info
)

__all__ = [
    'MoleculeDrawer',
    'create_interactive_canvas',
    'get_canvas_info',
    'get_atom_color',
    'calculate_distance', 
    'point_on_line',
    'get_common_smarts_patterns',
    'get_element_info'
]

__version__ = '1.0.0'
__author__ = 'Chemical Structure Editor Team'
__description__ = 'Drawing components for chemical structure editor'
