"""
Canvas 모듈 패키지 초기화
화학 구조 에디터의 캔버스 관련 기능들을 모아놓은 패키지
"""

from .canvas_styles import get_canvas_styles
from .canvas_template import get_canvas_template
from .element_manager import get_element_data
from .canvas_core import get_canvas_core_js
from .atom_bond_manager import get_atom_bond_manager
from .selection_drag_manager import get_selection_drag_manager
from .rendering_engine import get_rendering_engine
from .utilities_tests import get_utilities_and_tests

__all__ = [
    'get_canvas_styles',
    'get_canvas_template', 
    'get_element_data',
    'get_canvas_core_js',
    'get_atom_bond_manager',
    'get_selection_drag_manager',
    'get_rendering_engine',
    'get_utilities_and_tests'
]

__version__ = '1.0.0'
__author__ = 'Chemical Structure Editor Team'
__description__ = 'Canvas components for chemical structure drawing and editing'
