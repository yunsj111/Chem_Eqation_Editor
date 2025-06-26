import streamlit as st
import streamlit.components.v1 as components
import json

# Canvas 패키지에서 모듈들 import
from ..canvas import (
    get_canvas_styles,
    get_canvas_template,
    get_element_data,
    get_canvas_core_js,
    get_atom_bond_manager,
    get_selection_drag_manager,
    get_rendering_engine,
    get_utilities_and_tests
)

def create_interactive_canvas():
    """JavaScript 기반 인터랙티브 분자 캔버스 - 리팩토링된 버전"""
    
    # HTML/CSS/JavaScript 코드 조합
    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        {get_canvas_styles()}
    </head>
    <body>
        {get_canvas_template()}
        
        <script>
            {get_element_data()}
            {get_canvas_core_js()}
            {get_atom_bond_manager()}
            {get_selection_drag_manager()}
            {get_rendering_engine()}
            {get_utilities_and_tests()}
        </script>
    </body>
    </html>
    """
    
    # Streamlit에서 HTML 컴포넌트 표시
    components.html(html_code, height=800, scrolling=False)
    
    return {
        'status': 'Canvas created successfully',
        'features': [
            '드래그앤 드롭 결합 생성',
            '다중 선택 및 이동',
            '주기율표 UI',
            'SMARTS 패턴 지원',
            '골격 구조식 스타일',
            '모듈화된 코드 구조'
        ]
    }

def get_canvas_info():
    """캔버스 정보 반환 (디버깅용)"""
    return {
        'modules': [
            'canvas_styles.py - CSS 스타일',
            'canvas_template.py - HTML 템플릿',
            'element_manager.py - 원소 및 주기율표 관리',
            'canvas_core.py - 핵심 JavaScript 기능',
            'atom_bond_manager.py - 원자와 결합 관리',
            'selection_drag_manager.py - 선택과 드래그 기능',
            'rendering_engine.py - 렌더링 엔진',
            'utilities_tests.py - 유틸리티 및 테스트'
        ],
        'total_lines_reduced': 'from ~1660 to ~200 per module'
    }
