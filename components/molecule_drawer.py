import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd
from typing import List, Dict, Optional, Tuple, Any
import math
import json

from models.molecule import Molecule
from models.atom import Atom, AtomType, AtomQuery
from models.bond import Bond, BondType, BondQuery
from components.utils import get_atom_color, calculate_distance, point_on_line

class MoleculeDrawer:
    """분자 그리기 컴포넌트"""

    def __init__(self):
        self.molecule = Molecule()
        self.canvas_width = 800
        self.canvas_height = 600
        self.selected_atom = None
        self.temp_bond_start = None

    def initialize_session_state(self):
        """세션 상태 초기화"""
        if 'molecule_data' not in st.session_state:
            st.session_state.molecule_data = {
                'atoms': [],
                'bonds': [],
                'atom_counter': 0,
                'bond_counter': 0
            }

        if 'drawing_mode' not in st.session_state:
            st.session_state.drawing_mode = 'atom'

        if 'current_element' not in st.session_state:
            st.session_state.current_element = 'C'

        if 'current_bond_type' not in st.session_state:
            st.session_state.current_bond_type = BondType.SINGLE

        if 'selected_atom_id' not in st.session_state:
            st.session_state.selected_atom_id = None

        if 'temp_bond_start_id' not in st.session_state:
            st.session_state.temp_bond_start_id = None

    def create_drawing_interface(self):
        """그리기 인터페이스 생성"""
        self.initialize_session_state()

        # 상단 도구 모음 (ChemDraw 스타일)
        st.markdown("### 🧪 분자 구조 편집기")
        
        # 도구 모음을 컬럼으로 배치
        col1, col2, col3, col4, col5, col6 = st.columns([2, 2, 2, 2, 2, 2])
        
        with col1:
            # 편집 모드 (아이콘 스타일)
            if st.button("🔵 원자", help="원자 추가 모드"):
                st.session_state.drawing_mode = 'atom'
            if st.button("🔗 결합", help="결합 생성 모드"):
                st.session_state.drawing_mode = 'bond'
        
        with col2:
            if st.button("👆 선택", help="선택/이동 모드"):
                st.session_state.drawing_mode = 'select'
            if st.button("🗑️ 삭제", help="삭제 모드"):
                st.session_state.drawing_mode = 'delete'
        
        with col3:
            # 자주 사용하는 원소들
            st.write("**원소:**")
            element_col1, element_col2 = st.columns(2)
            with element_col1:
                if st.button("C", key="quick_C", use_container_width=True):
                    st.session_state.current_element = 'C'
                    st.session_state.drawing_mode = 'atom'
                if st.button("O", key="quick_O", use_container_width=True):
                    st.session_state.current_element = 'O'
                    st.session_state.drawing_mode = 'atom'
            with element_col2:
                if st.button("N", key="quick_N", use_container_width=True):
                    st.session_state.current_element = 'N'
                    st.session_state.drawing_mode = 'atom'
                if st.button("S", key="quick_S", use_container_width=True):
                    st.session_state.current_element = 'S'
                    st.session_state.drawing_mode = 'atom'
        
        with col4:
            # 결합 타입들
            st.write("**결합:**")
            if st.button("—", key="single_bond", help="단일 결합", use_container_width=True):
                st.session_state.current_bond_type = BondType.SINGLE
                st.session_state.drawing_mode = 'bond'
            if st.button("=", key="double_bond", help="이중 결합", use_container_width=True):
                st.session_state.current_bond_type = BondType.DOUBLE
                st.session_state.drawing_mode = 'bond'
        
        with col5:
            # 기본 구조들
            st.write("**구조:**")
            if st.button("⬟", key="benzene", help="벤젠", use_container_width=True):
                self.add_predefined_structure("c1ccccc1")
            if st.button("⬡", key="cyclohexane", help="사이클로헥산", use_container_width=True):
                self.add_predefined_structure("C1CCCCC1")
        
        with col6:
            # 액션 버튼들
            st.write("**액션:**")
            if st.button("🗑️ 전체삭제", use_container_width=True):
                self.clear_molecule()
            if st.button("↶ 실행취소", use_container_width=True):
                self.undo_last_action()

        # 현재 상태 표시 바
        st.markdown("---")
        status_col1, status_col2, status_col3, status_col4 = st.columns(4)
        
        with status_col1:
            mode_names = {
                'atom': '🔵 원자 추가',
                'bond': '🔗 결합 생성', 
                'select': '👆 선택/이동',
                'delete': '🗑️ 삭제'
            }
            st.write(f"**모드:** {mode_names.get(st.session_state.drawing_mode, st.session_state.drawing_mode)}")
        
        with status_col2:
            st.write(f"**원소:** {st.session_state.current_element}")
        
        with status_col3:
            bond_names = {
                BondType.SINGLE: "단일 (—)",
                BondType.DOUBLE: "이중 (=)",
                BondType.TRIPLE: "삼중 (≡)",
                BondType.AROMATIC: "방향족 (:)",
                BondType.ANY: "임의 (~)"
            }
            st.write(f"**결합:** {bond_names.get(st.session_state.current_bond_type, '단일')}")
        
        with status_col4:
            atom_count = len(st.session_state.molecule_data['atoms'])
            bond_count = len(st.session_state.molecule_data['bonds'])
            st.write(f"**분자:** 원자 {atom_count}개, 결합 {bond_count}개")

        # 메인 캔버스 (더 큰 크기로)
        canvas_result = self.create_canvas()
        
        # 하단 고급 옵션 (접을 수 있게)
        with st.expander("🔧 고급 옵션", expanded=False):
            adv_col1, adv_col2, adv_col3 = st.columns(3)
            
            with adv_col1:
                st.subheader("원소 선택")
                # 더 많은 원소들
                element_grid = st.columns(4)
                elements = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'H', '*', 'X']
                for i, elem in enumerate(elements):
                    with element_grid[i % 4]:
                        if st.button(elem, key=f"adv_elem_{elem}", use_container_width=True):
                            st.session_state.current_element = elem
                            st.session_state.drawing_mode = 'atom'
            
            with adv_col2:
                st.subheader("SMARTS 옵션")
                
                # 방향족/지방족
                atom_aromaticity = st.selectbox(
                    "방향족성:",
                    ["자동", "방향족", "지방족"],
                    key="atom_aromaticity"
                )
                
                # 전하
                charge = st.selectbox(
                    "전하:",
                    [0, +1, +2, -1, -2],
                    key="atom_charge"
                )
                
                # 결합 타입
                bond_type = st.selectbox(
                    "결합 타입:",
                    [BondType.SINGLE, BondType.DOUBLE, BondType.TRIPLE, 
                    BondType.AROMATIC, BondType.ANY],
                    format_func=lambda x: {
                        BondType.SINGLE: "단일 결합 (—)",
                        BondType.DOUBLE: "이중 결합 (=)",
                        BondType.TRIPLE: "삼중 결합 (≡)",
                        BondType.AROMATIC: "방향족 결합 (:)",
                        BondType.ANY: "임의 결합 (~)"
                    }[x],
                    key='current_bond_type'
                )
            
            with adv_col3:
                st.subheader("미리 정의된 구조")
                
                predefined_structures = {
                    "벤젠": "c1ccccc1",
                    "사이클로헥산": "C1CCCCC1",
                    "나프탈렌": "c1ccc2ccccc2c1",
                    "피리딘": "c1ccncc1",
                    "푸란": "c1ccoc1",
                    "이미다졸": "c1c[nH]cn1",
                    "페닐": "c1ccccc1",
                    "메틸": "C"
                }
                
                for name, smiles in predefined_structures.items():
                    if st.button(f"+ {name}", key=f"struct_{name}", use_container_width=True):
                        self.add_predefined_structure(smiles)

        return canvas_result

    def create_canvas(self):
        """그리기 캔버스 생성 (개선된 버전)"""
        # 현재 분자 데이터 로드
        self.load_molecule_from_session()

        # 더 큰 캔버스 크기
        self.canvas_width = 1000
        self.canvas_height = 700

        # Plotly 그래프 생성
        fig = go.Figure()

        # 격자 배경 (ChemDraw 스타일)
        for i in range(0, self.canvas_width, 50):
            fig.add_shape(
                type="line",
                x0=i, y0=0, x1=i, y1=self.canvas_height,
                line=dict(color="lightgray", width=0.5, dash="dot")
            )
        
        for i in range(0, self.canvas_height, 50):
            fig.add_shape(
                type="line",
                x0=0, y0=i, x1=self.canvas_width, y1=i,
                line=dict(color="lightgray", width=0.5, dash="dot")
            )

        # 캔버스 경계
        fig.add_shape(
            type="rect",
            x0=0, y0=0, x1=self.canvas_width, y1=self.canvas_height,
            fillcolor="white",
            line=dict(color="gray", width=2)
        )

        # 결합 그리기 (원자보다 먼저)
        self.draw_bonds(fig)

        # 원자 그리기
        self.draw_atoms(fig)

        # 임시 결합 표시 (결합 모드에서 첫 번째 원자 선택됨)
        if (st.session_state.drawing_mode == 'bond' and 
            st.session_state.temp_bond_start_id is not None):
            start_atom = self.get_atom_by_id(st.session_state.temp_bond_start_id)
            if start_atom:
                # 선택된 원자를 하이라이트
                fig.add_trace(go.Scatter(
                    x=[start_atom['x']],
                    y=[start_atom['y']],
                    mode='markers',
                    marker=dict(
                        size=30,
                        color='rgba(255,0,0,0.3)',
                        line=dict(color='red', width=3)
                    ),
                    showlegend=False,
                    hoverinfo='skip',
                    name='selected'
                ))

        # 디버깅: 현재 분자 데이터 확인
        if st.session_state.molecule_data['atoms']:
            st.write(f"**디버그**: {len(st.session_state.molecule_data['atoms'])}개의 원자가 있습니다.")
            for atom in st.session_state.molecule_data['atoms'][:3]:  # 처음 3개만 표시
                st.write(f"원자 {atom['id']}: {atom['element']} at ({atom['x']:.1f}, {atom['y']:.1f})")

        # 캔버스 설정
        fig.update_layout(
            width=self.canvas_width,
            height=self.canvas_height,
            xaxis=dict(
                range=[0, self.canvas_width], 
                showgrid=False, 
                showticklabels=False,
                zeroline=False,
                scaleanchor="y",
                scaleratio=1
            ),
            yaxis=dict(
                range=[0, self.canvas_height], 
                showgrid=False, 
                showticklabels=False,
                zeroline=False,
                autorange='reversed'  # Y축 뒤집기
            ),
            showlegend=False,
            margin=dict(l=0, r=0, t=0, b=0),
            plot_bgcolor='white',
            paper_bgcolor='white'
        )

        # 인터렉션 설정
        config = {
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['pan2d', 'lasso2d', 'autoScale2d'],
            'staticPlot': False
        }

        # Streamlit에서 Plotly 차트 표시
        st.plotly_chart(
            fig, 
            use_container_width=True, 
            config=config,
            key="molecule_canvas"
        )

        # 빠른 편집 패널
        st.markdown("### ⚡ 빠른 편집")
        
        quick_col1, quick_col2, quick_col3 = st.columns(3)
        
        with quick_col1:
            st.write("**위치 지정 추가**")
            x_pos = st.number_input("X", 0, self.canvas_width, 500, key="quick_x")
            y_pos = st.number_input("Y", 0, self.canvas_height, 350, key="quick_y")
            
            if st.button("📍 현재 모드로 추가", type="primary"):
                self.handle_canvas_click({'x': x_pos, 'y': y_pos})
        
        with quick_col2:
            st.write("**자동 배치**")
            if st.button("🎯 중앙에 원자"):
                self.handle_canvas_click({'x': 500, 'y': 350})
            if st.button("📐 격자에 맞춤"):
                self.snap_to_grid()
            if st.button("🔄 구조 정렬"):
                self.auto_arrange()
        
        with quick_col3:
            st.write("**파일 작업**")
            
            # 저장
            if st.button("💾 JSON 저장"):
                molecule_json = json.dumps(st.session_state.molecule_data, indent=2)
                st.download_button(
                    label="📥 다운로드",
                    data=molecule_json,
                    file_name="molecule_structure.json",
                    mime="application/json",
                    key="download_btn"
                )
            
            # 불러오기
            uploaded_file = st.file_uploader(
                "📁 구조 불러오기", 
                type=['json'], 
                key="quick_upload"
            )
            if uploaded_file:
                try:
                    molecule_data = json.load(uploaded_file)
                    st.session_state.molecule_data = molecule_data
                    st.success("구조를 불러왔습니다!")
                    st.rerun()
                except:
                    st.error("파일 형식이 올바르지 않습니다.")

        # 상세 정보 (축소 가능)
        with st.expander("📊 분자 정보", expanded=False):
            info_col1, info_col2 = st.columns(2)
            
            with info_col1:
                if st.session_state.molecule_data['atoms']:
                    st.write("**원자 목록:**")
                    atom_df = pd.DataFrame([
                        {
                            'ID': atom['id'],
                            '원소': atom['element'],
                            'X': f"{atom['x']:.0f}",
                            'Y': f"{atom['y']:.0f}"
                        }
                        for atom in st.session_state.molecule_data['atoms']
                    ])
                    st.dataframe(atom_df, use_container_width=True, height=200)
            
            with info_col2:
                if st.session_state.molecule_data['bonds']:
                    st.write("**결합 목록:**")
                    bond_df = pd.DataFrame([
                        {
                            'ID': bond['id'],
                            '원자1': bond['atom1_id'],
                            '원자2': bond['atom2_id'],
                            '타입': {1: "단일", 2: "이중", 3: "삼중", 4: "방향족", 5: "임의"}.get(bond['bond_type'], "알 수 없음")
                        }
                        for bond in st.session_state.molecule_data['bonds']
                    ])
                    st.dataframe(bond_df, use_container_width=True, height=200)

        return fig

    def draw_atoms(self, fig):
        """원자 그리기"""
        atoms = st.session_state.molecule_data['atoms']
        
        if not atoms:
            # 원자가 없으면 빈 trace 추가 (차트가 제대로 렌더링되도록)
            fig.add_trace(go.Scatter(
                x=[],
                y=[],
                mode='markers',
                showlegend=False
            ))
            return
        
        for atom_data in atoms:
            try:
                color = get_atom_color(atom_data['element'])

                # 원자 원 그리기
                fig.add_trace(go.Scatter(
                    x=[atom_data['x']],
                    y=[atom_data['y']],
                    mode='markers+text',
                    marker=dict(
                        size=25,
                        color=color,
                        line=dict(color='black', width=2)
                    ),
                    text=atom_data.get('display_text', atom_data['element']),
                    textfont=dict(
                        color='white' if color != '#FFFFFF' else 'black', 
                        size=14,
                        family="Arial"
                    ),
                    textposition='middle center',
                    hovertemplate=f"원자: {atom_data['element']}<br>ID: {atom_data['id']}<br>위치: ({atom_data['x']:.1f}, {atom_data['y']:.1f})<extra></extra>",
                    showlegend=False,
                    name=f"atom_{atom_data['id']}"
                ))
            except Exception as e:
                st.error(f"원자 그리기 오류 (ID: {atom_data['id']}): {e}")

    def draw_bonds(self, fig):
        """결합 그리기"""
        bonds = st.session_state.molecule_data['bonds']
        
        for bond_data in bonds:
            try:
                atom1 = self.get_atom_by_id(bond_data['atom1_id'])
                atom2 = self.get_atom_by_id(bond_data['atom2_id'])

                if atom1 and atom2:
                    bond_type = BondType(bond_data['bond_type'])

                    if bond_type == BondType.SINGLE:
                        self.draw_single_bond(fig, atom1, atom2)
                    elif bond_type == BondType.DOUBLE:
                        self.draw_double_bond(fig, atom1, atom2)
                    elif bond_type == BondType.TRIPLE:
                        self.draw_triple_bond(fig, atom1, atom2)
                    elif bond_type == BondType.AROMATIC:
                        self.draw_aromatic_bond(fig, atom1, atom2)
                    elif bond_type == BondType.ANY:
                        self.draw_any_bond(fig, atom1, atom2)
            except Exception as e:
                st.error(f"결합 그리기 오류 (ID: {bond_data['id']}): {e}")

    def snap_to_grid(self):
        """원자를 격자에 맞춤"""
        for atom in st.session_state.molecule_data['atoms']:
            atom['x'] = round(atom['x'] / 50) * 50
            atom['y'] = round(atom['y'] / 50) * 50
        st.success("격자에 맞춤 완료!")
        st.rerun()

    def auto_arrange(self):
        """구조 자동 정렬"""
        # 간단한 자동 정렬 (원자들을 중앙으로 이동)
        if st.session_state.molecule_data['atoms']:
            atoms = st.session_state.molecule_data['atoms']
            
            # 현재 중심점 계산
            avg_x = sum(atom['x'] for atom in atoms) / len(atoms)
            avg_y = sum(atom['y'] for atom in atoms) / len(atoms)
            
            # 캔버스 중심으로 이동
            center_x, center_y = self.canvas_width // 2, self.canvas_height // 2
            offset_x = center_x - avg_x
            offset_y = center_y - avg_y
            
            for atom in atoms:
                atom['x'] += offset_x
                atom['y'] += offset_y
            
            st.success("구조가 정렬되었습니다!")
            st.rerun()

    # def draw_atoms(self, fig):
    #     """원자 그리기"""
    #     for atom_data in st.session_state.molecule_data['atoms']:
    #         color = get_atom_color(atom_data['element'])

    #         # 원자 원 그리기
    #         fig.add_trace(go.Scatter(
    #             x=[atom_data['x']],
    #             y=[atom_data['y']],
    #             mode='markers+text',
    #             marker=dict(
    #                 size=20,
    #                 color=color,
    #                 line=dict(color='black', width=2)
    #             ),
    #             text=atom_data.get('display_text', atom_data['element']),
    #             textfont=dict(color='white' if color != '#FFFFFF' else 'black', size=12),
    #             textposition='middle center',
    #             hovertemplate=f"원자: {atom_data['element']}<br>ID: {atom_data['id']}<extra></extra>",
    #             showlegend=False
    #         ))

    # def draw_bonds(self, fig):
    #     """결합 그리기"""
    #     for bond_data in st.session_state.molecule_data['bonds']:
    #         atom1 = self.get_atom_by_id(bond_data['atom1_id'])
    #         atom2 = self.get_atom_by_id(bond_data['atom2_id'])

    #         if atom1 and atom2:
    #             bond_type = BondType(bond_data['bond_type'])

    #             if bond_type == BondType.SINGLE:
    #                 self.draw_single_bond(fig, atom1, atom2)
    #             elif bond_type == BondType.DOUBLE:
    #                 self.draw_double_bond(fig, atom1, atom2)
    #             elif bond_type == BondType.TRIPLE:
    #                 self.draw_triple_bond(fig, atom1, atom2)
    #             elif bond_type == BondType.AROMATIC:
    #                 self.draw_aromatic_bond(fig, atom1, atom2)
    #             elif bond_type == BondType.ANY:
    #                 self.draw_any_bond(fig, atom1, atom2)

    def draw_single_bond(self, fig, atom1, atom2):
        """단일 결합 그리기"""
        fig.add_trace(go.Scatter(
            x=[atom1['x'], atom2['x']],
            y=[atom1['y'], atom2['y']],
            mode='lines',
            line=dict(color='black', width=2),
            showlegend=False,
            hoverinfo='skip'
        ))

    def draw_double_bond(self, fig, atom1, atom2):
        """이중 결합 그리기"""
        # 두 평행선으로 이중 결합 표현
        dx = atom2['x'] - atom1['x']
        dy = atom2['y'] - atom1['y']
        length = math.sqrt(dx*dx + dy*dy)

        if length > 0:
            # 수직 벡터 계산
            ux, uy = -dy/length, dx/length
            offset = 3

            # 첫 번째 선
            fig.add_trace(go.Scatter(
                x=[atom1['x'] + ux*offset, atom2['x'] + ux*offset],
                y=[atom1['y'] + uy*offset, atom2['y'] + uy*offset],
                mode='lines',
                line=dict(color='black', width=2),
                showlegend=False,
                hoverinfo='skip'
            ))

            # 두 번째 선
            fig.add_trace(go.Scatter(
                x=[atom1['x'] - ux*offset, atom2['x'] - ux*offset],
                y=[atom1['y'] - uy*offset, atom2['y'] - uy*offset],
                mode='lines',
                line=dict(color='black', width=2),
                showlegend=False,
                hoverinfo='skip'
            ))

    def draw_triple_bond(self, fig, atom1, atom2):
        """삼중 결합 그리기"""
        # 중앙선
        fig.add_trace(go.Scatter(
            x=[atom1['x'], atom2['x']],
            y=[atom1['y'], atom2['y']],
            mode='lines',
            line=dict(color='black', width=2),
            showlegend=False,
            hoverinfo='skip'
        ))

        # 두 평행선
        dx = atom2['x'] - atom1['x']
        dy = atom2['y'] - atom1['y']
        length = math.sqrt(dx*dx + dy*dy)

        if length > 0:
            ux, uy = -dy/length, dx/length
            offset = 4

            # 위쪽 선
            fig.add_trace(go.Scatter(
                x=[atom1['x'] + ux*offset, atom2['x'] + ux*offset],
                y=[atom1['y'] + uy*offset, atom2['y'] + uy*offset],
                mode='lines',
                line=dict(color='black', width=2),
                showlegend=False,
                hoverinfo='skip'
            ))

            # 아래쪽 선
            fig.add_trace(go.Scatter(
                x=[atom1['x'] - ux*offset, atom2['x'] - ux*offset],
                y=[atom1['y'] - uy*offset, atom2['y'] - uy*offset],
                mode='lines',
                line=dict(color='black', width=2),
                showlegend=False,
                hoverinfo='skip'
            ))

    def draw_aromatic_bond(self, fig, atom1, atom2):
        """방향족 결합 그리기 (점선)"""
        fig.add_trace(go.Scatter(
            x=[atom1['x'], atom2['x']],
            y=[atom1['y'], atom2['y']],
            mode='lines',
            line=dict(color='black', width=2, dash='dot'),
            showlegend=False,
            hoverinfo='skip'
        ))

    def draw_any_bond(self, fig, atom1, atom2):
        """임의 결합 그리기 (물결선 효과)"""
        fig.add_trace(go.Scatter(
            x=[atom1['x'], atom2['x']],
            y=[atom1['y'], atom2['y']],
            mode='lines',
            line=dict(color='gray', width=3, dash='dashdot'),
            showlegend=False,
            hoverinfo='skip'
        ))

    def handle_canvas_click(self, point_data):
        """캔버스 클릭 이벤트 처리"""
        if 'x' not in point_data or 'y' not in point_data:
            return

        x, y = point_data['x'], point_data['y']
        mode = st.session_state.drawing_mode

        if mode == 'atom':
            self.add_atom_at(x, y)
        elif mode == 'bond':
            self.handle_bond_creation(x, y)
        elif mode == 'select':
            self.select_atom_at(x, y)
        elif mode == 'delete':
            self.delete_at(x, y)

    def add_atom_at(self, x: float, y: float):
        """지정된 위치에 원자 추가"""
        element = st.session_state.current_element

        # SMARTS 옵션 적용
        atom_data = {
            'id': st.session_state.molecule_data['atom_counter'],
            'x': x,
            'y': y,
            'element': element,
            'atom_type': 'wildcard' if element == '*' else 'element'
        }

        # SMARTS 쿼리 조건 추가
        if 'atom_aromaticity' in st.session_state:
            if st.session_state.atom_aromaticity == '방향족':
                atom_data['aromatic'] = True
            elif st.session_state.atom_aromaticity == '지방족':
                atom_data['aromatic'] = False

        if 'atom_charge' in st.session_state and st.session_state.atom_charge != 0:
            atom_data['charge'] = st.session_state.atom_charge

        if 'hydrogen_count' in st.session_state and st.session_state.hydrogen_count != '자동':
            atom_data['hydrogen_count'] = st.session_state.hydrogen_count

        if 'atom_negation' in st.session_state and st.session_state.atom_negation:
            atom_data['negation'] = True

        # 표시 텍스트 생성
        atom_data['display_text'] = self.generate_atom_display_text(atom_data)

        st.session_state.molecule_data['atoms'].append(atom_data)
        st.session_state.molecule_data['atom_counter'] += 1

        st.rerun()

    def generate_atom_display_text(self, atom_data):
        """원자 표시 텍스트 생성"""
        element = atom_data['element']

        if atom_data.get('negation'):
            return f"[!{element}]"
        elif atom_data.get('aromatic'):
            return element.lower()
        elif atom_data.get('charge'):
            charge = atom_data['charge']
            if charge > 0:
                return f"{element}+{charge}" if charge > 1 else f"{element}+"
            else:
                return f"{element}{charge}" if charge < -1 else f"{element}-"
        else:
            return element

    def handle_bond_creation(self, x: float, y: float):
        """결합 생성 처리"""
        clicked_atom = self.find_atom_at(x, y)

        if clicked_atom:
            if st.session_state.temp_bond_start_id is None:
                st.session_state.temp_bond_start_id = clicked_atom['id']
                st.success(f"첫 번째 원자 선택됨 (ID: {clicked_atom['id']})")
            else:
                if clicked_atom['id'] != st.session_state.temp_bond_start_id:
                    self.create_bond(st.session_state.temp_bond_start_id, clicked_atom['id'])
                    st.success("결합이 생성되었습니다!")
                st.session_state.temp_bond_start_id = None
                st.rerun()

    def create_bond(self, atom1_id: int, atom2_id: int):
        """결합 생성"""
        # 이미 존재하는 결합인지 확인
        for bond in st.session_state.molecule_data['bonds']:
            if ((bond['atom1_id'] == atom1_id and bond['atom2_id'] == atom2_id) or
                (bond['atom1_id'] == atom2_id and bond['atom2_id'] == atom1_id)):
                return

        bond_data = {
            'id': st.session_state.molecule_data['bond_counter'],
            'atom1_id': atom1_id,
            'atom2_id': atom2_id,
            'bond_type': st.session_state.current_bond_type.value
        }

        # 고리 조건 추가
        if 'bond_ring_condition' in st.session_state:
            if st.session_state.bond_ring_condition == '고리 내 (@)':
                bond_data['ring_membership'] = True
            elif st.session_state.bond_ring_condition == '고리 외 (!@)':
                bond_data['ring_membership'] = False

        st.session_state.molecule_data['bonds'].append(bond_data)
        st.session_state.molecule_data['bond_counter'] += 1

    def find_atom_at(self, x: float, y: float, radius: float = 20):
        """지정된 위치에서 원자 찾기"""
        for atom in st.session_state.molecule_data['atoms']:
            distance = calculate_distance(atom['x'], atom['y'], x, y)
            if distance <= radius:
                return atom
        return None

    def get_atom_by_id(self, atom_id: int):
        """ID로 원자 찾기"""
        for atom in st.session_state.molecule_data['atoms']:
            if atom['id'] == atom_id:
                return atom
        return None

    def select_atom_at(self, x: float, y: float):
        """원자 선택"""
        atom = self.find_atom_at(x, y)
        if atom:
            st.session_state.selected_atom_id = atom['id']
            st.info(f"원자 선택됨: {atom['element']} (ID: {atom['id']})")

    def delete_at(self, x: float, y: float):
        """지정된 위치에서 삭제"""
        # 원자 삭제 시도
        atom = self.find_atom_at(x, y)
        if atom:
            self.delete_atom(atom['id'])
            st.success("원자가 삭제되었습니다.")
            st.rerun()
            return

        # 결합 삭제 시도
        bond = self.find_bond_at(x, y)
        if bond:
            self.delete_bond(bond['id'])
            st.success("결합이 삭제되었습니다.")
            st.rerun()

    def find_bond_at(self, x: float, y: float):
        """지정된 위치에서 결합 찾기"""
        for bond in st.session_state.molecule_data['bonds']:
            atom1 = self.get_atom_by_id(bond['atom1_id'])
            atom2 = self.get_atom_by_id(bond['atom2_id'])

            if atom1 and atom2:
                if point_on_line(x, y, atom1['x'], atom1['y'], atom2['x'], atom2['y'], 10):
                    return bond
        return None

    def delete_atom(self, atom_id: int):
        """원자 삭제 (연결된 결합도 함께)"""
        # 연결된 결합 삭제
        st.session_state.molecule_data['bonds'] = [
            bond for bond in st.session_state.molecule_data['bonds']
            if bond['atom1_id'] != atom_id and bond['atom2_id'] != atom_id
        ]

        # 원자 삭제
        st.session_state.molecule_data['atoms'] = [
            atom for atom in st.session_state.molecule_data['atoms']
            if atom['id'] != atom_id
        ]

    def delete_bond(self, bond_id: int):
        """결합 삭제"""
        st.session_state.molecule_data['bonds'] = [
            bond for bond in st.session_state.molecule_data['bonds']
            if bond['id'] != bond_id
        ]

    def clear_molecule(self):
        """분자 전체 삭제"""
        st.session_state.molecule_data = {
            'atoms': [],
            'bonds': [],
            'atom_counter': 0,
            'bond_counter': 0
        }
        st.session_state.selected_atom_id = None
        st.session_state.temp_bond_start_id = None
        st.success("분자가 모두 삭제되었습니다.")
        st.rerun()

    def undo_last_action(self):
        """마지막 작업 취소"""
        if st.session_state.molecule_data['atoms']:
            last_atom = st.session_state.molecule_data['atoms'][-1]
            self.delete_atom(last_atom['id'])
            st.success("마지막 작업이 취소되었습니다.")
            st.rerun()

    def save_molecule(self):
        """분자 저장"""
        molecule_json = json.dumps(st.session_state.molecule_data, indent=2)
        st.download_button(
            label="💾 분자 구조 다운로드",
            data=molecule_json,
            file_name="molecule.json",
            mime="application/json"
        )

    def load_molecule(self):
        """분자 불러오기"""
        uploaded_file = st.file_uploader(
            "분자 구조 파일 업로드",
            type=['json'],
            key="molecule_upload"
        )

        if uploaded_file is not None:
            try:
                molecule_data = json.load(uploaded_file)
                st.session_state.molecule_data = molecule_data
                st.success("분자 구조가 성공적으로 불러와졌습니다!")
                st.rerun()
            except Exception as e:
                st.error(f"파일 로드 중 오류가 발생했습니다: {str(e)}")

    def load_molecule_from_session(self):
        """세션에서 분자 데이터 로드"""
        self.molecule.clear()

        # 원자 로드
        for atom_data in st.session_state.molecule_data['atoms']:
            atom = Atom(
                id=atom_data['id'],
                x=atom_data['x'],
                y=atom_data['y'],
                element=atom_data['element'],
                atom_type=AtomType.WILDCARD if atom_data['element'] == '*' else AtomType.ELEMENT
            )

            # SMARTS 쿼리 조건 설정
            query = AtomQuery(
                element=atom_data['element'],
                aromatic=atom_data.get('aromatic'),
                charge=atom_data.get('charge'),
                hydrogen_count=atom_data.get('hydrogen_count'),
                is_wildcard=atom_data['element'] == '*',
                negation=atom_data.get('negation', False)
            )
            atom.query = query
            self.molecule.atoms.append(atom)

        # 결합 로드
        for bond_data in st.session_state.molecule_data['bonds']:
            bond = Bond(
                id=bond_data['id'],
                atom1_id=bond_data['atom1_id'],
                atom2_id=bond_data['atom2_id'],
                bond_type=BondType(bond_data['bond_type'])
            )

            # SMARTS 쿼리 조건 설정
            query = BondQuery(
                bond_type=BondType(bond_data['bond_type']),
                ring_membership=bond_data.get('ring_membership')
            )
            bond.query = query
            self.molecule.bonds.append(bond)

    def add_predefined_structure(self, smiles: str):
        """미리 정의된 구조 추가"""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # 2D 좌표 생성
                AllChem.Compute2DCoords(mol)

                # 현재 분자에 추가 (기존 것 유지)
                conf = mol.GetConformer()
                atom_id_map = {}

                # 원자 추가
                for i in range(mol.GetNumAtoms()):
                    pos = conf.GetAtomPosition(i)
                    atom_data = {
                        'id': st.session_state.molecule_data['atom_counter'],
                        'x': pos.x * 50 + 400,  # 스케일링 및 위치 조정
                        'y': pos.y * 50 + 300,
                        'element': mol.GetAtomWithIdx(i).GetSymbol(),
                        'atom_type': 'aromatic' if mol.GetAtomWithIdx(i).GetIsAromatic() else 'element',
                        'aromatic': mol.GetAtomWithIdx(i).GetIsAromatic(),
                        'display_text': mol.GetAtomWithIdx(i).GetSymbol().lower() if mol.GetAtomWithIdx(i).GetIsAromatic() else mol.GetAtomWithIdx(i).GetSymbol()
                    }

                    atom_id_map[i] = st.session_state.molecule_data['atom_counter']
                    st.session_state.molecule_data['atoms'].append(atom_data)
                    st.session_state.molecule_data['atom_counter'] += 1

                # 결합 추가
                for bond in mol.GetBonds():
                    bond_type = BondType.SINGLE
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        bond_type = BondType.DOUBLE
                    elif bond.GetBondType() == Chem.BondType.TRIPLE:
                        bond_type = BondType.TRIPLE
                    elif bond.GetBondType() == Chem.BondType.AROMATIC:
                        bond_type = BondType.AROMATIC

                    bond_data = {
                        'id': st.session_state.molecule_data['bond_counter'],
                        'atom1_id': atom_id_map[bond.GetBeginAtomIdx()],
                        'atom2_id': atom_id_map[bond.GetEndAtomIdx()],
                        'bond_type': bond_type.value
                    }

                    st.session_state.molecule_data['bonds'].append(bond_data)
                    st.session_state.molecule_data['bond_counter'] += 1

                st.success(f"'{smiles}' 구조가 추가되었습니다!")
                st.rerun()

        except Exception as e:
            st.error(f"구조 추가 중 오류가 발생했습니다: {str(e)}")

    def get_molecule(self) -> Molecule:
        """현재 분자 객체 반환"""
        self.load_molecule_from_session()
        return self.molecule