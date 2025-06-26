import streamlit as st
import json
from typing import List, Dict, Optional, Tuple, Any
import sys
import os

# 프로젝트 루트를 Python 경로에 추가
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from models.molecule import Molecule
from models.atom import Atom, AtomType, AtomQuery
from models.bond import Bond, BondType, BondQuery
from .utils import get_atom_color, calculate_distance, point_on_line

class MoleculeDrawer:
    """분자 그리기 컴포넌트 - JavaScript 캔버스 전용"""

    def __init__(self):
        self.molecule = Molecule()
        self.canvas_width = 900
        self.canvas_height = 600

    def initialize_session_state(self):
        """세션 상태 초기화"""
        if 'molecule_data' not in st.session_state:
            st.session_state.molecule_data = {
                'atoms': [],
                'bonds': [],
                'atom_counter': 1,
                'bond_counter': 1
            }

        if 'drawing_mode' not in st.session_state:
            st.session_state.drawing_mode = 'atom'

        if 'current_element' not in st.session_state:
            st.session_state.current_element = 'C'

        if 'current_bond_type' not in st.session_state:
            st.session_state.current_bond_type = BondType.SINGLE

    def create_drawing_interface(self):
        """JavaScript 캔버스 기반 그리기 인터페이스"""
        self.initialize_session_state()

        # JavaScript 캔버스 컴포넌트 임포트 및 사용
        try:
            from components.drawing.js_molecule_canvas import create_interactive_canvas
            
            st.markdown("### 🧪 인터랙티브 분자 편집기")
            st.markdown("---")
            
            # JavaScript 캔버스 표시
            canvas_result = create_interactive_canvas()
            
            # 현재 분자 정보 표시
            col1, col2 = st.columns(2)
            
            with col1:
                atom_count = len(st.session_state.molecule_data.get('atoms', []))
                bond_count = len(st.session_state.molecule_data.get('bonds', []))
                st.metric("원자 수", atom_count)
                st.metric("결합 수", bond_count)
            
            with col2:
                # 분자 데이터 JSON 표시 (디버그용)
                if st.checkbox("🔍 분자 데이터 표시", key="show_molecule_data"):
                    if st.session_state.molecule_data['atoms'] or st.session_state.molecule_data['bonds']:
                        st.json(st.session_state.molecule_data)
                    else:
                        st.info("분자 데이터가 없습니다. 캔버스에 원자를 추가해보세요!")
            
            # 파일 작업 섹션
            with st.expander("💾 파일 작업", expanded=False):
                file_col1, file_col2 = st.columns(2)
                
                with file_col1:
                    st.subheader("저장")
                    if st.button("📁 JSON 내보내기", use_container_width=True):
                        if st.session_state.molecule_data['atoms']:
                            molecule_json = json.dumps(st.session_state.molecule_data, indent=2)
                            st.download_button(
                                label="💾 파일 다운로드",
                                data=molecule_json,
                                file_name="molecule_structure.json",
                                mime="application/json",
                                key="download_molecule"
                            )
                        else:
                            st.warning("저장할 분자 데이터가 없습니다.")
                
                with file_col2:
                    st.subheader("불러오기")
                    uploaded_file = st.file_uploader(
                        "JSON 파일 선택",
                        type=['json'],
                        key="upload_molecule"
                    )
                    
                    if uploaded_file is not None:
                        try:
                            molecule_data = json.load(uploaded_file)
                            # 데이터 유효성 검사
                            if self.validate_molecule_data(molecule_data):
                                st.session_state.molecule_data = molecule_data
                                st.success("분자 구조를 성공적으로 불러왔습니다!")
                                st.rerun()
                            else:
                                st.error("올바르지 않은 분자 데이터 형식입니다.")
                        except Exception as e:
                            st.error(f"파일 읽기 오류: {str(e)}")
            
            # 고급 분석 도구
            with st.expander("🔬 분자 분석", expanded=False):
                if st.session_state.molecule_data['atoms']:
                    analysis_col1, analysis_col2 = st.columns(2)
                    
                    with analysis_col1:
                        st.subheader("원소 분포")
                        element_counts = {}
                        for atom in st.session_state.molecule_data['atoms']:
                            element = atom['element']
                            element_counts[element] = element_counts.get(element, 0) + 1
                        
                        for element, count in element_counts.items():
                            st.write(f"**{element}**: {count}개")
                    
                    with analysis_col2:
                        st.subheader("결합 유형")
                        bond_type_names = {
                            1: "단일 결합",
                            2: "이중 결합", 
                            3: "삼중 결합",
                            4: "방향족 결합"
                        }
                        
                        bond_counts = {}
                        for bond in st.session_state.molecule_data['bonds']:
                            bond_type = bond['bond_type']
                            type_name = bond_type_names.get(bond_type, f"타입 {bond_type}")
                            bond_counts[type_name] = bond_counts.get(type_name, 0) + 1
                        
                        if bond_counts:
                            for bond_type, count in bond_counts.items():
                                st.write(f"**{bond_type}**: {count}개")
                        else:
                            st.write("결합이 없습니다.")
                else:
                    st.info("분석할 분자 데이터가 없습니다.")
            
            return canvas_result
            
        except ImportError as e:
            st.error("JavaScript 캔버스 컴포넌트를 불러올 수 없습니다.")
            st.error(f"오류 세부사항: {str(e)}")
            
            # 대체 인터페이스 제공
            st.warning("대신 기본 인터페이스를 사용합니다.")
            return self.create_fallback_interface()

    def create_fallback_interface(self):
        """JavaScript 캔버스를 사용할 수 없을 때의 대체 인터페이스"""
        st.markdown("### 📝 텍스트 기반 분자 편집기")
        
        # 간단한 분자 정보 입력
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("원자 추가")
            element = st.selectbox("원소", ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'H'])
            x_pos = st.number_input("X 좌표", 0, 900, 450)
            y_pos = st.number_input("Y 좌표", 0, 600, 300)
            
            if st.button("원자 추가"):
                self.add_atom_manual(x_pos, y_pos, element)
                st.success(f"{element} 원자가 추가되었습니다.")
                st.rerun()
        
        with col2:
            st.subheader("결합 추가")
            if st.session_state.molecule_data['atoms']:
                atom_options = [(atom['id'], f"{atom['element']} (ID: {atom['id']})") 
                              for atom in st.session_state.molecule_data['atoms']]
                
                atom1_id = st.selectbox("첫 번째 원자", 
                                      options=[opt[0] for opt in atom_options],
                                      format_func=lambda x: next(opt[1] for opt in atom_options if opt[0] == x))
                
                atom2_id = st.selectbox("두 번째 원자",
                                      options=[opt[0] for opt in atom_options if opt[0] != atom1_id],
                                      format_func=lambda x: next(opt[1] for opt in atom_options if opt[0] == x))
                
                bond_type = st.selectbox("결합 타입", [1, 2, 3, 4], 
                                       format_func=lambda x: {1: "단일", 2: "이중", 3: "삼중", 4: "방향족"}[x])
                
                if st.button("결합 추가"):
                    self.add_bond_manual(atom1_id, atom2_id, bond_type)
                    st.success("결합이 추가되었습니다.")
                    st.rerun()
            else:
                st.info("먼저 원자를 추가하세요.")
        
        # 현재 분자 상태 표시
        st.markdown("---")
        st.subheader("현재 분자 상태")
        
        if st.session_state.molecule_data['atoms'] or st.session_state.molecule_data['bonds']:
            col1, col2 = st.columns(2)
            
            with col1:
                if st.session_state.molecule_data['atoms']:
                    st.write("**원자 목록:**")
                    for atom in st.session_state.molecule_data['atoms']:
                        st.write(f"ID {atom['id']}: {atom['element']} at ({atom['x']}, {atom['y']})")
            
            with col2:
                if st.session_state.molecule_data['bonds']:
                    st.write("**결합 목록:**")
                    for bond in st.session_state.molecule_data['bonds']:
                        bond_name = {1: "단일", 2: "이중", 3: "삼중", 4: "방향족"}.get(bond['bond_type'], "알 수 없음")
                        st.write(f"ID {bond['id']}: {bond['atom1_id']} - {bond['atom2_id']} ({bond_name})")
        else:
            st.info("분자 데이터가 없습니다.")
        
        # 초기화 버튼
        if st.button("🗑️ 전체 삭제"):
            self.clear_molecule()
            st.success("분자가 삭제되었습니다.")
            st.rerun()

    def add_atom_manual(self, x: float, y: float, element: str):
        """수동으로 원자 추가"""
        atom_data = {
            'id': st.session_state.molecule_data['atom_counter'],
            'x': x,
            'y': y,
            'element': element,
            'atom_type': 'element',
            'display_text': element
        }
        
        st.session_state.molecule_data['atoms'].append(atom_data)
        st.session_state.molecule_data['atom_counter'] += 1

    def add_bond_manual(self, atom1_id: int, atom2_id: int, bond_type: int):
        """수동으로 결합 추가"""
        # 중복 결합 체크
        for bond in st.session_state.molecule_data['bonds']:
            if ((bond['atom1_id'] == atom1_id and bond['atom2_id'] == atom2_id) or
                (bond['atom1_id'] == atom2_id and bond['atom2_id'] == atom1_id)):
                st.warning("이미 존재하는 결합입니다.")
                return
        
        bond_data = {
            'id': st.session_state.molecule_data['bond_counter'],
            'atom1_id': atom1_id,
            'atom2_id': atom2_id,
            'bond_type': bond_type
        }
        
        st.session_state.molecule_data['bonds'].append(bond_data)
        st.session_state.molecule_data['bond_counter'] += 1

    def validate_molecule_data(self, data: dict) -> bool:
        """분자 데이터 유효성 검사"""
        try:
            # 필수 키 확인
            required_keys = ['atoms', 'bonds', 'atom_counter', 'bond_counter']
            if not all(key in data for key in required_keys):
                return False
            
            # 원자 데이터 검사
            for atom in data['atoms']:
                required_atom_keys = ['id', 'x', 'y', 'element']
                if not all(key in atom for key in required_atom_keys):
                    return False
            
            # 결합 데이터 검사
            for bond in data['bonds']:
                required_bond_keys = ['id', 'atom1_id', 'atom2_id', 'bond_type']
                if not all(key in bond for key in required_bond_keys):
                    return False
            
            return True
        except:
            return False

    def get_atom_by_id(self, atom_id: int):
        """ID로 원자 찾기"""
        for atom in st.session_state.molecule_data['atoms']:
            if atom['id'] == atom_id:
                return atom
        return None

    def clear_molecule(self):
        """분자 전체 삭제"""
        st.session_state.molecule_data = {
            'atoms': [],
            'bonds': [],
            'atom_counter': 1,
            'bond_counter': 1
        }

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

    def get_molecule(self) -> Molecule:
        """현재 분자 객체 반환"""
        self.load_molecule_from_session()
        return self.molecule

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