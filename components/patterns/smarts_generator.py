import streamlit as st
from typing import List, Dict, Optional, Tuple
import networkx as nx
from rdkit import Chem
import sys
import os

# 프로젝트 루트를 Python 경로에 추가
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from models.molecule import Molecule
from models.atom import Atom, AtomType, AtomQuery
from models.bond import Bond, BondType, BondQuery
from components.drawing.utils import validate_smarts

class SmartsGenerator:
    """SMARTS 패턴 생성기"""

    def __init__(self):
        pass

    def generate_smarts_from_molecule(self, molecule: Molecule) -> str:
        """분자 객체에서 SMARTS 패턴 생성"""
        if not molecule.atoms:
            return ""

        try:
            # NetworkX 그래프로 변환
            graph = molecule.to_networkx()

            # 연결 컴포넌트별로 처리
            components = list(nx.connected_components(graph))
            smarts_parts = []

            for component in components:
                if len(component) == 1:
                    # 단일 원자
                    atom_id = list(component)[0]
                    atom = molecule.get_atom_by_id(atom_id)
                    if atom:
                        smarts_parts.append(atom.to_smarts())
                else:
                    # 연결된 분자 부분
                    subgraph = graph.subgraph(component)
                    smarts = self._generate_smarts_from_subgraph(subgraph, molecule)
                    if smarts:
                        smarts_parts.append(smarts)

            # 여러 컴포넌트를 점(.)으로 연결
            return ".".join(smarts_parts)

        except Exception as e:
            st.error(f"SMARTS 생성 중 오류: {str(e)}")
            return ""

    def _generate_smarts_from_subgraph(self, subgraph: nx.Graph, molecule: Molecule) -> str:
        """서브그래프에서 SMARTS 생성"""
        if len(subgraph.nodes()) == 0:
            return ""

        # 시작 원자 선택 (도수가 가장 낮은 원자부터)
        start_atom = min(subgraph.nodes(), key=lambda x: subgraph.degree(x))

        # DFS로 SMARTS 문자열 구성
        visited = set()
        ring_closures = {}
        ring_counter = 1

        def dfs_smarts(atom_id: int, parent_id: Optional[int] = None) -> str:
            nonlocal ring_counter

            if atom_id in visited:
                # 고리 닫힘 처리
                if atom_id not in ring_closures:
                    ring_closures[atom_id] = ring_counter
                    ring_counter += 1
                return str(ring_closures[atom_id])

            visited.add(atom_id)
            atom = molecule.get_atom_by_id(atom_id)
            if not atom:
                return ""

            # 원자 SMARTS
            atom_smarts = atom.to_smarts()

            # 연결된 원자들 처리
            neighbors = list(subgraph.neighbors(atom_id))
            if parent_id:
                neighbors = [n for n in neighbors if n != parent_id]

            # 결합과 원자 추가
            branches = []
            for neighbor_id in neighbors:
                # 결합 SMARTS
                bond = self._find_bond_between_atoms(molecule, atom_id, neighbor_id)
                bond_smarts = bond.to_smarts() if bond else ""

                # 재귀적으로 이웃 원자 처리
                neighbor_smarts = dfs_smarts(neighbor_id, atom_id)

                if neighbor_smarts:
                    full_branch = bond_smarts + neighbor_smarts
                    branches.append(full_branch)

            # 분기 처리
            if len(branches) == 0:
                result = atom_smarts
            elif len(branches) == 1:
                result = atom_smarts + branches[0]
            else:
                # 여러 분기는 괄호로 묶기
                main_branch = branches[0]
                side_branches = [f"({branch})" for branch in branches[1:]]
                result = atom_smarts + main_branch + "".join(side_branches)

            # 고리 닫힘 번호 추가
            if atom_id in ring_closures:
                result += str(ring_closures[atom_id])

            return result

        return dfs_smarts(start_atom)

    def _find_bond_between_atoms(self, molecule: Molecule, atom1_id: int, atom2_id: int) -> Optional[Bond]:
        """두 원자 사이의 결합 찾기"""
        for bond in molecule.bonds:
            if ((bond.atom1_id == atom1_id and bond.atom2_id == atom2_id) or
                (bond.atom1_id == atom2_id and bond.atom2_id == atom1_id)):
                return bond
        return None

    def create_smarts_interface(self, molecule: Molecule):
        """SMARTS 생성 인터페이스"""
        st.header("🧬 SMARTS 패턴 생성")

        col1, col2 = st.columns([2, 1])

        with col1:
            # 자동 생성된 SMARTS
            smarts = self.generate_smarts_from_molecule(molecule)

            st.subheader("생성된 SMARTS 패턴")
            smarts_input = st.text_area(
                "SMARTS 패턴:",
                value=smarts,
                height=100,
                help="자동 생성된 SMARTS 패턴입니다. 직접 수정할 수 있습니다."
            )

            # 유효성 검사
            if smarts_input:
                is_valid, message = validate_smarts(smarts_input)
                if is_valid:
                    st.success(f"✅ {message}")
                else:
                    st.error(f"❌ {message}")

            # SMARTS 패턴 설명
            if smarts_input:
                st.subheader("패턴 분석")
                self._explain_smarts_pattern(smarts_input)

        with col2:
            st.subheader("SMARTS 도구")

            # 일반적인 SMARTS 패턴들
            st.write("**자주 사용되는 패턴:**")

            common_patterns = {
                "임의 원자": "*",
                "방향족 탄소": "c",
                "지방족 탄소": "C",
                "비수소 원자": "[!H]",
                "양전하": "[+]",
                "음전하": "[-]",
                "고리 내 원자": "[R]",
                "고리 외 원자": "[R0]",
                "1차 탄소": "[CH3]",
                "2차 탄소": "[CH2]",
                "3차 탄소": "[CH]",
                "4차 탄소": "[C]",
                "방향족 질소": "n",
                "지방족 질소": "N",
                "카르보닐 탄소": "[C]=[O]",
                "하이드록실": "[OH]",
                "할로겐": "[F,Cl,Br,I]"
            }

            for name, pattern in common_patterns.items():
                if st.button(f"{name}", key=f"pattern_{name}"):
                    st.session_state['selected_pattern'] = pattern
                    st.info(f"선택된 패턴: `{pattern}`")

            # 패턴 조합 도구
            st.write("**패턴 조합:**")

            if st.button("OR 조건 [A,B]"):
                st.info("예: [C,N] = 탄소 또는 질소")

            if st.button("NOT 조건 [!A]"):
                st.info("예: [!C] = 탄소가 아닌 원자")

            if st.button("AND 조건 [A&B]"):
                st.info("예: [C&R] = 고리 내 탄소")

            # 결합 패턴
            st.write("**결합 패턴:**")
            bond_patterns = {
                "단일결합": "-",
                "이중결합": "=",
                "삼중결합": "#",
                "방향족결합": ":",
                "임의결합": "~",
                "고리내결합": "@",
                "고리외결합": "!@"
            }

            for name, pattern in bond_patterns.items():
                if st.button(f"{name}", key=f"bond_{name}"):
                    st.info(f"결합 패턴: `{pattern}`")

        return smarts_input

    def _explain_smarts_pattern(self, smarts: str):
        """SMARTS 패턴 설명"""
        try:
            # RDKit으로 패턴 파싱
            mol = Chem.MolFromSmarts(smarts)
            if mol is None:
                st.warning("패턴을 분석할 수 없습니다.")
                return

            # 기본 정보
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()

            st.write(f"**원자 수:** {num_atoms}")
            st.write(f"**결합 수:** {num_bonds}")

            # 원자별 분석
            if num_atoms > 0:
                st.write("**원자 분석:**")
                for i, atom in enumerate(mol.GetAtoms()):
                    atom_info = []

                    # 원소
                    if atom.GetAtomicNum() != 0:
                        atom_info.append(f"원소: {atom.GetSymbol()}")
                    else:
                        atom_info.append("임의 원자")

                    # 방향족성
                    if atom.GetIsAromatic():
                        atom_info.append("방향족")

                    # 전하
                    if atom.GetFormalCharge() != 0:
                        atom_info.append(f"전하: {atom.GetFormalCharge():+d}")

                    # 수소 개수
                    if atom.GetTotalNumHs() > 0:
                        atom_info.append(f"수소: {atom.GetTotalNumHs()}")

                    st.write(f"  원자 {i+1}: {', '.join(atom_info)}")

            # 결합별 분석
            if num_bonds > 0:
                st.write("**결합 분석:**")
                for i, bond in enumerate(mol.GetBonds()):
                    bond_info = []

                    # 결합 타입
                    bond_type = bond.GetBondType()
                    if bond_type == Chem.BondType.SINGLE:
                        bond_info.append("단일결합")
                    elif bond_type == Chem.BondType.DOUBLE:
                        bond_info.append("이중결합")
                    elif bond_type == Chem.BondType.TRIPLE:
                        bond_info.append("삼중결합")
                    elif bond_type == Chem.BondType.AROMATIC:
                        bond_info.append("방향족결합")

                    # 연결된 원자
                    atom1_idx = bond.GetBeginAtomIdx()
                    atom2_idx = bond.GetEndAtomIdx()
                    bond_info.append(f"원자 {atom1_idx+1}-{atom2_idx+1}")

                    st.write(f"  결합 {i+1}: {', '.join(bond_info)}")

        except Exception as e:
            st.warning(f"패턴 분석 중 오류: {str(e)}")

    def create_pattern_library_interface(self):
        """패턴 라이브러리 인터페이스"""
        st.header("📚 SMARTS 패턴 라이브러리")

        # 카테고리별 패턴
        categories = {
            "기본 원자 패턴": {
                "임의 원자": "*",
                "탄소": "C",
                "질소": "N",
                "산소": "O",
                "황": "S",
                "인": "P",
                "할로겐": "[F,Cl,Br,I]",
                "비수소": "[!H]",
                "방향족 탄소": "c",
                "지방족 탄소": "C"
            },
            "전하 패턴": {
                "양전하": "[+]",
                "음전하": "[-]",
                "+1 전하": "[+1]",
                "-1 전하": "[-1]",
                "중성": "[+0]"
            },
            "수소 패턴": {
                "수소 없음": "[H0]",
                "수소 1개": "[H1]",
                "수소 2개": "[H2]",
                "수소 3개": "[H3]",
                "수소 있음": "[H]"
            },
            "고리 패턴": {
                "고리 내": "[R]",
                "고리 외": "[R0]",
                "5원환": "[R1]",
                "6원환": "[R2]",
                "방향족 고리": "[r6]"
            },
            "작용기 패턴": {
                "하이드록실": "[OH]",
                "카르보닐": "C=O",
                "카르복실": "C(=O)O",
                "에스테르": "C(=O)O[!H]",
                "아미드": "C(=O)N",
                "니트로": "[N+](=O)[O-]",
                "술폰산": "S(=O)(=O)O",
                "아미노": "N[!H]",
                "알데히드": "C(=O)[H]",
                "케톤": "C(=O)[!H]",
                "에테르": "O([!H])[!H]",
                "티올": "[SH]"
            },
            "결합 패턴": {
                "단일결합": "-",
                "이중결합": "=",
                "삼중결합": "#",
                "방향족결합": ":",
                "임의결합": "~",
                "고리내결합": "@",
                "고리외결합": "!@"
            }
        }

        # 카테고리 선택
        selected_category = st.selectbox(
            "카테고리 선택:",
            list(categories.keys())
        )

        # 선택된 카테고리의 패턴들 표시
        if selected_category:
            patterns = categories[selected_category]

            st.subheader(f"{selected_category}")

            # 패턴을 표 형태로 표시
            for name, pattern in patterns.items():
                col1, col2, col3 = st.columns([2, 2, 1])

                with col1:
                    st.write(f"**{name}**")

                with col2:
                    st.code(pattern)

                with col3:
                    if st.button("사용", key=f"use_{name}"):
                        st.session_state['selected_pattern'] = pattern
                        st.success(f"'{pattern}' 패턴이 선택되었습니다!")

        # 사용자 정의 패턴 저장
        st.subheader("사용자 정의 패턴")

        with st.expander("새 패턴 추가"):
            pattern_name = st.text_input("패턴 이름:")
            pattern_smarts = st.text_input("SMARTS 패턴:")
            pattern_description = st.text_area("설명:")

            if st.button("패턴 저장"):
                if pattern_name and pattern_smarts:
                    # 세션 상태에 사용자 패턴 저장
                    if 'user_patterns' not in st.session_state:
                        st.session_state.user_patterns = {}

                    st.session_state.user_patterns[pattern_name] = {
                        'smarts': pattern_smarts,
                        'description': pattern_description
                    }

                    st.success(f"패턴 '{pattern_name}'이 저장되었습니다!")
                else:
                    st.error("패턴 이름과 SMARTS를 모두 입력해주세요.")

        # 저장된 사용자 패턴 표시
        if 'user_patterns' in st.session_state and st.session_state.user_patterns:
            st.subheader("저장된 사용자 패턴")

            for name, data in st.session_state.user_patterns.items():
                with st.expander(f"📌 {name}"):
                    st.code(data['smarts'])
                    if data['description']:
                        st.write(data['description'])

                    col1, col2 = st.columns(2)
                    with col1:
                        if st.button("사용", key=f"user_use_{name}"):
                            st.session_state['selected_pattern'] = data['smarts']
                            st.success(f"'{data['smarts']}' 패턴이 선택되었습니다!")

                    with col2:
                        if st.button("삭제", key=f"user_delete_{name}"):
                            del st.session_state.user_patterns[name]
                            st.success(f"패턴 '{name}'이 삭제되었습니다!")
                            st.rerun()