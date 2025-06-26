import streamlit as st
from typing import List, Dict, Optional, Tuple, Set
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
from PIL import Image
import io
import base64

from models.molecule import Molecule
from components.utils import mol_to_image, image_to_base64, validate_smarts

class PatternMatcher:
    """SMARTS 패턴 매칭 기능"""

    def __init__(self):
        self.test_molecules = self._get_default_test_molecules()

    def _get_default_test_molecules(self) -> Dict[str, str]:
        """기본 테스트 분자들"""
        return {
            "아스피린": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "카페인": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "이부프로펜": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            "벤젠": "C1=CC=CC=C1",
            "톨루엔": "CC1=CC=CC=C1",
            "페놀": "C1=CC=C(C=C1)O",
            "아닐린": "C1=CC=C(C=C1)N",
            "아세트산": "CC(=O)O",
            "에탄올": "CCO",
            "메탄올": "CO",
            "아세톤": "CC(=O)C",
            "포름알데히드": "C=O",
            "글루코스": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
            "아데닌": "C1=NC(=C2C(=N1)N=CN2)N",
            "구아닌": "C1=NC2=C(N1)C(=O)NC(=N2)N",
            "시토신": "C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)O",
            "티민": "CC1=CN(C(=O)NC1=O)C2CC(C(O2)CO)O",
            "콜레스테롤": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
            "도파민": "C1=CC(=C(C=C1CCN)O)O",
            "세로토닌": "C1=CC2=C(C=C1O)C(=CN2)CCN",
            "아드레날린": "CNC(CC1=CC(=C(C=C1)O)O)C(=O)O"
        }

    def create_pattern_matching_interface(self):
        """패턴 매칭 인터페이스 생성"""
        st.header("🔍 SMARTS 패턴 매칭")

        # 패턴 입력
        col1, col2 = st.columns([2, 1])

        with col1:
            st.subheader("SMARTS 패턴")

            # 패턴 입력 방법 선택
            pattern_input_method = st.radio(
                "패턴 입력 방법:",
                ["직접 입력", "저장된 패턴에서 선택"],
                horizontal=True
            )

            if pattern_input_method == "직접 입력":
                smarts_pattern = st.text_input(
                    "SMARTS 패턴을 입력하세요:",
                    placeholder="예: c1ccccc1 (벤젠 고리)",
                    help="SMARTS 문법을 사용하여 패턴을 입력하세요."
                )
            else:
                # 저장된 패턴에서 선택
                if 'selected_pattern' in st.session_state:
                    default_pattern = st.session_state.selected_pattern
                else:
                    default_pattern = ""

                smarts_pattern = st.text_input(
                    "선택된 SMARTS 패턴:",
                    value=default_pattern,
                    help="패턴 라이브러리에서 선택된 패턴입니다."
                )

            # 패턴 유효성 검사
            if smarts_pattern:
                is_valid, message = validate_smarts(smarts_pattern)
                if is_valid:
                    st.success(f"✅ {message}")
                else:
                    st.error(f"❌ {message}")
                    return

        with col2:
            st.subheader("매칭 옵션")

            # 매칭 옵션들
            use_chirality = st.checkbox("입체화학 고려", value=False)
            use_aromaticity = st.checkbox("방향족성 고려", value=True)
            max_matches = st.slider("최대 매칭 수", 1, 100, 10)

            # 결과 표시 옵션
            show_molecule_images = st.checkbox("분자 구조 이미지 표시", value=True)
            show_match_details = st.checkbox("매칭 세부사항 표시", value=True)
            highlight_matches = st.checkbox("매칭 부분 하이라이트", value=True)

        if not smarts_pattern:
            st.info("SMARTS 패턴을 입력하세요.")
            return

        # 테스트할 분자 선택
        st.subheader("테스트 분자")

        # 분자 입력 방법 선택
        molecule_input_method = st.radio(
            "분자 입력 방법:",
            ["기본 분자에서 선택", "SMILES 직접 입력", "파일 업로드"],
            horizontal=True
        )

        test_molecules = []

        if molecule_input_method == "기본 분자에서 선택":
            selected_molecules = st.multiselect(
                "테스트할 분자를 선택하세요:",
                list(self.test_molecules.keys()),
                default=["벤젠", "톨루엔", "페놀", "아닐린"]
            )

            for mol_name in selected_molecules:
                test_molecules.append({
                    'name': mol_name,
                    'smiles': self.test_molecules[mol_name]
                })

        elif molecule_input_method == "SMILES 직접 입력":
            smiles_input = st.text_area(
                "SMILES를 한 줄에 하나씩 입력하세요:",
                placeholder="CCO\nCCCO\nC1=CC=CC=C1",
                height=150
            )

            if smiles_input:
                smiles_list = [s.strip() for s in smiles_input.split('\n') if s.strip()]
                for i, smiles in enumerate(smiles_list):
                    test_molecules.append({
                        'name': f"분자 {i+1}",
                        'smiles': smiles
                    })

        elif molecule_input_method == "파일 업로드":
            uploaded_file = st.file_uploader(
                "분자 파일 업로드 (CSV, TXT)",
                type=['csv', 'txt'],
                help="CSV: name,smiles 형식 또는 TXT: 한 줄에 하나씩 SMILES"
            )

            if uploaded_file:
                try:
                    if uploaded_file.name.endswith('.csv'):
                        df = pd.read_csv(uploaded_file)
                        if 'smiles' in df.columns:
                            name_col = 'name' if 'name' in df.columns else None
                            for i, row in df.iterrows():
                                name = row[name_col] if name_col else f"분자 {i+1}"
                                test_molecules.append({
                                    'name': str(name),
                                    'smiles': str(row['smiles'])
                                })
                    else:  # txt file
                        content = uploaded_file.read().decode('utf-8')
                        smiles_list = [s.strip() for s in content.split('\n') if s.strip()]
                        for i, smiles in enumerate(smiles_list):
                            test_molecules.append({
                                'name': f"분자 {i+1}",
                                'smiles': smiles
                            })

                    st.success(f"{len(test_molecules)}개 분자가 로드되었습니다.")

                except Exception as e:
                    st.error(f"파일 로드 중 오류: {str(e)}")

        if not test_molecules:
            st.info("테스트할 분자를 선택하거나 입력하세요.")
            return

        # 패턴 매칭 실행
        if st.button("🔍 패턴 매칭 실행", type="primary"):
            self.run_pattern_matching(
                smarts_pattern, 
                test_molecules,
                use_chirality=use_chirality,
                use_aromaticity=use_aromaticity,
                max_matches=max_matches,
                show_molecule_images=show_molecule_images,
                show_match_details=show_match_details,
                highlight_matches=highlight_matches
            )

    def run_pattern_matching(self, smarts_pattern: str, test_molecules: List[Dict], 
                           use_chirality: bool = False, use_aromaticity: bool = True,
                           max_matches: int = 10, show_molecule_images: bool = True,
                           show_match_details: bool = True, highlight_matches: bool = True):
        """패턴 매칭 실행"""

        st.subheader("📊 매칭 결과")

        # 진행률 표시
        progress_bar = st.progress(0)
        status_text = st.empty()

        # 패턴 분자 생성
        try:
            pattern_mol = Chem.MolFromSmarts(smarts_pattern)
            if pattern_mol is None:
                st.error("유효하지 않은 SMARTS 패턴입니다.")
                return
        except Exception as e:
            st.error(f"SMARTS 패턴 파싱 오류: {str(e)}")
            return

        # 결과 저장
        results = []
        matched_count = 0
        total_count = len(test_molecules)

        for i, mol_data in enumerate(test_molecules):
            status_text.text(f"처리 중: {mol_data['name']} ({i+1}/{total_count})")
            progress_bar.progress((i + 1) / total_count)

            try:
                # 분자 생성
                mol = Chem.MolFromSmiles(mol_data['smiles'])
                if mol is None:
                    results.append({
                        'name': mol_data['name'],
                        'smiles': mol_data['smiles'],
                        'matched': False,
                        'match_count': 0,
                        'error': 'Invalid SMILES'
                    })
                    continue

                # 패턴 매칭
                matches = mol.GetSubstructMatches(
                    pattern_mol,
                    useChirality=use_chirality,
                    useQueryQueryMatches=use_aromaticity,
                    maxMatches=max_matches
                )

                is_matched = len(matches) > 0
                if is_matched:
                    matched_count += 1

                result = {
                    'name': mol_data['name'],
                    'smiles': mol_data['smiles'],
                    'matched': is_matched,
                    'match_count': len(matches),
                    'matches': matches,
                    'mol': mol
                }

                # 분자 이미지 생성
                if show_molecule_images:
                    if highlight_matches and is_matched:
                        # 매칭된 부분 하이라이트
                        img = self._draw_molecule_with_highlights(mol, matches)
                    else:
                        img = mol_to_image(mol)
                    result['image'] = img

                results.append(result)

            except Exception as e:
                results.append({
                    'name': mol_data['name'],
                    'smiles': mol_data['smiles'],
                    'matched': False,
                    'match_count': 0,
                    'error': str(e)
                })

        # 진행률 바 제거
        progress_bar.empty()
        status_text.empty()

        # 결과 요약
        self._display_results_summary(results, matched_count, total_count, smarts_pattern)

        # 상세 결과 표시
        if show_match_details:
            self._display_detailed_results(results, show_molecule_images, highlight_matches)

        # 결과 다운로드 옵션
        self._create_download_options(results, smarts_pattern)

    def _draw_molecule_with_highlights(self, mol, matches):
        """매칭된 부분을 하이라이트하여 분자 그리기"""
        try:
            # 모든 매칭된 원자들을 하이라이트
            highlight_atoms = set()
            highlight_bonds = set()

            for match in matches:
                highlight_atoms.update(match)

                # 매칭된 원자들 사이의 결합도 하이라이트
                for i in range(len(match)):
                    for j in range(i + 1, len(match)):
                        bond = mol.GetBondBetweenAtoms(match[i], match[j])
                        if bond:
                            highlight_bonds.add(bond.GetIdx())

            # 하이라이트 색상 설정
            highlight_atom_colors = {atom_idx: (1.0, 0.0, 0.0) for atom_idx in highlight_atoms}
            highlight_bond_colors = {bond_idx: (1.0, 0.0, 0.0) for bond_idx in highlight_bonds}

            # 이미지 생성
            drawer = Draw.rdMolDraw2D.MolDraw2DCairo(300, 300)
            drawer.DrawMolecule(
                mol,
                highlightAtoms=list(highlight_atoms),
                highlightBonds=list(highlight_bonds),
                highlightAtomColors=highlight_atom_colors,
                highlightBondColors=highlight_bond_colors
            )
            drawer.FinishDrawing()

            # PIL 이미지로 변환
            img_data = drawer.GetDrawingText()
            img = Image.open(io.BytesIO(img_data))
            return img

        except Exception as e:
            # 오류 시 일반 이미지 반환
            return mol_to_image(mol)

    def _display_results_summary(self, results: List[Dict], matched_count: int, 
                                total_count: int, smarts_pattern: str):
        """결과 요약 표시"""
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("총 분자 수", total_count)

        with col2:
            st.metric("매칭된 분자", matched_count)

        with col3:
            match_rate = (matched_count / total_count * 100) if total_count > 0 else 0
            st.metric("매칭률", f"{match_rate:.1f}%")

        with col4:
            total_matches = sum(r.get('match_count', 0) for r in results)
            st.metric("총 매칭 수", total_matches)

        # 패턴 정보
        st.info(f"**검색 패턴:** `{smarts_pattern}`")

        # 매칭 통계 차트
        if len(results) > 1:
            self._create_matching_chart(results)

    def _create_matching_chart(self, results: List[Dict]):
        """매칭 통계 차트 생성"""
        import plotly.express as px
        import plotly.graph_objects as go

        # 매칭 여부별 분포
        matched_data = {
            'Status': ['매칭됨', '매칭되지 않음'],
            'Count': [
                len([r for r in results if r.get('matched', False)]),
                len([r for r in results if not r.get('matched', False)])
            ]
        }

        fig_pie = px.pie(
            values=matched_data['Count'],
            names=matched_data['Status'],
            title="매칭 분포",
            color_discrete_map={'매칭됨': '#00CC96', '매칭되지 않음': '#EF553B'}
        )

        col1, col2 = st.columns(2)

        with col1:
            st.plotly_chart(fig_pie, use_container_width=True)

        with col2:
            # 매칭 수별 분포
            match_counts = [r.get('match_count', 0) for r in results if r.get('matched', False)]
            if match_counts:
                fig_hist = px.histogram(
                    x=match_counts,
                    title="매칭 수 분포",
                    labels={'x': '매칭 수', 'y': '분자 수'},
                    nbins=max(10, len(set(match_counts)))
                )
                st.plotly_chart(fig_hist, use_container_width=True)

    def _display_detailed_results(self, results: List[Dict], show_images: bool, 
                                 highlight_matches: bool):
        """상세 결과 표시"""
        st.subheader("상세 매칭 결과")

        # 필터링 옵션
        col1, col2 = st.columns(2)

        with col1:
            show_filter = st.selectbox(
                "결과 필터:",
                ["전체", "매칭된 분자만", "매칭되지 않은 분자만"]
            )

        with col2:
            sort_by = st.selectbox(
                "정렬 기준:",
                ["이름", "매칭 수", "매칭 여부"]
            )

        # 결과 필터링 및 정렬
        filtered_results = results.copy()

        if show_filter == "매칭된 분자만":
            filtered_results = [r for r in filtered_results if r.get('matched', False)]
        elif show_filter == "매칭되지 않은 분자만":
            filtered_results = [r for r in filtered_results if not r.get('matched', False)]

        if sort_by == "매칭 수":
            filtered_results.sort(key=lambda x: x.get('match_count', 0), reverse=True)
        elif sort_by == "매칭 여부":
            filtered_results.sort(key=lambda x: x.get('matched', False), reverse=True)
        else:  # 이름
            filtered_results.sort(key=lambda x: x.get('name', ''))

        # 페이지네이션
        results_per_page = 10
        total_pages = (len(filtered_results) + results_per_page - 1) // results_per_page

        if total_pages > 1:
            page = st.selectbox(
                f"페이지 (총 {total_pages}페이지)",
                range(1, total_pages + 1)
            )
            start_idx = (page - 1) * results_per_page
            end_idx = start_idx + results_per_page
            page_results = filtered_results[start_idx:end_idx]
        else:
            page_results = filtered_results

        # 결과 표시
        for i, result in enumerate(page_results):
            with st.expander(
                f"{'✅' if result.get('matched', False) else '❌'} "
                f"{result['name']} "
                f"(매칭: {result.get('match_count', 0)}회)",
                expanded=result.get('matched', False)
            ):
                col1, col2 = st.columns([1, 2] if show_images else [1, 3])

                with col1:
                    if show_images and 'image' in result:
                        st.image(result['image'], caption=result['name'], width=200)

                with col2:
                    st.write(f"**분자명:** {result['name']}")
                    st.write(f"**SMILES:** `{result['smiles']}`")
                    st.write(f"**매칭 여부:** {'✅ 예' if result.get('matched', False) else '❌ 아니오'}")
                    st.write(f"**매칭 수:** {result.get('match_count', 0)}")

                    if 'error' in result:
                        st.error(f"오류: {result['error']}")

                    # 매칭 세부사항
                    if result.get('matched', False) and 'matches' in result:
                        with st.expander("매칭 세부사항"):
                            for j, match in enumerate(result['matches']):
                                st.write(f"**매칭 {j+1}:** 원자 인덱스 {list(match)}")

                    # 분자 특성 계산
                    if 'mol' in result:
                        self._display_molecule_properties(result['mol'])

    def _display_molecule_properties(self, mol):
        """분자 특성 표시"""
        try:
            with st.expander("분자 특성"):
                col1, col2 = st.columns(2)

                with col1:
                    st.write(f"**분자량:** {Descriptors.MolWt(mol):.2f}")
                    st.write(f"**LogP:** {Descriptors.MolLogP(mol):.2f}")
                    st.write(f"**원자 수:** {mol.GetNumAtoms()}")
                    st.write(f"**결합 수:** {mol.GetNumBonds()}")

                with col2:
                    st.write(f"**고리 수:** {Descriptors.RingCount(mol)}")
                    st.write(f"**방향족 고리 수:** {Descriptors.NumAromaticRings(mol)}")
                    st.write(f"**HBD:** {Descriptors.NumHDonors(mol)}")
                    st.write(f"**HBA:** {Descriptors.NumHAcceptors(mol)}")

        except Exception as e:
            st.write(f"특성 계산 오류: {str(e)}")

    def _create_download_options(self, results: List[Dict], smarts_pattern: str):
        """결과 다운로드 옵션"""
        st.subheader("📥 결과 다운로드")

        col1, col2, col3 = st.columns(3)

        with col1:
            # CSV 다운로드
            df_data = []
            for result in results:
                df_data.append({
                    'Name': result['name'],
                    'SMILES': result['smiles'],
                    'Matched': result.get('matched', False),
                    'Match_Count': result.get('match_count', 0),
                    'Error': result.get('error', '')
                })

            df = pd.DataFrame(df_data)
            csv = df.to_csv(index=False)

            st.download_button(
                label="📊 CSV 다운로드",
                data=csv,
                file_name=f"pattern_matching_results.csv",
                mime="text/csv"
            )

        with col2:
            # 매칭된 분자만 SMILES 다운로드
            matched_smiles = [
                result['smiles'] for result in results 
                if result.get('matched', False)
            ]

            if matched_smiles:
                smiles_text = '\n'.join(matched_smiles)
                st.download_button(
                    label="✅ 매칭된 SMILES",
                    data=smiles_text,
                    file_name=f"matched_molecules.smi",
                    mime="text/plain"
                )

        with col3:
            # 매칭되지 않은 분자 SMILES 다운로드
            unmatched_smiles = [
                result['smiles'] for result in results 
                if not result.get('matched', False)
            ]

            if unmatched_smiles:
                smiles_text = '\n'.join(unmatched_smiles)
                st.download_button(
                    label="❌ 매칭되지 않은 SMILES",
                    data=smiles_text,
                    file_name=f"unmatched_molecules.smi",
                    mime="text/plain"
                )

    def create_batch_matching_interface(self):
        """배치 매칭 인터페이스"""
        st.header("🔄 배치 패턴 매칭")
        st.write("여러 SMARTS 패턴을 한 번에 테스트할 수 있습니다.")

        # 패턴 입력
        st.subheader("SMARTS 패턴들")

        pattern_input_method = st.radio(
            "패턴 입력 방법:",
            ["직접 입력", "파일 업로드"],
            horizontal=True,
            key="batch_pattern_method"
        )

        patterns = []

        if pattern_input_method == "직접 입력":
            patterns_text = st.text_area(
                "SMARTS 패턴을 한 줄에 하나씩 입력하세요:",
                placeholder="c1ccccc1\nC=O\n[OH]\nN",
                height=150
            )

            if patterns_text:
                pattern_lines = [line.strip() for line in patterns_text.split('\n') if line.strip()]
                for i, pattern in enumerate(pattern_lines):
                    patterns.append({
                        'name': f"패턴 {i+1}",
                        'smarts': pattern
                    })

        else:  # 파일 업로드
            uploaded_file = st.file_uploader(
                "패턴 파일 업로드 (TXT 또는 CSV)",
                type=['txt', 'csv'],
                key="batch_pattern_file"
            )

            if uploaded_file:
                try:
                    if uploaded_file.name.endswith('.csv'):
                        df = pd.read_csv(uploaded_file)
                        if 'smarts' in df.columns:
                            name_col = 'name' if 'name' in df.columns else None
                            for i, row in df.iterrows():
                                name = row[name_col] if name_col else f"패턴 {i+1}"
                                patterns.append({
                                    'name': str(name),
                                    'smarts': str(row['smarts'])
                                })
                    else:  # txt file
                        content = uploaded_file.read().decode('utf-8')
                        pattern_lines = [line.strip() for line in content.split('\n') if line.strip()]
                        for i, pattern in enumerate(pattern_lines):
                            patterns.append({
                                'name': f"패턴 {i+1}",
                                'smarts': pattern
                            })

                except Exception as e:
                    st.error(f"파일 로드 중 오류: {str(e)}")

        if not patterns:
            st.info("SMARTS 패턴들을 입력하세요.")
            return

        st.write(f"**로드된 패턴 수:** {len(patterns)}")

        # 테스트 분자 선택 (간단화)
        st.subheader("테스트 분자")

        test_molecule_set = st.selectbox(
            "테스트 분자 세트:",
            ["기본 약물 분자", "기본 유기 분자", "모든 기본 분자"]
        )

        if test_molecule_set == "기본 약물 분자":
            test_mols = {k: v for k, v in self.test_molecules.items() 
                        if k in ["아스피린", "카페인", "이부프로펜", "도파민", "세로토닌", "아드레날린"]}
        elif test_molecule_set == "기본 유기 분자":
            test_mols = {k: v for k, v in self.test_molecules.items() 
                        if k in ["벤젠", "톨루엔", "페놀", "아닐린", "아세트산", "에탄올", "아세톤"]}
        else:
            test_mols = self.test_molecules

        test_molecules = [{'name': k, 'smiles': v} for k, v in test_mols.items()]

        # 배치 매칭 실행
        if st.button("🚀 배치 매칭 실행", type="primary"):
            self.run_batch_matching(patterns, test_molecules)

    def run_batch_matching(self, patterns: List[Dict], test_molecules: List[Dict]):
        """배치 매칭 실행"""
        st.subheader("📊 배치 매칭 결과")

        # 진행률 표시
        total_combinations = len(patterns) * len(test_molecules)
        progress_bar = st.progress(0)
        status_text = st.empty()

        # 결과 매트릭스 초기화
        results_matrix = []
        pattern_names = [p['name'] for p in patterns]
        molecule_names = [m['name'] for m in test_molecules]

        current_step = 0

        for i, pattern_data in enumerate(patterns):
            pattern_results = []

            try:
                pattern_mol = Chem.MolFromSmarts(pattern_data['smarts'])
                if pattern_mol is None:
                    pattern_results = [False] * len(test_molecules)
                    results_matrix.append(pattern_results)
                    continue
            except:
                pattern_results = [False] * len(test_molecules)
                results_matrix.append(pattern_results)
                continue

            for j, mol_data in enumerate(test_molecules):
                current_step += 1
                status_text.text(f"처리 중: {pattern_data['name']} vs {mol_data['name']} ({current_step}/{total_combinations})")
                progress_bar.progress(current_step / total_combinations)

                try:
                    mol = Chem.MolFromSmiles(mol_data['smiles'])
                    if mol is None:
                        pattern_results.append(False)
                        continue

                    matches = mol.GetSubstructMatches(pattern_mol)
                    pattern_results.append(len(matches) > 0)

                except:
                    pattern_results.append(False)

            results_matrix.append(pattern_results)

        # 진행률 바 제거
        progress_bar.empty()
        status_text.empty()

        # 결과 표시
        self._display_batch_results(results_matrix, pattern_names, molecule_names, patterns, test_molecules)

    def _display_batch_results(self, results_matrix: List[List[bool]], 
                              pattern_names: List[str], molecule_names: List[str],
                              patterns: List[Dict], test_molecules: List[Dict]):
        """배치 결과 표시"""

        # 결과 매트릭스를 DataFrame으로 변환
        df = pd.DataFrame(results_matrix, index=pattern_names, columns=molecule_names)

        # 요약 통계
        col1, col2, col3 = st.columns(3)

        with col1:
            total_matches = df.sum().sum()
            st.metric("총 매칭 수", int(total_matches))

        with col2:
            total_tests = len(patterns) * len(test_molecules)
            match_rate = (total_matches / total_tests * 100) if total_tests > 0 else 0
            st.metric("전체 매칭률", f"{match_rate:.1f}%")

        with col3:
            active_patterns = (df.sum(axis=1) > 0).sum()
            st.metric("활성 패턴 수", int(active_patterns))

        # 히트맵 생성
        import plotly.express as px

        # boolean을 숫자로 변환
        df_numeric = df.astype(int)

        fig = px.imshow(
            df_numeric,
            labels=dict(x="분자", y="SMARTS 패턴", color="매칭"),
            x=molecule_names,
            y=pattern_names,
            color_continuous_scale=['white', 'red'],
            aspect="auto"
        )

        fig.update_layout(
            title="SMARTS 패턴 매칭 히트맵",
            xaxis_title="테스트 분자",
            yaxis_title="SMARTS 패턴",
            height=max(400, len(patterns) * 30)
        )

        st.plotly_chart(fig, use_container_width=True)

        # 상세 결과 테이블
        st.subheader("상세 결과")

        # 패턴별 통계
        pattern_stats = []
        for i, pattern_name in enumerate(pattern_names):
            matches = df.iloc[i].sum()
            pattern_stats.append({
                'Pattern': pattern_name,
                'SMARTS': patterns[i]['smarts'],
                'Matches': int(matches),
                'Match_Rate': f"{matches/len(test_molecules)*100:.1f}%"
            })

        pattern_df = pd.DataFrame(pattern_stats)
        pattern_df = pattern_df.sort_values('Matches', ascending=False)

        st.write("**패턴별 매칭 통계:**")
        st.dataframe(pattern_df, use_container_width=True)

        # 분자별 통계
        molecule_stats = []
        for j, molecule_name in enumerate(molecule_names):
            matches = df.iloc[:, j].sum()
            molecule_stats.append({
                'Molecule': molecule_name,
                'SMILES': test_molecules[j]['smiles'],
                'Matches': int(matches),
                'Match_Rate': f"{matches/len(patterns)*100:.1f}%"
            })

        molecule_df = pd.DataFrame(molecule_stats)
        molecule_df = molecule_df.sort_values('Matches', ascending=False)

        st.write("**분자별 매칭 통계:**")
        st.dataframe(molecule_df, use_container_width=True)

        # 결과 다운로드
        st.subheader("결과 다운로드")

        col1, col2 = st.columns(2)

        with col1:
            # 매칭 매트릭스 다운로드
            csv_matrix = df.to_csv()
            st.download_button(
                label="📊 매칭 매트릭스 (CSV)",
                data=csv_matrix,
                file_name="batch_matching_matrix.csv",
                mime="text/csv"
            )

        with col2:
            # 통계 다운로드
            stats_data = {
                'Pattern_Stats': pattern_df.to_dict('records'),
                'Molecule_Stats': molecule_df.to_dict('records')
            }

            import json
            stats_json = json.dumps(stats_data, indent=2)

            st.download_button(
                label="📈 통계 데이터 (JSON)",
                data=stats_json,
                file_name="batch_matching_stats.json",
                mime="application/json"
            )