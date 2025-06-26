import streamlit as st
import sys
import os
from pathlib import Path

# 프로젝트 루트 디렉토리를 Python 패스에 추가
project_root = Path(__file__).parent
sys.path.append(str(project_root))

# 컴포넌트 임포트
from components.drawing import MoleculeDrawer, get_common_smarts_patterns, get_element_info
from components.patterns import SmartsGenerator, PatternMatcher

# 페이지 설정
st.set_page_config(
    page_title="SMARTS 패턴 편집기",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 사이드바 네비게이션
def create_sidebar():
    """사이드바 네비게이션 생성"""
    with st.sidebar:
        st.title("🧪 SMARTS 편집기")
        st.markdown("---")

        # 메인 네비게이션
        page = st.selectbox(
            "페이지 선택:",
            [
                "🏠 홈",
                "✏️ 분자 그리기",
                "🧬 SMARTS 생성",
                "🔍 패턴 매칭",
                "🔄 배치 매칭",
                "📚 패턴 라이브러리",
                "ℹ️ 도움말"
            ]
        )

        st.markdown("---")

        # 빠른 액세스
        st.subheader("빠른 액세스")

        # 자주 사용하는 패턴
        st.write("**자주 사용하는 패턴:**")
        common_patterns = get_common_smarts_patterns()

        for name, pattern in list(common_patterns.items())[:5]:
            if st.button(f"{pattern}", key=f"quick_{name}", help=name):
                st.session_state['selected_pattern'] = pattern
                st.success(f"'{pattern}' 선택됨")

        # 세션 정보
        st.markdown("---")
        st.subheader("세션 정보")

        if 'molecule_data' in st.session_state:
            atom_count = len(st.session_state.molecule_data.get('atoms', []))
            bond_count = len(st.session_state.molecule_data.get('bonds', []))
            st.write(f"**원자 수:** {atom_count}")
            st.write(f"**결합 수:** {bond_count}")

        if 'selected_pattern' in st.session_state:
            st.write(f"**선택된 패턴:** `{st.session_state.selected_pattern}`")

        # 세션 초기화
        if st.button("🔄 세션 초기화"):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.success("세션이 초기화되었습니다!")
            st.rerun()

    return page

def create_home_page():
    """홈 페이지"""
    st.title("🧪 SMARTS 패턴 편집기")
    st.markdown("---")

    # 소개
    st.markdown("""
    ## 환영합니다! 👋

    이 애플리케이션은 **SMARTS (SMiles ARbitrary Target Specification)** 패턴을 
    시각적으로 편집하고 테스트할 수 있는 통합 도구입니다.
    
    ### 🚀 새로운 JavaScript 기반 인터랙티브 편집기
    
    마우스 클릭만으로 분자를 그리고 편집할 수 있는 새로운 인터페이스를 제공합니다!
    """)

    # 주요 기능
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🎨 주요 기능")
        st.markdown("""
        - **🖱️ 인터랙티브 분자 그리기**: 마우스 클릭으로 직관적인 분자 편집
        - **🔵 원자 추가**: 원하는 위치에 클릭으로 원자 배치
        - **🔗 결합 생성**: 두 번의 클릭으로 원자 간 결합 생성
        - **⬟ 미리 정의된 구조**: 벤젠, 사이클로헥산 등 원클릭 추가
        - **🧬 SMARTS 생성**: 그린 분자에서 자동으로 SMARTS 패턴 생성
        - **🔍 패턴 매칭**: SMARTS 패턴으로 분자 검색
        - **💾 파일 관리**: JSON 형식으로 저장/불러오기
        """)

    with col2:
        st.subheader("🚀 시작하기")
        st.markdown("""
        1. **✏️ 분자 그리기** 페이지로 이동
        2. 캔버스에서 원하는 위치를 클릭하여 원자 추가
        3. 결합 모드로 전환하여 원자들을 연결
        4. **🧬 SMARTS 생성** 페이지에서 패턴 확인
        5. **🔍 패턴 매칭**에서 다른 분자들과 비교
        
        ### 🎯 편집 모드
        - **🔵 원자 모드**: 클릭으로 원자 추가
        - **🔗 결합 모드**: 두 원자를 클릭하여 결합 생성
        - **👆 선택 모드**: 원자 선택 및 이동
        - **🗑️ 삭제 모드**: 클릭하여 원자/결합 삭제
        """)

    # 인터페이스 미리보기
    st.markdown("---")
    st.subheader("🖥️ 인터페이스 미리보기")
    
    # 현재 세션 상태 확인
    if 'molecule_data' in st.session_state and st.session_state.molecule_data['atoms']:
        st.success("✅ 현재 세션에 분자 데이터가 있습니다!")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            atom_count = len(st.session_state.molecule_data['atoms'])
            st.metric("현재 원자 수", atom_count)
        
        with col2:
            bond_count = len(st.session_state.molecule_data['bonds'])
            st.metric("현재 결합 수", bond_count)
        
        with col3:
            elements = set(atom['element'] for atom in st.session_state.molecule_data['atoms'])
            st.metric("원소 종류", len(elements))
            if elements:
                st.write(f"**원소들**: {', '.join(sorted(elements))}")
    else:
        st.info("🎨 아직 분자가 그려지지 않았습니다. '분자 그리기' 페이지에서 시작해보세요!")

    # 기능 소개
    st.markdown("---")
    st.subheader("🔧 고급 기능")

    tab1, tab2, tab3 = st.tabs(["기본 편집", "구조 라이브러리", "분석 도구"])

    with tab1:
        st.markdown("""
        ### 기본 편집 도구
        
        **원소 지원:**
        - 주요 원소: C, N, O, S, P
        - 할로겐: F, Cl, Br, I  
        - 기타: H, 와일드카드 (*)
        
        **결합 타입:**
        - 단일 결합 (—)
        - 이중 결합 (=)
        - 삼중 결합 (≡)
        - 방향족 결합 (:)
        
        **편집 기능:**
        - 격자 스냅 (50px 단위)
        - 실시간 원자/결합 카운터
        - 실행 취소 기능
        """)

    with tab2:
        st.markdown("""
        ### 미리 정의된 구조
        
        **기본 고리:**
        - ⬟ 벤젠 (방향족)
        - ⬡ 사이클로헥산 (지방족)
        
        **복잡한 구조:**
        - 나프탈렌
        - 피리딘
        - 푸란
        - 이미다졸
        
        **작용기:**
        - 페닐기
        - 메틸기
        - 카르복실기
        """)

    with tab3:
        st.markdown("""
        ### 분석 및 검증 도구
        
        **분자 정보:**
        - 원소별 개수 통계
        - 결합 타입별 분포
        - 분자식 자동 계산
        
        **SMARTS 검증:**
        - 패턴 문법 체크
        - 구조 일치 확인
        - 부분 구조 하이라이트
        
        **파일 관리:**
        - JSON 형식 저장/불러오기
        - 구조 데이터 내보내기
        - 세션 백업/복원
        """)

    # 업데이트 정보
    st.markdown("---")
    st.subheader("📝 최근 업데이트")

    with st.expander("버전 2.0.0 - JavaScript 기반 리뉴얼", expanded=True):
        st.markdown("""
        ### 🆕 새로운 기능
        - ✅ **JavaScript 기반 인터랙티브 캔버스**: 실시간 편집 가능
        - ✅ **마우스 클릭 인터랙션**: 직관적인 분자 그리기
        - ✅ **실시간 미리보기**: 변경사항 즉시 반영
        - ✅ **향상된 성능**: 빠른 렌더링과 반응성
        - ✅ **격자 스냅 기능**: 정확한 원자 배치
        - ✅ **다양한 결합 타입**: 시각적으로 구분되는 결합 표현
        
        ### 🔧 개선사항
        - 불필요한 Plotly/Matplotlib 의존성 제거
        - 코드 최적화 및 정리
        - 세션 상태 관리 개선
        - 오류 처리 강화
        """)

def create_help_page():
    """도움말 페이지"""
    st.title("ℹ️ 도움말 - JavaScript 캔버스 사용법")
    st.markdown("---")

    # JavaScript 캔버스 사용법 중심으로 업데이트
    st.subheader("🖱️ 인터랙티브 캔버스 사용법")

    tab1, tab2, tab3, tab4 = st.tabs(["기본 조작", "고급 기능", "단축키", "문제 해결"])

    with tab1:
        st.markdown("""
        ### 🔵 원자 추가
        1. 상단 툴바에서 **🔵 원자** 버튼 클릭
        2. 원하는 원소 버튼 클릭 (C, N, O, S, P)
        3. 캔버스에서 원하는 위치 클릭
        4. 원자가 자동으로 격자에 맞춰 배치됩니다
        
        ### 🔗 결합 생성
        1. **🔗 결합** 버튼 클릭
        2. 결합 타입 선택 (—, =, ≡, :)
        3. 첫 번째 원자 클릭 (빨간색으로 선택됨)
        4. 두 번째 원자 클릭
        5. 결합이 자동으로 생성됩니다
        
        ### 👆 선택 및 이동
        1. **👆 선택** 버튼 클릭
        2. 원자를 클릭하여 선택
        3. 드래그하여 이동 (개발 예정)
        
        ### 🗑️ 삭제
        1. **🗑️ 삭제** 버튼 클릭
        2. 삭제할 원자나 결합 클릭
        3. 연결된 결합도 함께 삭제됩니다
        """)

    with tab2:
        st.markdown("""
        ### ⬟ 미리 정의된 구조
        - **⬟ 벤젠**: 방향족 6각 고리 자동 생성
        - **⬡ 사이클로헥산**: 지방족 6각 고리 자동 생성
        - 중앙 위치 (450, 300)에 배치됨
        
        ### 🎨 시각적 기능
        - **격자 배경**: 50px 간격 격자로 정확한 배치
        - **원소별 색상**: 각 원소마다 고유한 색상
        - **결합 스타일**: 
          - 단일: 실선 (—)
          - 이중: 평행선 (=)  
          - 삼중: 3개 평행선 (≡)
          - 방향족: 점선 (:)
        
        ### 📊 실시간 정보
        - 현재 모드 표시
        - 원자/결합 개수 실시간 업데이트
        - 선택된 원자 하이라이트
        """)

    with tab3:
        st.markdown("""
        ### ⌨️ 빠른 조작법
        
        **원소 빠른 선택:**
        - C, N, O, S, P 버튼으로 원소 변경과 동시에 원자 모드 활성화
        
        **결합 빠른 생성:**
        - —, =, ≡, : 버튼으로 결합 타입 변경과 동시에 결합 모드 활성화
        
        **빠른 초기화:**
        - **🗑️ 전체삭제** 버튼으로 모든 원자와 결합 삭제
        
        **구조 빠른 추가:**
        - **⬟ 벤젠**, **⬡ 사이클로헥산** 버튼으로 즉시 구조 추가
        
        **데이터 관리:**
        - **💾 내보내기** 버튼으로 콘솔에 JSON 데이터 출력
        """)

    with tab4:
        st.markdown("""
        ### 🔧 일반적인 문제들

        **1. 캔버스가 로드되지 않을 때**
        - 브라우저 새로고침 (F5)
        - 다른 브라우저에서 시도
        - JavaScript가 활성화되어 있는지 확인
        - 브라우저 콘솔에서 오류 메시지 확인

        **2. 원자가 추가되지 않을 때**
        - 원자 모드가 활성화되어 있는지 확인
        - 캔버스 내부를 클릭했는지 확인
        - 브라우저 콘솔에서 오류 메시지 확인

        **3. 결합이 생성되지 않을 때**
        - 결합 모드가 활성화되어 있는지 확인
        - 두 개의 서로 다른 원자를 클릭했는지 확인
        - 이미 존재하는 결합은 중복 생성되지 않음

        **4. 데이터가 Streamlit에 전달되지 않을 때**
        - 세션 상태 확인
        - 페이지 새로고침 후 재시도
        - 분자 데이터 표시 체크박스로 데이터 확인

        **5. 성능 문제**
        - 너무 많은 원자 (100개 이상) 사용 시 느려질 수 있음
        - 브라우저 탭을 너무 많이 열지 말 것
        - 필요 없는 분자는 삭제하여 성능 향상
        """)

    # 기술적 세부사항
    st.markdown("---")
    st.subheader("🔬 기술적 세부사항")

    with st.expander("JavaScript 캔버스 아키텍처", expanded=False):
        st.markdown("""
        ### 구현 세부사항
        
        **렌더링 엔진:**
        - SVG 기반 벡터 그래픽
        - 실시간 DOM 조작
        - 이벤트 기반 인터랙션
        
        **데이터 구조:**
        ```javascript
        atoms = [{
            id: 1,
            x: 450,
            y: 300,
            element: 'C'
        }]
        
        bonds = [{
            id: 1,
            atom1Id: 1,
            atom2Id: 2,
            type: 'single'
        }]
        ```
        
        **Streamlit 통합:**
        - `streamlit.components.v1.html()` 사용
        - `window.parent.postMessage()` 를 통한 데이터 전송
        - 세션 상태와 실시간 동기화
        """)

def main():
    """메인 함수"""
    # 사이드바 네비게이션
    page = create_sidebar()

    # 페이지별 라우팅
    if page == "🏠 홈":
        create_home_page()

    elif page == "✏️ 분자 그리기":
        st.title("✏️ 인터랙티브 분자 그리기")
        st.markdown("---")
        
        # JavaScript 기반 인터페이스 사용
        drawer = MoleculeDrawer()
        drawer.create_drawing_interface()

    elif page == "🧬 SMARTS 생성":
        st.title("🧬 SMARTS 생성")
        st.markdown("---")

        # 현재 분자 로드
        drawer = MoleculeDrawer()
        drawer.initialize_session_state()
        
        # 분자 데이터 확인
        if not st.session_state.molecule_data['atoms']:
            st.warning("⚠️ 먼저 분자를 그려주세요!")
            st.info("**✏️ 분자 그리기** 페이지에서 JavaScript 캔버스를 사용하여 분자를 그린 후 다시 시도하세요.")
            
            # 빠른 시작 버튼들
            st.markdown("### 🚀 빠른 시작")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                if st.button("⬟ 벤젠으로 시작", type="primary"):
                    # 벤젠 구조 직접 추가
                    drawer.add_predefined_structure("c1ccccc1")
                    st.success("벤젠 구조가 추가되었습니다!")
                    st.rerun()
            
            with col2:
                if st.button("🧪 에탄올로 시작"):
                    # 에탄올 구조 직접 추가  
                    drawer.add_predefined_structure("CCO")
                    st.success("에탄올 구조가 추가되었습니다!")
                    st.rerun()
            
            with col3:
                if st.button("💊 아스피린으로 시작"):
                    # 아스피린 구조 직접 추가
                    drawer.add_predefined_structure("CC(=O)OC1=CC=CC=C1C(=O)O")
                    st.success("아스피린 구조가 추가되었습니다!")
                    st.rerun()
                    
        else:
            generator = SmartsGenerator()
            # 세션 상태의 분자 데이터를 직접 전달
            smarts = generator.create_smarts_interface(st.session_state.molecule_data)

            # 생성된 SMARTS를 세션에 저장
            if smarts:
                st.session_state['selected_pattern'] = smarts

    elif page == "🔍 패턴 매칭":
        st.title("🔍 패턴 매칭")
        st.markdown("---")

        matcher = PatternMatcher()
        matcher.create_pattern_matching_interface()

    elif page == "🔄 배치 매칭":
        st.title("🔄 배치 매칭")
        st.markdown("---")

        matcher = PatternMatcher()
        matcher.create_batch_matching_interface()

    elif page == "📚 패턴 라이브러리":
        st.title("📚 패턴 라이브러리")
        st.markdown("---")

        generator = SmartsGenerator()
        generator.create_pattern_library_interface()

    elif page == "ℹ️ 도움말":
        create_help_page()

    # 푸터
    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center; color: gray; padding: 20px;'>
            🧪 SMARTS 패턴 편집기 v2.0.0 | 
            JavaScript 기반 인터랙티브 캔버스 | 
            Made with ❤️ using Streamlit & JavaScript
        </div>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()