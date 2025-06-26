import streamlit as st
import sys
import os
from pathlib import Path

# 프로젝트 루트 디렉토리를 Python 패스에 추가
project_root = Path(__file__).parent
sys.path.append(str(project_root))

# 컴포넌트 임포트
from components.molecule_drawer import MoleculeDrawer
from components.smarts_generator import SmartsGenerator
from components.pattern_matcher import PatternMatcher
from components.utils import get_common_smarts_patterns, get_element_info

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
    """)

    # 주요 기능
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🎨 주요 기능")
        st.markdown("""
        - **분자 그리기**: 직관적인 GUI로 분자 구조 편집
        - **SMARTS 생성**: 그린 분자에서 자동으로 SMARTS 패턴 생성
        - **패턴 매칭**: SMARTS 패턴으로 분자 검색
        - **배치 처리**: 여러 패턴을 한 번에 테스트
        - **패턴 라이브러리**: 자주 사용하는 패턴 관리
        """)

    with col2:
        st.subheader("🚀 시작하기")
        st.markdown("""
        1. **분자 그리기** 페이지에서 분자를 그려보세요
        2. **SMARTS 생성** 페이지에서 패턴을 확인하세요
        3. **패턴 매칭** 페이지에서 다른 분자들과 비교하세요
        4. **패턴 라이브러리**에서 유용한 패턴들을 찾아보세요
        """)

    # 예시 섹션
    st.markdown("---")
    st.subheader("📋 사용 예시")

    tab1, tab2, tab3 = st.tabs(["기본 패턴", "작용기 검색", "약물 분석"])

    with tab1:
        st.markdown("""
        ### 기본 SMARTS 패턴 예시

        | 패턴 | 설명 | 예시 |
        |------|------|------|
        | `c1ccccc1` | 벤젠 고리 | 벤젠, 톨루엔, 페놀 |
        | `C=O` | 카르보닐 | 알데히드, 케톤 |
        | `[OH]` | 하이드록실 | 알코올, 페놀 |
        | `N` | 질소 원자 | 아민, 아미드 |
        | `*` | 임의 원자 | 모든 원자 |
        """)

    with tab2:
        st.markdown("""
        ### 작용기 검색 예시

        **카르복실산 검색:**
        ```
        패턴: C(=O)O
        매칭: 아세트산, 벤조산, 아미노산 등
        ```

        **에스테르 검색:**
        ```
        패턴: C(=O)O[!H]
        매칭: 아스피린, 지방산 에스테르 등
        ```

        **아미드 검색:**
        ```
        패턴: C(=O)N
        매칭: 단백질, 나일론 등
        ```
        """)

    with tab3:
        st.markdown("""
        ### 약물 분석 예시

        **NSAIDs (비스테로이드성 소염제) 패턴:**
        ```
        패턴: c1ccccc1C(=O)O
        매칭: 아스피린, 이부프로펜 등
        ```

        **베타 차단제 패턴:**
        ```
        패턴: c1ccccc1OCCNC
        매칭: 프로프라놀롤, 메토프롤롤 등
        ```
        """)

    # 통계 정보
    st.markdown("---")
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("지원 원소", "20+")

    with col2:
        st.metric("기본 패턴", "50+")

    with col3:
        st.metric("테스트 분자", "20+")

    with col4:
        st.metric("결합 타입", "5")

    # 최근 업데이트
    st.markdown("---")
    st.subheader("📝 최근 업데이트")

    with st.expander("버전 1.0.0 - 2024년"):
        st.markdown("""
        - ✅ 분자 그리기 기능 추가
        - ✅ SMARTS 패턴 자동 생성
        - ✅ 패턴 매칭 기능
        - ✅ 배치 처리 기능
        - ✅ 패턴 라이브러리
        - ✅ 파일 입출력 지원
        """)

def create_help_page():
    """도움말 페이지"""
    st.title("ℹ️ 도움말")
    st.markdown("---")

    # 목차
    st.subheader("📋 목차")
    st.markdown("""
    1. [SMARTS 문법 기초](#smarts-문법-기초)
    2. [분자 그리기 가이드](#분자-그리기-가이드)
    3. [패턴 매칭 가이드](#패턴-매칭-가이드)
    4. [자주 묻는 질문](#자주-묻는-질문)
    5. [문제 해결](#문제-해결)
    """)

    # SMARTS 문법
    st.markdown("---")
    st.subheader("📚 SMARTS 문법 기초")

    tab1, tab2, tab3, tab4 = st.tabs(["원자", "결합", "논리 연산", "고급"])

    with tab1:
        st.markdown("""
        ### 원자 표현

        | 표현 | 의미 | 예시 |
        |------|------|------|
        | `C` | 지방족 탄소 | 메탄, 에탄 |
        | `c` | 방향족 탄소 | 벤젠, 톨루엔 |
        | `N` | 질소 | 아민, 아미드 |
        | `O` | 산소 | 알코올, 에테르 |
        | `*` | 임의 원자 | 모든 원자 |
        | `[C]` | 명시적 탄소 | 특정 조건의 탄소 |
        | `[!C]` | 탄소가 아닌 원자 | 헤테로 원자 |
        | `[C,N]` | 탄소 또는 질소 | OR 조건 |
        """)

    with tab2:
        st.markdown("""
        ### 결합 표현

        | 표현 | 의미 | 예시 |
        |------|------|------|
        | `-` | 단일 결합 | C-C |
        | `=` | 이중 결합 | C=O |
        | `#` | 삼중 결합 | C#N |
        | `:` | 방향족 결합 | 벤젠 고리 |
        | `~` | 임의 결합 | 모든 결합 |
        | `@` | 고리 내 결합 | 고리 구조 |
        | `!@` | 고리 외 결합 | 선형 구조 |
        """)

    with tab3:
        st.markdown("""
        ### 논리 연산

        | 연산 | 의미 | 예시 |
        |------|------|------|
        | `,` | OR | `[C,N]` (탄소 또는 질소) |
        | `&` | AND | `[C&R]` (고리 내 탄소) |
        | `!` | NOT | `[!C]` (탄소가 아님) |
        | `;` | 낮은 우선순위 AND | 복잡한 조건 |
        """)

    with tab4:
        st.markdown("""
        ### 고급 표현

        | 표현 | 의미 | 예시 |
        |------|------|------|
        | `[H]` | 수소 개수 | `[CH3]` (메틸) |
        | `[+]` | 양전하 | `[N+]` (양이온) |
        | `[-]` | 음전하 | `[O-]` (음이온) |
        | `[R]` | 고리 내 원자 | 고리 구조 |
        | `[R0]` | 고리 외 원자 | 선형 구조 |
        | `[D2]` | 결합 차수 2 | 2개 결합 |
        | `[v4]` | 원자가 4 | 4가 원자 |
        """)

    # 분자 그리기 가이드
    st.markdown("---")
    st.subheader("🎨 분자 그리기 가이드")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        ### 기본 조작

        1. **원자 추가**
           - 원자 모드 선택
           - 원소 선택 (C, N, O 등)
           - 캔버스 클릭

        2. **결합 생성**
           - 결합 모드 선택
           - 결합 타입 선택
           - 첫 번째 원자 클릭
           - 두 번째 원자 클릭

        3. **원자 이동**
           - 선택 모드 활성화
           - 원자를 드래그
        """)

    with col2:
        st.markdown("""
        ### 고급 기능

        1. **SMARTS 옵션**
           - 방향족/지방족 설정
           - 전하 설정
           - 수소 개수 설정
           - 부정 조건 설정

        2. **미리 정의된 구조**
           - 벤젠, 사이클로헥산 등
           - 원클릭으로 추가

        3. **파일 관리**
           - JSON 형식으로 저장
           - 구조 불러오기
        """)

    # 패턴 매칭 가이드
    st.markdown("---")
    st.subheader("🔍 패턴 매칭 가이드")

    st.markdown("""
    ### 매칭 과정

    1. **패턴 준비**
       - SMARTS 패턴 입력 또는 생성
       - 패턴 유효성 확인

    2. **분자 준비**
       - 기본 분자에서 선택
       - SMILES 직접 입력
       - 파일에서 불러오기

    3. **매칭 실행**
       - 옵션 설정 (입체화학, 방향족성 등)
       - 매칭 실행
       - 결과 확인

    4. **결과 분석**
       - 매칭된 분자 확인
       - 하이라이트 표시
       - 통계 정보 확인
    """)

    # FAQ
    st.markdown("---")
    st.subheader("❓ 자주 묻는 질문")

    with st.expander("Q: SMARTS와 SMILES의 차이점은?"):
        st.markdown("""
        **SMILES**는 특정 분자를 나타내는 문자열이고,
        **SMARTS**는 분자 패턴을 나타내는 문자열입니다.

        - SMILES: `CCO` (에탄올)
        - SMARTS: `C-C-O` (알코올 패턴)
        """)

    with st.expander("Q: 방향족 원자는 어떻게 표현하나요?"):
        st.markdown("""
        소문자를 사용합니다:
        - `c`: 방향족 탄소
        - `n`: 방향족 질소
        - `o`: 방향족 산소

        예시: `c1ccccc1` (벤젠)
        """)

    with st.expander("Q: 고리 구조는 어떻게 그리나요?"):
        st.markdown("""
        1. 고리의 원자들을 순서대로 그립니다
        2. 마지막 원자와 첫 번째 원자를 연결합니다
        3. 자동으로 고리 번호가 할당됩니다

        또는 미리 정의된 구조를 사용하세요.
        """)

    with st.expander("Q: 패턴이 매칭되지 않는 이유는?"):
        st.markdown("""
        가능한 원인:
        1. 패턴 문법 오류
        2. 방향족성 불일치
        3. 입체화학 차이
        4. 전하 상태 차이

        해결 방법:
        1. 패턴 유효성 확인
        2. 매칭 옵션 조정
        3. 더 일반적인 패턴 사용
        """)

    # 문제 해결
    st.markdown("---")
    st.subheader("🔧 문제 해결")

    st.markdown("""
    ### 일반적인 문제들

    **1. 분자가 그려지지 않을 때**
    - 브라우저 새로고침
    - 다른 브라우저 시도
    - 캐시 삭제

    **2. SMARTS 생성이 안 될 때**
    - 원자가 연결되어 있는지 확인
    - 유효한 분자 구조인지 확인
    - 세션 초기화 후 재시도

    **3. 매칭 결과가 예상과 다를 때**
    - 패턴 문법 재확인
    - 매칭 옵션 조정
    - 더 구체적이거나 일반적인 패턴 시도

    **4. 성능이 느릴 때**
    - 분자 수 줄이기
    - 패턴 단순화
    - 배치 크기 조정
    """)

    # 연락처
    st.markdown("---")
    st.subheader("📞 지원")

    st.markdown("""
    추가 도움이 필요하시면:

    - 📧 이메일: support@example.com
    - 🐛 버그 리포트: GitHub Issues
    - 💡 기능 제안: GitHub Discussions
    - 📖 문서: 온라인 매뉴얼
    """)

def main():
    """메인 함수"""
    # 사이드바 네비게이션
    page = create_sidebar()

    # 페이지별 라우팅
    if page == "🏠 홈":
        create_home_page()

    elif page == "✏️ 분자 그리기":
        st.title("✏️ 분자 그리기")
        st.markdown("---")

        drawer = MoleculeDrawer()
        drawer.create_drawing_interface()

    elif page == "🧬 SMARTS 생성":
        st.title("🧬 SMARTS 생성")
        st.markdown("---")

        # 현재 분자 로드
        drawer = MoleculeDrawer()
        drawer.initialize_session_state()
        molecule = drawer.get_molecule()

        if not molecule.atoms:
            st.warning("⚠️ 먼저 분자를 그려주세요!")
            st.info("분자 그리기 페이지로 이동하여 분자를 그린 후 다시 시도하세요.")
        else:
            generator = SmartsGenerator()
            smarts = generator.create_smarts_interface(molecule)

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
            🧪 SMARTS 패턴 편집기 v1.0.0 | 
            Made with ❤️ using Streamlit & RDKit
        </div>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()