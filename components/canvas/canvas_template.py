"""
화학 구조 에디터 HTML 템플릿 정의
"""

def get_canvas_template():
    """캔버스 HTML 템플릿 반환"""
    return """
        <div class="toolbar">
            <div class="tool-group">
                <label>모드:</label>
                <button class="tool-btn active" id="atomBtn" onclick="setMode('atom')">원자</button>
                <button class="tool-btn" id="bondBtn" onclick="setMode('bond')">결합</button>
                <button class="tool-btn" id="selectBtn" onclick="setMode('select')">선택</button>
                <button class="tool-btn" id="deleteBtn" onclick="setMode('delete')">삭제</button>
            </div>
            
            <div class="tool-group">
                <label>원소:</label>
                <div class="selected-element-display">
                    <span>선택된 원소: </span>
                    <span id="selectedElement" class="element-badge">C</span>
                </div>
                <button class="tool-btn" id="periodicBtn" onclick="togglePeriodicTable()">주기율표</button>
            </div>
            
            <div class="tool-group">
                <label>결합:</label>
                <button class="tool-btn active" id="singleBtn" onclick="setBondType('single')">단일</button>
                <button class="tool-btn" id="doubleBtn" onclick="setBondType('double')">이중</button>
                <button class="tool-btn" id="tripleBtn" onclick="setBondType('triple')">삼중</button>
                <button class="tool-btn" id="aromaticBtn" onclick="setBondType('aromatic')">방향족</button>
            </div>
            
            <div class="tool-group">
                <div class="smarts-dropdown">
                    <button class="tool-btn" onclick="toggleSmartsDropdown()">SMARTS</button>
                    <div class="smarts-content" id="smartsDropdown">
                        <div class="smarts-item" onclick="setSmartsPattern('[!H]')">[!H] - 수소가 아닌 원자</div>
                        <div class="smarts-item" onclick="setSmartsPattern('[C&H3]')">[C&H3] - 3개 수소를 가진 탄소</div>
                        <div class="smarts-item" onclick="setSmartsPattern('[N+]')">[N+] - 양전하 질소</div>
                        <div class="smarts-item" onclick="setSmartsPattern('[O-]')">[O-] - 음전하 산소</div>
                        <div class="smarts-item" onclick="setSmartsPattern('[R]')">[R] - 고리 원자</div>
                        <div class="smarts-item" onclick="setSmartsPattern('[!R]')">[!R] - 비고리 원자</div>
                        <div class="smarts-item" onclick="setSmartsPattern('*')">* - 모든 원자</div>
                    </div>
                </div>
                <input type="text" class="smarts-input" id="customSmarts" placeholder="사용자 정의 SMARTS" 
                       onkeypress="if(event.key==='Enter') setSmartsPattern(this.value)">
            </div>
            
            <div class="tool-group">
                <button class="tool-btn" onclick="clearCanvas()">전체삭제</button>
                <button class="tool-btn" onclick="exportMolecule()">내보내기</button>
                <button class="tool-btn" onclick="loadTestMolecule()">테스트</button>
            </div>
        </div>
        
        <div class="main-content">
            <div class="canvas-container">
                <svg id="moleculeCanvas" width="800" height="600"></svg>
            </div>
            
            <div class="periodic-table" id="periodicTable">
                <div class="periodic-header">
                    <h3>주기율표</h3>
                    <button class="close-btn" onclick="togglePeriodicTable()">×</button>
                </div>
                <div class="periodic-grid" id="periodicGrid">
                    <!-- 주기율표 원소들이 JavaScript로 생성됩니다 -->
                </div>
            </div>
        </div>
        
        <div class="status-bar" id="statusBar">
            모드: 원자 | 원소: C | 결합: 단일 | 선택된 항목: 없음
        </div>
    """
