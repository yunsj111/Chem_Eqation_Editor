"""
화학 구조 에디터 JavaScript 코어 기능
"""

def get_canvas_core_js():
    """캔버스 코어 JavaScript 기능 반환"""
    return """
        // 전역 변수들
        let currentMode = 'atom';
        let currentElement = 'C';
        let currentBondType = 'single';
        let atoms = {};
        let bonds = {};
        let selectedAtoms = new Set();
        let selectedBonds = new Set();
        let nextAtomId = 1;
        let nextBondId = 1;
        let canvas = null;
        
        // 선택 관련 변수
        let isSelecting = false;
        let selectionStart = { x: 0, y: 0 };
        let selectionRect = null;
        
        // 드래그 관련 변수
        let isDragging = false;
        let dragStart = { x: 0, y: 0 };
        let draggedAtoms = new Set();
        
        // 결합 드래그 관련 변수
        let isBondDragging = false;
        let bondDragStart = null;
        let bondDragLine = null;
        let bondSourceAtom = null;
        
        // 첫 번째 결합 원자 (결합 모드용)
        let firstBondAtom = null;
        
        // 초기화 함수
        function initCanvas() {
            canvas = document.getElementById('moleculeCanvas');
            if (!canvas) return;
            
            setupEventListeners();
            updateStatusBar();
            
            // 주기율표 초기화 및 표시
            generatePeriodicTable();
            document.getElementById('periodicTable').classList.add('show');
            setCurrentElement('C'); // 기본 원소 설정
        }
        
        // 이벤트 리스너 설정
        function setupEventListeners() {
            canvas.addEventListener('mousedown', handleMouseDown);
            canvas.addEventListener('mousemove', handleMouseMove);
            canvas.addEventListener('mouseup', handleMouseUp);
            canvas.addEventListener('contextmenu', e => e.preventDefault());
            
            // 클릭 외부에서 드롭다운 닫기
            document.addEventListener('click', function(event) {
                if (!event.target.matches('.tool-btn')) {
                    closeAllDropdowns();
                }
            });
        }
        
        // 모든 드롭다운 닫기
        function closeAllDropdowns() {
            document.querySelectorAll('.smarts-content').forEach(dropdown => {
                dropdown.classList.remove('show');
            });
            
            if (!document.getElementById('periodicTable').contains(event.target)) {
                document.getElementById('periodicTable').classList.remove('show');
            }
        }
        
        // 마우스 이벤트 핸들러들
        function handleMouseDown(event) {
            event.preventDefault();
            const rect = canvas.getBoundingClientRect();
            const x = event.clientX - rect.left;
            const y = event.clientY - rect.top;
            
            const clickedAtom = getAtomAt(x, y);
            
            if (currentMode === 'atom') {
                if (clickedAtom) {
                    changeAtomElement(clickedAtom, currentElement);
                } else {
                    createAtom(x, y, currentElement);
                }
            } else if (currentMode === 'bond') {
                if (clickedAtom) {
                    // 원자를 클릭한 경우 - 기존 결합 모드 유지
                    handleBondMode(clickedAtom);
                } else {
                    // 빈 공간을 클릭한 경우 - 드래그 모드 시작
                    startBondDrag(x, y);
                }
            } else if (currentMode === 'select') {
                if (clickedAtom) {
                    handleAtomSelection(clickedAtom, event.ctrlKey || event.metaKey);
                    if (selectedAtoms.has(clickedAtom)) {
                        startDrag(x, y);
                    }
                } else {
                    if (!event.ctrlKey && !event.metaKey) {
                        clearSelection();
                    }
                    startAreaSelection(x, y);
                }
            } else if (currentMode === 'delete') {
                if (clickedAtom) {
                    deleteAtom(clickedAtom);
                } else {
                    const clickedBond = getBondAt(x, y);
                    if (clickedBond) {
                        deleteBond(clickedBond);
                    }
                }
            }
        }
        
        function handleMouseMove(event) {
            const rect = canvas.getBoundingClientRect();
            const x = event.clientX - rect.left;
            const y = event.clientY - rect.top;
            
            if (isSelecting) {
                updateAreaSelection(x, y);
            } else if (isDragging) {
                updateDrag(x, y);
            } else if (isBondDragging) {
                updateBondDrag(x, y);
            }
        }
        
        function handleMouseUp(event) {
            const rect = canvas.getBoundingClientRect();
            const x = event.clientX - rect.left;
            const y = event.clientY - rect.top;
            
            if (isSelecting) {
                finishAreaSelection();
            } else if (isDragging) {
                finishDrag();
            } else if (isBondDragging) {
                finishBondDrag(x, y);
            }
        }
        
        // 모드 설정 함수들
        function setMode(mode) {
            currentMode = mode;
            firstBondAtom = null;
            clearSelection();
            
            // 버튼 스타일 업데이트
            document.querySelectorAll('.tool-btn').forEach(btn => {
                btn.classList.remove('active');
            });
            document.getElementById(mode + 'Btn').classList.add('active');
            
            updateStatusBar();
        }
        
        function setCurrentElement(element) {
            currentElement = element;
            
            // 선택된 원소 표시 업데이트
            const selectedElementDisplay = document.getElementById('selectedElement');
            if (selectedElementDisplay) {
                selectedElementDisplay.textContent = element;
                selectedElementDisplay.style.backgroundColor = getElementColor(element);
            }
            
            updateStatusBar();
        }
        
        function setBondType(bondType) {
            currentBondType = bondType;
            
            // 버튼 스타일 업데이트
            document.querySelectorAll('#singleBtn, #doubleBtn, #tripleBtn, #aromaticBtn').forEach(btn => {
                btn.classList.remove('active');
            });
            document.getElementById(bondType + 'Btn').classList.add('active');
            
            updateStatusBar();
        }
        
        // 상태바 업데이트
        function updateStatusBar() {
            const selectedCount = selectedAtoms.size + selectedBonds.size;
            const selectionText = selectedCount > 0 ? `${selectedCount}개 선택됨` : '없음';
            
            document.getElementById('statusBar').textContent = 
                `모드: ${currentMode} | 원소: ${currentElement} | 결합: ${currentBondType} | 선택된 항목: ${selectionText}`;
        }
        
        // 좌표에서 원자 찾기
        function getAtomAt(x, y) {
            for (const [atomId, atom] of Object.entries(atoms)) {
                const distance = Math.sqrt((x - atom.x) ** 2 + (y - atom.y) ** 2);
                const radius = getElementRadius(atom.element);
                if (distance <= radius + 5) {
                    return atomId;
                }
            }
            return null;
        }
        
        // 좌표에서 결합 찾기
        function getBondAt(x, y) {
            for (const [bondId, bond] of Object.entries(bonds)) {
                const atom1 = atoms[bond.atom1];
                const atom2 = atoms[bond.atom2];
                if (!atom1 || !atom2) continue;
                
                const distance = distanceToLine(x, y, atom1.x, atom1.y, atom2.x, atom2.y);
                if (distance <= 5) {
                    return bondId;
                }
            }
            return null;
        }
        
        // 점과 선분 사이의 거리 계산
        function distanceToLine(px, py, x1, y1, x2, y2) {
            const A = px - x1;
            const B = py - y1;
            const C = x2 - x1;
            const D = y2 - y1;
            
            const dot = A * C + B * D;
            const lenSq = C * C + D * D;
            
            if (lenSq === 0) return Math.sqrt(A * A + B * B);
            
            let param = dot / lenSq;
            param = Math.max(0, Math.min(1, param));
            
            const xx = x1 + param * C;
            const yy = y1 + param * D;
            
            const dx = px - xx;
            const dy = py - yy;
            
            return Math.sqrt(dx * dx + dy * dy);
        }
        
        // 유틸리티 함수들
        function clearCanvas() {
            atoms = {};
            bonds = {};
            selectedAtoms.clear();
            selectedBonds.clear();
            nextAtomId = 1;
            nextBondId = 1;
            firstBondAtom = null;
            redrawCanvas();
            updateStatusBar();
        }
        
        // 초기화 실행
        document.addEventListener('DOMContentLoaded', initCanvas);
    """
