"""
선택과 드래그 기능 관리
"""

def get_selection_drag_manager():
    """선택과 드래그 JavaScript 함수들 반환"""
    return """
        // 원자 선택 처리
        function handleAtomSelection(atomId, multiSelect = false) {
            if (!multiSelect) {
                clearSelection();
            }
            
            if (selectedAtoms.has(atomId)) {
                selectedAtoms.delete(atomId);
            } else {
                selectedAtoms.add(atomId);
            }
            
            updateSelectionVisuals();
            updateStatusBar();
        }
        
        // 선택 영역 시작
        function startAreaSelection(x, y) {
            isSelecting = true;
            selectionStart = { x, y };
            
            // 선택 사각형 생성
            selectionRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            selectionRect.setAttribute('class', 'selection-area');
            selectionRect.setAttribute('x', x);
            selectionRect.setAttribute('y', y);
            selectionRect.setAttribute('width', 0);
            selectionRect.setAttribute('height', 0);
            canvas.appendChild(selectionRect);
        }
        
        // 선택 영역 업데이트
        function updateAreaSelection(x, y) {
            if (!selectionRect) return;
            
            const width = Math.abs(x - selectionStart.x);
            const height = Math.abs(y - selectionStart.y);
            const rectX = Math.min(x, selectionStart.x);
            const rectY = Math.min(y, selectionStart.y);
            
            selectionRect.setAttribute('x', rectX);
            selectionRect.setAttribute('y', rectY);
            selectionRect.setAttribute('width', width);
            selectionRect.setAttribute('height', height);
        }
        
        // 선택 영역 완료
        function finishAreaSelection() {
            if (!selectionRect) return;
            
            const rect = selectionRect.getBoundingClientRect();
            const canvasRect = canvas.getBoundingClientRect();
            
            // 캔버스 좌표로 변환
            const selectionBounds = {
                left: parseFloat(selectionRect.getAttribute('x')),
                top: parseFloat(selectionRect.getAttribute('y')),
                right: parseFloat(selectionRect.getAttribute('x')) + parseFloat(selectionRect.getAttribute('width')),
                bottom: parseFloat(selectionRect.getAttribute('y')) + parseFloat(selectionRect.getAttribute('height'))
            };
            
            // 영역 내의 원자들 선택
            for (const [atomId, atom] of Object.entries(atoms)) {
                if (atom.x >= selectionBounds.left && atom.x <= selectionBounds.right &&
                    atom.y >= selectionBounds.top && atom.y <= selectionBounds.bottom) {
                    selectedAtoms.add(atomId);
                }
            }
            
            // 선택 사각형 제거
            canvas.removeChild(selectionRect);
            selectionRect = null;
            isSelecting = false;
            
            updateSelectionVisuals();
            updateStatusBar();
        }
        
        // 드래그 시작
        function startDrag(x, y) {
            isDragging = true;
            dragStart = { x, y };
            draggedAtoms = new Set(selectedAtoms);
        }
        
        // 드래그 업데이트
        function updateDrag(x, y) {
            if (!isDragging) return;
            
            const deltaX = x - dragStart.x;
            const deltaY = y - dragStart.y;
            
            // 선택된 원자들 이동
            draggedAtoms.forEach(atomId => {
                if (atoms[atomId]) {
                    const originalX = atoms[atomId].x - deltaX;
                    const originalY = atoms[atomId].y - deltaY;
                    atoms[atomId].x = originalX + deltaX;
                    atoms[atomId].y = originalY + deltaY;
                }
            });
            
            redrawCanvas();
        }
        
        // 드래그 완료
        function finishDrag() {
            if (!isDragging) return;
            
            isDragging = false;
            draggedAtoms.clear();
            
            // 연결된 결합들 업데이트
            selectedAtoms.forEach(atomId => {
                updateConnectedBonds(atomId);
            });
            
            redrawCanvas();
        }
        
        // 선택 해제
        function clearSelection() {
            selectedAtoms.clear();
            selectedBonds.clear();
            updateSelectionVisuals();
            updateStatusBar();
        }
        
        // 선택 시각화 업데이트
        function updateSelectionVisuals() {
            // 모든 원자의 선택 스타일 제거
            document.querySelectorAll('.selected-atom').forEach(element => {
                element.classList.remove('selected-atom');
            });
            
            // 모든 결합의 선택 스타일 제거
            document.querySelectorAll('.selected-bond').forEach(element => {
                element.classList.remove('selected-bond');
            });
            
            // 선택된 원자들에 스타일 적용
            selectedAtoms.forEach(atomId => {
                const atomElement = document.getElementById(atomId);
                if (atomElement) {
                    atomElement.classList.add('selected-atom');
                }
            });
            
            // 선택된 결합들에 스타일 적용
            selectedBonds.forEach(bondId => {
                const bondElement = document.getElementById(bondId);
                if (bondElement) {
                    bondElement.classList.add('selected-bond');
                }
            });
        }
        
        // 키보드 이벤트 처리
        function setupKeyboardEvents() {
            document.addEventListener('keydown', function(event) {
                if (event.key === 'Delete' || event.key === 'Backspace') {
                    deleteSelected();
                } else if (event.ctrlKey || event.metaKey) {
                    if (event.key === 'a' || event.key === 'A') {
                        event.preventDefault();
                        selectAll();
                    } else if (event.key === 'c' || event.key === 'C') {
                        event.preventDefault();
                        copySelected();
                    } else if (event.key === 'v' || event.key === 'V') {
                        event.preventDefault();
                        pasteSelected();
                    }
                }
            });
        }
        
        // 선택된 항목들 삭제
        function deleteSelected() {
            // 선택된 결합 삭제
            selectedBonds.forEach(bondId => deleteBond(bondId));
            
            // 선택된 원자 삭제
            selectedAtoms.forEach(atomId => deleteAtom(atomId));
            
            clearSelection();
        }
        
        // 모든 항목 선택
        function selectAll() {
            selectedAtoms = new Set(Object.keys(atoms));
            selectedBonds = new Set(Object.keys(bonds));
            updateSelectionVisuals();
            updateStatusBar();
        }
        
        // 복사 기능 (향후 구현)
        function copySelected() {
            // TODO: 선택된 원자와 결합을 클립보드에 복사
            console.log('복사 기능은 아직 구현되지 않았습니다.');
        }
        
        // 붙여넣기 기능 (향후 구현)
        function pasteSelected() {
            // TODO: 클립보드에서 원자와 결합을 붙여넣기
            console.log('붙여넣기 기능은 아직 구현되지 않았습니다.');
        }
        
        // 키보드 이벤트 초기화
        document.addEventListener('DOMContentLoaded', setupKeyboardEvents);
    """
