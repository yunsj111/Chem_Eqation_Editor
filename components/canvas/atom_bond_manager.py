"""
원자와 결합 관리 기능
"""

def get_atom_bond_manager():
    """원자와 결합 관리 JavaScript 함수들 반환"""
    return """
        // 원자 생성
        function createAtom(x, y, element) {
            const atomId = 'atom_' + nextAtomId++;
            atoms[atomId] = {
                id: atomId,
                x: x,
                y: y,
                element: element,
                charge: 0,
                bonds: []
            };
            
            redrawCanvas();
            return atomId;
        }
        
        // 원자 원소 변경
        function changeAtomElement(atomId, newElement) {
            if (atoms[atomId]) {
                atoms[atomId].element = newElement;
                redrawCanvas();
            }
        }
        
        // 원자 삭제
        function deleteAtom(atomId) {
            if (!atoms[atomId]) return;
            
            // 연결된 결합들 삭제
            const connectedBonds = atoms[atomId].bonds.slice();
            connectedBonds.forEach(bondId => deleteBond(bondId));
            
            delete atoms[atomId];
            selectedAtoms.delete(atomId);
            
            redrawCanvas();
            updateStatusBar();
        }
        
        // 결합 생성
        function createBond(atom1Id, atom2Id, bondType = currentBondType) {
            if (!atoms[atom1Id] || !atoms[atom2Id] || atom1Id === atom2Id) return null;
            
            // 이미 존재하는 결합 확인
            for (const [bondId, bond] of Object.entries(bonds)) {
                if ((bond.atom1 === atom1Id && bond.atom2 === atom2Id) ||
                    (bond.atom1 === atom2Id && bond.atom2 === atom1Id)) {
                    // 기존 결합의 타입 변경
                    bond.type = bondType;
                    redrawCanvas();
                    return bondId;
                }
            }
            
            const bondId = 'bond_' + nextBondId++;
            bonds[bondId] = {
                id: bondId,
                atom1: atom1Id,
                atom2: atom2Id,
                type: bondType
            };
            
            // 원자에 결합 정보 추가
            atoms[atom1Id].bonds.push(bondId);
            atoms[atom2Id].bonds.push(bondId);
            
            redrawCanvas();
            return bondId;
        }
        
        // 결합 삭제
        function deleteBond(bondId) {
            if (!bonds[bondId]) return;
            
            const bond = bonds[bondId];
            
            // 원자에서 결합 정보 제거
            if (atoms[bond.atom1]) {
                atoms[bond.atom1].bonds = atoms[bond.atom1].bonds.filter(id => id !== bondId);
            }
            if (atoms[bond.atom2]) {
                atoms[bond.atom2].bonds = atoms[bond.atom2].bonds.filter(id => id !== bondId);
            }
            
            delete bonds[bondId];
            selectedBonds.delete(bondId);
            
            redrawCanvas();
            updateStatusBar();
        }
        
        // 결합 모드 처리
        function handleBondMode(atomId) {
            if (!firstBondAtom) {
                firstBondAtom = atomId;
                highlightAtom(atomId);
            } else if (firstBondAtom !== atomId) {
                createBond(firstBondAtom, atomId);
                unhighlightAtom(firstBondAtom);
                firstBondAtom = null;
            } else {
                unhighlightAtom(firstBondAtom);
                firstBondAtom = null;
            }
        }
        
        // 원자 하이라이트
        function highlightAtom(atomId) {
            const atomElement = document.getElementById(atomId);
            if (atomElement) {
                atomElement.style.stroke = '#007bff';
                atomElement.style.strokeWidth = '3';
            }
        }
        
        // 원자 하이라이트 해제
        function unhighlightAtom(atomId) {
            const atomElement = document.getElementById(atomId);
            if (atomElement) {
                atomElement.style.stroke = getElementColor(atoms[atomId].element);
                atomElement.style.strokeWidth = '2';
            }
        }
        
        // 결합 드래그 시작 (개선된 버전)
        function startBondDrag(x, y) {
            // 클릭 지점에서 가장 가까운 원자 찾기
            let closestAtom = null;
            let minDistance = Infinity;
            
            for (const [atomId, atom] of Object.entries(atoms)) {
                const distance = Math.sqrt((x - atom.x) ** 2 + (y - atom.y) ** 2);
                if (distance < minDistance && distance <= 50) { // 범위를 늘림
                    minDistance = distance;
                    closestAtom = atomId;
                }
            }
            
            if (closestAtom) {
                isBondDragging = true;
                bondSourceAtom = closestAtom;
                bondDragStart = { x: atoms[closestAtom].x, y: atoms[closestAtom].y };
                
                // 소스 원자 하이라이트
                highlightAtom(closestAtom);
                
                // 드래그 라인 생성
                bondDragLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                bondDragLine.setAttribute('class', 'bond-drag-line');
                bondDragLine.setAttribute('x1', atoms[closestAtom].x);
                bondDragLine.setAttribute('y1', atoms[closestAtom].y);
                bondDragLine.setAttribute('x2', atoms[closestAtom].x);
                bondDragLine.setAttribute('y2', atoms[closestAtom].y);
                canvas.appendChild(bondDragLine);
                
                // 상태바 업데이트
                document.getElementById('statusBar').textContent = 
                    `결합 생성 중... 대상 원자로 드래그하세요. (${currentBondType} 결합)`;
            }
        }
        
        // 결합 드래그 업데이트
        function updateBondDrag(x, y) {
            if (bondDragLine) {
                bondDragLine.setAttribute('x2', x);
                bondDragLine.setAttribute('y2', y);
            }
        }
        
        // 결합 드래그 완료 (개선된 버전)
        function finishBondDrag(x, y) {
            if (!isBondDragging) return;
            
            // 드래그 라인 제거
            if (bondDragLine) {
                canvas.removeChild(bondDragLine);
                bondDragLine = null;
            }
            
            // 소스 원자 하이라이트 해제
            if (bondSourceAtom) {
                unhighlightAtom(bondSourceAtom);
            }
            
            // 대상 원자 찾기
            const targetAtom = getAtomAt(x, y);
            
            if (targetAtom && targetAtom !== bondSourceAtom) {
                // 결합 생성
                const newBondId = createBond(bondSourceAtom, targetAtom);
                if (newBondId) {
                    document.getElementById('statusBar').textContent = 
                        `${currentBondType} 결합이 생성되었습니다.`;
                }
            } else if (targetAtom === bondSourceAtom) {
                document.getElementById('statusBar').textContent = 
                    '동일한 원자에는 결합을 생성할 수 없습니다.';
            } else {
                // 빈 공간에 드롭한 경우 - 새 원자 생성 후 결합
                if (Math.sqrt((x - atoms[bondSourceAtom].x) ** 2 + (y - atoms[bondSourceAtom].y) ** 2) > 20) {
                    const newAtomId = createAtom(x, y, currentElement);
                    createBond(bondSourceAtom, newAtomId);
                    document.getElementById('statusBar').textContent = 
                        `새 원자(${currentElement})와 ${currentBondType} 결합이 생성되었습니다.`;
                }
            }
            
            // 상태 초기화
            isBondDragging = false;
            bondSourceAtom = null;
            bondDragStart = null;
            
            // 상태바 원상복구
            setTimeout(() => {
                updateStatusBar();
            }, 2000);
        }
        
        // 원자 이동
        function moveAtom(atomId, newX, newY) {
            if (atoms[atomId]) {
                atoms[atomId].x = newX;
                atoms[atomId].y = newY;
            }
        }
        
        // 여러 원자 이동
        function moveAtoms(atomIds, deltaX, deltaY) {
            atomIds.forEach(atomId => {
                if (atoms[atomId]) {
                    atoms[atomId].x += deltaX;
                    atoms[atomId].y += deltaY;
                }
            });
        }
        
        // 원자와 연결된 결합들 업데이트
        function updateConnectedBonds(atomId) {
            if (!atoms[atomId]) return;
            
            atoms[atomId].bonds.forEach(bondId => {
                if (bonds[bondId]) {
                    // 결합 시각적 업데이트는 redrawCanvas에서 처리
                }
            });
        }
    """
