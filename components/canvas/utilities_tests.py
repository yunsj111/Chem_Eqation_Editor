"""
유틸리티 및 테스트 기능
"""

def get_utilities_and_tests():
    """유틸리티 및 테스트 JavaScript 함수들 반환"""
    return """
        // 분자 내보내기
        function exportMolecule() {
            const moleculeData = {
                atoms: atoms,
                bonds: bonds,
                metadata: {
                    createdAt: new Date().toISOString(),
                    atomCount: Object.keys(atoms).length,
                    bondCount: Object.keys(bonds).length
                }
            };
            
            const jsonString = JSON.stringify(moleculeData, null, 2);
            
            // 다운로드 링크 생성
            const blob = new Blob([jsonString], { type: 'application/json' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'molecule.json';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
            
            console.log('분자 데이터:', moleculeData);
        }
        
        // 분자 가져오기
        function importMolecule(moleculeData) {
            try {
                atoms = moleculeData.atoms || {};
                bonds = moleculeData.bonds || {};
                
                // ID 카운터 업데이트
                nextAtomId = Math.max(...Object.keys(atoms).map(id => parseInt(id.split('_')[1]) + 1), 1);
                nextBondId = Math.max(...Object.keys(bonds).map(id => parseInt(id.split('_')[1]) + 1), 1);
                
                clearSelection();
                redrawCanvas();
                updateStatusBar();
                
                console.log('분자를 성공적으로 가져왔습니다.');
            } catch (error) {
                console.error('분자 가져오기 오류:', error);
                alert('분자 데이터를 가져오는 중 오류가 발생했습니다.');
            }
        }
        
        // 테스트 분자 로드
        function loadTestMolecule() {
            clearCanvas();
            
            // 다양한 원소와 결합 타입을 포함한 테스트 분자
            const testAtoms = {
                'atom_1': { id: 'atom_1', x: 100, y: 200, element: 'C', charge: 0, bonds: ['bond_1', 'bond_2'] },
                'atom_2': { id: 'atom_2', x: 150, y: 150, element: 'N', charge: 0, bonds: ['bond_1', 'bond_3'] },
                'atom_3': { id: 'atom_3', x: 200, y: 200, element: 'O', charge: 0, bonds: ['bond_2', 'bond_4'] },
                'atom_4': { id: 'atom_4', x: 250, y: 150, element: 'S', charge: 0, bonds: ['bond_3', 'bond_5'] },
                'atom_5': { id: 'atom_5', x: 300, y: 200, element: 'P', charge: 0, bonds: ['bond_4', 'bond_6'] },
                'atom_6': { id: 'atom_6', x: 350, y: 150, element: 'H', charge: 0, bonds: ['bond_5'] },
                'atom_7': { id: 'atom_7', x: 175, y: 250, element: 'F', charge: 0, bonds: ['bond_6'] },
                'atom_8': { id: 'atom_8', x: 225, y: 100, element: 'Cl', charge: 0, bonds: [] },
                'atom_9': { id: 'atom_9', x: 275, y: 250, element: '[!H]', charge: 0, bonds: [] },
                'atom_10': { id: 'atom_10', x: 325, y: 100, element: '*', charge: 0, bonds: [] }
            };
            
            const testBonds = {
                'bond_1': { id: 'bond_1', atom1: 'atom_1', atom2: 'atom_2', type: 'single' },
                'bond_2': { id: 'bond_2', atom1: 'atom_1', atom2: 'atom_3', type: 'double' },
                'bond_3': { id: 'bond_3', atom1: 'atom_2', atom2: 'atom_4', type: 'triple' },
                'bond_4': { id: 'bond_4', atom1: 'atom_3', atom2: 'atom_5', type: 'aromatic' },
                'bond_5': { id: 'bond_5', atom1: 'atom_4', atom2: 'atom_6', type: 'single' },
                'bond_6': { id: 'bond_6', atom1: 'atom_5', atom2: 'atom_7', type: 'double' }
            };
            
            atoms = testAtoms;
            bonds = testBonds;
            nextAtomId = 11;
            nextBondId = 7;
            
            redrawCanvas();
            updateStatusBar();
        }
        
        // 분자 검증
        function validateMolecule() {
            const issues = [];
            
            // 원자 검증
            for (const [atomId, atom] of Object.entries(atoms)) {
                if (!atom.x || !atom.y || !atom.element) {
                    issues.push(`원자 ${atomId}에 필수 속성이 누락되었습니다.`);
                }
                
                // 결합 참조 검증
                atom.bonds.forEach(bondId => {
                    if (!bonds[bondId]) {
                        issues.push(`원자 ${atomId}가 존재하지 않는 결합 ${bondId}를 참조합니다.`);
                    }
                });
            }
            
            // 결합 검증
            for (const [bondId, bond] of Object.entries(bonds)) {
                if (!atoms[bond.atom1] || !atoms[bond.atom2]) {
                    issues.push(`결합 ${bondId}가 존재하지 않는 원자를 참조합니다.`);
                }
                
                if (bond.atom1 === bond.atom2) {
                    issues.push(`결합 ${bondId}가 동일한 원자를 연결합니다.`);
                }
            }
            
            return issues;
        }
        
        // 분자 통계
        function getMoleculeStats() {
            const elementCounts = {};
            const bondTypeCounts = {};
            
            // 원소별 개수
            for (const atom of Object.values(atoms)) {
                elementCounts[atom.element] = (elementCounts[atom.element] || 0) + 1;
            }
            
            // 결합 타입별 개수
            for (const bond of Object.values(bonds)) {
                bondTypeCounts[bond.type] = (bondTypeCounts[bond.type] || 0) + 1;
            }
            
            return {
                totalAtoms: Object.keys(atoms).length,
                totalBonds: Object.keys(bonds).length,
                elementCounts,
                bondTypeCounts,
                selectedAtoms: selectedAtoms.size,
                selectedBonds: selectedBonds.size
            };
        }
        
        // 화면 중앙으로 분자 이동
        function centerMolecule() {
            if (Object.keys(atoms).length === 0) return;
            
            const atomPositions = Object.values(atoms);
            const minX = Math.min(...atomPositions.map(a => a.x));
            const maxX = Math.max(...atomPositions.map(a => a.x));
            const minY = Math.min(...atomPositions.map(a => a.y));
            const maxY = Math.max(...atomPositions.map(a => a.y));
            
            const centerX = (minX + maxX) / 2;
            const centerY = (minY + maxY) / 2;
            
            const canvasRect = canvas.getBoundingClientRect();
            const targetX = canvasRect.width / 2;
            const targetY = canvasRect.height / 2;
            
            const deltaX = targetX - centerX;
            const deltaY = targetY - centerY;
            
            // 모든 원자 이동
            for (const atom of Object.values(atoms)) {
                atom.x += deltaX;
                atom.y += deltaY;
            }
            
            redrawCanvas();
        }
        
        // 캔버스 크기 조정
        function resizeCanvas(width, height) {
            canvas.setAttribute('width', width);
            canvas.setAttribute('height', height);
            redrawCanvas();
        }
        
        // 확대/축소 (향후 구현)
        function zoomIn() {
            console.log('확대 기능은 아직 구현되지 않았습니다.');
        }
        
        function zoomOut() {
            console.log('축소 기능은 아직 구현되지 않았습니다.');
        }
        
        function resetZoom() {
            console.log('확대/축소 초기화 기능은 아직 구현되지 않았습니다.');
        }
        
        // 실행 취소/다시 실행 (향후 구현)
        function undo() {
            console.log('실행 취소 기능은 아직 구현되지 않았습니다.');
        }
        
        function redo() {
            console.log('다시 실행 기능은 아직 구현되지 않았습니다.');
        }
        
        // 디버그 정보 출력
        function debugInfo() {
            console.log('=== 분자 에디터 디버그 정보 ===');
            console.log('현재 모드:', currentMode);
            console.log('현재 원소:', currentElement);
            console.log('현재 결합 타입:', currentBondType);
            console.log('원자 수:', Object.keys(atoms).length);
            console.log('결합 수:', Object.keys(bonds).length);
            console.log('선택된 원자:', Array.from(selectedAtoms));
            console.log('선택된 결합:', Array.from(selectedBonds));
            console.log('분자 통계:', getMoleculeStats());
            console.log('검증 결과:', validateMolecule());
        }
        
        // 전역 함수로 디버그 정보 노출
        window.debugMoleculeEditor = debugInfo;
    """
