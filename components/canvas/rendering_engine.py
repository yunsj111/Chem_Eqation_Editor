"""
렌더링 엔진 - 원자와 결합 그리기
"""

def get_rendering_engine():
    """렌더링 JavaScript 함수들 반환"""
    return """
        // 캔버스 다시 그리기
        function redrawCanvas() {
            if (!canvas) return;
            
            // 기존 내용 지우기 (선택 영역과 드래그 라인 제외)
            const elementsToRemove = canvas.querySelectorAll('.atom-circle, .bond-line, .atom-label');
            elementsToRemove.forEach(element => element.remove());
            
            // 결합 먼저 그리기 (원자 아래에 표시되도록)
            for (const [bondId, bond] of Object.entries(bonds)) {
                drawBond(bond);
            }
            
            // 원자 그리기
            for (const [atomId, atom] of Object.entries(atoms)) {
                drawAtom(atom);
            }
            
            // 선택 시각화 업데이트
            updateSelectionVisuals();
        }
        
        // 원자 그리기
        function drawAtom(atom) {
            const element = atom.element;
            const color = getElementColor(element);
            const radius = getElementRadius(element);
            
            if (element === 'C') {
                // 탄소는 작은 점으로 표시
                drawCarbonDot(atom, color, radius);
            } else if (element === 'H') {
                // 수소는 작은 회색 점으로 표시
                drawHydrogenDot(atom, color, radius);
            } else if (element.startsWith('[') && element.endsWith(']')) {
                // SMARTS 패턴은 사각형으로 표시
                drawSmartsPattern(atom, color, radius);
            } else if (element === '*') {
                // 와일드카드는 다이아몬드로 표시
                drawWildcard(atom, color, radius);
            } else {
                // 헤테로원자는 색깔 있는 원과 라벨로 표시
                drawHeteroatom(atom, color, radius);
            }
        }
        
        // 탄소 점 그리기
        function drawCarbonDot(atom, color, radius) {
            const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
            circle.setAttribute('id', atom.id);
            circle.setAttribute('class', 'atom-circle carbon-dot');
            circle.setAttribute('cx', atom.x);
            circle.setAttribute('cy', atom.y);
            circle.setAttribute('r', radius / 2);
            circle.setAttribute('fill', color);
            circle.setAttribute('stroke', '#666666');
            circle.setAttribute('stroke-width', '1');
            canvas.appendChild(circle);
        }
        
        // 수소 점 그리기
        function drawHydrogenDot(atom, color, radius) {
            const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
            circle.setAttribute('id', atom.id);
            circle.setAttribute('class', 'atom-circle hydrogen-dot');
            circle.setAttribute('cx', atom.x);
            circle.setAttribute('cy', atom.y);
            circle.setAttribute('r', radius / 2);
            circle.setAttribute('fill', color);
            circle.setAttribute('stroke', '#999999');
            circle.setAttribute('stroke-width', '1');
            canvas.appendChild(circle);
        }
        
        // SMARTS 패턴 그리기
        function drawSmartsPattern(atom, color, radius) {
            const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('id', atom.id);
            rect.setAttribute('class', 'atom-circle');
            rect.setAttribute('x', atom.x - radius);
            rect.setAttribute('y', atom.y - radius);
            rect.setAttribute('width', radius * 2);
            rect.setAttribute('height', radius * 2);
            rect.setAttribute('fill', color);
            rect.setAttribute('stroke', '#333');
            rect.setAttribute('stroke-width', '2');
            rect.setAttribute('rx', '3');
            canvas.appendChild(rect);
            
            // SMARTS 패턴 텍스트
            const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('class', 'atom-label');
            text.setAttribute('x', atom.x);
            text.setAttribute('y', atom.y + 4);
            text.setAttribute('text-anchor', 'middle');
            text.setAttribute('font-family', 'monospace');
            text.setAttribute('font-size', '10');
            text.setAttribute('fill', 'white');
            text.textContent = atom.element.length > 6 ? atom.element.substring(0, 4) + '...' : atom.element;
            canvas.appendChild(text);
        }
        
        // 와일드카드 그리기
        function drawWildcard(atom, color, radius) {
            const diamond = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
            const points = [
                [atom.x, atom.y - radius],
                [atom.x + radius, atom.y],
                [atom.x, atom.y + radius],
                [atom.x - radius, atom.y]
            ].map(point => point.join(',')).join(' ');
            
            diamond.setAttribute('id', atom.id);
            diamond.setAttribute('class', 'atom-circle');
            diamond.setAttribute('points', points);
            diamond.setAttribute('fill', color);
            diamond.setAttribute('stroke', '#333');
            diamond.setAttribute('stroke-width', '2');
            canvas.appendChild(diamond);
            
            // 와일드카드 텍스트
            const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('class', 'atom-label');
            text.setAttribute('x', atom.x);
            text.setAttribute('y', atom.y + 4);
            text.setAttribute('text-anchor', 'middle');
            text.setAttribute('font-family', 'Arial, sans-serif');
            text.setAttribute('font-size', '12');
            text.setAttribute('font-weight', 'bold');
            text.setAttribute('fill', 'white');
            text.textContent = '*';
            canvas.appendChild(text);
        }
        
        // 헤테로원자 그리기
        function drawHeteroatom(atom, color, radius) {
            const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
            circle.setAttribute('id', atom.id);
            circle.setAttribute('class', 'atom-circle heteroatom-circle');
            circle.setAttribute('cx', atom.x);
            circle.setAttribute('cy', atom.y);
            circle.setAttribute('r', radius);
            circle.setAttribute('fill', color);
            circle.setAttribute('stroke', '#333');
            circle.setAttribute('stroke-width', '2');
            canvas.appendChild(circle);
            
            // 원소 기호 텍스트
            const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('class', 'atom-label');
            text.setAttribute('x', atom.x);
            text.setAttribute('y', atom.y + 4);
            text.setAttribute('text-anchor', 'middle');
            text.setAttribute('font-family', 'Arial, sans-serif');
            text.setAttribute('font-size', '14');
            text.setAttribute('font-weight', 'bold');
            text.setAttribute('fill', 'white');
            text.textContent = atom.element;
            canvas.appendChild(text);
        }
        
        // 결합 그리기
        function drawBond(bond) {
            const atom1 = atoms[bond.atom1];
            const atom2 = atoms[bond.atom2];
            
            if (!atom1 || !atom2) return;
            
            switch (bond.type) {
                case 'single':
                    drawSingleBond(bond, atom1, atom2);
                    break;
                case 'double':
                    drawDoubleBond(bond, atom1, atom2);
                    break;
                case 'triple':
                    drawTripleBond(bond, atom1, atom2);
                    break;
                case 'aromatic':
                    drawAromaticBond(bond, atom1, atom2);
                    break;
            }
        }
        
        // 단일 결합 그리기
        function drawSingleBond(bond, atom1, atom2) {
            const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line.setAttribute('id', bond.id);
            line.setAttribute('class', 'bond-line');
            line.setAttribute('x1', atom1.x);
            line.setAttribute('y1', atom1.y);
            line.setAttribute('x2', atom2.x);
            line.setAttribute('y2', atom2.y);
            line.setAttribute('stroke', '#333');
            line.setAttribute('stroke-width', '2');
            canvas.appendChild(line);
        }
        
        // 이중 결합 그리기
        function drawDoubleBond(bond, atom1, atom2) {
            const dx = atom2.x - atom1.x;
            const dy = atom2.y - atom1.y;
            const length = Math.sqrt(dx * dx + dy * dy);
            const unitX = dx / length;
            const unitY = dy / length;
            const perpX = -unitY * 3; // 선 간격
            const perpY = unitX * 3;
            
            // 첫 번째 선
            const line1 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line1.setAttribute('id', bond.id + '_1');
            line1.setAttribute('class', 'bond-line');
            line1.setAttribute('x1', atom1.x + perpX);
            line1.setAttribute('y1', atom1.y + perpY);
            line1.setAttribute('x2', atom2.x + perpX);
            line1.setAttribute('y2', atom2.y + perpY);
            line1.setAttribute('stroke', '#333');
            line1.setAttribute('stroke-width', '2');
            canvas.appendChild(line1);
            
            // 두 번째 선
            const line2 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line2.setAttribute('id', bond.id + '_2');
            line2.setAttribute('class', 'bond-line');
            line2.setAttribute('x1', atom1.x - perpX);
            line2.setAttribute('y1', atom1.y - perpY);
            line2.setAttribute('x2', atom2.x - perpX);
            line2.setAttribute('y2', atom2.y - perpY);
            line2.setAttribute('stroke', '#333');
            line2.setAttribute('stroke-width', '2');
            canvas.appendChild(line2);
        }
        
        // 삼중 결합 그리기
        function drawTripleBond(bond, atom1, atom2) {
            const dx = atom2.x - atom1.x;
            const dy = atom2.y - atom1.y;
            const length = Math.sqrt(dx * dx + dy * dy);
            const unitX = dx / length;
            const unitY = dy / length;
            const perpX = -unitY * 4; // 선 간격
            const perpY = unitX * 4;
            
            // 중앙 선
            const line1 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line1.setAttribute('id', bond.id + '_center');
            line1.setAttribute('class', 'bond-line');
            line1.setAttribute('x1', atom1.x);
            line1.setAttribute('y1', atom1.y);
            line1.setAttribute('x2', atom2.x);
            line1.setAttribute('y2', atom2.y);
            line1.setAttribute('stroke', '#333');
            line1.setAttribute('stroke-width', '2');
            canvas.appendChild(line1);
            
            // 위쪽 선
            const line2 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line2.setAttribute('id', bond.id + '_top');
            line2.setAttribute('class', 'bond-line');
            line2.setAttribute('x1', atom1.x + perpX);
            line2.setAttribute('y1', atom1.y + perpY);
            line2.setAttribute('x2', atom2.x + perpX);
            line2.setAttribute('y2', atom2.y + perpY);
            line2.setAttribute('stroke', '#333');
            line2.setAttribute('stroke-width', '2');
            canvas.appendChild(line2);
            
            // 아래쪽 선
            const line3 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line3.setAttribute('id', bond.id + '_bottom');
            line3.setAttribute('class', 'bond-line');
            line3.setAttribute('x1', atom1.x - perpX);
            line3.setAttribute('y1', atom1.y - perpY);
            line3.setAttribute('x2', atom2.x - perpX);
            line3.setAttribute('y2', atom2.y - perpY);
            line3.setAttribute('stroke', '#333');
            line3.setAttribute('stroke-width', '2');
            canvas.appendChild(line3);
        }
        
        // 방향족 결합 그리기
        function drawAromaticBond(bond, atom1, atom2) {
            const dx = atom2.x - atom1.x;
            const dy = atom2.y - atom1.y;
            const length = Math.sqrt(dx * dx + dy * dy);
            const unitX = dx / length;
            const unitY = dy / length;
            const perpX = -unitY * 3;
            const perpY = unitX * 3;
            
            // 실선
            const line1 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line1.setAttribute('id', bond.id + '_solid');
            line1.setAttribute('class', 'bond-line');
            line1.setAttribute('x1', atom1.x + perpX);
            line1.setAttribute('y1', atom1.y + perpY);
            line1.setAttribute('x2', atom2.x + perpX);
            line1.setAttribute('y2', atom2.y + perpY);
            line1.setAttribute('stroke', '#333');
            line1.setAttribute('stroke-width', '2');
            canvas.appendChild(line1);
            
            // 점선
            const line2 = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line2.setAttribute('id', bond.id + '_dashed');
            line2.setAttribute('class', 'bond-line');
            line2.setAttribute('x1', atom1.x - perpX);
            line2.setAttribute('y1', atom1.y - perpY);
            line2.setAttribute('x2', atom2.x - perpX);
            line2.setAttribute('y2', atom2.y - perpY);
            line2.setAttribute('stroke', '#333');
            line2.setAttribute('stroke-width', '2');
            line2.setAttribute('stroke-dasharray', '5,5');
            canvas.appendChild(line2);
        }
    """
