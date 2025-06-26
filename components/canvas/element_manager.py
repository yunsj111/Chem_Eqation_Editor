"""
원소 관리 및 주기율표 관련 기능
"""

def get_element_data():
    """원소 데이터와 주기율표 관련 JavaScript 함수들 반환"""
    return """
        // 원소 데이터와 색상 정의
        const elementData = {
            'H': { name: 'Hydrogen', color: '#FFFFFF', number: 1, row: 1, col: 1, type: 'nonmetal' },
            'He': { name: 'Helium', color: '#D9FFFF', number: 2, row: 1, col: 18, type: 'noble-gas' },
            'Li': { name: 'Lithium', color: '#CC80FF', number: 3, row: 2, col: 1, type: 'alkali' },
            'Be': { name: 'Beryllium', color: '#C2FF00', number: 4, row: 2, col: 2, type: 'alkaline' },
            'B': { name: 'Boron', color: '#FFB5B5', number: 5, row: 2, col: 13, type: 'metalloid' },
            'C': { name: 'Carbon', color: '#909090', number: 6, row: 2, col: 14, type: 'nonmetal' },
            'N': { name: 'Nitrogen', color: '#3050F8', number: 7, row: 2, col: 15, type: 'nonmetal' },
            'O': { name: 'Oxygen', color: '#FF0D0D', number: 8, row: 2, col: 16, type: 'nonmetal' },
            'F': { name: 'Fluorine', color: '#90E050', number: 9, row: 2, col: 17, type: 'halogen' },
            'Ne': { name: 'Neon', color: '#B3E3F5', number: 10, row: 2, col: 18, type: 'noble-gas' },
            'Na': { name: 'Sodium', color: '#AB5CF2', number: 11, row: 3, col: 1, type: 'alkali' },
            'Mg': { name: 'Magnesium', color: '#8AFF00', number: 12, row: 3, col: 2, type: 'alkaline' },
            'Al': { name: 'Aluminum', color: '#BFA6A6', number: 13, row: 3, col: 13, type: 'metal' },
            'Si': { name: 'Silicon', color: '#F0C8A0', number: 14, row: 3, col: 14, type: 'metalloid' },
            'P': { name: 'Phosphorus', color: '#FF8000', number: 15, row: 3, col: 15, type: 'nonmetal' },
            'S': { name: 'Sulfur', color: '#FFFF30', number: 16, row: 3, col: 16, type: 'nonmetal' },
            'Cl': { name: 'Chlorine', color: '#1FF01F', number: 17, row: 3, col: 17, type: 'halogen' },
            'Ar': { name: 'Argon', color: '#80D1E3', number: 18, row: 3, col: 18, type: 'noble-gas' },
            'K': { name: 'Potassium', color: '#8F40D4', number: 19, row: 4, col: 1, type: 'alkali' },
            'Ca': { name: 'Calcium', color: '#3DFF00', number: 20, row: 4, col: 2, type: 'alkaline' },
            'Sc': { name: 'Scandium', color: '#E6E6E6', number: 21, row: 4, col: 3, type: 'transition' },
            'Ti': { name: 'Titanium', color: '#BFC2C7', number: 22, row: 4, col: 4, type: 'transition' },
            'V': { name: 'Vanadium', color: '#A6A6AB', number: 23, row: 4, col: 5, type: 'transition' },
            'Cr': { name: 'Chromium', color: '#8A99C7', number: 24, row: 4, col: 6, type: 'transition' },
            'Mn': { name: 'Manganese', color: '#9C7AC7', number: 25, row: 4, col: 7, type: 'transition' },
            'Fe': { name: 'Iron', color: '#E06633', number: 26, row: 4, col: 8, type: 'transition' },
            'Co': { name: 'Cobalt', color: '#F090A0', number: 27, row: 4, col: 9, type: 'transition' },
            'Ni': { name: 'Nickel', color: '#50D050', number: 28, row: 4, col: 10, type: 'transition' },
            'Cu': { name: 'Copper', color: '#C88033', number: 29, row: 4, col: 11, type: 'transition' },
            'Zn': { name: 'Zinc', color: '#7D80B0', number: 30, row: 4, col: 12, type: 'transition' },
            'Ga': { name: 'Gallium', color: '#C28F8F', number: 31, row: 4, col: 13, type: 'metal' },
            'Ge': { name: 'Germanium', color: '#668F8F', number: 32, row: 4, col: 14, type: 'metalloid' },
            'As': { name: 'Arsenic', color: '#BD80E3', number: 33, row: 4, col: 15, type: 'metalloid' },
            'Se': { name: 'Selenium', color: '#FFA100', number: 34, row: 4, col: 16, type: 'nonmetal' },
            'Br': { name: 'Bromine', color: '#A62929', number: 35, row: 4, col: 17, type: 'halogen' },
            'Kr': { name: 'Krypton', color: '#5CB8D1', number: 36, row: 4, col: 18, type: 'noble-gas' },
            'I': { name: 'Iodine', color: '#940094', number: 53, row: 5, col: 17, type: 'halogen' }
        };
        
        // SMARTS 패턴 데이터
        const smartsPatterns = {
            '[!H]': '수소가 아닌 원자',
            '[C&H3]': '3개 수소를 가진 탄소',
            '[N+]': '양전하 질소',
            '[O-]': '음전하 산소',
            '[R]': '고리 원자',
            '[!R]': '비고리 원자',
            '*': '모든 원자'
        };
        
        // 주기율표 생성 함수
        function generatePeriodicTable() {
            const grid = document.getElementById('periodicGrid');
            if (!grid) return;
            
            grid.innerHTML = '';
            
            // 7행 18열 그리드 생성
            for (let row = 1; row <= 7; row++) {
                for (let col = 1; col <= 18; col++) {
                    const element = findElementByPosition(row, col);
                    const cell = document.createElement('div');
                    
                    if (element) {
                        cell.className = `pt-element ${element.type}`;
                        cell.innerHTML = element.symbol;
                        cell.onclick = () => selectElementFromTable(element.symbol);
                        cell.title = `${element.name} (${element.symbol})`;
                        
                        // 현재 선택된 원소 표시
                        if (element.symbol === currentElement) {
                            cell.classList.add('selected');
                        }
                    } else {
                        cell.className = 'pt-element pt-empty';
                    }
                    
                    grid.appendChild(cell);
                }
            }
        }
        
        // 위치로 원소 찾기
        function findElementByPosition(row, col) {
            for (const [symbol, data] of Object.entries(elementData)) {
                if (data.row === row && data.col === col) {
                    return { symbol, ...data };
                }
            }
            return null;
        }
        
        // 주기율표에서 원소 선택
        function selectElementFromTable(element) {
            // 이전 선택 제거
            const prevSelected = document.querySelector('.pt-element.selected');
            if (prevSelected) {
                prevSelected.classList.remove('selected');
            }
            
            // 새로운 선택 추가
            const elements = document.querySelectorAll('.pt-element');
            elements.forEach(el => {
                if (el.title && el.title.includes(`(${element})`)) {
                    el.classList.add('selected');
                }
            });
            
            // 현재 원소 설정
            setCurrentElement(element);
        }
        
        // 주기율표 토글
        function togglePeriodicTable() {
            const table = document.getElementById('periodicTable');
            if (table.classList.contains('show')) {
                table.classList.remove('show');
            } else {
                table.classList.add('show');
                generatePeriodicTable();
            }
        }
        
        // SMARTS 드롭다운 토글
        function toggleSmartsDropdown() {
            const dropdown = document.getElementById('smartsDropdown');
            dropdown.classList.toggle('show');
        }
        
        // SMARTS 패턴 설정
        function setSmartsPattern(pattern) {
            currentElement = pattern;
            document.getElementById('customSmarts').value = pattern;
            toggleSmartsDropdown();
            updateStatusBar();
        }
        
        // 원소 색상 가져오기
        function getElementColor(element) {
            if (element.startsWith('[') && element.endsWith(']'))
                return '#FF6B6B'; // SMARTS 패턴 색상
            if (element === '*')
                return '#4ECDC4'; // 와일드카드 색상
            return elementData[element]?.color || '#909090';
        }
        
        // 원소 반지름 가져오기
        function getElementRadius(element) {
            if (element === 'H') return 6;
            if (element === 'C') return 8;
            if (element.startsWith('[') || element === '*') return 12;
            return 10;
        }
    """
