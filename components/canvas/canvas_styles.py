"""
화학 구조 에디터 CSS 스타일 정의
"""

def get_canvas_styles():
    """캔버스와 UI 요소들의 CSS 스타일 반환"""
    return """
        <style>
            body {
                margin: 0;
                padding: 10px;
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background: #f8f9fa;
            }
            
            .main-content {
                display: flex;
                gap: 20px;
                margin-top: 10px;
            }
            
            .canvas-container {
                flex: 1;
            }
            
            .selected-element-display {
                display: flex;
                align-items: center;
                gap: 8px;
                margin-right: 10px;
            }
            
            .element-badge {
                background: #007bff;
                color: white;
                padding: 4px 8px;
                border-radius: 4px;
                font-weight: bold;
                min-width: 20px;
                text-align: center;
            }
            
            #moleculeCanvas {
                border: 2px solid #ccc;
                border-radius: 8px;
                cursor: crosshair;
                background: white;
                background-image: 
                    repeating-linear-gradient(0deg, transparent, transparent 49px, #f0f0f0 50px),
                    repeating-linear-gradient(90deg, transparent, transparent 49px, #f0f0f0 50px);
            }
            
            .toolbar {
                margin: 10px 0;
                padding: 10px;
                background: #f8f9fa;
                border-radius: 5px;
                display: flex;
                gap: 10px;
                flex-wrap: wrap;
                align-items: center;
            }
            
            .tool-group {
                display: flex;
                gap: 5px;
                align-items: center;
                border-right: 1px solid #ddd;
                padding-right: 10px;
            }
            
            .tool-btn {
                padding: 8px 12px;
                border: 1px solid #ddd;
                background: white;
                cursor: pointer;
                border-radius: 4px;
                font-size: 14px;
                transition: all 0.2s;
            }
            
            .tool-btn:hover {
                background: #e9ecef;
            }
            
            .tool-btn.active {
                background: #007bff;
                color: white;
                border-color: #007bff;
            }
            
            .tool-group label {
                margin-right: 5px;
                font-weight: bold;
                font-size: 14px;
                color: #333;
            }
            
            .element-btn {
                width: 40px;
                height: 40px;
                border-radius: 50%;
                font-weight: bold;
                display: flex;
                align-items: center;
                justify-content: center;
            }
            
            .status-bar {
                margin-top: 10px;
                padding: 5px 10px;
                background: #e9ecef;
                border-radius: 4px;
                font-size: 12px;
                color: #495057;
            }
            
            .atom-circle {
                cursor: pointer;
            }
            
            .carbon-dot {
                stroke: #666666;
                stroke-width: 1;
                opacity: 0.8;
            }
            
            .carbon-dot:hover {
                opacity: 1;
                stroke-width: 2;
            }
            
            .hydrogen-dot {
                stroke: #999999;
                stroke-width: 1;
                opacity: 0.6;
            }
            
            .heteroatom-circle {
                cursor: pointer;
                transition: all 0.2s;
            }
            
            .heteroatom-circle:hover {
                stroke-width: 3;
            }
            
            .selection-area {
                fill: rgba(0, 123, 255, 0.1);
                stroke: #007bff;
                stroke-width: 2;
                stroke-dasharray: 5,5;
                pointer-events: none;
            }
            
            .selected-atom {
                stroke: #007bff !important;
                stroke-width: 3 !important;
                filter: drop-shadow(0 0 8px rgba(0, 123, 255, 0.5));
            }
            
            .selected-bond {
                stroke: #007bff !important;
                stroke-width: 4 !important;
                filter: drop-shadow(0 0 8px rgba(0, 123, 255, 0.5));
            }
            
            .periodic-table {
                display: none;
                width: 400px;
                background: white;
                border: 2px solid #ddd;
                border-radius: 8px;
                padding: 15px;
                box-shadow: 0 4px 20px rgba(0,0,0,0.15);
                height: fit-content;
                max-height: 600px;
                overflow-y: auto;
            }
            
            .periodic-table.show {
                display: block;
            }
            
            .periodic-header {
                display: flex;
                justify-content: space-between;
                align-items: center;
                margin-bottom: 15px;
                padding-bottom: 10px;
                border-bottom: 1px solid #eee;
            }
            
            .periodic-header h3 {
                margin: 0;
                color: #495057;
                font-size: 18px;
            }
            
            .close-btn {
                background: #dc3545;
                color: white;
                border: none;
                width: 30px;
                height: 30px;
                border-radius: 50%;
                cursor: pointer;
                display: flex;
                align-items: center;
                justify-content: center;
                font-size: 18px;
                font-weight: bold;
            }
            
            .close-btn:hover {
                background: #c82333;
            }
            
            .periodic-grid {
                display: grid;
                grid-template-columns: repeat(18, 1fr);
                gap: 2px;
                margin-bottom: 10px;
            }
            
            .pt-element {
                width: 18px;
                height: 18px;
                border: 1px solid #ddd;
                border-radius: 3px;
                display: flex;
                align-items: center;
                justify-content: center;
                font-size: 10px;
                font-weight: bold;
                cursor: pointer;
                transition: all 0.2s;
                background: white;
                position: relative;
            }
            
            .pt-element:hover {
                border-color: #007bff;
                background: #e7f3ff;
                transform: scale(1.2);
                z-index: 10;
                box-shadow: 0 2px 8px rgba(0,123,255,0.3);
            }
            
            .pt-element.selected {
                border-color: #007bff;
                background: #007bff;
                color: white;
                transform: scale(1.1);
            }
            
            .pt-element.nonmetal { background-color: #ffffcc; }
            .pt-element.noble-gas { background-color: #ffccff; }
            .pt-element.alkali { background-color: #ffcccc; }
            .pt-element.alkaline { background-color: #ffddcc; }
            .pt-element.metalloid { background-color: #ccffcc; }
            .pt-element.halogen { background-color: #ccffff; }
            .pt-element.metal { background-color: #ffeecc; }
            .pt-element.transition { background-color: #ddddff; }
            .pt-element.lanthanide { background-color: #ffddff; }
            .pt-element.actinide { background-color: #ddffdd; }
            
            .pt-empty {
                background: transparent;
                border: none;
                cursor: default;
            }
            
            .pt-element .symbol {
                font-size: 14px;
                line-height: 1;
            }
            
            .pt-element .number {
                position: absolute;
                top: 2px;
                left: 2px;
                font-size: 8px;
                color: #666;
            }
            
            .smarts-dropdown {
                position: relative;
                display: inline-block;
            }
            
            .smarts-content {
                display: none;
                position: absolute;
                background-color: white;
                min-width: 200px;
                box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
                z-index: 1001;
                border-radius: 4px;
                border: 1px solid #ddd;
                max-height: 300px;
                overflow-y: auto;
            }
            
            .smarts-content.show {
                display: block;
            }
            
            .smarts-item {
                color: black;
                padding: 8px 12px;
                text-decoration: none;
                display: block;
                cursor: pointer;
                border-bottom: 1px solid #eee;
                font-family: monospace;
                font-size: 12px;
            }
            
            .smarts-item:hover {
                background-color: #f1f1f1;
            }
            
            .smarts-input {
                width: 150px;
                padding: 6px;
                border: 1px solid #ddd;
                border-radius: 4px;
                font-family: monospace;
                font-size: 12px;
            }
            
            .bond-drag-line {
                stroke: #007bff;
                stroke-width: 2;
                stroke-dasharray: 5,5;
                pointer-events: none;
                opacity: 0.7;
            }
            
            .drag-preview {
                fill: rgba(0, 123, 255, 0.3);
                stroke: #007bff;
                stroke-width: 2;
                stroke-dasharray: 3,3;
                pointer-events: none;
            }
        </style>
    """
