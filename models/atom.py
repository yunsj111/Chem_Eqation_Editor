from dataclasses import dataclass
from typing import Optional, List, Dict, Any
from enum import Enum

class AtomType(Enum):
    """원자 타입 정의"""
    ELEMENT = "element"      # 특정 원소 (C, N, O, etc.)
    WILDCARD = "wildcard"    # * (any atom)
    AROMATIC = "aromatic"    # 방향족 원자 (c, n, o, etc.)
    ALIPHATIC = "aliphatic"  # 지방족 원자
    CUSTOM = "custom"        # 사용자 정의

@dataclass
class AtomQuery:
    """SMARTS 원자 쿼리 조건"""
    element: Optional[str] = None
    aromatic: Optional[bool] = None
    charge: Optional[int] = None
    hydrogen_count: Optional[int] = None
    degree: Optional[int] = None
    valence: Optional[int] = None
    ring_membership: Optional[bool] = None
    ring_size: Optional[int] = None
    is_wildcard: bool = False
    negation: bool = False  # [!C] 같은 부정 조건

    def to_smarts(self) -> str:
        """SMARTS 문자열로 변환"""
        if self.is_wildcard:
            return "*"

        if self.element and not any([self.aromatic, self.charge, self.hydrogen_count, 
                                   self.degree, self.valence, self.ring_membership]):
            return self.element.lower() if self.aromatic else self.element

        # 복잡한 쿼리는 대괄호로 감싸기
        parts = []

        if self.negation:
            parts.append("!")

        if self.element:
            parts.append(self.element.lower() if self.aromatic else self.element)

        if self.charge is not None:
            if self.charge > 0:
                parts.append(f"+{self.charge}" if self.charge > 1 else "+")
            elif self.charge < 0:
                parts.append(f"{self.charge}" if self.charge < -1 else "-")

        if self.hydrogen_count is not None:
            parts.append(f"H{self.hydrogen_count}" if self.hydrogen_count > 1 else "H")

        if self.degree is not None:
            parts.append(f"D{self.degree}")

        if self.valence is not None:
            parts.append(f"v{self.valence}")

        if self.ring_membership is not None:
            if self.ring_membership:
                if self.ring_size:
                    parts.append(f"R{self.ring_size}")
                else:
                    parts.append("R")
            else:
                parts.append("R0")

        if len(parts) == 1 and not self.negation:
            return parts[0]
        else:
            return f"[{''.join(parts)}]"

@dataclass
class Atom:
    """원자 클래스"""
    id: int
    x: float
    y: float
    atom_type: AtomType = AtomType.ELEMENT
    element: str = "C"
    query: Optional[AtomQuery] = None

    def __post_init__(self):
        if self.query is None:
            self.query = AtomQuery(element=self.element)

    def to_smarts(self) -> str:
        """SMARTS 표현으로 변환"""
        return self.query.to_smarts()

    def get_display_text(self) -> str:
        """화면 표시용 텍스트"""
        if self.atom_type == AtomType.WILDCARD:
            return "*"
        elif self.query and self.query.negation:
            return f"[!{self.element}]"
        elif self.atom_type == AtomType.AROMATIC:
            return self.element.lower()
        else:
            return self.element