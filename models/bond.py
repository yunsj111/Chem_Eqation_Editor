from dataclasses import dataclass
from typing import Optional
from enum import Enum

class BondType(Enum):
    """결합 타입"""
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4
    ANY = 5        # ~
    RING = 6       # @
    NOT_RING = 7   # !@

@dataclass
class BondQuery:
    """SMARTS 결합 쿼리 조건"""
    bond_type: BondType = BondType.SINGLE
    ring_membership: Optional[bool] = None
    negation: bool = False

    def to_smarts(self) -> str:
        """SMARTS 문자열로 변환"""
        if self.bond_type == BondType.SINGLE and not self.ring_membership:
            return ""  # 단일결합은 생략 가능

        symbol = ""
        if self.bond_type == BondType.SINGLE:
            symbol = "-"
        elif self.bond_type == BondType.DOUBLE:
            symbol = "="
        elif self.bond_type == BondType.TRIPLE:
            symbol = "#"
        elif self.bond_type == BondType.AROMATIC:
            symbol = ":"
        elif self.bond_type == BondType.ANY:
            symbol = "~"

        if self.ring_membership is not None:
            if self.ring_membership:
                symbol += "@"
            else:
                symbol += "!@"

        return symbol

@dataclass
class Bond:
    """결합 클래스"""
    id: int
    atom1_id: int
    atom2_id: int
    bond_type: BondType = BondType.SINGLE
    query: Optional[BondQuery] = None

    def __post_init__(self):
        if self.query is None:
            self.query = BondQuery(bond_type=self.bond_type)

    def to_smarts(self) -> str:
        """SMARTS 표현으로 변환"""
        return self.query.to_smarts()

    def get_display_width(self) -> int:
        """화면 표시용 선 두께"""
        if self.bond_type == BondType.SINGLE:
            return 2
        elif self.bond_type == BondType.DOUBLE:
            return 3
        elif self.bond_type == BondType.TRIPLE:
            return 4
        elif self.bond_type == BondType.AROMATIC:
            return 2
        else:
            return 2