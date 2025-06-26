from typing import List, Dict, Optional, Tuple
from .atom import Atom, AtomType
from .bond import Bond, BondType
import networkx as nx

class Molecule:
    """분자 클래스"""

    def __init__(self):
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self.atom_counter = 0
        self.bond_counter = 0

    def add_atom(self, x: float, y: float, element: str = "C", 
                 atom_type: AtomType = AtomType.ELEMENT) -> Atom:
        """원자 추가"""
        atom = Atom(
            id=self.atom_counter,
            x=x,
            y=y,
            element=element,
            atom_type=atom_type
        )
        self.atoms.append(atom)
        self.atom_counter += 1
        return atom

    def add_bond(self, atom1_id: int, atom2_id: int, 
                 bond_type: BondType = BondType.SINGLE) -> Bond:
        """결합 추가"""
        # 이미 존재하는 결합인지 확인
        for bond in self.bonds:
            if ((bond.atom1_id == atom1_id and bond.atom2_id == atom2_id) or
                (bond.atom1_id == atom2_id and bond.atom2_id == atom1_id)):
                return bond

        bond = Bond(
            id=self.bond_counter,
            atom1_id=atom1_id,
            atom2_id=atom2_id,
            bond_type=bond_type
        )
        self.bonds.append(bond)
        self.bond_counter += 1
        return bond

    def remove_atom(self, atom_id: int):
        """원자 제거 (연결된 결합도 함께 제거)"""
        self.atoms = [a for a in self.atoms if a.id != atom_id]
        self.bonds = [b for b in self.bonds 
                     if b.atom1_id != atom_id and b.atom2_id != atom_id]

    def remove_bond(self, bond_id: int):
        """결합 제거"""
        self.bonds = [b for b in self.bonds if b.id != bond_id]

    def get_atom_by_id(self, atom_id: int) -> Optional[Atom]:
        """ID로 원자 찾기"""
        for atom in self.atoms:
            if atom.id == atom_id:
                return atom
        return None

    def get_bond_by_id(self, bond_id: int) -> Optional[Bond]:
        """ID로 결합 찾기"""
        for bond in self.bonds:
            if bond.id == bond_id:
                return bond
        return None

    def find_atom_at(self, x: float, y: float, radius: float = 15) -> Optional[Atom]:
        """좌표에서 원자 찾기"""
        for atom in self.atoms:
            distance = ((atom.x - x) ** 2 + (atom.y - y) ** 2) ** 0.5
            if distance <= radius:
                return atom
        return None

    def get_connected_atoms(self, atom_id: int) -> List[int]:
        """연결된 원자 ID 목록 반환"""
        connected = []
        for bond in self.bonds:
            if bond.atom1_id == atom_id:
                connected.append(bond.atom2_id)
            elif bond.atom2_id == atom_id:
                connected.append(bond.atom1_id)
        return connected

    def to_networkx(self) -> nx.Graph:
        """NetworkX 그래프로 변환"""
        G = nx.Graph()

        # 노드 추가
        for atom in self.atoms:
            G.add_node(atom.id, 
                      element=atom.element,
                      atom_type=atom.atom_type,
                      query=atom.query)

        # 엣지 추가
        for bond in self.bonds:
            G.add_edge(bond.atom1_id, bond.atom2_id,
                      bond_type=bond.bond_type,
                      query=bond.query)

        return G

    def clear(self):
        """모든 원자와 결합 제거"""
        self.atoms.clear()
        self.bonds.clear()
        self.atom_counter = 0
        self.bond_counter = 0

    def get_bounds(self) -> Tuple[float, float, float, float]:
        """분자의 경계 좌표 반환 (min_x, min_y, max_x, max_y)"""
        if not self.atoms:
            return 0, 0, 0, 0

        min_x = min(atom.x for atom in self.atoms)
        max_x = max(atom.x for atom in self.atoms)
        min_y = min(atom.y for atom in self.atoms)
        max_y = max(atom.y for atom in self.atoms)

        return min_x, min_y, max_x, max_y