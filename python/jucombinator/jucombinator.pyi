class Bond:
    def __init__(self, atom1_idx: int, atom2_idx: int, order: int) -> None: ...
    @property
    def atom1_idx(self) -> int: ...
    @property
    def atom2_idx(self) -> int: ...
    @property
    def order(self) -> int: ...

class AromaticBond:
    def __init__(self, atom1_idx: int, atom2_idx: int) -> None: ...
    @property
    def atom1_idx(self) -> int: ...
    @property
    def atom2_idx(self) -> int: ...

class Molecule:
    def __init__(
        self,
        atomic_numbers: list[int],
        num_implicit_hs: list[int],
        bonds: list[Bond],
        aromatic_bonds: list[AromaticBond],
    ) -> None: ...

class Substituent:
    def __init__(
        self,
        atomic_numbers: list[int],
        bonds: list[Bond],
    ) -> None: ...

class SubstitutedMolecule:
    @property
    def atomic_numbers(self) -> list[int]: ...
    @property
    def bonds(self) -> list[Bond]: ...
    @property
    def aromatic_bonds(self) -> list[AromaticBond]: ...

def substitute_1(
    skeleton: Molecule,
    substituents: list[Substituent],
) -> list[SubstitutedMolecule]: ...
def substitute(
    skeleton: Molecule,
    substituents: list[Substituent],
    n: int,
) -> list[SubstitutedMolecule]: ...
