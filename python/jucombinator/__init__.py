import rdkit.Chem as rdkit

from jucombinator import jucombinator
from jucombinator.jucombinator import (
    AromaticBond,
    Bond,
    Molecule,
    Substituent,
    SubstitutedMolecule,
)


def _smiles_to_molecule(smiles: str) -> Molecule:
    molecule = rdkit.MolFromSmiles(smiles)
    atomic_numbers = []
    num_implicit_hs = []
    for atom in molecule.GetAtoms():
        atomic_numbers.append(atom.GetAtomicNum())
        num_implicit_hs.append(atom.GetNumImplicitHs())
    bonds = []
    aromatic_bonds = []
    for bond in molecule.GetBonds():
        match bond.GetBondType():
            case (
                rdkit.BondType.SINGLE
                | rdkit.BondType.DOUBLE
                | rdkit.BondType.TRIPLE
                | rdkit.BondType.QUADRUPLE
                | rdkit.BondType.QUINTUPLE
                | rdkit.BondType.HEXTUPLE
            ):
                bonds.append(
                    Bond(
                        bond.GetBeginAtomIdx(),
                        bond.GetEndAtomIdx(),
                        int(bond.GetBondTypeAsDouble()),
                    )
                )
            case rdkit.BondType.AROMATIC:
                aromatic_bonds.append(
                    AromaticBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                )
            case _:
                raise RuntimeError(
                    f"unsupported bond type: {bond.GetBondType()}"
                )

    return Molecule(atomic_numbers, num_implicit_hs, bonds, aromatic_bonds)


def _smiles_to_substituent(smiles: str) -> Substituent:
    molecule = rdkit.MolFromSmiles(smiles)
    atomic_numbers = []
    for atom in molecule.GetAtoms():
        atomic_numbers.append(atom.GetAtomicNum())
    bonds = []
    for bond in molecule.GetBonds():
        match bond.GetBondType():
            case (
                rdkit.BondType.SINGLE
                | rdkit.BondType.DOUBLE
                | rdkit.BondType.TRIPLE
                | rdkit.BondType.QUADRUPLE
                | rdkit.BondType.QUINTUPLE
                | rdkit.BondType.HEXTUPLE
            ):
                bonds.append(
                    Bond(
                        bond.GetBeginAtomIdx(),
                        bond.GetEndAtomIdx(),
                        int(bond.GetBondTypeAsDouble()),
                    )
                )
            case _:
                raise RuntimeError(
                    f"unsupported bond type: {bond.GetBondType()}"
                )

    return Substituent(atomic_numbers, bonds)


def _substituted_molecule_to_smiles(molecule: SubstitutedMolecule) -> str:
    result = rdkit.EditableMol(rdkit.Mol())
    for atomic_number in molecule.atomic_numbers:
        result.AddAtom(rdkit.Atom(atomic_number))
    for bond in molecule.bonds:
        result.AddBond(
            bond.atom1_idx, bond.atom2_idx, rdkit.BondType(bond.order)
        )
    for aromatic_bond in molecule.aromatic_bonds:
        result.AddBond(
            aromatic_bond.atom1_idx,
            aromatic_bond.atom2_idx,
            rdkit.BondType.AROMATIC,
        )
    return rdkit.MolToSmiles(result.GetMol(), canonical=True)


def substitute_1(skeleton: str, substituents: list[str]) -> set[str]:
    skeleton_ = _smiles_to_molecule(skeleton)
    substituents_ = list(map(_smiles_to_substituent, substituents))
    return set(
        map(
            _substituted_molecule_to_smiles,
            jucombinator.substitute_1(skeleton_, substituents_),
        )
    )


def substitute(skeleton: str, substituents: list[str], n: int) -> set[str]:
    skeleton_ = _smiles_to_molecule(skeleton)
    substituents_ = list(map(_smiles_to_substituent, substituents))
    return set(
        map(
            _substituted_molecule_to_smiles,
            jucombinator.substitute(skeleton_, substituents_, n),
        )
    )
