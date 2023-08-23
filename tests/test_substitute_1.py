import rdkit.Chem as rdkit
from jucombinator import (
    Bond,
    Molecule,
    Substituent,
    SubstitutedMolecule,
    substitute_1,
)


def smiles_to_molecule(smiles: str) -> Molecule:
    molecule = rdkit.MolFromSmiles(smiles)
    atomic_numbers = []
    num_implicit_hs = []
    for atom in molecule.GetAtoms():
        atomic_numbers.append(atom.GetAtomicNum())
        num_implicit_hs.append(atom.GetNumImplicitHs())
    bonds = []
    for bond in molecule.GetBonds():
        bonds.append(
            Bond(
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond.GetBondTypeAsDouble(),
            )
        )

    return Molecule(atomic_numbers, num_implicit_hs, bonds)


def smiles_to_substituent(smiles: str) -> Substituent:
    molecule = rdkit.MolFromSmiles(smiles)
    atomic_numbers = []
    for atom in molecule.GetAtoms():
        atomic_numbers.append(atom.GetAtomicNum())
    bonds = []
    for bond in molecule.GetBonds():
        bonds.append(
            Bond(
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond.GetBondTypeAsDouble(),
            )
        )

    return Substituent(atomic_numbers, bonds)


def substituted_molecule_to_smiles(molecule: SubstitutedMolecule) -> str:
    result = rdkit.EditableMol(rdkit.Mol())
    for atomic_number in molecule.atomic_numbers:
        result.AddAtom(rdkit.Atom(atomic_number))
    for bond in molecule.bonds:
        result.AddBond(
            bond.atom1_idx, bond.atom2_idx, rdkit.BondType(bond.order)
        )
    return rdkit.MolToSmiles(result.GetMol())


def to_canonical_smiles(smiles: str) -> str:
    return rdkit.MolToSmiles(rdkit.MolFromSmiles(smiles), canonical=True)


def s1(skeleton: str, substituents: list[str]) -> list[str]:
    skeleton_ = smiles_to_molecule(skeleton)
    substituents_ = list(map(smiles_to_substituent, substituents))
    return list(
        map(
            substituted_molecule_to_smiles,
            substitute_1(skeleton_, substituents_),
        )
    )


def test_substitute_1() -> None:
    result = sorted(map(to_canonical_smiles, s1("CCC", ["Br", "NO"])))
    expected = sorted(
        map(
            to_canonical_smiles,
            [
                "C(Br)CC",
                "CC(Br)C",
                "CCC(Br)",
                "C(NO)CC",
                "CC(NO)C",
                "CCC(NO)",
            ],
        )
    )
    assert result == expected
