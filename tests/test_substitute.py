import rdkit.Chem as rdkit
from jucombinator import substitute


def to_canonical_smiles(smiles: str) -> str:
    return rdkit.MolToSmiles(rdkit.MolFromSmiles(smiles), canonical=True)


def test_substitute_1() -> None:
    result = set(map(to_canonical_smiles, substitute("CCC", ["Br", "NO"], 1)))
    expected = set(
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


def test_substitute_1_aromatic() -> None:
    result = set(
        map(to_canonical_smiles, substitute("Cc1cccc(C)c1", ["Br", "NO"], 1))
    )
    expected = set(
        map(
            to_canonical_smiles,
            [
                "C(Br)c1cccc(C)c1",
                "Cc1c(Br)ccc(C)c1",
                "Cc1cc(Br)cc(C)c1",
                "Cc1ccc(Br)c(C)c1",
                "Cc1cccc(CBr)c1",
                "Cc1cccc(C)c1(Br)",
                "C(NO)c1cccc(C)c1",
                "Cc1c(NO)ccc(C)c1",
                "Cc1cc(NO)cc(C)c1",
                "Cc1ccc(NO)c(C)c1",
                "Cc1cccc(CNO)c1",
                "Cc1cccc(C)c1(NO)",
            ],
        )
    )
    assert result == expected
