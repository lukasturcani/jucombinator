use pyo3::prelude::*;

#[derive(Clone, Copy)]
#[pyclass]
struct Bond {
    #[pyo3(get)]
    atom1_idx: u16,
    #[pyo3(get)]
    atom2_idx: u16,
    #[pyo3(get)]
    order: f32,
}

#[pymethods]
impl Bond {
    #[new]
    fn new(atom1_idx: u16, atom2_idx: u16, order: f32) -> Self {
        Self {
            atom1_idx,
            atom2_idx,
            order,
        }
    }
}

#[pyclass]
struct Molecule {
    atomic_numbers: Vec<u8>,
    num_implicit_hs: Vec<u8>,
    bonds: Vec<Bond>,
}

#[pymethods]
impl Molecule {
    #[new]
    fn new(atomic_numbers: Vec<u8>, num_implicit_hs: Vec<u8>, bonds: Vec<Bond>) -> Self {
        Self {
            atomic_numbers,
            num_implicit_hs,
            bonds,
        }
    }
}

#[pyclass]
struct Substituent {
    atomic_numbers: Vec<u8>,
    bonds: Vec<Bond>,
}

#[pymethods]
impl Substituent {
    #[new]
    fn new(atomic_numbers: Vec<u8>, bonds: Vec<Bond>) -> Self {
        Self {
            atomic_numbers,
            bonds,
        }
    }
}

#[pyclass]
struct SubstitutedMolecule {
    #[pyo3(get)]
    atomic_numbers: Vec<u8>,
    #[pyo3(get)]
    bonds: Vec<Bond>,
}

fn substitute_at_index(
    skeleton: &Molecule,
    substitutent: &Substituent,
    id: usize,
) -> SubstitutedMolecule {
    let mut atomic_numbers = skeleton.atomic_numbers.clone();
    atomic_numbers.extend(substitutent.atomic_numbers.iter());
    let mut bonds = skeleton.bonds.clone();
    let num_atoms = skeleton.atomic_numbers.len() as u16;
    bonds.extend(substitutent.bonds.iter().map(|bond| Bond {
        atom1_idx: bond.atom1_idx + num_atoms,
        atom2_idx: bond.atom2_idx + num_atoms,
        ..*bond
    }));
    bonds.push(Bond {
        atom1_idx: id as u16,
        atom2_idx: num_atoms as u16,
        order: 1.0,
    });
    SubstitutedMolecule {
        atomic_numbers,
        bonds,
    }
}

#[pyfunction]
fn substitute_1(
    skeleton: &Molecule,
    substituents: Vec<PyRef<Substituent>>,
) -> Vec<SubstitutedMolecule> {
    substituents
        .iter()
        .flat_map(|substituent| {
            skeleton
                .atomic_numbers
                .iter()
                .zip(skeleton.num_implicit_hs.iter())
                .enumerate()
                .filter(|(_, (&atomic_number, &num_implicit_hs))| {
                    atomic_number == 6 && num_implicit_hs > 0
                })
                .map(|(id, _)| substitute_at_index(skeleton, substituent, id))
        })
        .collect()
}

/// A Python module implemented in Rust.
#[pymodule]
fn jucombinator(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(substitute_1, m)?)?;
    m.add_class::<Bond>()?;
    m.add_class::<Molecule>()?;
    m.add_class::<Substituent>()?;
    m.add_class::<SubstitutedMolecule>()?;
    Ok(())
}
