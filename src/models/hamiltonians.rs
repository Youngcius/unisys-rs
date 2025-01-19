
use std::collections::HashMap;

use super::paulis::BSF;

// class HamiltonianModel:
//     def __init__(self, paulis, coeffs, *args, **kwargs):
//         self.paulis = list(paulis)
//         self.coeffs = np.array(coeffs)
//         self.num_qubits = len(paulis[0])

#[derive(Debug)]
pub struct Hamiltonian {
    paulis: Vec<String>,
    coeffs: Vec<f64>,
    num_qubits: usize,
}

impl Hamiltonian {
    pub fn new(paulis: Vec<String>, coeffs: Vec<f64>) -> Self {
        let num_qubits = paulis[0].len();
        for pauli in paulis.iter() {
            for c in pauli.chars() {
                assert!(c == 'I' || c == 'X' || c == 'Y' || c == 'Z');
            }
        }
        assert!(paulis.len() == coeffs.len());
        Hamiltonian {
            paulis,
            coeffs,
            num_qubits,
        }
    }

/*
    def norm(self) -> float:
        """Return the norm of the Hamiltonian, i.e., sum of all spectral norms of each iterm"""
        norm = 0
        for p, c in zip(self.paulis, self.coeffs):
            norm += np.abs(c) * linalg.norm(qi.Pauli(p).to_matrix(), 2)
        return norm

*/

    // pub fn norm(&self) -> f64 {


    // }

    pub fn to_bsf(&self) -> BSF {
        BSF::new(self.paulis.clone(), Some(self.coeffs.clone()), None)
    }

    pub fn group_paulis_and_coeffs(&self) {
        panic!("Not implemented yet");
    }

}
