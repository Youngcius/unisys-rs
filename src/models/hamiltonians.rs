use super::paulis;
use super::paulis::BSF;
use crate::basic::circuits::Circuit;
use crate::phoenix::{
    grouping,
    simplification::simplify_bsf,
    utils::{config_to_circuit, SimpItem},
};
use ndarray::{Array1, Array2};
use ndarray_linalg::c64;

#[derive(Debug)]
pub struct Hamiltonian {
    paulis: Vec<String>,
    coeffs: Vec<f64>, // in Hamiltonian, coeffs is Vec<f64>; in BSF, coeffs is Option<Array1<f64>>
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

    pub fn norm(&self) -> f64 {
        self.coeffs.iter().map(|c| c.abs()).sum()
    }

    pub fn normalize(&self) -> Self {
        let norm = self.norm();
        let coeffs = self.coeffs.iter().map(|c| c / norm).collect();
        Hamiltonian::new(self.paulis.clone(), coeffs)
    }

    pub fn to_bsf(&self) -> BSF {
        BSF::new(
            self.paulis.clone(),
            Some(Array1::from_vec(self.coeffs.clone())),
            None,
        )
    }

    pub fn to_matrix(&self) -> Array2<c64> {
        let dim: usize = 1 << self.num_qubits;
        let mut result = Array2::<c64>::zeros((dim, dim));
        for (pauli, coeff) in self.paulis.iter().zip(self.coeffs.iter()) {
            result = result + paulis::pauli_to_matrix(pauli) * coeff.clone();
        }
        result
    }

    pub fn unitary_evolution(&self) -> Array2<c64> {
        panic!("Not implemented yet")
    }

    pub fn generate_circuit(&self) -> Circuit {
        // Currently only support 1-order Trotterization
        BSF::new(
            self.paulis.clone(),
            Some(Array1::from_vec(self.coeffs.clone())),
            None,
        )
        .as_cnot_circuit()
    }

    pub fn reconfigure(&self) -> Vec<SimpItem> {
        let groups = self.group_paulis_and_coeffs();
        let mut config = Vec::new();

        for (idx, (paulis, coeffs)) in groups.iter() {
            let mut res = Vec::new();
            let bsf = BSF::new(paulis.clone(), Some(Array1::from_vec(coeffs.clone())), None);
            let (bsf_, cliffords_with_locals) = simplify_bsf(&bsf);
            res.push(SimpItem::BSFItem(bsf_));
            for (cliff, local_bsf) in cliffords_with_locals.iter().rev() {
                res.insert(0, SimpItem::UCGItem(cliff.clone()));
                res.push(SimpItem::UCGItem(cliff.clone()));
                if local_bsf.num_paulis() > 0 {
                    res.push(SimpItem::BSFItem(local_bsf.clone()));
                }
            }
            config.extend(res);
        }

        config
    }

    pub fn reconfigure_and_generate_circuit(&self) -> Circuit {
        config_to_circuit(self.reconfigure())
    }

    pub fn group_paulis_and_coeffs(&self) -> Vec<(Vec<usize>, (Vec<String>, Vec<f64>))> {
        grouping::group_paulis_and_coeffs(&self.paulis, &self.coeffs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn large_paulis() -> Vec<String> {
        [
            "ZYZZXZZYIZ",
            "IYZYXYZIYX",
            "ZXIYIZYXXX",
            "XYIIZXIZXI",
            "XYXIYIXXXZ",
            "YYZZIZYXXI",
            "ZYXIZZZZYZ",
            "XIIIYXZXZY",
            "IZXZZIXIIZ",
            "IYIZYIYZXX",
            "YXZYYIXYYX",
            "XIYZIIIYYI",
            "YZIIIIXZZY",
            "IZZIYZXZIZ",
            "YYZZIXIYIY",
            "XYXYZZIIXY",
            "YYXIYYIXZI",
            "YZXYXZXIXI",
            "IZXZYXXIXX",
            "IXXYXXXZIZ",
        ]
        .iter()
        .map(|&x| x.to_string())
        .collect()
    }

    fn heisenberg_hamiltonian() -> Hamiltonian {
        let paulis = vec![
            "XXIII".to_string(),
            "YYIII".to_string(),
            "ZZIII".to_string(),
            "IXXII".to_string(),
            "IYYII".to_string(),
            "IZZII".to_string(),
            "IIXXI".to_string(),
            "IIYYI".to_string(),
            "IIZZI".to_string(),
            "IIIXX".to_string(),
            "IIIYY".to_string(),
            "IIIZZ".to_string(),
            "ZIIII".to_string(),
            "IZIII".to_string(),
            "IIZII".to_string(),
            "IIIZI".to_string(),
            "IIIIZ".to_string(),
        ];
        let coeffs = vec![
            0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0,
        ];

        Hamiltonian::new(paulis, coeffs)
    }

    #[test]
    fn test_heisenberg_circuit() {
        let ham = heisenberg_hamiltonian();
        let circ = ham.generate_circuit();
        for op in circ.ops.iter() {
            println!("{}", op);
        }
        let formatted_matrix = circ.unitary().row(0).mapv(|x| format!("{:.3}", x));
        println!("{}", formatted_matrix);
    }

    #[test]
    fn test_heisenberg_circuit_reconfigure() {
        let ham = heisenberg_hamiltonian();
        let config = ham.reconfigure();
        for item in config.iter() {
            println!("{}", item);
        }
    }

    #[test]
    fn test_large_paulis_reconfigure() {
        let start = std::time::Instant::now();
        let paulis = large_paulis();
        let coeffs = vec![0.1; paulis.len()];
        let ham = Hamiltonian::new(paulis, coeffs);
        let config = ham.reconfigure();
        for item in config.iter() {
            println!("{}", item);
        }
        println!("Elapsed time: {:?}", start.elapsed());
        println!("config length: {}", config.len());

        let circ = ham.reconfigure_and_generate_circuit();

        println!();
        println!("circ.num_ops: {}", circ.num_ops());
        println!("circ.num_nonlocal_ops: {}", circ.num_nonlocal_ops());
        println!("circ.depth: {}", circ.depth());
        println!("circ.depth_nonlocal: {}", circ.depth_nonlocal());
        println!("{:?}", circ.gate_stats());
        println!("Elapsed time: {:?}", start.elapsed());
    }
}
