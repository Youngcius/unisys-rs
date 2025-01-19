use crate::basic::circuits::Circuit;
use crate::basic::gates::Gate;
use crate::basic::operations::Operation;
use crate::{i, r};
use itertools::izip;
use ndarray::{array, s, Array1, Array2};
use ndarray_linalg::c64;
use std::sync::LazyLock;
use std::vec;

pub const I: LazyLock<Array2<c64>> =
    LazyLock::new(|| array![[r!(1.0), r!(0.0)], [r!(0.0), r!(1.0)]]);
pub const X: LazyLock<Array2<c64>> =
    LazyLock::new(|| array![[r!(0.0), r!(1.0)], [r!(1.0), r!(0.0)]]);
pub const Y: LazyLock<Array2<c64>> =
    LazyLock::new(|| array![[r!(0.0), -i!(1.0)], [i!(1.0), r!(0.0)]]);
pub const Z: LazyLock<Array2<c64>> =
    LazyLock::new(|| array![[r!(1.0), r!(0.0)], [r!(0.0), r!(-1.0)]]);

pub fn pauli_to_bsf_vec(pauli: &String) -> Array1<i8> {
    // e.g. 'XXX' --> [1, 1, 1, 0, 0, 0]
    // e.g. 'XXZ' --> [1, 1, 0, 0, 0, 1]
    // e.g. 'XXY' --> [1, 1, 1, 0, 0, 1]
    let pauli = pauli.to_uppercase();
    let mut vec = Array1::<i8>::zeros(pauli.len() * 2);
    for (i, c) in pauli.chars().enumerate() {
        match c {
            'I' => {
                vec[i] = 0;
                vec[i + pauli.len()] = 0;
            }
            'X' => {
                vec[i] = 1;
                vec[i + pauli.len()] = 0;
            }
            'Y' => {
                vec[i] = 1;
                vec[i + pauli.len()] = 1;
            }
            'Z' => {
                vec[i] = 0;
                vec[i + pauli.len()] = 1;
            }
            _ => panic!("Invalid Pauli string"),
        }
    }
    vec
}

pub fn bsf_vec_to_pauli(vec: &Array1<i8>) -> String {
    let len = vec.len();
    assert!(len % 2 == 0, "Input vector length must be even.");
    let half_len = len / 2;

    // 左右两部分
    let left = &vec.slice(s![..half_len]);
    let right = &vec.slice(s![half_len..]);

    // 映射规则: (left, right) -> 'I', 'Z', 'X', 'Y'
    let pauli_map = ['I', 'Z', 'X', 'Y']; // 对应 (0,0)->I, (0,1)->Z, (1,0)->X, (1,1)->Y,

    // 初始化结果字符串，预分配容量
    let mut pauli = String::with_capacity(half_len);

    // 遍历每个位置，生成对应的 Pauli 字符
    for (l, r) in left.iter().zip(right.iter()) {
        let index = ((*l as usize) << 1) | (*r as usize); // 计算索引
        pauli.push(pauli_map[index]); // 添加字符
    }

    pauli
}

pub fn paulis_to_bsf_mat(paulis: &Vec<String>) -> Array2<i8> {
    let mut mat = Array2::<i8>::zeros((paulis.len(), paulis[0].len() * 2));
    for (i, pauli) in paulis.iter().enumerate() {
        mat.slice_mut(s![i, ..]).assign(&pauli_to_bsf_vec(pauli));
    }
    mat
}

pub fn bsf_mat_to_paulis(mat: &Array2<i8>) -> Vec<String> {
    let mut paulis = Vec::new();
    for i in 0..mat.nrows() {
        paulis.push(bsf_vec_to_pauli(&mat.slice(s![i, ..]).to_owned()));
    }
    paulis
}

// binary simplectic form of Pauli strings
pub struct BSF {
    pub mat: Array2<i8>,
    pub coeffs: Option<Vec<f64>>,
    pub signs: Option<Vec<i8>>,
}

impl BSF {
    pub fn new(paulis: Vec<String>, coeffs: Option<Vec<f64>>, signs: Option<Vec<i8>>) -> Self {
        let n_paulis = paulis.len();
        if let Some(coeffs) = &coeffs {
            if coeffs.len() != n_paulis {
                panic!("The number of coefficients must be equal to the number of Paulis");
            }
        }
        if let Some(signs) = &signs {
            if signs.len() != n_paulis {
                panic!("The number of signs must be equal to the number of Paulis");
            }
        }
        BSF {
            mat: paulis_to_bsf_mat(&paulis),
            coeffs,
            signs,
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.mat.ncols() / 2
    }

    pub fn num_paulis(&self) -> usize {
        self.mat.nrows()
    }

    pub fn paulis(&self) -> Vec<String> {
        bsf_mat_to_paulis(&self.mat)
    }

    pub fn x(&self) -> Array2<i8> {
        // return the X-part of BSF tableau
        self.mat.slice(s![.., 0..self.num_qubits()]).to_owned()
    }

    pub fn z(&self) -> Array2<i8> {
        // return the Z-part of BSF tableau
        self.mat.slice(s![.., self.num_qubits()..]).to_owned()
    }

    pub fn with_ops(&self) -> Array2<i8> {
        // return the XOR of X-part and the Z-part
        self.x() | &self.z()
    }

    pub fn reverse(&self) -> Self {
        // reverse orders of tableau
        let mut mat = self.mat.clone();
        mat.slice_mut(s![.., ..]).invert_axis(ndarray::Axis(0));

        // reverse coefficients
        let coeffs = if let Some(coeffs) = &self.coeffs {
            let mut reversed_coeffs = coeffs.clone();
            reversed_coeffs.reverse();
            Some(reversed_coeffs)
        } else {
            None
        };
        // reverse signs
        let signs = if let Some(signs) = &self.signs {
            let mut reversed_signs = signs.clone();
            reversed_signs.reverse();
            Some(reversed_signs)
        } else {
            None
        };

        BSF { mat, coeffs, signs }
    }

    pub fn as_cnot_circuit(&self) -> Circuit {
        let mut circ = Circuit::new();

        let paulis = self.paulis();

        println!("paulis: {:?}", paulis);

        let coeffs = self
            .coeffs
            .clone()
            .unwrap_or_else(|| vec![1.0; self.num_paulis()]);
        let signs = self
            .signs
            .clone()
            .unwrap_or_else(|| vec![0; self.num_paulis()]);

        println!("coeffs: {:?}", coeffs);
        println!("signs: {:?}", signs);

        for (pauli, coeff, sign) in izip!(paulis, coeffs, signs) {
            let mut theta = coeff * 2.0_f64;
            if sign % 2 == 1 {
                theta *= -1.0;
            }

            // indices = np.where(np.array(list(paulistr)) != 'I')[0].tolist()

            let indices = pauli
                .chars()
                .enumerate()
                .filter(|(_, c)| *c != 'I')
                .map(|(i, _)| i)
                .collect::<Vec<_>>();
            if indices.len() == 0 {
                continue;
            } else if indices.len() == 1 {
                match pauli.chars().nth(indices[0]).unwrap() {
                    'X' => {
                        // circ.append(Operation::new(Gate::rx(theta), indices, None));
                        circ.append(Operation::new(Gate::h(), indices.clone(), None));
                        circ.append(Operation::new(Gate::rz(theta), indices.clone(), None));
                        circ.append(Operation::new(Gate::h(), indices, None));
                    }
                    'Y' => {
                        // circ.append(Operation::new(Gate::ry(theta), indices, None));
                        circ.append(Operation::new(Gate::sdg(), indices.clone(), None));
                        circ.append(Operation::new(Gate::h(), indices.clone(), None));
                        circ.append(Operation::new(Gate::s(), indices.clone(), None));
                        circ.append(Operation::new(Gate::rz(theta), indices.clone(), None));
                        circ.append(Operation::new(Gate::sdg(), indices.clone(), None));
                        circ.append(Operation::new(Gate::h(), indices.clone(), None));
                        circ.append(Operation::new(Gate::s(), indices, None));
                    }
                    'Z' => {
                        circ.append(Operation::new(Gate::rz(theta), indices, None));
                    }
                    _ => panic!("Invalid Pauli string"),
                }
            } else if indices.len() == 2 {
                // single-qubit Clifford (left)
                for &idx in indices.iter() {
                    match pauli.chars().nth(idx).unwrap() {
                        'X' => {
                            circ.append(Operation::new(Gate::h(), vec![idx], None));
                        }
                        'Y' => {
                            circ.append(Operation::new(Gate::sdg(), vec![idx], None));
                            circ.append(Operation::new(Gate::h(), vec![idx], None));
                            circ.append(Operation::new(Gate::s(), vec![idx], None));
                        }
                        'Z' => {}
                        _ => panic!("Invalid Pauli string"),
                    }
                }

                // CNOT tree (left)
                for i in 0..indices.len() - 1 {
                    circ.append(Operation::new(
                        Gate::x(),
                        vec![indices[i + 1]],
                        Some(vec![indices[i]]),
                    ));
                }

                // Rz rotation
                circ.append(Operation::new(
                    Gate::rz(theta),
                    vec![*indices.last().unwrap()],
                    None,
                ));

                // CNOT tree (right)
                for i in 0..indices.len() - 1 {
                    circ.append(Operation::new(
                        Gate::x(),
                        vec![indices[i + 1]],
                        Some(vec![indices[i]]),
                    ));
                }
                // single-qubit Clifford (right)
                for &idx in indices.iter() {
                    match pauli.chars().nth(idx).unwrap() {
                        'X' => {
                            circ.append(Operation::new(Gate::h(), vec![idx], None));
                        }
                        'Y' => {
                            circ.append(Operation::new(Gate::sdg(), vec![idx], None));
                            circ.append(Operation::new(Gate::h(), vec![idx], None));
                            circ.append(Operation::new(Gate::s(), vec![idx], None));
                        }
                        'Z' => {}
                        _ => panic!("Invalid Pauli string"),
                    }
                }
            } else {
                panic!("Only support the case of weight-1 and weight-2 Pauli strings");
            }
        }

        circ
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constr_bsf() {
        let paulis = vec!["XXX".to_string(), "XYY".to_string(), "XZZ".to_string()];

        let bsf = BSF::new(paulis, None, None);

        assert_eq!(bsf.num_qubits(), 3);
        assert_eq!(bsf.num_paulis(), 3);

        println!("Tableau:");
        println!("{}", bsf.mat);
        println!("X-part:");
        println!("{}", bsf.x());
        println!("Z-part:");
        println!("{}", bsf.z());
        println!("With ops:");
        println!("{}", bsf.with_ops());
        println!("Paulis:");
        println!("{:?}", bsf.paulis());
        assert_eq!(bsf.paulis(), vec!["XXX", "XYY", "XZZ"]);
    }

    #[test]
    fn test_bsf_to_circuit() {
        let paulis = vec!["XXI".to_string(), "IYY".to_string(), "ZIZ".to_string()];
        let coeffs = vec![1.0, 0.5, 0.25];
        let signs = vec![0, 1, 0];

        let bsf = BSF::new(paulis, Some(coeffs), Some(signs));

        let circ = bsf.as_cnot_circuit();
        println!("{}", circ.qasm());
    }

    #[test]
    fn test_pauli_matrices() {
        println!("sigmai:\n{}", *I);
        println!("sigmax:\n{}", *X);
        println!("sigmay:\n{}", *Y);
        println!("sigmaz:\n{}", *Z);
    }
}
