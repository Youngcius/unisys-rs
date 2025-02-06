use crate::basic::circuits::Circuit;
use crate::basic::gates::{Clifford1Q, Clifford2Q, Gate};
use crate::basic::matrices::Kronecker;
use crate::basic::operations::Operation;
use crate::{i, r};
use itertools::izip;
use ndarray::{array, s, Array1, Array2, Axis, Zip};
use ndarray_linalg::c64;
use std::collections::HashMap;
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
pub const II: LazyLock<Array2<c64>> = LazyLock::new(|| I.clone().kron(&I.clone()));
pub const IX: LazyLock<Array2<c64>> = LazyLock::new(|| I.clone().kron(&X.clone()));
pub const IY: LazyLock<Array2<c64>> = LazyLock::new(|| I.clone().kron(&Y.clone()));
pub const IZ: LazyLock<Array2<c64>> = LazyLock::new(|| I.clone().kron(&Z.clone()));
pub const XI: LazyLock<Array2<c64>> = LazyLock::new(|| X.clone().kron(&I.clone()));
pub const XX: LazyLock<Array2<c64>> = LazyLock::new(|| X.clone().kron(&X.clone()));
pub const XY: LazyLock<Array2<c64>> = LazyLock::new(|| X.clone().kron(&Y.clone()));
pub const XZ: LazyLock<Array2<c64>> = LazyLock::new(|| X.clone().kron(&Z.clone()));
pub const YI: LazyLock<Array2<c64>> = LazyLock::new(|| Y.clone().kron(&I.clone()));
pub const YX: LazyLock<Array2<c64>> = LazyLock::new(|| Y.clone().kron(&X.clone()));
pub const YY: LazyLock<Array2<c64>> = LazyLock::new(|| Y.clone().kron(&Y.clone()));
pub const YZ: LazyLock<Array2<c64>> = LazyLock::new(|| Y.clone().kron(&Z.clone()));
pub const ZI: LazyLock<Array2<c64>> = LazyLock::new(|| Z.clone().kron(&I.clone()));
pub const ZX: LazyLock<Array2<c64>> = LazyLock::new(|| Z.clone().kron(&X.clone()));
pub const ZY: LazyLock<Array2<c64>> = LazyLock::new(|| Z.clone().kron(&Y.clone()));
pub const ZZ: LazyLock<Array2<c64>> = LazyLock::new(|| Z.clone().kron(&Z.clone()));

pub const SinglePaulis: LazyLock<HashMap<String, Array2<c64>>> = LazyLock::new(|| {
    let mut m = HashMap::new();
    m.insert("I".to_string(), I.clone());
    m.insert("X".to_string(), X.clone());
    m.insert("Y".to_string(), Y.clone());
    m.insert("Z".to_string(), Z.clone());
    m
});

pub const DoublePaulis: LazyLock<HashMap<String, Array2<c64>>> = LazyLock::new(|| {
    let mut m = HashMap::new();
    m.insert("II".to_string(), II.clone());
    m.insert("IX".to_string(), IX.clone());
    m.insert("IY".to_string(), IY.clone());
    m.insert("IZ".to_string(), IZ.clone());
    m.insert("XI".to_string(), XI.clone());
    m.insert("XX".to_string(), XX.clone());
    m.insert("XY".to_string(), XY.clone());
    m.insert("XZ".to_string(), XZ.clone());
    m.insert("YI".to_string(), YI.clone());
    m.insert("YX".to_string(), YX.clone());
    m.insert("YY".to_string(), YY.clone());
    m.insert("YZ".to_string(), YZ.clone());
    m.insert("ZI".to_string(), ZI.clone());
    m.insert("ZX".to_string(), ZX.clone());
    m.insert("ZY".to_string(), ZY.clone());
    m.insert("ZZ".to_string(), ZZ.clone());
    m
});

pub fn nontrivial_indices(pauli: &String) -> Vec<usize> {
    pauli
        .chars()
        .enumerate()
        .filter(|(_, c)| *c != 'I')
        .map(|(i, _)| i)
        .collect()
}

pub fn pauli_to_matrix(pauli: &String) -> Array2<c64> {
    let pauli = pauli.to_uppercase();
    pauli.chars().fold(Array2::eye(1), |acc, c| match c {
        'I' => acc.kron(&I.clone()),
        'X' => acc.kron(&X.clone()),
        'Y' => acc.kron(&Y.clone()),
        'Z' => acc.kron(&Z.clone()),
        _ => panic!("Invalid Pauli string"),
    })
}

pub fn pauli_to_bsf_vec(pauli: &String) -> Array1<u8> {
    // e.g. 'XXX' --> [1, 1, 1, 0, 0, 0]
    // e.g. 'XXZ' --> [1, 1, 0, 0, 0, 1]
    // e.g. 'XXY' --> [1, 1, 1, 0, 0, 1]
    let pauli = pauli.to_uppercase();
    let mut vec = Array1::<u8>::zeros(pauli.len() * 2);
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

pub fn bsf_vec_to_pauli(vec: &Array1<u8>) -> String {
    let len = vec.len();
    assert!(len % 2 == 0, "Input vector length must be even.");
    let half_len = len / 2;

    // left and right halves of the vector
    let left = &vec.slice(s![..half_len]);
    let right = &vec.slice(s![half_len..]);

    // mapping rules: (left, right) -> 'I', 'Z', 'X', 'Y'
    let pauli_map = ['I', 'Z', 'X', 'Y']; // (0,0)->I, (0,1)->Z, (1,0)->X, (1,1)->Y,

    // initialize the result string with pre-allocated capacity
    let mut pauli = String::with_capacity(half_len);

    // traverse each position and generate the corresponding Pauli character
    for (l, r) in left.iter().zip(right.iter()) {
        let index = ((*l as usize) << 1) | (*r as usize);
        pauli.push(pauli_map[index]);
    }

    pauli
}

pub fn paulis_to_bsf_mat(paulis: &Vec<String>) -> Array2<u8> {
    let mut mat = Array2::<u8>::zeros((paulis.len(), paulis[0].len() * 2));
    for (i, pauli) in paulis.iter().enumerate() {
        mat.slice_mut(s![i, ..]).assign(&pauli_to_bsf_vec(pauli));
    }
    mat
}

pub fn bsf_mat_to_paulis(mat: &Array2<u8>) -> Vec<String> {
    let mut paulis = Vec::new();
    for i in 0..mat.nrows() {
        paulis.push(bsf_vec_to_pauli(&mat.slice(s![i, ..]).to_owned()));
    }
    paulis
}

#[derive(Clone, Debug)]
pub struct BSF {
    pub mat: Array2<u8>,
    pub coeffs: Option<Array1<f64>>,
    pub signs: Option<Array1<u8>>,
}

impl BSF {
    pub fn new(
        paulis: Vec<String>,
        coeffs: Option<Array1<f64>>,
        signs: Option<Array1<u8>>,
    ) -> Self {
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

    pub fn x(&self) -> Array2<u8> {
        // return the X-part of BSF tableau
        self.mat.slice(s![.., 0..self.num_qubits()]).to_owned()
    }

    pub fn z(&self) -> Array2<u8> {
        // return the Z-part of BSF tableau
        self.mat.slice(s![.., self.num_qubits()..]).to_owned()
    }

    pub fn weights(&self) -> Array2<usize> {
        let x_weights = self.x().sum_axis(Axis(1)).mapv(|w| w as usize);
        let z_weights = self.z().sum_axis(Axis(1)).mapv(|w| w as usize);
        ndarray::stack(Axis(1), &[x_weights.view(), z_weights.view()]).unwrap()
    }

    pub fn num_nonlocal_paulis(&self) -> usize {
        self.with_ops()
            .sum_axis(Axis(1))
            .mapv(|x| x > 1)
            .iter()
            .filter(|&&x| x)
            .count()
    }

    pub fn num_local_paulis(&self) -> usize {
        self.with_ops()
            .sum_axis(Axis(1))
            .mapv(|x| x <= 1)
            .iter()
            .filter(|&&x| x)
            .count()
    }

    pub fn which_nonlocal_paulis(&self) -> Vec<usize> {
        self.with_ops()
            .sum_axis(Axis(1))
            .indexed_iter()
            .filter_map(|(i, &x)| if x > 1 { Some(i) } else { None })
            .collect()
    }

    pub fn which_local_paulis(&self) -> Vec<usize> {
        self.with_ops()
            .sum_axis(Axis(1))
            .indexed_iter()
            .filter_map(|(i, &x)| if x <= 1 { Some(i) } else { None })
            .collect()
    }

    pub fn total_weight(&self) -> usize {
        if self.num_paulis() == 0 {
            return 0;
        }

        if self
            .paulis()
            .iter()
            .all(|p| p.chars().filter(|&c| c == 'I').count() == self.num_qubits())
        {
            return 0;
        }

        if self.num_nonlocal_paulis() == 0 {
            return 1;
        }

        // let selected_ops = self.with_ops().select(Axis(0), &self.which_nonlocal_paulis());
        // let mut combined = selected_ops.row(0).to_owned();
        // for row in selected_ops.rows().into_iter().skip(1) {
        //     combined = combined.mapv2(|a, b| a | b, &row);
        // }

        // use ndarray "fold" iteration to rewrite the above code
        let combined = self
            .with_ops()
            .select(Axis(0), &self.which_nonlocal_paulis())
            .fold_axis(Axis(0), 0, |acc, &x| acc | x);

        combined.sum() as usize
    }

    pub fn with_ops(&self) -> Array2<u8> {
        // return the XOR of X-part and the Z-part
        self.x() | &self.z()
    }

    pub fn qubits_with_ops(&self) -> Vec<usize> {
        self.with_ops()
            .sum_axis(Axis(0))
            .indexed_iter()
            .filter_map(|(i, &x)| if x > 0 { Some(i) } else { None })
            .collect()
    }

    pub fn qubits_with_nonlocal_ops(&self) -> Vec<usize> {
        let mut qubits: indexmap::IndexSet<usize> = indexmap::IndexSet::new();
        for pauli in self.paulis() {
            if self.num_qubits() - pauli.chars().filter(|&c| c == 'I').count() > 1 {
                qubits.extend(
                    pauli
                        .chars()
                        .enumerate()
                        .filter(|(_, c)| *c != 'I')
                        .map(|(i, _)| i),
                );
            }
        }
        let mut qubits: Vec<usize> = qubits.into_iter().collect();
        qubits.sort();
        qubits
    }

    pub fn apply_h(&self, idx: usize) -> Self {
        let num_qubits = self.num_qubits();
        let mut bsf = self.clone();

        let (mut col_x, mut col_z) = bsf
            .mat
            .multi_slice_mut((s![.., idx], s![.., idx + num_qubits]));
        Zip::from(&mut col_x).and(&mut col_z).for_each(|x, z| {
            std::mem::swap(x, z);
        });

        let x_col = self.mat.column(idx);
        let z_col = self.mat.column(idx + num_qubits);
        // if let Some(signs) = &mut bsf.signs {
        //     signs
        //         .iter_mut()
        //         .zip(x_col.iter().zip(z_col.iter()))
        //         .for_each(|(sign, (&x, &z))| {
        //             *sign ^= x & z; // 异或赋值
        //         });
        // }
        // Now signs is modified from  Option<Vec<u8>> to Option<Array1<u8>>
        if let Some(signs) = &mut bsf.signs {
            Zip::from(signs)
                .and(&x_col)
                .and(&z_col)
                .for_each(|sign, &x, &z| {
                    *sign ^= x & z;
                });
        }

        bsf
    }

    pub fn apply_s(&self, idx: usize) -> Self {
        let num_qubits = self.num_qubits();
        let mut bsf = self.clone();

        let x_col = self.mat.column(idx);
        let z_col_mut = bsf.mat.column_mut(idx + num_qubits);

        // element-wise z = z ^ x
        Zip::from(z_col_mut).and(&x_col).for_each(|z, &x| {
            *z ^= x;
        });

        // update signs
        if let Some(signs) = &mut bsf.signs {
            let z_col = self.mat.column(idx + num_qubits); // 取原始 z 列进行符号计算
            Zip::from(signs)
                .and(&x_col)
                .and(&z_col)
                .for_each(|sign, &x, &z| {
                    *sign ^= x & z;
                });
        }

        bsf
    }

    pub fn apply_sdg(&self, idx: usize) -> Self {
        let num_qubits = self.num_qubits();
        let mut bsf = self.clone();

        let x_col = self.mat.column(idx);
        let z_col_mut = bsf.mat.column_mut(idx + num_qubits);

        // element-wise z = z ^ x
        Zip::from(z_col_mut).and(&x_col).for_each(|z, &x| {
            *z ^= x;
        });

        // update signs
        if let Some(signs) = &mut bsf.signs {
            let z_col = self.mat.column(idx + num_qubits); // 取原始 z 列进行符号计算
            Zip::from(signs)
                .and(&x_col)
                .and(&z_col)
                .for_each(|sign, &x, &z| {
                    *sign ^= ((x == 1) & (z == 0)) as u8;
                });
        }

        bsf
    }

    pub fn apply_cx(&self, ctrl: usize, targ: usize) -> Self {
        let num_qubits = self.num_qubits();
        let mut bsf = self.clone();

        // update mat
        let x_ctrl_col = self.mat.column(ctrl);
        let z_targ_col = self.mat.column(targ + num_qubits);

        let x_targ_col_mut = bsf.mat.column_mut(targ);
        Zip::from(x_targ_col_mut)
            .and(&x_ctrl_col)
            .for_each(|targ_x_mut, &ctrl_x| {
                *targ_x_mut ^= ctrl_x;
            });

        let z_ctrl_col_mut = bsf.mat.column_mut(ctrl + num_qubits);
        Zip::from(z_ctrl_col_mut)
            .and(&z_targ_col)
            .for_each(|ctrl_z_mut, &targ_z| {
                *ctrl_z_mut ^= targ_z;
            });

        // update signs
        if let Some(signs) = &mut bsf.signs {
            let x_ctrl_col = self.mat.column(ctrl);
            let x_targ_col = self.mat.column(targ);
            let z_ctrl_col = self.mat.column(ctrl + num_qubits);
            let z_targ_col = self.mat.column(targ + num_qubits);

            Zip::from(signs)
                .and(&x_ctrl_col)
                .and(&x_targ_col)
                .and(&z_ctrl_col)
                .and(&z_targ_col)
                .for_each(|sign, &x_ctrl, &x_targ, &z_ctrl, &z_targ| {
                    *sign ^= (((x_ctrl == 1) & (x_targ == 0) & (z_ctrl == 0) & (z_targ == 1))
                        as u8)
                        | (((x_ctrl == 1) & (x_targ == 1) & (z_ctrl == 1) & (z_targ == 1)) as u8);
                });
        }

        bsf
    }

    pub fn apply_cz(&self, ctrl: usize, targ: usize) -> Self {
        let num_qubits = self.num_qubits();
        let mut bsf = self.clone();

        // update mat
        let x_ctrl_col = self.mat.column(ctrl);
        let x_targ_col = self.mat.column(targ);
        let z_ctrl_col = self.mat.column(ctrl + num_qubits);
        let z_targ_col = self.mat.column(targ + num_qubits);

        let z_ctrl_col_mut = bsf.mat.column_mut(ctrl + num_qubits);

        Zip::from(z_ctrl_col_mut)
            .and(&x_targ_col)
            .for_each(|ctrl_z, &targ_x| {
                *ctrl_z ^= targ_x;
            });

        let z_targ_col_mut = bsf.mat.column_mut(targ + num_qubits);
        Zip::from(z_targ_col_mut)
            .and(&x_ctrl_col)
            .for_each(|targ_z, &ctrl_x| {
                *targ_z ^= ctrl_x;
            });

        // update signs
        if let Some(signs) = &mut bsf.signs {
            Zip::from(signs)
                .and(&x_ctrl_col)
                .and(&x_targ_col)
                .and(&z_ctrl_col)
                .and(&z_targ_col)
                .for_each(|sign, &x_ctrl, &x_targ, &z_ctrl, &z_targ| {
                    *sign ^= (((x_ctrl == 1) & (x_targ == 1) & (z_ctrl == 0) & (z_targ == 1))
                        as u8)
                        | (((x_ctrl == 1) & (x_targ == 1) & (z_ctrl == 1) & (z_targ == 0)) as u8);
                });
        }

        bsf
    }

    pub fn apply_clifford2q(&self, cliff: &Clifford2Q, ctrl: usize, targ: usize) -> Self {
        let (wrt0_cz, wrt1_cz) = cliff.wrt_cz();

        let mut bsf = self.clone();

        // e.g., Y = S H S† Z S H S†, wrt0_cz = {'S†', 'H', 'S'}
        //      {SDG, H, S, Z, SDG, H, S}
        for cliff1q in wrt0_cz.iter() {
            match cliff1q {
                Clifford1Q::H => bsf = bsf.apply_h(ctrl),
                Clifford1Q::S => bsf = bsf.apply_s(ctrl),
                Clifford1Q::SDG => bsf = bsf.apply_sdg(ctrl),
            }
        }

        for cliff1q in wrt1_cz.iter() {
            match cliff1q {
                Clifford1Q::H => bsf = bsf.apply_h(targ),
                Clifford1Q::S => bsf = bsf.apply_s(targ),
                Clifford1Q::SDG => bsf = bsf.apply_sdg(targ),
            }
        }

        bsf = bsf.apply_cz(ctrl, targ);

        for cliff1q in wrt0_cz.iter().rev() {
            match cliff1q {
                Clifford1Q::H => bsf = bsf.apply_h(ctrl),
                Clifford1Q::S => bsf = bsf.apply_s(ctrl),
                Clifford1Q::SDG => bsf = bsf.apply_sdg(ctrl),
            }
        }

        for cliff1q in wrt1_cz.iter().rev() {
            match cliff1q {
                Clifford1Q::H => bsf = bsf.apply_h(targ),
                Clifford1Q::S => bsf = bsf.apply_s(targ),
                Clifford1Q::SDG => bsf = bsf.apply_sdg(targ),
            }
        }

        bsf
    }

    pub fn reverse(&self) -> Self {
        // reverse orders of tableau
        let mut mat = self.mat.clone();
        mat.slice_mut(s![.., ..]).invert_axis(Axis(0));

        // reverse coefficients
        let coeffs = if let Some(coeffs) = &self.coeffs {
            let mut reversed_coeffs = coeffs.clone();
            reversed_coeffs.invert_axis(Axis(0));
            Some(reversed_coeffs)
        } else {
            None
        };
        // reverse signs
        let signs = if let Some(signs) = &self.signs {
            let mut reversed_signs = signs.clone();
            reversed_signs.invert_axis(Axis(0));
            Some(reversed_signs)
        } else {
            None
        };

        BSF { mat, coeffs, signs }
    }

    pub fn pop_local_paulis(&mut self) -> BSF {
        let which_local_paulis = self.which_local_paulis();
        let which_nonlocal_paulis = self.which_nonlocal_paulis();

        // compute popped local_bsf
        let mut local_bsf = self.clone();
        local_bsf.mat = local_bsf.mat.select(Axis(0), &which_local_paulis);
        if let Some(coeffs) = &mut local_bsf.coeffs {
            *coeffs = coeffs
                .iter()
                .enumerate()
                .filter_map(|(i, _)| {
                    if which_local_paulis.contains(&i) {
                        Some(coeffs[i])
                    } else {
                        None
                    }
                })
                .collect();
        }
        if let Some(signs) = &mut local_bsf.signs {
            *signs = signs
                .iter()
                .enumerate()
                .filter_map(|(i, _)| {
                    if which_local_paulis.contains(&i) {
                        Some(signs[i])
                    } else {
                        None
                    }
                })
                .collect();
        }

        // compute modified self
        self.mat = self.mat.select(Axis(0), &which_nonlocal_paulis);
        if let Some(coeffs) = &mut self.coeffs {
            *coeffs = coeffs
                .iter()
                .enumerate()
                .filter_map(|(i, _)| {
                    if which_nonlocal_paulis.contains(&i) {
                        Some(coeffs[i])
                    } else {
                        None
                    }
                })
                .collect();
        }
        if let Some(signs) = &mut self.signs {
            *signs = signs
                .iter()
                .enumerate()
                .filter_map(|(i, _)| {
                    if which_nonlocal_paulis.contains(&i) {
                        Some(signs[i])
                    } else {
                        None
                    }
                })
                .collect();
        }

        local_bsf
    }

    pub fn as_cnot_circuit(&self) -> Circuit {
        let mut circ = Circuit::new();

        let paulis = self.paulis();

        // println!("paulis: {:?}", paulis);

        let coeffs = self
            .coeffs
            .clone()
            .unwrap_or_else(|| Array1::from_elem(self.num_paulis(), 1.0));
        let signs = self
            .signs
            .clone()
            .unwrap_or_else(|| Array1::from_elem(self.num_paulis(), 0));

        // println!("coeffs: {:?}", coeffs);
        // println!("signs: {:?}", signs);

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
                        // herein we use Y = S H S† Z S H S†
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
        let coeffs = Array1::from_vec(vec![1.0, 0.5, 0.25]);
        let signs = Array1::from_vec(vec![0, 1, 0]);

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

    #[test]
    fn test_pauli_to_matrix() {
        let pauli = "XYZ".to_string();
        let matrix = pauli_to_matrix(&pauli);
        println!("{}", matrix);
    }

    #[test]
    fn test_nontrivial_indices() {
        let pauli = "XYZIZIX".to_string();
        let indices = nontrivial_indices(&pauli);
        println!("{:?}", indices);
    }

    #[test]
    fn test_apply_cz() {
        // E.g., ["II", "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX", "ZY", "ZZ"]
        //   --> ["II", 'ZX', 'ZY', 'IZ', 'XZ', 'YY', '-YX', 'XI', 'YZ', '-XY', 'XX', 'YI', 'ZI', 'IX', 'IY', 'ZZ']
        let paulis: Vec<String> = [
            "II", "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX",
            "ZY", "ZZ",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let signs: Array1<u8> = Array1::from_vec(vec![0; 16]);

        let paulis_desired: Vec<String> = [
            "II", "ZX", "ZY", "IZ", "XZ", "YY", "YX", "XI", "YZ", "XY", "XX", "YI", "ZI", "IX",
            "IY", "ZZ",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let signs_desired: Array1<u8> =
            Array1::from_vec(vec![0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0]);

        let bsf = BSF::new(paulis, None, Some(signs));
        let bsf_cz = bsf.apply_cz(0, 1);
        println!("Tableau:");
        println!("{}", bsf_cz.mat);
        println!("Paulis:");
        println!("{:?}", bsf_cz.paulis());
        assert!(bsf_cz.paulis() == paulis_desired);
        println!("Signs:");
        println!("{:?}", bsf_cz.signs);
        assert!(bsf_cz.signs.unwrap() == signs_desired);
    }

    #[test]
    fn test_apply_cx() {
        // E.g., ["II", "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX", "ZY", "ZZ"]
        //   --> ['II', 'IX', 'ZY', 'ZZ', 'XX', 'XI', 'YZ', '-YY', 'YX', 'YI', '-XZ', 'XY', 'ZI', 'ZX', 'IY', 'IZ']
        let paulis = [
            "II", "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX",
            "ZY", "ZZ",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let signs = Array1::from_vec(vec![0; 16]);

        let paulis_desired: Vec<String> = [
            "II", "IX", "ZY", "ZZ", "XX", "XI", "YZ", "YY", "YX", "YI", "XZ", "XY", "ZI", "ZX",
            "IY", "IZ",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let signs_desired: Array1<u8> =
            Array1::from_vec(vec![0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0]);

        let bsf = BSF::new(paulis, None, Some(signs));
        let bsf_cx = bsf.apply_cx(0, 1);
        println!("Tableau:");
        println!("{}", bsf_cx.mat);
        println!("Paulis:");
        println!("{:?}", bsf_cx.paulis());
        assert!(bsf_cx.paulis() == paulis_desired);
        println!("Signs:");
        println!("{:?}", bsf_cx.signs);
        assert!(bsf_cx.signs.unwrap() == signs_desired);
    }
}
