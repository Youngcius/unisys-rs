// Extended matrix traits for ndarray
use super::gates::Gate;
use crate::utils::functions::is_power_of_two;
use crate::{c, i, r};
use ndarray::{array, s, Array};
use ndarray::{Array1, Array2, ArrayView2};
use ndarray_linalg::c64;
use ndarray_rand::rand_distr::Uniform;
use ndarray_rand::RandomExt;
use std::f64::consts::PI;

pub trait Real {
    fn real(&self) -> Array2<f64>;
}

impl Real for Array2<c64> {
    fn real(&self) -> Array2<f64> {
        self.mapv(|i| i.re)
    }
}

impl Real for ArrayView2<'_, c64> {
    fn real(&self) -> Array2<f64> {
        self.mapv(|i| i.re)
    }
}

pub trait Imag {
    fn imag(&self) -> Array2<f64>;
}

impl Imag for Array2<c64> {
    fn imag(&self) -> Array2<f64> {
        self.mapv(|i| i.im)
    }
}

impl Imag for ArrayView2<'_, c64> {
    fn imag(&self) -> Array2<f64> {
        self.mapv(|i| i.im)
    }
}

// conjugation of a complex matrix
pub trait Conj {
    fn conj(&self) -> Array2<c64>;
}

impl Conj for Array2<c64> {
    fn conj(&self) -> Array2<c64> {
        self.mapv(|i| i.conj())
    }
}

impl Conj for ArrayView2<'_, c64> {
    fn conj(&self) -> Array2<c64> {
        self.mapv(|i| i.conj())
    }
}

// dagger (hermitian) of a complex matrix
pub trait Dagger {
    fn dagger(&self) -> Array2<c64>;
}

impl Dagger for Array2<c64> {
    fn dagger(&self) -> Array2<c64> {
        self.mapv(|i| i.conj()).reversed_axes()
    }
}

impl<'a> Dagger for ArrayView2<'a, c64> {
    fn dagger(&self) -> Array2<c64> {
        self.mapv(|i| i.conj()).reversed_axes()
    }
}

pub trait Kronecker {
    fn kron(&self, other: &Array2<c64>) -> Array2<c64>;
}

impl Kronecker for Array2<c64> {
    fn kron(&self, other: &Array2<c64>) -> Array2<c64> {
        kronecker_product(self, other)
    }
}

// Kronecker product of two complex matrices
fn kronecker_product(a: &Array2<c64>, b: &Array2<c64>) -> Array2<c64> {
    let (a_rows, a_cols) = a.dim();
    let (b_rows, b_cols) = b.dim();

    // Resulting matrix dimensions
    let result_rows = a_rows * b_rows;
    let result_cols = a_cols * b_cols;

    // Initialize the resulting matrix with zeros
    let mut result = Array2::zeros((result_rows, result_cols));

    // Compute the Kronecker product
    for i in 0..a_rows {
        for j in 0..a_cols {
            let a_value = a[[i, j]];

            // Compute the block corresponding to a[i, j] * B
            let sub_matrix = b.mapv(|b_value| a_value * b_value);

            // Place the block into the resulting matrix
            let row_start = i * b_rows;
            let col_start = j * b_cols;

            result
                .slice_mut(s![
                    row_start..row_start + b_rows,
                    col_start..col_start + b_cols
                ])
                .assign(&sub_matrix);
        }
    }

    result
}

pub fn random_hermitian(d: usize) -> Array2<c64> {
    let real_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
    let imag_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
    let mat: Array2<c64> = real_part.mapv(|re| r!(re)) + imag_part.mapv(|im| i!(im));
    (mat.clone() + mat.dagger()) / 2.0
}

// import numpy as np

// # 自定义采样器，用 sin(theta) 的概率密度函数
// def sin_sampler(size):
//     """
//     根据概率密度函数 f(theta) = 0.5 * sin(theta) 进行采样
//     使用逆变换采样法生成随机变量
//     """
//     # 计算累积分布函数 (CDF) 的反函数
//     # CDF: F(theta) = -0.5 * cos(theta) + 0.5
//     # Inverse CDF: F^(-1)(u) = arccos(1 - 2 * u)
//     uniform_samples = np.random.uniform(size=size)
//     theta_samples = np.arccos(1 - 2 * uniform_samples)
//     return theta_samples

fn sin_sampler(size: usize) -> Array1<f64> {
    let uniform_samples = Array::random(size, Uniform::new(0.0, 1.0));
    let theta_samples = uniform_samples.mapv(|u| (1.0_f64 - 2.0_f64 * u).acos());
    theta_samples
}

pub fn random_su2() -> Array2<c64> {
    // Sampling SU(2) based on Haar random measure
    // Reference: https://pennylane.ai/qml/demos/tutorial_haar_measure
    let angles = Array::random((2,), Uniform::new(0.0, 2.0 * PI));
    let (phi, lambda) = (angles[0], angles[1]);
    let theta = sin_sampler(1)[0];
    Gate::u3(theta, phi, lambda).data
}

// pub fn random_unitary(d: usize) -> Array2<c64> {
//     // let real_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
//     // let imag_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
//     // let mat: Array2<c64> = real_part.mapv(|re| r!(re)) + imag_part.mapv(|im| i!(im));

//     // let (q, _r) = mat.qr().unwrap(); // error in ".qr()": Method `qr` not found in the current scope for type `Array2<Complex64>` [E0599]
//     // q

//     // ! Alternatively, we use another construction method
//     // ! Reference: https://pennylane.ai/qml/demos/tutorial_haar_measure
//     let mut mat = Array2::zeros((d, d));
//     for i in 0..d {
//         for j in 0..d {
//             mat[[i, j]] = random_su2()[[0, 0]];
//         }
//     }
//     mat
// }

// pub fn random_su4() -> Array2<c64> {
//     let coord = Array::random((3,), Uniform::new(0.0, 1.0));
//     let can = Gate::can(coord[0], coord[1], coord[2]).data;

//     let a1_params = Array::random((3,), Uniform::new(0.0, 2.0 * PI));
//     let a2_params = Array::random((3,), Uniform::new(0.0, 2.0 * PI));
//     let b1_params = Array::random((3,), Uniform::new(0.0, 2.0 * PI));
//     let b2_params = Array::random((3,), Uniform::new(0.0, 2.0 * PI));

//     let a1 = Gate::u3(a1_params[0], a1_params[1], a1_params[2]).data;
//     let a2 = Gate::u3(a2_params[0], a2_params[1], a2_params[2]).data;
//     let b1 = Gate::u3(b1_params[0], b1_params[1], b1_params[2]).data;
//     let b2 = Gate::u3(b2_params[0], b2_params[1], b2_params[2]).data;

//     kronecker_product(&a1, &a2)
//         .dot(&can)
//         .dot(&kronecker_product(&b1, &b2))
// }

pub fn controlled_unitary_matrix(u: &Array2<c64>, num_ctrl: usize) -> Array2<c64> {
    let proj_0 = Array2::from_diag(&array![r!(1.0), r!(0.0)]);
    let proj_1 = Array2::from_diag(&array![r!(0.0), r!(1.0)]);

    let mut u = u.clone();
    for _ in 0..num_ctrl {
        let ident = Array2::eye(1 << (u.dim().0 as f64).log2() as usize);
        u = proj_0.kron(&ident) + proj_1.kron(&u);
    }
    u
}

pub fn tensor_1_slot(u: &Array2<c64>, n: usize, tq: usize) -> Array2<c64> {
    if tq >= n {
        panic!("the qubit index is out of range");
    }

    let mut ops = vec![Array2::eye(2); n];
    ops[tq] = u.clone();
    ops.iter().fold(Array2::eye(1), |acc, x| acc.kron(x))
}

pub fn tensor_slots(u: &Array2<c64>, n: usize, indices: &[usize]) -> Array2<c64> {
    let (d1, d2) = u.dim();
    if d1 != d2 {
        panic!("U must be a square matrix!");
    }

    if !is_power_of_two(d1) {
        panic!(
            "Dimension of input matrix need to be power of 2, but got {}",
            d1
        );
    }
    let m = (d1 as f64).log2() as usize;
    if indices.len() != m {
        panic!("Length of indices does not match log2(dim(U))");
    }
    if indices.iter().max().unwrap_or(&0) >= &n {
        panic!("Some qubit index is out of range");
    }

    if m == 1 {
        tensor_1_slot(u, n, indices[0])
    } else {
        let mut arr_list: Vec<Array2<c64>> = Vec::with_capacity(n);
        arr_list.push(u.clone());
        for _ in 0..(n - m) {
            arr_list.push(Array2::eye(2));
        }

        let mut iter = arr_list.into_iter();
        let first = iter.next().unwrap();
        let res = iter.fold(first, |acc, x| kronecker_product(&acc, &x));

        // Now res is a matrix of dimension (2^n, 2^n), reshape it 2n times (i.e. shape is [2;2*n])
        let full_shape = vec![2; 2 * n];
        let res = res
            .into_shape_with_order(full_shape)
            .expect("reshape error");

        let mut idx = vec![-1; n];
        for (i, &k) in indices.iter().enumerate() {
            idx[k] = i as isize;
        }

        let mut fill_range = m;
        for x in idx.iter_mut() {
            if *x < 0 {
                *x = fill_range as isize;
                fill_range += 1;
            }
        }
        let idx_latter: Vec<isize> = idx.iter().map(|&i| i + n as isize).collect();

        // idx_combined = idx + idx_latter
        let mut idx_combined = idx.clone();
        idx_combined.extend(idx_latter);

        // Convert isize to usize and check non-negative
        let idx_usize: Vec<usize> = idx_combined.into_iter().map(|x| x as usize).collect();

        // Permute the axes by permuted_axes
        let permuted = res.permuted_axes(idx_usize);

        // extract all elements from permuted
        let data: Vec<c64> = permuted.iter().cloned().collect();
        Array2::from_shape_vec((1 << n, 1 << n), data).expect("from_shape_vec error")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::gates;
    use crate::utils::ops::allclose;
    use ndarray::{array, Array};
    use ndarray_linalg::assert;
    // use ndarray_rand::rand_distr::num_traits::zero;

    #[test]
    fn test_kronecker_product() {
        // Define two complex matrices
        let a = array![[c!(1.0, 0.0), c!(2.0, 1.0)], [c!(0.0, -1.0), c!(1.0, 1.0)]];

        let b = array![[c!(0.0, 1.0), c!(1.0, 0.0)], [c!(1.0, 1.0), c!(0.0, 0.0)]];

        // Desired result
        let desired = array![
            [c!(0.0, 1.0), c!(1.0, 0.0), c!(-1.0, 2.0), c!(2.0, 1.0)],
            [c!(1.0, 1.0), c!(0.0, 0.0), c!(1.0, 3.0), c!(0.0, 0.0)],
            [c!(1.0, 0.0), c!(0.0, -1.0), c!(-1.0, 1.0), c!(1.0, 1.0)],
            [c!(1.0, -1.0), c!(0.0, 0.0), c!(0.0, 2.0), c!(0.0, 0.0)]
        ];

        // Compute their Kronecker product
        let result = kronecker_product(&a, &b);

        assert!(allclose(&result, &desired));
    }

    #[test]
    fn test_tensor_1_slot() {
        let u = gates::Gate::x().data;
        let n = 3;
        let tq = 1;
        let zero = r!(0.0);
        let one = r!(1.0);
        let res = tensor_1_slot(&u, n, tq);
        let desired = array![
            [zero, zero, one, zero, zero, zero, zero, zero],
            [zero, zero, zero, one, zero, zero, zero, zero],
            [one, zero, zero, zero, zero, zero, zero, zero],
            [zero, one, zero, zero, zero, zero, zero, zero],
            [zero, zero, zero, zero, zero, zero, one, zero],
            [zero, zero, zero, zero, zero, zero, zero, one],
            [zero, zero, zero, zero, one, zero, zero, zero],
            [zero, zero, zero, zero, zero, one, zero, zero]
        ];
        assert!(allclose(&res, &desired));
    }

    #[test]
    fn test_tensor_slots() {
        let u = gates::Gate::swap().data;
        let n = 3;
        let indices = vec![0, 2];
        let zero = r!(0.0);
        let one = r!(1.0);
        let res = tensor_slots(&u, n, &indices);
        println!("{}", res.real());
        let desired = array![
            [one, zero, zero, zero, zero, zero, zero, zero],
            [zero, zero, zero, zero, one, zero, zero, zero],
            [zero, zero, one, zero, zero, zero, zero, zero],
            [zero, zero, zero, zero, zero, zero, one, zero],
            [zero, one, zero, zero, zero, zero, zero, zero],
            [zero, zero, zero, zero, zero, one, zero, zero],
            [zero, zero, zero, one, zero, zero, zero, zero],
            [zero, zero, zero, zero, zero, zero, zero, one]
        ];
        println!();
        println!("{}", desired.real());
        assert!(allclose(&res, &desired));
        assert_eq!(res, desired);
    }

    #[test]
    fn test_random_hermitian() {
        let d = 8;
        let mat = random_hermitian(d);
        assert!(allclose(&mat, &mat.dagger()));
    }

    // #[test]
    // fn test_random_su4() {
    //     let mat = random_su4();
    //     let id = Array2::eye(4);
    //     let prod = mat.dot(&mat.dagger());
    //     assert!(allclose(&prod, &id));
    // }

    // #[test]
    // fn test_random_unitary() {
    //     let d = 10;
    //     let mat = random_unitary(d);
    //     let id = Array2::eye(d);
    //     let prod = mat.dot(&mat.dagger());
    //     assert!(allclose(&prod, &id));
    // }

    #[test]
    fn test_controlled_gate() {
        let zero = r!(0.0);
        let one = r!(1.0);
        let u = array![
            [one, zero, zero, zero],
            [zero, one, zero, zero],
            [zero, zero, zero, one],
            [zero, zero, one, zero]
        ];
        let cx = controlled_unitary_matrix(&Gate::x().data, 1);

        println!("{}", cx);
        assert!(allclose(&cx, &u));
    }
}
