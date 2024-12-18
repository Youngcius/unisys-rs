// Extended matrix traits for ndarray
use crate::utils::functions::is_power_of_two;
use crate::{c, i, r};
use ndarray::s;
use ndarray::{Array2, ArrayView2};
use ndarray_linalg::c64;

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

// TODO: Implement controlled_matrix function
// pub fn controlled_matrix(
//     matrix: &Array2<c64>,
//     control_qubit: usize,
//     target_qubit: usize,
//     num_qubits: usize,
// ) -> Array2<c64> {
//     let dim = 1 << num_qubits;
//     let mut result = Array2::zeros((dim, dim));

//     for i in 0..dim {
//         let control_bit = (i >> control_qubit) & 1;
//         let target_bit = (i >> target_qubit) & 1;

//         let index = if control_bit == 1 { 1 } else { 0 };
//         let j = i ^ (index << target_qubit);

//         result[[i, j]] = matrix[[target_bit, target_bit]];
//     }

//     result
// }

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
    use ndarray::array;

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
    }
}
