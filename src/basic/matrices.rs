// Extended matrix traits for ndarray
use crate::{c, i, r};
use ndarray::{array, s, Array};
use ndarray::{Array1, Array2, ArrayView2};
use ndarray_linalg::c64;
use ndarray_rand::rand_distr::Uniform;
use ndarray_rand::RandomExt;

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
pub fn kronecker_product(a: &Array2<c64>, b: &Array2<c64>) -> Array2<c64> {
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
}
