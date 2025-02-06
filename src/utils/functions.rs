use crate::basic::matrices::Dagger;
use ndarray::Array2;
use ndarray_linalg::c64;
use ndarray_linalg::svd::SVD;
use ndarray_linalg::{Norm, Trace};

pub fn is_power_of_two(num: usize) -> bool {
    (num & (num - 1) == 0) && num != 0
}

pub fn infidelity(u: &Array2<c64>, v: &Array2<c64>) -> f64 {
    let d = u.dim().0;
    1.0 - u.dagger().dot(v).trace().unwrap().norm() / d as f64
}

pub fn matrix_spec_norm(data: &Array2<c64>) -> f64 {
    data.svd(false, false).unwrap().1[0]
}

pub fn spectral_distance(u: &Array2<c64>, v: &Array2<c64>) -> f64 {
    if u.dim() != v.dim() {
        panic!("u and v must have the same shape.");
    }
    matrix_spec_norm(&(u - v)) // by default, .norm() is alias for .norm_l2() in ndarray-linalg
}

pub fn average_case_error(u: &Array2<c64>, v: &Array2<c64>) -> f64 {
    if u.dim() != v.dim() {
        panic!("u and v must have the same shape.");
    }
    let d = u.dim().0;
    (u - v).norm() / (d as f64).sqrt() // by default, .norm() is Frobenius norm of matrix in ndarray-linalg
}
