use crate::basic::gates::Gate;
use crate::basic::matrices::kronecker_product;
use crate::basic::matrices::{Dagger, Kronecker};
use crate::utils::functions::is_power_of_two;
use crate::{i, r};
use ndarray::{array, Array, Array1, Array2, ArrayBase};
use ndarray_linalg::c64;
use ndarray_linalg::QR;
use ndarray_rand::rand_distr::Uniform;
use ndarray_rand::RandomExt;
use std::f64::consts::PI;

pub fn allclose(data1: &Array2<c64>, data2: &Array2<c64>) -> bool {
    let diff = data1 - data2;
    let diff = diff.mapv(|x| x.norm());
    let tol = 1e-8;
    diff.iter().all(|&x| x < tol)
}

// use ndarray Array2<c64> to compute its matrix trace:
// pub fn trace(data: &Array2<c64>) -> c64 {
//     data.diag().sum()
// }

pub fn match_global_phase(a: Array2<c64>, b: Array2<c64>) -> (Array2<c64>, Array2<c64>) {
    // If the shapes are different or either is empty, return copies
    if a.shape() != b.shape() || a.len() == 0 {
        return (a.clone(), b.clone());
    }

    // Find the index k of the element with the largest magnitude in b
    let mut max_mag = 0.0;
    let mut k = (0, 0);
    let shape = b.shape();
    for i in 0..shape[0] {
        for j in 0..shape[1] {
            let mag = b[(i, j)].norm();
            if mag > max_mag {
                max_mag = mag;
                k = (i, j);
            }
        }
    }

    // Calculate the multiplier required to zero the phase of a complex number
    fn dephase(v: c64) -> c64 {
        let re = v.re;
        let im = v.im;
        if im == 0.0 {
            // 避免浮点误差: 实数轴上，如果 r<0 则乘 -1，否则乘 +1
            if re < 0.0 {
                r!(-1.0)
            } else {
                r!(1.0)
            }
        } else if re == 0.0 {
            // English: If i < 0, multiply by 1j; otherwise, multiply by -1j
            if im < 0.0 {
                i!(1.0)
            } else {
                i!(-1.0)
            }
        } else {
            // Multiply by e^{-i * arg(v)} to zero the phase
            let angle = im.atan2(re); // i.e., arg(v)
            i!(-angle).exp()
        }
    }

    // Calculate the multiplier required to zero the phase at index k
    let phase_factor_a = dephase(a[k]);
    let phase_factor_b = dephase(b[k]);

    // Apply the phase correction to the entire matrices a and b
    let a_prime = a.mapv(|val| val * phase_factor_a);
    let b_prime = b.mapv(|val| val * phase_factor_b);

    (a_prime, b_prime)
}

pub fn equiv_unitary(data1: &Array2<c64>, data2: &Array2<c64>) -> bool {
    let (data1, data2) = match_global_phase(data1.clone(), data2.clone());
    allclose(&data1, &data2)
}

pub fn random_hermitian(d: usize) -> Array2<c64> {
    let real_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
    let imag_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
    let mat: Array2<c64> = real_part.mapv(|re| r!(re)) + imag_part.mapv(|im| i!(im));
    (mat.clone() + mat.dagger()) / 2.0
}

fn sin_sampler(size: usize) -> Array1<f64> {
    // Customize the sampler using the probability density function of sin (theta)
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

pub fn random_su4() -> Array2<c64> {
    random_unitary(4)
}

pub fn random_unitary(d: usize) -> Array2<c64> {
    let real_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
    let imag_part: Array2<f64> = Array::random((d, d), Uniform::new(0.0, 1.0));
    let mat: Array2<c64> = real_part.mapv(|re| r!(re)) + imag_part.mapv(|im| i!(im));
    let (q, _r) = mat.qr().unwrap();
    q
}

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

// def kak_coefficients(U: np.ndarray) -> np.ndarray:
//     r"""Sort of different from weyl_coordinates(), herein we cohere to such a tetrahedron with coordinates:
//         (x, y, z) ~ e^{i \frac{\pi}{2}(x XX + y YY + z ZZ)}

//     0 ≤ abs(z2) ≤ y2 ≤ x2 ≤ π/4
//     """
//     if U.shape != (4, 4):
//         raise ValueError('U should be a 4*4 matrix')
//     d = U.shape[0]
//     if not np.allclose(U @ U.conj().T, np.identity(d)):
//         raise ValueError('U is not unitary')
//     return replace_close_to_zero_with_zero(cirq.kak_decomposition(U).interaction_coefficients)

// def weyl_coordinate(U: np.ndarray) -> np.ndarray:
//     r"""
//     Reduce the rotation angles to coord of Weyl chamber.

//     Herein we cohere to such a tetrahedron with coord:
//         (x, y, z) ~ e^{- i \frac{\pi}{2}(x XX + y YY + z ZZ)}
//         where (x, y, z) ∈ {1/2 ≥ x ≥ y ≥ z ≥ 0} ∪ {1/2 ≥ (1-x) ≥ y ≥ z ≥ 0}

//     According to the series of symmetry characteristics of Weyl gate.
//     """
//     coord = - 2 / np.pi * kak_coefficients(U)
//     x, y, z = replace_close_to_zero_with_zero(coord)  # there must be 1/2 >= -x >= -y >= |z|

//     if x < 0 and y < 0:
//         x, y = -x, -y
//     if x < 0 and y == 0:
//         x = -x
//     if x == 0 and y < 0:
//         y = -y

//     if z < 0:
//         if np.allclose(x, 0.5):
//             z = -z
//         else:
//             x, y, z = 1 - x, y, -z

//     if z == 0 and 0.5 < x < 1:
//         x = 1 - x

//     return np.array([x, y, z])

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::gates;
    use crate::basic::matrices::Real;
    use ndarray::Array;
    use ndarray_rand::rand_distr::Uniform;
    use ndarray_rand::RandomExt;

    #[test]
    fn test_equiv_unitary() {
        let angles = Array::random((3,), Uniform::new(0.0, 3.0));
        let theta1 = angles[0];
        let theta2 = angles[1];
        let theta3 = angles[2];
        let u3 = gates::Gate::u3(theta1, theta2, theta3);
        let mat1 = u3.data;
        let phase = i!(1.234).exp();
        let mat2 = &mat1 * phase;

        println!("mat1: {}", mat1);
        println!("mat2: {}", mat2);
        assert!(equiv_unitary(&mat1, &mat2));
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

    #[test]
    fn test_random_su4() {
        let mat = random_su4();
        let id = Array2::eye(4);
        let prod = mat.dot(&mat.dagger());
        assert!(allclose(&prod, &id));
    }

    #[test]
    fn test_random_unitary() {
        // SU(2)
        let mat = random_su2();
        let id = Array2::eye(2);
        let prod = mat.dot(&mat.dagger());
        assert!(allclose(&prod, &id));

        // SU(4)
        let mat = random_su4();
        let id = Array2::eye(4);
        let prod = mat.dot(&mat.dagger());
        assert!(allclose(&prod, &id));

        // SU(10)
        let d = 10;
        let mat = random_unitary(d);
        let id = Array2::eye(d);
        let prod = mat.dot(&mat.dagger());
        assert!(allclose(&prod, &id));
    }

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
